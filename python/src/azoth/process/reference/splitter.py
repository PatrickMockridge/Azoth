"""``process.splitter`` - the splitter's kernel.

Spec: ``specs/models/process/splitter.toml``

The Python twin of ``crates/azoth-process/src/models/splitter.rs`` over
``crates/azoth-process/src/kernels/splitter.rs``, written to mirror them.

# Nothing changes but the flow

A splitter adds no heat and no work and does not flash. Every outlet therefore carries
the feed's composition, pressure, temperature and molar enthalpy, and the fractions are
the whole of the arithmetic. That is also why the outlet enthalpy is not something an
oracle can be asked for directly - see the divergence below.

# NeqSim's ``Splitter.run`` does not conserve enthalpy

It reaches its branch state through ``clone()``, then ``init(0)`` - which leaves a
``SystemPrEos`` holding *two* phase objects, each with the full inventory - then
``addComponent(int, double)``, which mutates every phase in ``phaseArray`` while
``setTotalNumberOfMolesRaw`` counts the change once. The ``TPflash`` after that converges
the active phase and leaves the stale one in place, and ``getEnthalpy()`` sums what it
finds there. Measured on the captured fluid: -20689.86 J/mol in, -93092.87 and -28142.87
out, weighting to -47627.87. This diverges deliberately and the spec says so.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SplitterResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.process.reference import capacity

#: `Splitter.designPressureDrop`, the class's `0.05` bar, in pascals.
DEFAULT_DESIGN_PRESSURE_DROP_PA = 5000.0
#: `Splitter.maxDesignVelocity`, m/s.
DEFAULT_MAX_DESIGN_VELOCITY_M_S = 30.0
#: `CapacityConstraint.MAX_UTILIZATION`.
MAX_UTILIZATION = 9.99


def splitter(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    split_factors: list[float],
    design_pressure_drop: Q | None = None,
    max_design_velocity: Q | None = None,
    header_diameter: Q | None = None,
    capacity_analysis_enabled: bool | None = None,
) -> SplitterResult:
    """Split a stream into several with the same state, in proportion to ``split_factors``.

    Args:
        components: the fluid's substances, by name.
        feed_n: the inlet's molar flow.
        feed_z: the inlet composition.
        feed_p: the inlet pressure.
        feed_t: the inlet temperature.
        split_factors: the fraction of the feed to each outlet. Normalised, so only the
            ratios matter and the vector's length is the outlet count.
        design_pressure_drop: ``setDesignPressureDrop``'s limit; absent is the class's
            0.05 bar, used when a setter initializes the map.
        max_design_velocity: ``setMaxDesignVelocity``'s limit; absent is the class's 30 m/s.
        header_diameter: the mechanical design's header diameter; absent is 0.0, unsized.
        capacity_analysis_enabled: ``splitterCapacityAnalysisEnabled``; absent is false.
            It gates initialization only - the capacity getters ignore it.

    Returns:
        Every outlet's record and the capacity map, with the per-constraint fields absent
        when that constraint is not in the map.

    Raises:
        InvalidInputError: where there is no outlet, or the factors do not sum to a
            positive value.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = splitter(
        ...     ["n-butane", "n-pentane"],
        ...     q(1.0, "mol/s"),
        ...     [0.6, 0.4],
        ...     q(10.0, "bar"),
        ...     q(300.0, "K"),
        ...     [0.3, 0.7],
        ... )
        >>> [round(n.to("mol/s").magnitude, 12) for n in r.products_n]
        [0.3, 0.7]
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "feed_n", feed_n)
    t = input_to_si(spec, "feed_t", feed_t)
    p = input_to_si(spec, "feed_p", feed_p)
    design_pressure_drop_si = (
        None
        if design_pressure_drop is None
        else input_to_si(spec, "design_pressure_drop", design_pressure_drop)
    )
    max_design_velocity_si = (
        None
        if max_design_velocity is None
        else input_to_si(spec, "max_design_velocity", max_design_velocity)
    )
    header_diameter_si = (
        None if header_diameter is None else input_to_si(spec, "header_diameter", header_diameter)
    )

    apply_checks(checks.on_input, {"feed_t": t}.get, warnings)

    if not split_factors:
        raise InvalidInputError("split_factors", "a splitter needs at least one outlet")

    total = sum(split_factors)
    if total <= 0.0:
        raise InvalidInputError("split_factors", "the split fractions must sum to a positive value")

    states = _route(components, t, p, feed_z, split_factors)
    capacity_map = _capacity(
        components,
        t,
        p,
        feed_z,
        n,
        design_pressure_drop_si,
        max_design_velocity_si,
        header_diameter_si,
        capacity_analysis_enabled is True,
    )
    z = tuple(feed_z)
    return SplitterResult(
        products_n=tuple(from_si(n * fraction, "mol/s") for fraction in states.fractions),
        products_z=tuple(z for _ in states.fractions),
        products_p=tuple(feed_p for _ in states.fractions),
        products_t=tuple(feed_t for _ in states.fractions),
        products_h=tuple(from_si(states.h_in, "J/mol") for _ in states.fractions),
        pressure_drop_value=_quantity(capacity_map.pressure_drop_value, "Pa"),
        pressure_drop_design=_quantity(capacity_map.pressure_drop_design, "Pa"),
        pressure_drop_utilization=capacity_map.pressure_drop_utilization,
        velocity_value=_quantity(capacity_map.velocity_value, "m/s"),
        velocity_design=_quantity(capacity_map.velocity_design, "m/s"),
        velocity_utilization=capacity_map.velocity_utilization,
        max_utilization=capacity_map.max_utilization,
        capacity_exceeded=capacity_map.capacity_exceeded,
        bottleneck_constraint=capacity_map.bottleneck_constraint,
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class SplitterCapacity:
    """The class's capacity map and its three summaries, in SI.

    The six per-constraint fields are ``None`` when that constraint is not in the map, which
    is the class's lazy initialization and not an absent answer.
    """

    pressure_drop_value: float | None
    pressure_drop_design: float | None
    pressure_drop_utilization: float | None
    velocity_value: float | None
    velocity_design: float | None
    velocity_utilization: float | None
    max_utilization: float
    capacity_exceeded: bool
    bottleneck_constraint: str | None


def _capacity(
    components: list[str],
    t: float,
    p: float,
    feed_z: list[float],
    n: float,
    design_pressure_drop: float | None,
    max_design_velocity: float | None,
    header_diameter: float | None,
    analysis_enabled: bool,
) -> SplitterCapacity:
    """``Splitter.splitterCapacityConstraints``, evaluated on the feed.

    **The map is lazy.** A supplied limit is a setter that was called and initializes the map;
    the enable flag initializes it only when it would otherwise be empty. A non-positive limit
    suppresses its own constraint but still initializes.
    """
    if design_pressure_drop is None and max_design_velocity is None and not analysis_enabled:
        return SplitterCapacity(None, None, None, None, None, None, 0.0, False, None)

    pressure_drop_design = (
        DEFAULT_DESIGN_PRESSURE_DROP_PA if design_pressure_drop is None else design_pressure_drop
    )
    velocity_design = (
        DEFAULT_MAX_DESIGN_VELOCITY_M_S if max_design_velocity is None else max_design_velocity
    )

    pressure_drop = None
    if pressure_drop_design > 0.0:
        pressure_drop = (0.0, pressure_drop_design, _utilization(0.0, pressure_drop_design))

    velocity = None
    if velocity_design > 0.0:
        diameter = 0.0 if header_diameter is None else header_diameter
        if diameter > 0.0:
            volume_flow = capacity.system_volumetric_flow(
                components, t, p, feed_z, n, from_feed=True
            )
            area = math.pi * diameter * diameter / 4.0
            value = volume_flow / area
        else:
            value = 0.0
        velocity = (value, velocity_design, _utilization(value, velocity_design))

    # `getBottleneckConstraint` keeps a constraint only when its utilization is strictly
    # above the running maximum, which starts at zero - so the always-zero pressure drop
    # never wins and a zero-utilization map has no bottleneck.
    max_utilization = 0.0
    bottleneck_constraint = None
    if pressure_drop is not None and pressure_drop[2] > max_utilization:
        max_utilization = pressure_drop[2]
        bottleneck_constraint = "pressureDrop"
    if velocity is not None and velocity[2] > max_utilization:
        max_utilization = velocity[2]
        bottleneck_constraint = "velocity"

    return SplitterCapacity(
        pressure_drop_value=None if pressure_drop is None else pressure_drop[0],
        pressure_drop_design=None if pressure_drop is None else pressure_drop[1],
        pressure_drop_utilization=None if pressure_drop is None else pressure_drop[2],
        velocity_value=None if velocity is None else velocity[0],
        velocity_design=None if velocity is None else velocity[1],
        velocity_utilization=None if velocity is None else velocity[2],
        max_utilization=max_utilization,
        capacity_exceeded=(pressure_drop is not None and pressure_drop[2] > 1.0)
        or (velocity is not None and velocity[2] > 1.0),
        bottleneck_constraint=bottleneck_constraint,
    )


def _utilization(value: float, design: float) -> float:
    """`CapacityConstraint.getUtilization` for one of the splitter's maximum constraints."""
    if design <= 0.0:
        return 0.0
    return min(value / design, MAX_UTILIZATION)


def _quantity(value: float | None, unit: str) -> Q | None:
    """A result field's optional SI magnitude, rebuilt in the spec's unit."""
    return None if value is None else from_si(value, unit)


@dataclass(frozen=True, slots=True)
class SplitStates:
    """The feed's own state and the share each outlet takes of it, in SI.

    The refusals are not here: they are about the *declared* factors, and
    :func:`splitter` states them before the fluid is resolved. What is here is what the
    split is - one enthalpy, and the normalised fractions that divide the flow.
    """

    h_in: float
    fractions: tuple[float, ...]


def _route(
    components: list[str],
    t: float,
    p: float,
    feed_z: list[float],
    split_factors: list[float],
) -> SplitStates:
    """The splitter's intermediates, in the order it computes them, in SI.

    **One arithmetic, two consumers.** :func:`splitter` builds the result from this and
    :mod:`azoth.process.layers` reports the same numbers against a NeqSim capture, so a
    layer the harness compares is a layer the model actually took. The feed's enthalpy is
    the whole interior: a split is a flow division and nothing else, which is exactly why
    the branch enthalpies the capture carries are the row where the two libraries part.
    """
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    # The feed's own state, which every outlet carries: `h` is a function of `(T, P, z)`
    # and none of the three moves.
    h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, feed_z)
    total = sum(split_factors)
    return SplitStates(h_in=h_in, fractions=tuple(factor / total for factor in split_factors))


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.splitter")
