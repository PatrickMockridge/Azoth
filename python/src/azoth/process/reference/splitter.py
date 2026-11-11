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
#: `SplitterMechanicalDesign.designVelocity`, m/s.
DEFAULT_DESIGN_VELOCITY_M_S = 15.0
#: `MechanicalDesign.maxOperationPressure`, the class's `100` bara, in pascals.
DEFAULT_MAX_OPERATION_PRESSURE_PA = 100.0e5
#: `MechanicalDesign.tensileStrength`, the class's `483` MPa, in pascals.
DEFAULT_TENSILE_STRENGTH_PA = 483.0e6
#: `MechanicalDesign.jointEfficiency`: read by the wall formula and then discarded.
DEFAULT_JOINT_EFFICIENCY = 1.0
#: `SplitterMechanicalDesign.calculateWeights`' steel density, kg/m3.
STEEL_DENSITY_KG_M3 = 7850.0
#: The class's wall-thickness floor, mm.
MIN_WALL_THICKNESS_MM = 3.0
#: The class's branch-spacing floor, m.
MIN_BRANCH_SPACING_M = 0.5
#: The class's assumed branch length, m.
BRANCH_LENGTH_M = 0.5
#: `SplitterMechanicalDesign.calcDesign`'s splitting loss coefficient.
SPLITTING_LOSS_COEFFICIENT = 0.5
#: The B31.3 temperature coefficient the class hard-sets.
WALL_Y_COEFFICIENT = 0.4
#: The class's standard pipe sizes, in metres, from 1/2" to 48".
STANDARD_PIPE_SIZES_M = (
    0.0127,
    0.0254,
    0.0381,
    0.0508,
    0.0635,
    0.0762,
    0.1016,
    0.1270,
    0.1524,
    0.2032,
    0.2540,
    0.3048,
    0.3556,
    0.4064,
    0.4572,
    0.5080,
    0.6096,
    0.7620,
    0.9144,
    1.0668,
    1.2192,
)


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
    mechanical_design: bool | None = None,
    design_velocity: Q | None = None,
    max_operation_pressure: Q | None = None,
    tensile_strength: Q | None = None,
    joint_efficiency: float | None = None,
    corrosion_allowance: Q | None = None,
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
        mechanical_design: run ``SplitterMechanicalDesign.calcDesign``; absent is false. The
            five design parameters below are refused without it, and a stated
            ``header_diameter`` is refused with it.
        design_velocity: the design's sizing velocity; absent is the class's 15 m/s.
        max_operation_pressure: the wall formula's pressure; absent is 100 bara.
        tensile_strength: the wall formula's tensile strength; absent is 483 MPa.
        joint_efficiency: read by the wall formula and discarded; absent is 1.0.
        corrosion_allowance: the wall thickness's corrosion allowance; absent is 0 mm.

    Returns:
        Every outlet's record, the capacity map and, with ``mechanical_design``, the sized
        header, branch, wall, manifold volume, pressure drop and weights. The per-constraint
        and design fields are absent when they were not asked for.

    Raises:
        InvalidInputError: where there is no outlet, the factors do not sum to a positive
            value, a design input is stated without ``mechanical_design``, a stated
            ``header_diameter`` meets the flag, or a design division has no answer.

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
    design_velocity_si = (
        None if design_velocity is None else input_to_si(spec, "design_velocity", design_velocity)
    )
    max_operation_pressure_si = (
        None
        if max_operation_pressure is None
        else input_to_si(spec, "max_operation_pressure", max_operation_pressure)
    )
    tensile_strength_si = (
        None
        if tensile_strength is None
        else input_to_si(spec, "tensile_strength", tensile_strength)
    )
    corrosion_allowance_si = (
        None
        if corrosion_allowance is None
        else input_to_si(spec, "corrosion_allowance", corrosion_allowance)
    )

    apply_checks(
        checks.on_input,
        {"feed_t": t, "design_velocity": design_velocity_si}.get,
        warnings,
    )

    if not split_factors:
        raise InvalidInputError("split_factors", "a splitter needs at least one outlet")

    total = sum(split_factors)
    if total <= 0.0:
        raise InvalidInputError("split_factors", "the split fractions must sum to a positive value")

    states = _route(components, t, p, feed_z, split_factors)
    design = _design(
        components,
        t,
        p,
        feed_z,
        n,
        split_factors,
        header_diameter_si,
        mechanical_design is True,
        design_velocity_si,
        max_operation_pressure_si,
        tensile_strength_si,
        joint_efficiency,
        corrosion_allowance_si,
    )
    # With the design on, the capacity map reads the sized diameter live; with it off, the
    # stated input stands.
    effective_header_diameter = header_diameter_si if design is None else design.header_diameter
    capacity_map = _capacity(
        components,
        t,
        p,
        feed_z,
        n,
        design_pressure_drop_si,
        max_design_velocity_si,
        effective_header_diameter,
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
        design_header_diameter=(None if design is None else _quantity(design.header_diameter, "m")),
        design_branch_diameter=(None if design is None else _quantity(design.branch_diameter, "m")),
        design_header_length=(None if design is None else _quantity(design.header_length, "m")),
        design_wall_thickness=(None if design is None else _quantity(design.wall_thickness, "m")),
        manifold_volume=(None if design is None else _quantity(design.manifold_volume, "m**3")),
        design_total_pressure_drop=(
            None if design is None else _quantity(design.total_pressure_drop, "Pa")
        ),
        header_weight=(None if design is None else _quantity(design.header_weight, "kg")),
        branch_weight=(None if design is None else _quantity(design.branch_weight, "kg")),
        total_weight=(None if design is None else _quantity(design.total_weight, "kg")),
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


@dataclass(frozen=True, slots=True)
class SplitterDesign:
    """`SplitterMechanicalDesign.calcDesign`'s answer, in SI (m, m3, kg, Pa)."""

    header_diameter: float
    branch_diameter: float
    header_length: float
    wall_thickness: float
    manifold_volume: float
    total_pressure_drop: float
    header_weight: float
    branch_weight: float
    total_weight: float


def _design(
    components: list[str],
    t: float,
    p: float,
    feed_z: list[float],
    n: float,
    split_factors: list[float],
    header_diameter: float | None,
    mechanical_design: bool,
    design_velocity: float | None,
    max_operation_pressure: float | None,
    tensile_strength: float | None,
    joint_efficiency: float | None,
    corrosion_allowance: float | None,
) -> SplitterDesign | None:
    """`SplitterMechanicalDesign.calcDesign`, with the class's lazy call and its guards.

    **The class does not call this from `run`.** The flag is the caller's explicit request;
    with it off, any stated design input is refused by name. The formula keeps the class's
    units: `P` in bar, `S` = 40 % of the tensile strength in MPa, and `joint_efficiency`
    read and discarded.
    """
    if not mechanical_design:
        for name, stated in (
            ("design_velocity", design_velocity is not None),
            ("max_operation_pressure", max_operation_pressure is not None),
            ("tensile_strength", tensile_strength is not None),
            ("joint_efficiency", joint_efficiency is not None),
            ("corrosion_allowance", corrosion_allowance is not None),
        ):
            if stated:
                raise InvalidInputError(name, f"`{name}` was stated without `mechanical_design`")
        return None
    if header_diameter is not None:
        raise InvalidInputError(
            "header_diameter",
            "`mechanical_design` sizes the header; state one or the other, not both",
        )

    total = sum(split_factors)
    if not split_factors or total <= 0.0:
        raise InvalidInputError(
            "split_factors", "a splitter needs at least one outlet whose factors sum above zero"
        )

    design_velocity_value = (
        DEFAULT_DESIGN_VELOCITY_M_S if design_velocity is None else design_velocity
    )
    max_operation_pressure_value = (
        DEFAULT_MAX_OPERATION_PRESSURE_PA
        if max_operation_pressure is None
        else max_operation_pressure
    )
    tensile_strength_mpa = (
        DEFAULT_TENSILE_STRENGTH_PA if tensile_strength is None else tensile_strength
    ) / 1.0e6
    joint_efficiency_value = (
        DEFAULT_JOINT_EFFICIENCY if joint_efficiency is None else joint_efficiency
    )
    corrosion_allowance_mm = (0.0 if corrosion_allowance is None else corrosion_allowance) * 1000.0

    from azoth.process.reference import capacity as capacity_module

    volume_flow = capacity_module.system_volumetric_flow(
        components, t, p, feed_z, n, from_feed=True
    )
    if not math.isfinite(volume_flow) or volume_flow <= 0.0:
        # `calcDesign` returns before it divides or sizes anything; the class's fields stay at
        # their initial zeros, and they are answers, not absent ones.
        return SplitterDesign(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    if not design_velocity_value > 0.0 or not math.isfinite(design_velocity_value):
        raise InvalidInputError(
            "design_velocity",
            "the header area divides by the design velocity, so it must be positive",
        )

    design_pressure_bar = max_operation_pressure_value / 1.0e5 * 1.1
    allowable_stress_mpa = tensile_strength_mpa * 0.4
    denominator = 2.0 * allowable_stress_mpa + 2.0 * WALL_Y_COEFFICIENT * design_pressure_bar
    if not denominator > 0.0 or not math.isfinite(denominator):
        raise InvalidInputError(
            "tensile_strength",
            "the wall-thickness denominator is zero when the tensile strength and the design "
            "pressure are both zero",
        )
    # Read and discarded, as the class does.
    _ = joint_efficiency_value

    branches = float(len(split_factors))
    max_ratio = max(factor / total for factor in split_factors)
    header_area = volume_flow / design_velocity_value
    header_diameter_value = _round_to_standard_pipe_size(math.sqrt(4.0 * header_area / math.pi))
    branch_area = volume_flow * max_ratio / design_velocity_value
    branch_diameter = _round_to_standard_pipe_size(math.sqrt(4.0 * branch_area / math.pi))

    branch_spacing = max(2.0 * branch_diameter, MIN_BRANCH_SPACING_M)
    header_length = (branches + 1.0) * branch_spacing
    wall_thickness_mm = max(
        design_pressure_bar * header_diameter_value * 1000.0 / denominator + corrosion_allowance_mm,
        MIN_WALL_THICKNESS_MM,
    )

    header_section = math.pi * (header_diameter_value / 2.0) ** 2
    manifold_volume = header_section * header_length
    velocity = volume_flow / header_section
    density = capacity_module.system_mass_density(components, t, p, feed_z, from_feed=True)
    total_pressure_drop = SPLITTING_LOSS_COEFFICIENT * 0.5 * density * velocity * velocity

    header_outer = header_diameter_value + 2.0 * wall_thickness_mm / 1000.0
    header_weight = (
        math.pi
        / 4.0
        * (header_outer**2 - header_diameter_value**2)
        * header_length
        * STEEL_DENSITY_KG_M3
    )
    branch_outer = branch_diameter + 2.0 * wall_thickness_mm / 1000.0
    branch_weight = (
        math.pi
        / 4.0
        * (branch_outer**2 - branch_diameter**2)
        * BRANCH_LENGTH_M
        * STEEL_DENSITY_KG_M3
        * branches
    )

    return SplitterDesign(
        header_diameter=header_diameter_value,
        branch_diameter=branch_diameter,
        header_length=header_length,
        wall_thickness=wall_thickness_mm / 1000.0,
        manifold_volume=manifold_volume,
        total_pressure_drop=total_pressure_drop,
        header_weight=header_weight,
        branch_weight=branch_weight,
        total_weight=header_weight + branch_weight,
    )


def _round_to_standard_pipe_size(diameter: float) -> float:
    """`SplitterMechanicalDesign.roundToStandardPipeSize`: the first size at or above."""
    for size in STANDARD_PIPE_SIZES_M:
        if size >= diameter:
            return size
    return STANDARD_PIPE_SIZES_M[-1]


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
