"""``process.anti_surge_recycle`` - the steady recycle that holds a machine on its control line.

Spec: ``specs/models/process/anti_surge_recycle.toml``

The Python twin of ``crates/azoth-process/src/models/anti_surge_recycle.rs`` over
``crates/azoth-process/src/kernels/anti_surge_recycle.rs``, written to mirror them.

# A loop around a machine, not a machine

``AntiSurgeRecycleCalculator.solve`` runs the compressor once to find the natural operating
point, compares its inlet volume flow with the control line, and - where the natural point is
short - wires a cooled clone of the suction through a mixer and converges the recycle mass.
It adds no compressor physics: the surge flow and the distance to surge are the machine's own
read-outs through ``process.compressor``'s chart branch, which is run at every mixed state.

# The update is on the recycle mass, through the mixed inlet's volume per mass

Each pass divides the volume-flow error by the mixed inlet's volume per unit mass, applies the
damping, and floors the result at zero. **A converged pass breaks before updating**, so the
reported recycle mass is the one the machine last ran at; an exhausted cap reports the last
update instead, which the machine was not run at.

# The inactive branch is not a degenerate loop

A natural inlet at or above the control line reports ``recycle_active`` false, zero recycle,
zero passes and the natural distance to surge, with ``converged`` true - the class's own
``Result`` constructor for that branch.
"""

from __future__ import annotations

from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import AntiSurgeRecycleResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.process.reference import capacity
from azoth.process.reference.compressor import _feed_state
from azoth.process.reference.compressor import compressor as compressor_run
from azoth.process.reference.mixer import _route

#: ``Compressor.speed``'s initialiser, rpm.
DEFAULT_SPEED = 3000.0
#: ``AntiSurgeRecycleCalculator.surgeControlMargin``'s default.
DEFAULT_MARGIN = 0.05
#: ``recycleCoolerTemperature``'s default: 40 C.
DEFAULT_COOLER_TEMPERATURE = 313.15
#: ``AntiSurgeRecycleCalculator.damping``'s default.
DEFAULT_DAMPING = 0.8
#: ``AntiSurgeRecycleCalculator.tolerance``'s default, as a fraction of the control flow.
DEFAULT_TOLERANCE = 0.005
#: ``AntiSurgeRecycleCalculator.maxIterations``'s default.
DEFAULT_MAX_ITERATIONS = 20


def anti_surge_recycle(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    curve: list[dict[str, Q]],
    surge_curve: list[dict[str, Q]],
    speed: float | None = None,
    surge_control_margin: float | None = None,
    recycle_cooler_temperature: Q | None = None,
    damping: float | None = None,
    tolerance: float | None = None,
    max_iterations: int | None = None,
) -> AntiSurgeRecycleResult:
    """Converge the cooled recycle that holds a charted compressor on its control line.

    Args:
        components: the fluid's substances, by name.
        feed_n: the natural suction's molar flow.
        feed_z: the natural suction's composition.
        feed_p: the suction pressure.
        feed_t: the suction temperature.
        curve: the chart's measured points, in the compressor's own units.
        surge_curve: the surge boundary's measured points, in the same units.
        speed: the machine's speed, rpm; absent takes the class's 3000.
        surge_control_margin: the control line as a fraction above the surge flow; absent is
            0.05.
        recycle_cooler_temperature: the recycle clone's temperature; absent is 313.15 K.
        damping: the fixed-point relaxation; absent is 0.8.
        tolerance: the convergence gate, as a fraction of the control flow; absent is 0.005.
        max_iterations: the fixed-point cap; absent is 20.

    Returns:
        Whether a recycle was opened and converged, how many passes it took, and the recycle
        mass, inlet, surge and control flows and distance to surge at the final machine state.

    Raises:
        InvalidInputError: where a feed is not a composition, the machine reports no surge
            flow, or a bound the spec declares is violated.

    See :func:`azoth.process.reference.compressor.compressor` for the machine.
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    if len(feed_z) != len(components):
        raise InvalidInputError(
            "feed_z", f"a feed for {len(components)} components has {len(feed_z)} entries"
        )
    n = input_to_si(spec, "feed_n", feed_n)
    p = input_to_si(spec, "feed_p", feed_p)
    t = input_to_si(spec, "feed_t", feed_t)
    cooler = (
        DEFAULT_COOLER_TEMPERATURE
        if recycle_cooler_temperature is None
        else input_to_si(spec, "recycle_cooler_temperature", recycle_cooler_temperature)
    )
    rate = DEFAULT_SPEED if speed is None else float(speed)
    margin = DEFAULT_MARGIN if surge_control_margin is None else float(surge_control_margin)
    relax = DEFAULT_DAMPING if damping is None else float(damping)
    gate = DEFAULT_TOLERANCE if tolerance is None else float(tolerance)
    cap = DEFAULT_MAX_ITERATIONS if max_iterations is None else int(max_iterations)

    apply_checks(
        checks.on_input,
        {
            "feed_n": n,
            "feed_p": p,
            "feed_t": t,
            "speed": speed,
            "surge_control_margin": surge_control_margin,
            "recycle_cooler_temperature": cooler,
            "damping": damping,
            "tolerance": tolerance,
            "max_iterations": max_iterations,
        }.get,
        warnings,
    )

    mixture, _ = _components.mixture_of(components, eos="pr")
    molar_mass = capacity._molar_mass(mixture, list(feed_z))
    feed_volume = _volume_flow(components, n, feed_z, t, p)

    natural = _machine(components, n, feed_z, p, t, curve, surge_curve, rate, margin)
    surge = _m3hr(natural.surge_flow_rate, "surge_flow_rate")
    control = surge * (1.0 + margin)
    distance = _required(natural.distance_to_surge, "distance_to_surge")

    if not feed_volume < control:
        return AntiSurgeRecycleResult(
            recycle_active=False,
            converged=True,
            iterations=0,
            recycle_mass_flow=from_si(0.0, "kg/h"),
            inlet_volume_flow=from_si(feed_volume / 3600.0, "m**3/h"),
            surge_flow=from_si(surge / 3600.0, "m**3/h"),
            control_flow=from_si(control / 3600.0, "m**3/h"),
            distance_to_surge=distance,
            warnings=tuple(warnings),
        )

    # ``feedMass`` in the class: the natural suction's mass flow, kg/hr.
    feed_mass = n * molar_mass * 3600.0
    recycle_mass = max(
        (control - feed_volume) / max(feed_volume, 1.0e-6) * feed_mass,
        0.0,
    )
    inlet_volume = feed_volume
    converged = False
    iterations = 0
    for index in range(cap):
        iterations = index + 1
        recycle_n = recycle_mass / 3600.0 / molar_mass
        mixed = _route(
            components,
            [n, recycle_n],
            [list(feed_z), list(feed_z)],
            [p, p],
            [t, cooler],
            None,
        )
        mixed_t = mixed.product_t.to("K").magnitude
        mixed_p = mixed.pressure
        inlet_volume = _volume_flow(components, mixed.n_total, list(mixed.z), mixed_t, mixed_p)
        machine = _machine(
            components,
            mixed.n_total,
            list(mixed.z),
            mixed_p,
            mixed_t,
            curve,
            surge_curve,
            rate,
            margin,
        )
        surge = _m3hr(machine.surge_flow_rate, "surge_flow_rate")
        control = surge * (1.0 + margin)
        distance = _required(machine.distance_to_surge, "distance_to_surge")
        error = control - inlet_volume
        if abs(error) < gate * control:
            converged = True
            break
        total_mass = mixed.n_total * molar_mass * 3600.0
        volume_per_mass = inlet_volume / max(total_mass, 1.0e-6)
        recycle_mass = max(
            recycle_mass + relax * error / max(volume_per_mass, 1.0e-9),
            0.0,
        )

    return AntiSurgeRecycleResult(
        recycle_active=True,
        converged=converged,
        iterations=iterations,
        recycle_mass_flow=from_si(recycle_mass / 3600.0, "kg/h"),
        inlet_volume_flow=from_si(inlet_volume / 3600.0, "m**3/h"),
        surge_flow=from_si(surge / 3600.0, "m**3/h"),
        control_flow=from_si(control / 3600.0, "m**3/h"),
        distance_to_surge=distance,
        warnings=tuple(warnings),
    )


def _machine(
    components: list[str],
    n: float,
    z: list[float],
    p: float,
    t: float,
    curve: list[dict[str, Q]],
    surge_curve: list[dict[str, Q]],
    speed: float,
    margin: float,
) -> Any:
    """The charted machine at one inlet state, through ``process.compressor``.

    **The run-time recycle block stays off**, which is what the calculator does: it owns the
    recycle. The isentropic efficiency is stated only because the compressor's boundary
    requires one; its chart branch never reads it.
    """
    return compressor_run(
        components,
        from_si(n, "mol/s"),
        list(z),
        from_si(p, "Pa"),
        from_si(t, "K"),
        isentropic_efficiency=1.0,
        curve=curve,
        speed=speed,
        surge_curve=surge_curve,
        surge_control_margin=margin,
    )


def _volume_flow(components: list[str], n: float, z: list[float], t: float, p: float) -> float:
    """The stream's actual volumetric flow, m3/hr, at the feed's own root.

    The same reading ``process.compressor``'s chart branch scales its reduced flow from -
    ``getFlowRate("m3/hr")`` on the machine's inlet - so the loop compares one convention.
    """
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    _, _, molar_volume = _feed_state(mixture, ideal_gas, t, p, list(z))
    return n * molar_volume * 3600.0


def _m3hr(value: Q | None, field: str) -> float:
    """One compressor flow read-out in m3/hr, the unit the loop compares in."""
    if value is None:
        raise InvalidInputError(
            field, "the machine reports no surge flow, and a control line needs one"
        )
    return float(value.to("m**3/h").magnitude)


def _required(value: float | None, field: str) -> float:
    """A read-out the surge curve's own existence guarantees, refused rather than omitted."""
    if value is None:
        raise InvalidInputError(field, "the machine reports no distance to surge")
    return float(value)


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.anti_surge_recycle")
