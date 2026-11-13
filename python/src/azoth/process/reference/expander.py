"""``process.expander`` - the isentropic route, multiplying by the efficiency.

Spec: ``specs/models/process/expander.toml``

The Python twin of ``crates/azoth-process/src/models/expander.rs`` over
``crates/azoth-process/src/kernels/expander.rs``, written to mirror them.

# The step is isentropic, and the inlet's entropy is derived

NeqSim reads ``getEntropy()`` off the inlet system, sets the outlet pressure, runs a
``PSflash`` on that entropy and takes the enthalpy it lands on as the reversible outlet. The
record has no entropy field - ``s = f(T, P, z)``, the same flash that gives ``h`` - so a port
that agrees with the captured efficiency-of-one rows has computed it rather than carried it.

# The efficiency multiplies the isentropic step

``h_out = h_in + (h_isentropic - h_in) * eta``, so an efficiency below one recovers
less work than the reversible step. That is ``Compressor.run``'s rule with the efficiency on
the other side of the division, and it is the same physical rule: an expansion's difference is
negative, so dividing would make the machine beat the reversible one.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ExpanderResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve
from azoth.eos.reference.ps_flash import entropy_at
from azoth.eos.reference.ps_flash import ps_flash as ps_flash_solve
from azoth.process.reference.expander_mechanical import design as design_expander


def expander(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_pressure: Q,
    isentropic_efficiency: float,
    use_polytropic_calc: bool | None = None,
    polytropic_efficiency: float | None = None,
    number_of_steps: int | None = None,
    rated_recovered_power: Q | None = None,
    shaft_power: Q | None = None,
    outlet_temperature: Q | None = None,
    mechanical_design: bool | None = None,
) -> ExpanderResult:
    """Drop a stream's pressure along an isentrope, or in polytropic steps.

    Args:
        components: the fluid's substances, by name.
        inlet_n: the inlet's molar flow.
        inlet_z: the inlet composition.
        inlet_p: the inlet pressure.
        inlet_t: the inlet temperature.
        outlet_pressure: the pressure the machine takes the stream to.
        isentropic_efficiency: the machine's isentropic efficiency, in ``(0, 1]``.
        use_polytropic_calc: whether the stepped branch runs instead.
        polytropic_efficiency: the stepped branch's efficiency, in ``(0, 1]``.
        number_of_steps: the stepped branch's pressure increments.
        rated_recovered_power: ``setRatedRecoveredPower``; absent is 0, which creates no
            constraint, and so does a non-positive value.
        shaft_power: the total shaft power the run applies instead of an efficiency.
        outlet_temperature: the outlet temperature the efficiency is solved for.
        mechanical_design: run ``ExpanderMechanicalDesign.calcDesign`` and publish its
            API 617 design; absent is false.

    Returns:
        The outlet's record, the recovered power (negative) and the recovered-power map.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = expander(
        ...     ["methane", "n-butane"],
        ...     1.0,
        ...     q([0.9, 0.1]),
        ...     q(60.0, "bar"),
        ...     q(320.0, "K"),
        ...     q(30.0, "bar"),
        ...     0.75,
        ... )
        >>> round(r.outlet_h.to("J/mol").magnitude, 3)
        -618.933
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "inlet_n", inlet_n)
    t = input_to_si(spec, "inlet_t", inlet_t)
    p = input_to_si(spec, "inlet_p", inlet_p)
    p_out = input_to_si(spec, "outlet_pressure", outlet_pressure)
    rated = (
        None
        if rated_recovered_power is None
        else input_to_si(spec, "rated_recovered_power", rated_recovered_power)
    )
    shaft = None if shaft_power is None else input_to_si(spec, "shaft_power", shaft_power)
    target = (
        None
        if outlet_temperature is None
        else input_to_si(spec, "outlet_temperature", outlet_temperature)
    )

    apply_checks(
        checks.on_input,
        {
            "inlet_t": t,
            "isentropic_efficiency": isentropic_efficiency,
            "polytropic_efficiency": polytropic_efficiency,
        }.get,
        warnings,
    )

    if not 0.0 < isentropic_efficiency <= 1.0:
        raise InvalidInputError(
            "isentropic_efficiency",
            f"the machine's isentropic efficiency must be in (0, 1], and "
            f"{isentropic_efficiency} is not",
        )

    if shaft is not None and target is not None:
        raise InvalidInputError(
            "shaft_power",
            "a stated power and an outlet temperature are two modes; the class's "
            "`setPower` and `solveEfficiency` are one at a time",
        )
    if (shaft is not None or target is not None) and use_polytropic_calc:
        raise InvalidInputError(
            "shaft_power",
            "the stated-power and outlet-temperature modes are the isentropic branch's",
        )
    for name, value in (("shaft_power", shaft), ("outlet_temperature", target)):
        if value is not None and not math.isfinite(value):
            raise InvalidInputError(name, f"{value} is not a finite number")

    efficiency = 1.0 if polytropic_efficiency is None else float(polytropic_efficiency)
    reported: float | None = None
    if shaft is not None:
        h_in, h_out, outlet_t, reported = _at_power(components, t, p, inlet_z, p_out, n, shaft)
    elif target is not None:
        h_in, h_out, outlet_t, reported = _solve_efficiency(
            components, t, p, inlet_z, p_out, target, isentropic_efficiency
        )
    elif use_polytropic_calc:
        h_in, h_out, outlet_t = _polytropic_states(
            components,
            inlet_z,
            t,
            p,
            p_out,
            efficiency,
            5 if number_of_steps is None else int(number_of_steps),
        )
    else:
        states = _states(components, t, p, inlet_z, p_out, isentropic_efficiency)
        h_in, h_out, outlet_t = states.h_in, states.h_out, states.outlet_t

    power = n * (h_out - h_in)
    design = (
        None
        if not mechanical_design
        else design_expander(
            components,
            n,
            inlet_z,
            t,
            p,
            outlet_t.to("K").magnitude,
            power,
        )
    )
    capacity = _recovered_power_capacity(power, rated)
    return ExpanderResult(
        outlet_n=from_si(n, "mol/s"),
        outlet_z=tuple(inlet_z),
        outlet_p=outlet_pressure,
        outlet_t=outlet_t,
        outlet_h=from_si(h_out, "J/mol"),
        power=from_si(power, "W"),
        recovered_power_utilization=capacity.recovered_power_utilization,
        max_utilization=capacity.max_utilization,
        capacity_exceeded=capacity.capacity_exceeded,
        hard_limit_exceeded=capacity.hard_limit_exceeded,
        bottleneck_constraint=capacity.bottleneck_constraint,
        reported_isentropic_efficiency=reported,
        design_expander_type=None if design is None else design.expander_type,
        design_load_type=None if design is None else design.load_type,
        design_number_of_stages=None if design is None else float(design.number_of_stages),
        design_wheel_diameter=(None if design is None else from_si(design.wheel_diameter, "m")),
        design_tip_speed=None if design is None else from_si(design.tip_speed, "m/s"),
        design_shaft_diameter=(None if design is None else from_si(design.shaft_diameter, "m")),
        design_rated_speed=None if design is None else design.rated_speed,
        design_max_continuous_speed=None if design is None else design.max_continuous_speed,
        design_trip_speed=None if design is None else design.trip_speed,
        design_first_critical_speed=None if design is None else design.first_critical_speed,
        design_casing_pressure=(None if design is None else from_si(design.casing_pressure, "Pa")),
        design_casing_temperature=(
            None if design is None else from_si(design.casing_temperature, "K")
        ),
        design_casing_wall_thickness=(
            None if design is None else from_si(design.casing_wall_thickness, "m")
        ),
        design_bearing_span=None if design is None else from_si(design.bearing_span, "m"),
        design_wheel_weight=None if design is None else from_si(design.wheel_weight, "kg"),
        design_shaft_weight=None if design is None else from_si(design.shaft_weight, "kg"),
        design_generator_weight=(
            None if design is None else from_si(design.generator_weight, "kg")
        ),
        design_seal_type=None if design is None else design.seal_type,
        design_bearing_type=None if design is None else design.bearing_type,
        design_inner_diameter=(None if design is None else from_si(design.inner_diameter, "m")),
        design_outer_diameter=(None if design is None else from_si(design.outer_diameter, "m")),
        design_tantan_length=None if design is None else from_si(design.tantan_length, "m"),
        design_weight_vessel=None if design is None else from_si(design.weight_vessel, "kg"),
        design_weight_internals=(
            None if design is None else from_si(design.weight_internals, "kg")
        ),
        design_weight_piping=None if design is None else from_si(design.weight_piping, "kg"),
        design_weight_electro_instrument=(
            None if design is None else from_si(design.weight_electro_instrument, "kg")
        ),
        design_weight_structural_steel=(
            None if design is None else from_si(design.weight_structural_steel, "kg")
        ),
        design_weight_total=None if design is None else from_si(design.weight_total, "kg"),
        design_module_length=(None if design is None else from_si(design.module_length, "m")),
        design_module_width=None if design is None else from_si(design.module_width, "m"),
        design_module_height=None if design is None else from_si(design.module_height, "m"),
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class ExpanderStates:
    """The machine's intermediates, in the order it computes them, in SI."""

    #: The inlet's molar enthalpy and entropy, J/mol and J/(mol*K).
    h_in: float
    s_in: float
    #: The reversible outlet: the temperature at the outlet pressure on the inlet's entropy,
    #: and its molar enthalpy.
    t_isentropic: Q
    h_isentropic: float
    #: The isentropic step, and the step the efficiency actually asks for.
    dh_isentropic: float
    dh_actual: float
    #: The outlet: its molar enthalpy, its temperature, and its molar entropy - which rises
    #: across an irreversible step and is unchanged across a reversible one.
    h_out: float
    outlet_t: Q
    s_out: float


def _states(
    components: list[str],
    t: float,
    p: float,
    inlet_z: list[float],
    p_out: float,
    isentropic_efficiency: float,
) -> ExpanderStates:
    """The machine's intermediates, in SI, in the order it computes them."""
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")

    h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, inlet_z)
    entropy, _ = entropy_at(mixture, ideal_gas, t, p, inlet_z)

    isentropic = ps_flash_solve(
        mixture, ideal_gas, from_si(p_out, "Pa"), from_si(entropy, "J/(mol*K)"), inlet_z
    )
    h_isentropic, _ = enthalpy_at(
        mixture, ideal_gas, isentropic.T.to("K").magnitude, p_out, inlet_z
    )

    dh_isentropic = h_isentropic - h_in
    dh_actual = dh_isentropic * isentropic_efficiency
    h_out = h_in + dh_actual
    outlet_t = ph_flash_solve(
        mixture, ideal_gas, from_si(p_out, "Pa"), from_si(h_out, "J/mol"), inlet_z
    )

    s_out, _ = entropy_at(mixture, ideal_gas, outlet_t.T.to("K").magnitude, p_out, inlet_z)

    return ExpanderStates(
        h_in=h_in,
        s_in=entropy,
        t_isentropic=isentropic.T,
        h_isentropic=h_isentropic,
        dh_isentropic=dh_isentropic,
        dh_actual=dh_actual,
        h_out=h_out,
        outlet_t=outlet_t.T,
        s_out=s_out,
    )


def _at_power(
    components: list[str],
    t: float,
    p: float,
    inlet_z: list[float],
    p_out: float,
    n: float,
    power: float,
) -> tuple[float, float, Q, float | None]:
    """The stated-power route: the enthalpy step is `power/n`, and the efficiency divides it."""
    if not n > 0.0 or not math.isfinite(n):
        raise InvalidInputError(
            "shaft_power",
            f"the stated power divides by the molar flow, and {n} mol/s is not positive",
        )
    reference = _states(components, t, p, inlet_z, p_out, 1.0)
    dh = power / n
    h_out = reference.h_in + dh
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    outlet_t = ph_flash_solve(
        mixture, ideal_gas, from_si(p_out, "Pa"), from_si(h_out, "J/mol"), inlet_z
    ).T
    efficiency = None if reference.dh_isentropic == 0.0 else dh / reference.dh_isentropic
    return reference.h_in, h_out, outlet_t, efficiency


def _solve_efficiency(
    components: list[str],
    t: float,
    p: float,
    inlet_z: list[float],
    p_out: float,
    target: float,
    start_efficiency: float,
) -> tuple[float, float, Q, float]:
    """`solveEfficiency`: the class's secant, its first step and its clamps."""
    efficiency = (
        min(max(start_efficiency, 1.0e-6), 1.0) if math.isfinite(start_efficiency) else 1.0e-6
    )
    states = _states(components, t, p, inlet_z, p_out, efficiency)
    function_old = 0.0
    delta_poly = 100.0
    iteration = 0
    while True:
        iteration += 1
        function = states.outlet_t.to("K").magnitude - target
        slope = (
            (function - function_old) / delta_poly
            if delta_poly != 0.0
            and abs(delta_poly) > 1.0e-12
            and math.isfinite(delta_poly)
            and math.isfinite(function)
            and math.isfinite(function_old)
            else math.nan
        )
        candidate = (
            efficiency + (0.01 if function > 0.0 else -0.01)
            if not math.isfinite(slope) or abs(slope) < 1.0e-12
            else efficiency - function / slope
        )
        if iteration <= 1:
            candidate = efficiency + 0.01
        if not math.isfinite(candidate):
            candidate = efficiency + (0.01 if function > 0.0 else -0.01)
        old = efficiency
        efficiency = min(max(candidate, 1.0e-6), 1.0)
        delta_poly = efficiency - old
        if abs(delta_poly) < 1.0e-12:
            delta_poly = 1.0e-12 if function > 0.0 else -1.0e-12
        function_old = function
        states = _states(components, t, p, inlet_z, p_out, efficiency)
        residual = states.outlet_t.to("K").magnitude - target
        if not (abs(residual) > 0.01 or iteration < 3) or iteration >= 50:
            return states.h_in, states.h_out, states.outlet_t, efficiency


@dataclass(frozen=True, slots=True)
class ExpanderCapacity:
    """The `recoveredPower` constraint and its summaries, in SI."""

    recovered_power_utilization: float | None
    max_utilization: float
    capacity_exceeded: bool
    hard_limit_exceeded: bool
    bottleneck_constraint: str | None


def _recovered_power_capacity(power: float, rated: float | None) -> ExpanderCapacity:
    """`Expander`'s recovered-power constraint, evaluated on the shaft power.

    The parent's speed, surge and stonewall constraints are disabled without a chart, so this
    is the only enabled entry. The value is `abs(power)`, and a non-finite power answers 0.0,
    as the class's supplier does. A non-positive rating creates no constraint at all.
    """
    if rated is None or not rated > 0.0:
        return ExpanderCapacity(None, 0.0, False, False, None)
    recovered = abs(power) if math.isfinite(power) else 0.0
    utilization = min(recovered / rated, 9.99)
    return ExpanderCapacity(
        recovered_power_utilization=utilization,
        max_utilization=utilization if utilization > 0.0 else 0.0,
        capacity_exceeded=utilization > 1.0,
        hard_limit_exceeded=recovered > 1.1 * rated,
        bottleneck_constraint="recoveredPower" if utilization > 0.0 else None,
    )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.expander")


def _polytropic_states(
    components: list[str],
    inlet_z: list[float],
    t: float,
    p: float,
    p_out: float,
    polytropic_efficiency: float,
    number_of_steps: int,
) -> tuple[float, float, Q]:
    """The stepped branch: ``Expander.run``'s ``numbersteps`` increments of ``dp``.

    Each step takes the current entropy and enthalpy, flashes isentropically to the next
    pressure, and multiplies the step by the efficiency before a ``PHflash`` lands the state.
    The count is floored at one, as the class's ``Math.max(1, ...)`` does.
    """
    if not 0.0 < polytropic_efficiency <= 1.0:
        raise InvalidInputError(
            "polytropic_efficiency",
            "the expander's polytropic efficiency must be in (0, 1]",
        )
    steps = max(1, number_of_steps)
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    increment = (p_out - p) / steps
    state_t, state_p = t, p
    h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, inlet_z)
    h = h_in
    for _ in range(steps):
        entropy, _ = entropy_at(mixture, ideal_gas, state_t, state_p, inlet_z)
        next_p = state_p + increment
        isentropic = ps_flash_solve(
            mixture,
            ideal_gas,
            from_si(next_p, "Pa"),
            from_si(entropy, "J/(mol*K)"),
            inlet_z,
        )
        h_isentropic, _ = enthalpy_at(
            mixture, ideal_gas, isentropic.T.to("K").magnitude, next_p, inlet_z
        )
        h = h + (h_isentropic - h) * polytropic_efficiency
        flashed = ph_flash_solve(
            mixture, ideal_gas, from_si(next_p, "Pa"), from_si(h, "J/mol"), inlet_z
        )
        state_t = flashed.T.to("K").magnitude
        state_p = next_p
    return h_in, h, flashed.T
