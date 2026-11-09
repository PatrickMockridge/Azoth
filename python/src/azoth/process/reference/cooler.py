"""``process.cooler`` - ``Heater.run`` reached through ``Cooler``.

Spec: ``specs/models/process/cooler.toml``

The Python twin of ``crates/azoth-process/src/models/cooler.rs`` over
``crates/azoth-process/src/kernels/cooler.rs``, which is ``kernels/heater.rs``.

# The delegation is the port

``Cooler extends Heater`` and overrides ``runTransient``, the mechanical design and a
handful of getters - **not ``run``**, which is the whole of the steady-state arithmetic. Its
own 258 lines are a bounded dynamic temperature-control model that reaches the outlet only
through ``runTransient``.

So the arithmetic here is :func:`azoth.process.reference.heater._states`, called rather than
restated. Writing it out again would be a second implementation of one machine's physics, and
the two could then disagree - which is the failure the two-kernel rule exists to catch, and
not the failure a deliberate copy would avoid.

# What that claim rests on

The probe runs the same six rows through ``Heater`` and through ``Cooler``, and
``captures/process_heater.tsv`` and ``captures/process_cooler.tsv`` are byte-identical.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import CoolerResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.process.reference.heater import _states


def cooler(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_temperature: Q | None = None,
    duty: Q | None = None,
    pressure_drop: Q | None = None,
    time_step: Q | None = None,
    dynamic_temperature_control_enabled: bool | None = None,
    design_mass_flow: Q | None = None,
    design_inlet_temperature: Q | None = None,
    design_outlet_temperature: Q | None = None,
    cooling_medium_temperature: Q | None = None,
    design_valve_opening: float | None = None,
    actuator_time_constant: Q | None = None,
    thermal_time_constant: Q | None = None,
    cooling_valve_opening: float | None = None,
    requested_valve_opening: float | None = None,
    dynamic_outlet_temperature: Q | None = None,
) -> CoolerResult:
    """Cool a stream to a stated temperature, or by a stated duty.

    Args:
        components: the fluid's substances, by name.
        inlet_n: the inlet's molar flow.
        inlet_z: the inlet composition.
        inlet_p: the inlet pressure.
        inlet_t: the inlet temperature.
        outlet_temperature: the temperature to reach, or ``None`` for the enthalpy route.
        duty: the heat removed, as a negative number, or ``None`` to state a temperature.
        pressure_drop: the pressure loss across the cooler, or ``None`` for none.
        time_step: one ``runTransient`` step, or ``None`` for the steady run.
        dynamic_temperature_control_enabled: the configured dynamic branch's gate.
        design_mass_flow: the design point's process mass flow.
        design_inlet_temperature: the design point's inlet temperature.
        design_outlet_temperature: the design point's outlet temperature.
        cooling_medium_temperature: the cooling medium's supply temperature.
        design_valve_opening: the valve opening at the design point, percent.
        actuator_time_constant: the valve actuator's first-order constant.
        thermal_time_constant: the process side's first-order constant.
        cooling_valve_opening: the carried valve opening, percent.
        requested_valve_opening: the controller's response for this step, percent.
        dynamic_outlet_temperature: the carried outlet temperature.

    Returns:
        The outlet's record and the duty it moved.

    Raises:
        InvalidInputError: if a temperature and a duty are both given, if the pressure drop
            leaves a non-positive pressure, or if a duty is asked of a stream carrying no
            flow.
        OutOfRangeError: for a temperature the spec bounds.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = cooler(
        ...     ["methane", "n-butane"],
        ...     1.0,
        ...     q([0.9, 0.1]),
        ...     q(30.0, "bar"),
        ...     q(320.0, "K"),
        ...     duty=q(-5000.0, "W"),
        ... )
        >>> round(r.outlet_t.to("K").magnitude, 3)
        248.068
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "inlet_n", inlet_n)
    t = input_to_si(spec, "inlet_t", inlet_t)
    p = input_to_si(spec, "inlet_p", inlet_p)
    # The Rust boundary refuses a composition that does not pair with the component list;
    # the transient mass-flow sum must not be the first thing to notice it.
    if len(inlet_z) != len(components):
        raise InvalidInputError(
            "z",
            f"a feed for {len(components)} components has {len(inlet_z)} entries",
        )
    t_out = (
        None
        if outlet_temperature is None
        else input_to_si(spec, "outlet_temperature", outlet_temperature)
    )
    q = None if duty is None else input_to_si(spec, "duty", duty)
    drop = None if pressure_drop is None else input_to_si(spec, "pressure_drop", pressure_drop)
    dt = None if time_step is None else input_to_si(spec, "time_step", time_step)
    design_mass = (
        None
        if design_mass_flow is None
        else input_to_si(spec, "design_mass_flow", design_mass_flow)
    )
    design_in = (
        None
        if design_inlet_temperature is None
        else input_to_si(spec, "design_inlet_temperature", design_inlet_temperature)
    )
    design_out = (
        None
        if design_outlet_temperature is None
        else input_to_si(spec, "design_outlet_temperature", design_outlet_temperature)
    )
    medium = (
        None
        if cooling_medium_temperature is None
        else input_to_si(spec, "cooling_medium_temperature", cooling_medium_temperature)
    )
    # **A dimensionless input is the number it is**, as the public wrapper's `float`
    # annotation says; `input_to_si` expects a quantity and refuses a bare number.
    design_opening_si = design_valve_opening
    tau_a = (
        None
        if actuator_time_constant is None
        else input_to_si(spec, "actuator_time_constant", actuator_time_constant)
    )
    tau_t = (
        None
        if thermal_time_constant is None
        else input_to_si(spec, "thermal_time_constant", thermal_time_constant)
    )
    carried_opening_si = cooling_valve_opening
    requested_si = requested_valve_opening
    dynamic_out = (
        None
        if dynamic_outlet_temperature is None
        else input_to_si(spec, "dynamic_outlet_temperature", dynamic_outlet_temperature)
    )

    apply_checks(
        checks.on_input,
        {
            "inlet_t": t,
            "outlet_temperature": t_out,
            "time_step": dt,
            "design_mass_flow": design_mass,
            "design_inlet_temperature": design_in,
            "design_outlet_temperature": design_out,
            "cooling_medium_temperature": medium,
            "design_valve_opening": design_opening_si,
            "actuator_time_constant": tau_a,
            "thermal_time_constant": tau_t,
            "cooling_valve_opening": carried_opening_si,
            "requested_valve_opening": requested_si,
            "dynamic_outlet_temperature": dynamic_out,
        }.get,
        warnings,
    )

    if t_out is not None and q is not None:
        raise InvalidInputError(
            "duty",
            "a heater holds a temperature specification or a duty specification and not "
            "both: `Heater.setOutletTemperature` clears the duty flag and `setDuty` clears "
            "the temperature flag, so which one the class runs is the order they were set "
            "in and not a function of these inputs",
        )

    transient_requested = dt is not None
    if not transient_requested:
        for name, stated in (
            (
                "dynamic_temperature_control_enabled",
                dynamic_temperature_control_enabled is not None,
            ),
            ("design_mass_flow", design_mass is not None),
            ("design_inlet_temperature", design_in is not None),
            ("design_outlet_temperature", design_out is not None),
            ("cooling_medium_temperature", medium is not None),
            ("design_valve_opening", design_opening_si is not None),
            ("actuator_time_constant", tau_a is not None),
            ("thermal_time_constant", tau_t is not None),
            ("cooling_valve_opening", carried_opening_si is not None),
            ("requested_valve_opening", requested_si is not None),
            ("dynamic_outlet_temperature", dynamic_out is not None),
        ):
            if stated:
                raise InvalidInputError(name, "a transient input is stated without `time_step`")
    for name, value in (
        ("time_step", dt),
        ("design_mass_flow", design_mass),
        ("design_inlet_temperature", design_in),
        ("design_outlet_temperature", design_out),
        ("cooling_medium_temperature", medium),
        ("design_valve_opening", design_opening_si),
        ("actuator_time_constant", tau_a),
        ("thermal_time_constant", tau_t),
        ("cooling_valve_opening", carried_opening_si),
        ("requested_valve_opening", requested_si),
    ):
        if value is not None and not math.isfinite(value):
            raise InvalidInputError(name, f"a finite value is required, and {value} is not")
    if dynamic_out is not None and math.isinf(dynamic_out):
        raise InvalidInputError(
            "dynamic_outlet_temperature",
            f"a carried temperature is finite or NaN, and {dynamic_out} is neither",
        )

    active = transient_requested and bool(dynamic_temperature_control_enabled)
    if active:
        for name, stated in (
            ("design_mass_flow", design_mass is not None),
            ("design_inlet_temperature", design_in is not None),
            ("design_outlet_temperature", design_out is not None),
            ("cooling_medium_temperature", medium is not None),
        ):
            if not stated:
                raise InvalidInputError(
                    name,
                    "the configured branch needs the design point "
                    "`configureDynamicTemperatureControl` takes",
                )
        for name, stated in (
            ("outlet_temperature", t_out is not None),
            ("duty", q is not None),
        ):
            if stated:
                raise InvalidInputError(
                    name,
                    "the dynamic branch sets the outlet itself, so a steady specification "
                    "is ambiguous",
                )

    if active:
        if not math.isfinite(n) or n <= 0.0:
            raise InvalidInputError(
                "inlet_n",
                f"a unit operation needs a positive molar flow, and {n} mol/s is not",
            )
        dt_value = 0.0 if dt is None else dt
        design_in_value = 0.0 if design_in is None else design_in
        design_out_value = 0.0 if design_out is None else design_out
        medium_value = 0.0 if medium is None else medium
        design_opening = 50.0 if design_opening_si is None else design_opening_si
        opening = design_opening if carried_opening_si is None else carried_opening_si
        opening = min(100.0, max(0.0, opening))
        requested = opening if requested_si is None else min(100.0, max(0.0, requested_si))
        tau_actuator = 10.0 if tau_a is None else tau_a
        tau_thermal = 20.0 if tau_t is None else tau_t

        design_ntu = -math.log((design_out_value - medium_value) / (design_in_value - medium_value))
        actuator = dt_value / (tau_actuator + dt_value) if tau_actuator > 0.0 else 1.0
        opening += actuator * (requested - opening)

        mixture, _ = _components.mixture_of(components, eos="pr")
        molar_mass = sum(
            (component.molar_mass.to("kg/mol").magnitude if component.molar_mass else 0.0)
            * fraction
            for component, fraction in zip(mixture.components, inlet_z, strict=True)
        )
        mass_flow = max(n * molar_mass * 3600.0, 1.0e-12)
        design_mass_h = (0.0 if design_mass is None else design_mass) * 3600.0
        ntu = min(
            20.0,
            max(0.0, design_ntu * opening / design_opening * design_mass_h / mass_flow),
        )
        equilibrium = medium_value + (t - medium_value) * math.exp(-ntu) if t > medium_value else t
        # **The absent state is the design outlet**, which is the state
        # `configureDynamicTemperatureControl` writes before the first step.
        previous = (
            design_out_value if dynamic_out is None or math.isnan(dynamic_out) else dynamic_out
        )
        thermal = dt_value / (tau_thermal + dt_value) if tau_thermal > 0.0 else 1.0
        outlet_t = previous + thermal * (equilibrium - previous)
        outlet_t = min(t, max(min(medium_value, t), outlet_t))

        states = _states(components, n, inlet_z, p, t, outlet_t, None, drop)
        return CoolerResult(
            outlet_n=from_si(n, "mol/s"),
            outlet_z=tuple(inlet_z),
            outlet_p=from_si(states.p_out, "Pa"),
            outlet_t=states.outlet_t,
            outlet_h=from_si(states.h_out, "J/mol"),
            outlet_duty=from_si(states.duty, "W"),
            cooling_valve_opening=opening,
            equilibrium_outlet_temperature=from_si(equilibrium, "K"),
            warnings=tuple(warnings),
        )

    states = _states(components, n, inlet_z, p, t, t_out, q, drop)

    return CoolerResult(
        outlet_n=from_si(n, "mol/s"),
        outlet_z=tuple(inlet_z),
        outlet_p=from_si(states.p_out, "Pa"),
        outlet_t=states.outlet_t,
        outlet_h=from_si(states.h_out, "J/mol"),
        outlet_duty=from_si(states.duty, "W"),
        cooling_valve_opening=None,
        equilibrium_outlet_temperature=None,
        warnings=tuple(warnings),
    )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.cooler")
