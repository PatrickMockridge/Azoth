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

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import CoolerResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
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
    t_out = (
        None
        if outlet_temperature is None
        else input_to_si(spec, "outlet_temperature", outlet_temperature)
    )
    q = None if duty is None else input_to_si(spec, "duty", duty)
    drop = None if pressure_drop is None else input_to_si(spec, "pressure_drop", pressure_drop)

    apply_checks(checks.on_input, {"inlet_t": t, "outlet_temperature": t_out}.get, warnings)

    if t_out is not None and q is not None:
        raise InvalidInputError(
            "duty",
            "a heater holds a temperature specification or a duty specification and not "
            "both: `Heater.setOutletTemperature` clears the duty flag and `setDuty` clears "
            "the temperature flag, so which one the class runs is the order they were set "
            "in and not a function of these inputs",
        )

    states = _states(components, n, inlet_z, p, t, t_out, q, drop)

    return CoolerResult(
        outlet_n=from_si(n, "mol/s"),
        outlet_z=tuple(inlet_z),
        outlet_p=from_si(states.p_out, "Pa"),
        outlet_t=states.outlet_t,
        outlet_h=from_si(states.h_out, "J/mol"),
        outlet_duty=from_si(states.duty, "W"),
        warnings=tuple(warnings),
    )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.cooler")
