"""``process.heater`` - the heater's kernel.

Spec: ``specs/models/process/heater.toml``

The Python twin of ``crates/azoth-process/src/models/heater.rs`` over
``crates/azoth-process/src/kernels/heater.rs``, written to mirror them.

# The duty is the state's, and it is in W against an enthalpy in J/mol

``Heater.run`` adds ``getDuty()`` (W) to ``system.getEnthalpy()`` (a total, J) directly, and
the two are the same dimension only while a process fluid's total moles equal its mol/s rate.
Here the conversion is written out - ``duty / inlet_n`` is the shift in J/mol - and the duty
*reported* is ``inlet_n * (h_out - h_in)``, which is ``run``'s own recomputation after the
flash rather than whatever the caller set.

# A temperature and a duty together are refused

``setOutletTemperature`` sets the temperature flag and clears the duty flag; ``setDuty`` does
the reverse. The class's answer to both is therefore *the order they were called in*, which
these inputs cannot express, so the pair is refused rather than resolved.

# With neither, the drop is isothermal

``run``'s else branch is ``T_in + dT`` with ``dT`` defaulting to zero, and ``dT`` is not a
palette parameter - so a heater that specifies nothing still moves the pressure.
"""

from __future__ import annotations

from dataclasses import dataclass

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HeaterResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve


def heater(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_temperature: Q | None = None,
    duty: Q | None = None,
    pressure_drop: Q | None = None,
) -> HeaterResult:
    """Heat or cool a stream to a stated temperature, or by a stated duty.

    Args:
        components: the fluid's substances, by name.
        inlet_n: the inlet's molar flow.
        inlet_z: the inlet composition.
        inlet_p: the inlet pressure.
        inlet_t: the inlet temperature.
        outlet_temperature: the temperature to reach, or ``None`` for the enthalpy route.
        duty: the heat to add, or ``None`` to state a temperature instead. Negative removes
            heat.
        pressure_drop: the pressure loss across the heater, or ``None`` for none.

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
        >>> r = heater(
        ...     ["methane", "n-butane"],
        ...     1.0,
        ...     q([0.9, 0.1]),
        ...     q(30.0, "bar"),
        ...     q(320.0, "K"),
        ...     duty=q(5000.0, "W"),
        ... )
        >>> round(r.outlet_t.to("K").magnitude, 3)
        420.458
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

    apply_checks(
        checks.on_input,
        {"inlet_t": t, "outlet_temperature": t_out}.get,
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

    states = _states(components, n, inlet_z, p, t, t_out, q, drop)

    return HeaterResult(
        outlet_n=from_si(n, "mol/s"),
        outlet_z=tuple(inlet_z),
        outlet_p=from_si(states.p_out, "Pa"),
        outlet_t=states.outlet_t,
        outlet_h=from_si(states.h_out, "J/mol"),
        outlet_duty=from_si(states.duty, "W"),
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class HeaterStates:
    """The heater's intermediates, in the order it computes them, in SI."""

    #: The inlet's molar enthalpy, J/mol.
    h_in: float
    #: The outlet pressure: the inlet's less the drop.
    p_out: float
    #: The outlet molar enthalpy, J/mol.
    h_out: float
    #: The outlet temperature.
    outlet_t: Q
    #: The duty moved, W.
    duty: float


def _states(
    components: list[str],
    inlet_n: float,
    inlet_z: list[float],
    inlet_p: float,
    inlet_t: float,
    outlet_temperature: float | None,
    duty: float | None,
    pressure_drop: float | None,
) -> HeaterStates:
    """The heater's intermediates, in SI, in the order it computes them.

    **One arithmetic, two consumers**: :func:`heater` builds the result from this and
    :mod:`azoth.process.layers` reports the same numbers against a NeqSim capture.
    """
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")

    h_in, _ = enthalpy_at(mixture, ideal_gas, inlet_t, inlet_p, inlet_z)

    p_out = inlet_p - (0.0 if pressure_drop is None else pressure_drop)
    if p_out <= 0.0:
        raise InvalidInputError(
            "pressure_drop",
            f"a pressure drop of {0.0 if pressure_drop is None else pressure_drop} Pa takes "
            f"a {inlet_p} Pa inlet to {p_out} Pa, which is not a pressure",
        )

    if outlet_temperature is not None:
        # The state at the stated temperature. `Heater.run` sets the temperature and runs a
        # `TPflash`, which is what `enthalpy_at` computes the enthalpy of.
        h_out, _ = enthalpy_at(mixture, ideal_gas, outlet_temperature, p_out, inlet_z)
        outlet_t = from_si(outlet_temperature, "K")
    elif duty is not None:
        if inlet_n == 0.0:
            raise InvalidInputError(
                "duty",
                "a duty is moved per mole of flow, and this stream carries none",
            )
        h_out = h_in + duty / inlet_n
        outlet_t = ph_flash_solve(
            mixture, ideal_gas, from_si(p_out, "Pa"), from_si(h_out, "J/mol"), inlet_z
        ).T
    else:
        # `run`'s else branch: `T_in + dT`, with `dT` zero. **An isothermal drop and not a
        # throttling**, which is the distinction an earlier version of this file got wrong
        # and the cross-implementation test caught: the outlet holds the inlet's
        # temperature, and its enthalpy is the *state's* at the outlet pressure, which for a
        # real fluid is not the inlet's. Measured on the captured fluid, two bar of drop at
        # 320 K moves it `49.6` J/mol - and NeqSim reports that movement as the duty.
        h_out, _ = enthalpy_at(mixture, ideal_gas, inlet_t, p_out, inlet_z)
        outlet_t = from_si(inlet_t, "K")

    return HeaterStates(
        h_in=h_in,
        p_out=p_out,
        h_out=h_out,
        outlet_t=outlet_t,
        duty=inlet_n * (h_out - h_in),
    )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.heater")
