"""``process.heat_exchanger`` - the exchanger's kernel.

Spec: ``specs/models/process/heat_exchanger.toml``

The Python twin of ``crates/azoth-process/src/models/heat_exchanger.rs`` over
``crates/azoth-process/src/kernels/heat_exchanger.rs``, written to mirror them.

# The rating, not a duty

``HeatExchanger.run``'s default branch specifies ``UAvalue`` and runs an effectiveness-NTU
rating: each side's capacity is estimated by flashing it at the **other side's inlet
temperature** and dividing the total enthalpy change by that temperature change, then
``NTU = UA / C_min`` and ``calcThermalEffectivenes(NTU, C_min / C_max)`` scales that swing.
``duty`` is an *output* of this class - ``energyInput`` appears nowhere in it - so a kernel
that moved a caller's duty was a different unit operation under the same name.

``C`` is a difference of *totals* over a difference of temperatures, which is W/K; that is
what makes the NTU dimensionless. The absolute value matters: the hot side's seeded swing
is negative because it cools, and taking it at face value would make ``C_min`` and the NTU
negative, which is an exponential in the wrong direction.

At an effectiveness of one ``run`` leaves the capacity-limiting side at its seeded state -
"this side at the other's inlet temperature". That is the same *enthalpy* the arithmetic
gives, so the branch is an optimisation rather than a second answer, and this reaches it
through the flash either way.

# The other mode pins an outlet

``hot_outlet_temperature`` or ``cold_outlet_temperature`` pins that side at its inlet
pressure and energy-balances the other against it. **NeqSim's own ``runSpecifiedStream``
does not reach the state it names**: it clones the stream's already-flashed fluid, sets the
temperature and flashes again, and the enthalpy it lands on is one NeqSim's own ``PHflash``
places at another temperature. This reaches the pinned temperature directly, which is why
the case is pinned to the capture's ``reference_pin_*`` rows rather than to its
``out_temperature_*`` ones.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HeatExchangerResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve

#: The arrangements ``calcThermalEffectivenes`` names, in azoth's own spelling.
FLOW_ARRANGEMENTS = ("counterflow", "parallelflow", "shell_and_tube")


def heat_exchanger(
    hot_components: list[str],
    hot_in_n: Q,
    hot_in_z: list[float],
    hot_in_p: Q,
    hot_in_t: Q,
    cold_components: list[str],
    cold_in_n: Q,
    cold_in_z: list[float],
    cold_in_p: Q,
    cold_in_t: Q,
    ua: Q | None = None,
    flow_arrangement: str = "counterflow",
    hot_outlet_temperature: Q | None = None,
    cold_outlet_temperature: Q | None = None,
) -> HeatExchangerResult:
    """Exchange heat between two streams.

    Args:
        hot_components: the hot side's substances, by name.
        hot_in_n: the hot inlet's molar flow.
        hot_in_z: the hot inlet's composition.
        hot_in_p: the hot inlet's pressure.
        hot_in_t: the hot inlet's temperature.
        cold_components: the cold side's substances, by name.
        cold_in_n: the cold inlet's molar flow.
        cold_in_z: the cold inlet's composition.
        cold_in_p: the cold inlet's pressure.
        cold_in_t: the cold inlet's temperature.
        ua: the overall conductance, for the rating. Unread when an outlet is pinned.
        flow_arrangement: ``counterflow``, ``parallelflow`` or ``shell_and_tube``.
        hot_outlet_temperature: pin the hot outlet at this temperature.
        cold_outlet_temperature: pin the cold outlet at this temperature.

    Returns:
        Both outlets' records.

    Raises:
        InvalidInputError: for the mode and arrangement refusals the spec lists.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = heat_exchanger(
        ...     ["methane", "n-butane"],
        ...     q(1.0, "mol/s"),
        ...     [0.8, 0.2],
        ...     q(20.0, "bar"),
        ...     q(400.0, "K"),
        ...     ["n-butane", "n-pentane"],
        ...     q(1.0, "mol/s"),
        ...     [0.5, 0.5],
        ...     q(5.0, "bar"),
        ...     q(300.0, "K"),
        ...     q(100.0, "W/K"),
        ...     "counterflow",
        ... )
        >>> round(r.hot_out_t.to("K").magnitude, 4)
        313.7571
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    ua_si = None if ua is None else input_to_si(spec, "ua", ua)
    hot_t_si = input_to_si(spec, "hot_in_t", hot_in_t)
    cold_t_si = input_to_si(spec, "cold_in_t", cold_in_t)
    # Every declared bound's quantity is fed here: a closure that supplies some and not
    # the rest emits a `RANGE_CHECK_SKIPPED` warning for a check that could have run.
    apply_checks(
        checks.on_input,
        {"ua": ua_si, "hot_in_t": hot_t_si, "cold_in_t": cold_t_si}.get,
        warnings,
    )

    if hot_outlet_temperature is not None and cold_outlet_temperature is not None:
        raise InvalidInputError(
            "hot_outlet_temperature",
            "an exchanger pins one outlet at most: with both stated the other side has "
            "nothing left to solve for",
        )
    pinned = hot_outlet_temperature is not None or cold_outlet_temperature is not None
    if not pinned and ua_si is None:
        raise InvalidInputError(
            "ua",
            "the rating needs an overall conductance. NeqSim's `UAvalue` defaults to "
            "500 W/K, and this refuses to guess it: a default that only changes the "
            "answer is not one a caller can see",
        )
    if flow_arrangement not in FLOW_ARRANGEMENTS:
        raise InvalidInputError(
            "flow_arrangement",
            f"{flow_arrangement} is not one of {list(FLOW_ARRANGEMENTS)}. NeqSim's `run` "
            f"falls through to the counterflow relation for an arrangement it does not "
            f"know, which would make a misspelling a plausible answer",
        )

    hot_mixture, hot_gas = _components.mixture_of(hot_components, eos="pr")
    cold_mixture, cold_gas = _components.mixture_of(cold_components, eos="pr")
    hot_n = input_to_si(spec, "hot_in_n", hot_in_n)
    cold_n = input_to_si(spec, "cold_in_n", cold_in_n)
    hot_p = input_to_si(spec, "hot_in_p", hot_in_p)
    cold_p = input_to_si(spec, "cold_in_p", cold_in_p)
    hot_h, _ = enthalpy_at(hot_mixture, hot_gas, hot_t_si, hot_p, hot_in_z)
    cold_h, _ = enthalpy_at(cold_mixture, cold_gas, cold_t_si, cold_p, cold_in_z)

    if pinned:
        hot_pinned = hot_outlet_temperature is not None
        pinned_t = input_to_si(
            spec,
            "hot_outlet_temperature" if hot_pinned else "cold_outlet_temperature",
            hot_outlet_temperature if hot_pinned else cold_outlet_temperature,
        )
        if hot_pinned:
            specified, _ = enthalpy_at(hot_mixture, hot_gas, pinned_t, hot_p, hot_in_z)
            duty = hot_n * (specified - hot_h)
            hot_out_h, cold_out_h = specified, cold_h - duty / cold_n
        else:
            specified, _ = enthalpy_at(cold_mixture, cold_gas, pinned_t, cold_p, cold_in_z)
            duty = cold_n * (specified - cold_h)
            cold_out_h, hot_out_h = specified, hot_h - duty / hot_n
    else:
        span = abs(hot_t_si - cold_t_si)
        # The same guard the Rust uses, so the two refuse the same inputs rather than
        # agreeing by accident on the cases that do not reach it.
        if span <= 2.220446049250313e-16 * max(abs(hot_t_si), abs(cold_t_si)):
            raise InvalidInputError(
                "hot_in_t",
                f"both inlets are at {hot_t_si} K, so there is no driving force for the "
                f"rating to size against",
            )
        assert ua_si is not None  # the mode check above guarantees it
        hot_seeded, _ = enthalpy_at(hot_mixture, hot_gas, cold_t_si, hot_p, hot_in_z)
        cold_seeded, _ = enthalpy_at(cold_mixture, cold_gas, hot_t_si, cold_p, cold_in_z)
        hot_swing = hot_n * (hot_seeded - hot_h)
        cold_swing = cold_n * (cold_seeded - cold_h)
        c_min = min(abs(hot_swing), abs(cold_swing)) / span
        c_max = max(abs(hot_swing), abs(cold_swing)) / span
        effectiveness = _effectiveness(ua_si / c_min, c_min / c_max, flow_arrangement)
        # The side whose seeded swing is the larger is the one `run` keeps, and the other
        # is energy-balanced against it.
        if abs(cold_swing) > abs(hot_swing):
            duty = effectiveness * hot_swing
            hot_out_h, cold_out_h = hot_h + duty / hot_n, cold_h - duty / cold_n
        else:
            duty = effectiveness * cold_swing
            cold_out_h, hot_out_h = cold_h + duty / cold_n, hot_h - duty / hot_n

    return HeatExchangerResult(
        hot_out_n=from_si(hot_n, "mol/s"),
        hot_out_z=tuple(hot_in_z),
        hot_out_p=hot_in_p,
        hot_out_t=_temperature(hot_mixture, hot_gas, hot_p, hot_out_h, hot_in_z),
        hot_out_h=from_si(hot_out_h, "J/mol"),
        cold_out_n=from_si(cold_n, "mol/s"),
        cold_out_z=tuple(cold_in_z),
        cold_out_p=cold_in_p,
        cold_out_t=_temperature(cold_mixture, cold_gas, cold_p, cold_out_h, cold_in_z),
        cold_out_h=from_si(cold_out_h, "J/mol"),
        warnings=tuple(warnings),
    )


def _temperature(
    mixture: object, ideal_gas: object, pressure: float, h: float, z: list[float]
) -> Q:
    """The temperature a side reaches at a shifted enthalpy."""
    solved = ph_flash_solve(
        mixture,  # type: ignore[arg-type]
        ideal_gas,  # type: ignore[arg-type]
        from_si(pressure, "Pa"),
        from_si(h, "J/mol"),
        z,
    )
    return solved.T


def _effectiveness(ntu: float, capacity_ratio: float, arrangement: str) -> float:
    """The effectiveness ``calcThermalEffectivenes`` carries, for the three arrangements."""
    if capacity_ratio == 0.0:
        return 1.0 - math.exp(-ntu)
    if arrangement == "parallelflow":
        return (1.0 - math.exp(-ntu * (1.0 + capacity_ratio))) / (1.0 + capacity_ratio)
    if arrangement == "counterflow" and capacity_ratio == 1.0:
        return ntu / (1.0 + ntu)
    exp = math.exp(-ntu * (1.0 - capacity_ratio))
    return (1.0 - exp) / (1.0 - capacity_ratio * exp)


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.heat_exchanger")
