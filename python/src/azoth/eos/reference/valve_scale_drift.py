"""``eos.valve_scale_drift`` - NeqSim's ``ValveScaleDrift``.

A uniform radial deposit on a circular trim port: the open-area fraction, the fouling
fraction and the effective Kv at an elapsed time, with the time to plug the port and the
time for the valve opening to pin at 100 percent.

Spec: ``specs/models/eos/valve_scale_drift.toml``.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ValveScaleDriftResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.valve_scale_drift"

#: The fouling fraction's cap, as the class's `getFoulingFraction` states it.
FOULING_CAP = 0.999_999


def valve_scale_drift(
    port_diameter: Q,
    growth_rate: Q,
    elapsed_time: Q,
    clean_opening_percent: float,
    valve_kv: float,
) -> ValveScaleDriftResult:
    """The trim fouling a deposit growth rate produces at an elapsed time.

    Args:
        port_diameter: the trim port or seat diameter.
        growth_rate: the radial deposit growth rate.
        elapsed_time: the elapsed simulated time.
        clean_opening_percent: the opening the clean valve held, 0..100.
        valve_kv: the clean valve's flow coefficient.

    Returns:
        The thickness, open area, fouling, effective Kv and the two times.

    Raises:
        OutOfRangeError: if an input is outside the spec's hard range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = valve_scale_drift(
        ...     q(0.05, "m"), q(3.168808781402895e-11, "m/s"), q(2_592_000.0, "s"),
        ...     50.0, 100.0)
        >>> abs(r.open_area_fraction - 0.9934399521016661) < 1e-12
        True
    """
    from azoth._models_gen import model as _model_for

    spec: dict[str, Any] = _model_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "port_diameter": input_to_si(spec, "port_diameter", port_diameter),
        "growth_rate": input_to_si(spec, "growth_rate", growth_rate),
        "elapsed_time": input_to_si(spec, "elapsed_time", elapsed_time),
        "clean_opening_percent": clean_opening_percent,
        "valve_kv": valve_kv,
    }
    apply_checks(checks.on_input, values.get, warnings)

    growth = max(0.0, values["growth_rate"])
    elapsed = max(0.0, values["elapsed_time"])
    port = values["port_diameter"]
    thickness = growth * elapsed

    if port <= 0.0:
        open_area = 1.0
    else:
        open_diameter = port - 2.0 * thickness
        open_area = 0.0 if open_diameter <= 0.0 else (open_diameter / port) ** 2

    fouling = min(FOULING_CAP, max(0.0, 1.0 - open_area))

    if growth <= 0.0 or port <= 0.0:
        time_to_plug = math.inf
    else:
        time_to_plug = max(0.0, port / 2.0 - thickness) / growth

    opening_percent = 100.0 if open_area <= 0.0 else min(100.0, clean_opening_percent / open_area)

    if growth <= 0.0 or port <= 0.0:
        time_to_pin = math.inf
    elif clean_opening_percent >= 100.0:
        time_to_pin = 0.0
    elif clean_opening_percent <= 0.0:
        time_to_pin = math.inf
    else:
        pin_thickness = port / 2.0 * (1.0 - math.sqrt(clean_opening_percent / 100.0))
        remaining = pin_thickness - thickness
        time_to_pin = 0.0 if remaining <= 0.0 else remaining / growth

    return ValveScaleDriftResult(
        deposit_thickness=from_si(thickness, "m"),
        elapsed_time=from_si(elapsed, "s"),
        open_area_fraction=open_area,
        fouling_fraction=fouling,
        effective_kv=valve_kv * open_area,
        time_to_plug=from_si(time_to_plug, "s"),
        opening_percent=opening_percent,
        time_to_pin=from_si(time_to_pin, "s"),
        warnings=tuple(warnings),
    )
