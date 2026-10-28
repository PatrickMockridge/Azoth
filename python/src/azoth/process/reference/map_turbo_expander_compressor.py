"""``process.map_turbo_expander_compressor`` - the map-based single shaft.

Spec: ``specs/models/process/map_turbo_expander_compressor.toml``

The Python twin of ``crates/azoth-process/src/models/map_turbo_expander_compressor.rs`` over
``crates/azoth-process/src/kernels/map_turbo_expander.rs``, written to mirror them.

# Two modes, one identity

The expander's power, less the mechanical loss and the bearing loss, is what the compressor may
consume. In ``specified_pressures`` the compressor runs at its own setting and the residual is
reported; in ``balanced_speed`` the shaft speed is **bisected** over the compressor map's own
measured range until the compressor's power meets the available power - pegged at the lower
bound with ``under_power_surge`` and at the upper with ``over_power_max_speed``.

# This is the second shaft

``process.turbo_expander_compressor`` is NeqSim's older `TurboExpanderCompressor`; this one
drives a real brake compressor on a ``CompressorChart`` and its speed range is the map's.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import MapTurboExpanderResult
from azoth.core.units import Q, from_si, input_records_to_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.process.reference import compressor_chart
from azoth.process.reference.compressor import compressor as _compressor
from azoth.process.reference.expander import expander as _expander

MODEL_ID = "process.map_turbo_expander_compressor"

DEFAULT_MECHANICAL_EFFICIENCY = 0.98
DEFAULT_BEARING_LOSS = 0.0
DEFAULT_COMPRESSOR_SPEED = 3000.0
BISECTION_ITERATIONS = 80
SPEED_TOLERANCE = 1.0e-3


def map_turbo_expander_compressor(
    expander_components: list[str],
    compressor_components: list[str],
    expander_n: Q,
    expander_z: list[float],
    expander_p: Q,
    expander_t: Q,
    compressor_n: Q,
    compressor_z: list[float],
    compressor_p: Q,
    compressor_t: Q,
    expander_out_pressure: Q,
    expander_isentropic_efficiency: float,
    shaft_mode: str | None = None,
    mechanical_efficiency: float | None = None,
    bearing_loss: Q | None = None,
    compressor_outlet_pressure: Q | None = None,
    compressor_isentropic_efficiency: float | None = None,
    compressor_speed: float | None = None,
    min_shaft_speed: float | None = None,
    max_shaft_speed: float | None = None,
    curve: list[dict[str, Q]] | None = None,
    chart_type: str | None = None,
    chart_conditions: list[float] | None = None,
    impeller_diameter: Q | None = None,
) -> MapTurboExpanderResult:
    """Balance a map-driven compressor against an expander on one shaft.

    Args:
        expander_components: the expander feed's substances.
        compressor_components: the compressor feed's substances.
        expander_n: the expander feed's molar flow.
        expander_z: the expander feed's composition.
        expander_p: the expander feed's pressure.
        expander_t: the expander feed's temperature.
        compressor_n: the compressor feed's molar flow.
        compressor_z: the compressor feed's composition.
        compressor_p: the compressor feed's pressure.
        compressor_t: the compressor feed's temperature.
        expander_out_pressure: the expander's stated outlet pressure.
        expander_isentropic_efficiency: the expander's isentropic efficiency.
        shaft_mode: ``specified_pressures`` (absent) or ``balanced_speed``.
        mechanical_efficiency: the share of the expander's power that reaches the shaft.
        bearing_loss: the bearing loss, W.
        compressor_outlet_pressure: the isentropic route's outlet pressure.
        compressor_isentropic_efficiency: the isentropic route's efficiency.
        compressor_speed: the speed the map is read at in the specified mode, rpm.
        min_shaft_speed: the bisection's lower bound, rpm.
        max_shaft_speed: the bisection's upper bound, rpm.
        curve: the compressor's measured map.
        chart_type: which map class reads it.
        chart_conditions: the Khader map's chart conditions.
        impeller_diameter: the Khader map's impeller diameter.

    Returns:
        Both outlets, the shaft speed, both powers, the residual and the operating status.

    Raises:
        InvalidInputError: for a mechanical efficiency outside ``(0, 1]``, a negative bearing
            loss, a balanced run with no map or a map-less run with no outlet pressure.

    See :func:`azoth.process.map_turbo_expander_compressor`.
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    mode = "specified_pressures" if shaft_mode is None else shaft_mode
    if mode not in ("specified_pressures", "balanced_speed"):
        raise InvalidInputError(
            "shaft_mode", f"'{mode}' is not one of the two declared shaft modes"
        )
    mech = DEFAULT_MECHANICAL_EFFICIENCY if mechanical_efficiency is None else mechanical_efficiency
    loss = (
        DEFAULT_BEARING_LOSS
        if bearing_loss is None
        else input_to_si(spec, "bearing_loss", bearing_loss)
    )
    apply_checks(
        checks.on_input,
        {
            "expander_isentropic_efficiency": expander_isentropic_efficiency,
            "mechanical_efficiency": mech,
            "expander_t": input_to_si(spec, "expander_t", expander_t),
            "expander_p": input_to_si(spec, "expander_p", expander_p),
            "compressor_t": input_to_si(spec, "compressor_t", compressor_t),
            "compressor_p": input_to_si(spec, "compressor_p", compressor_p),
        }.get,
        warnings,
    )
    if not 0.0 < mech <= 1.0:
        raise InvalidInputError(
            "mechanical_efficiency",
            f"the mechanical efficiency must be in (0, 1], and {mech} is not",
        )
    if loss < 0.0 or math.isnan(loss):
        raise InvalidInputError(
            "bearing_loss", f"the bearing loss is {loss} W, and it cannot be negative"
        )

    rows: list[dict[str, float]] = (
        [] if curve is None else input_records_to_si(spec, "curve", curve)
    )
    points = [
        compressor_chart.CurvePoint(
            speed=float(row["speed"]),
            flow=float(row["flow"]) * 3600.0,
            head=float(row["head"]),
            efficiency=float(row["efficiency"]),
        )
        for row in rows
    ]

    exp_n = input_to_si(spec, "expander_n", expander_n)
    expander_outlet = _expander(
        expander_components,
        expander_n,
        expander_z,
        expander_p,
        expander_t,
        expander_out_pressure,
        expander_isentropic_efficiency,
    )
    mixture, ideal_gas = _components.mixture_of(expander_components, eos="pr")
    h_in, _ = enthalpy_at(
        mixture,
        ideal_gas,
        input_to_si(spec, "expander_t", expander_t),
        input_to_si(spec, "expander_p", expander_p),
        list(expander_z),
    )
    expander_power = exp_n * (expander_outlet.outlet_h.to("J/mol").magnitude - h_in)
    available = max(abs(expander_power) * mech - loss, 0.0)

    speed_default = (
        DEFAULT_COMPRESSOR_SPEED if compressor_speed is None else float(compressor_speed)
    )

    def power_at(speed: float) -> tuple[float, Any]:
        result = _compressor(
            components=compressor_components,
            inlet_n=compressor_n,
            inlet_z=compressor_z,
            inlet_p=compressor_p,
            inlet_t=compressor_t,
            isentropic_efficiency=(
                0.75
                if compressor_isentropic_efficiency is None
                else compressor_isentropic_efficiency
            ),
            outlet_pressure=None if points else compressor_outlet_pressure,
            curve=curve,
            speed=speed,
            chart_type=chart_type,
            chart_conditions=chart_conditions,
            impeller_diameter=impeller_diameter,
        )
        return result.power.to("W").magnitude, result

    if mode == "specified_pressures":
        consumed, compressor_outlet = power_at(speed_default)
        shaft_speed = speed_default
        speed_limited = False
        status = "balanced"
    else:
        if not points:
            raise InvalidInputError(
                "curve",
                "the balanced-speed mode bisects the compressor map, and this machine states none",
            )
        speeds = [point.speed for point in points]
        low = min(speeds) if min_shaft_speed is None else float(min_shaft_speed)
        high = max(speeds) if max_shaft_speed is None else float(max_shaft_speed)
        if not high > low:
            raise InvalidInputError(
                "curve", f"the map's speed range is {low} to {high} rpm, and a bisection needs one"
            )
        low_power, low_outlet = power_at(low)
        if low_power - available > 0.0:
            shaft_speed, speed_limited, status = low, True, "under_power_surge"
            consumed, compressor_outlet = low_power, low_outlet
        else:
            high_power, high_outlet = power_at(high)
            if high_power - available < 0.0:
                shaft_speed, speed_limited, status = high, True, "over_power_max_speed"
                consumed, compressor_outlet = high_power, high_outlet
            else:
                a, b = low, high
                for _ in range(BISECTION_ITERATIONS):
                    if abs(b - a) < SPEED_TOLERANCE:
                        break
                    mid = 0.5 * (a + b)
                    power, _ = power_at(mid)
                    if power - available > 0.0:
                        b = mid
                    else:
                        a = mid
                shaft_speed = 0.5 * (a + b)
                consumed, compressor_outlet = power_at(shaft_speed)
                speed_limited = False
                status = "balanced"

    return MapTurboExpanderResult(
        expander_outlet_n=from_si(exp_n, "mol/s"),
        expander_outlet_z=tuple(expander_z),
        expander_outlet_p=from_si(
            input_to_si(spec, "expander_out_pressure", expander_out_pressure), "Pa"
        ),
        expander_outlet_t=expander_outlet.outlet_t,
        expander_outlet_h=expander_outlet.outlet_h,
        compressor_outlet_n=compressor_outlet.outlet_n,
        compressor_outlet_z=compressor_outlet.outlet_z,
        compressor_outlet_p=compressor_outlet.outlet_p,
        compressor_outlet_t=compressor_outlet.outlet_t,
        compressor_outlet_h=compressor_outlet.outlet_h,
        shaft_speed=shaft_speed,
        available_shaft_power=from_si(available, "W"),
        consumed_compressor_power=from_si(consumed, "W"),
        power_balance_residual=from_si(available - consumed, "W"),
        speed_limited=speed_limited,
        operating_status=status,
        warnings=tuple(warnings),
    )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model(MODEL_ID)
