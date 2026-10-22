"""``eos.flowline_scale_profile`` - NeqSim's ``FlowlineScaleProfile``.

A water analysis walked along a linear temperature and pressure profile: ``N + 1`` points from
the inlet fraction zero to the outlet fraction one, each evaluated by the scale predictor.

Spec: ``specs/models/eos/flowline_scale_profile.toml``.
"""

from __future__ import annotations

from typing import Any

from azoth.core.range import apply_checks, checks_for
from azoth.core.record_outputs_gen import FlowlineScaleSegment
from azoth.core.result import FlowlineScaleProfileResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference.scale_prediction import scale_prediction as _predictor

CALC_ID = "eos.flowline_scale_profile"

#: The class's calculator default pH, which it never overwrites from `setWaterChemistry`.
DEFAULT_PH = 6.5


def flowline_scale_profile(
    calcium: Q,
    barium: Q,
    strontium: Q,
    iron: Q,
    magnesium: Q,
    sodium: Q,
    bicarbonate: Q,
    sulphate: Q,
    total_dissolved_solids: Q,
    co2_partial_pressure: Q,
    auto_ph: bool,
    inlet_temperature: Q,
    inlet_pressure: Q,
    outlet_temperature: Q,
    outlet_pressure: Q,
    segments: int,
) -> FlowlineScaleProfileResult:
    """Mineral-scale saturation indices along a flowline's profile.

    Args:
        calcium: calcium ion concentration.
        barium: barium ion concentration.
        strontium: strontium ion concentration.
        iron: ferrous iron concentration.
        magnesium: magnesium ion concentration.
        sodium: sodium ion concentration.
        bicarbonate: bicarbonate concentration.
        sulphate: sulphate concentration.
        total_dissolved_solids: TDS, which sets each point's ionic strength.
        co2_partial_pressure: CO2 partial pressure, read by the auto-pH path.
        auto_ph: whether each point estimates its pH with the class's Henry constant.
        inlet_temperature: the wellhead's absolute temperature.
        inlet_pressure: the wellhead's absolute pressure.
        outlet_temperature: the host's absolute temperature.
        outlet_pressure: the host's absolute pressure.
        segments: the number of segments; the grid is one point more, clamped to two.

    Returns:
        The points and the largest index each mineral reaches.

    Raises:
        OutOfRangeError: if an input is outside its declared range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = flowline_scale_profile(
        ...     q(0.4, "kg/m**3"), q(0.01, "kg/m**3"), q(0.005, "kg/m**3"),
        ...     q(0.0, "kg/m**3"), q(0.15, "kg/m**3"), q(0.01, "kg/m**3"),
        ...     q(0.15, "kg/m**3"), q(0.01, "kg/m**3"), q(40.0, "kg/m**3"),
        ...     q(2.0e5, "Pa"), False, q(353.15, "K"), q(2.0e7, "Pa"),
        ...     q(293.15, "K"), q(5.0e6, "Pa"), 10)
        >>> len(r.segments)
        11
        >>> abs(r.max_baso4 - 1.048081964483009) < 1e-9
        True
    """
    from azoth._models_gen import model as _model_for

    spec: dict[str, Any] = _model_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        name: input_to_si(spec, name, value)
        for name, value in (
            ("calcium", calcium),
            ("barium", barium),
            ("strontium", strontium),
            ("iron", iron),
            ("magnesium", magnesium),
            ("sodium", sodium),
            ("bicarbonate", bicarbonate),
            ("sulphate", sulphate),
            ("total_dissolved_solids", total_dissolved_solids),
            ("co2_partial_pressure", co2_partial_pressure),
            ("inlet_temperature", inlet_temperature),
            ("inlet_pressure", inlet_pressure),
            ("outlet_temperature", outlet_temperature),
            ("outlet_pressure", outlet_pressure),
        )
    }
    values["segments"] = segments
    apply_checks(checks.on_input, values.get, warnings)

    steps = max(2, int(segments))
    rows: list[FlowlineScaleSegment] = []
    maxima: list[float | None] = [None] * 5
    for index in range(steps + 1):
        fraction = index / steps
        temperature_c = (values["inlet_temperature"] - 273.15) + fraction * (
            (values["outlet_temperature"] - 273.15) - (values["inlet_temperature"] - 273.15)
        )
        pressure_bara = values["inlet_pressure"] / 1.0e5 + fraction * (
            values["outlet_pressure"] / 1.0e5 - values["inlet_pressure"] / 1.0e5
        )
        point = _predictor(
            calcium,
            barium,
            strontium,
            iron,
            magnesium,
            sodium,
            bicarbonate,
            sulphate,
            total_dissolved_solids,
            from_si(temperature_c + 273.15, "K"),
            from_si(pressure_bara * 1.0e5, "Pa"),
            co2_partial_pressure,
            DEFAULT_PH,
            auto_ph,
            False,
        )
        indices = (
            point.si_caco3,
            point.si_baso4,
            point.si_srso4,
            point.si_caso4,
            point.si_feco3,
        )
        for slot, value in enumerate(indices):
            current = maxima[slot]
            if value is not None and (current is None or value > current):
                maxima[slot] = value
        rows.append(
            FlowlineScaleSegment(
                distance_fraction=fraction,
                temperature=from_si(temperature_c + 273.15, "K"),
                pressure=from_si(pressure_bara * 1.0e5, "Pa"),
                si_caco3=indices[0],
                si_baso4=indices[1],
                si_srso4=indices[2],
                si_caso4=indices[3],
                si_feco3=indices[4],
            )
        )

    return FlowlineScaleProfileResult(
        segments=tuple(rows),
        max_caco3=maxima[0],
        max_baso4=maxima[1],
        max_srso4=maxima[2],
        max_caso4=maxima[3],
        max_feco3=maxima[4],
        warnings=tuple(warnings),
    )
