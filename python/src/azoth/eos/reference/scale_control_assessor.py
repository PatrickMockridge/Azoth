"""``eos.scale_control_assessor`` - NeqSim's ``ScaleControlAssessor``.

The scale predictor's five indices with a threshold inhibitor's efficiency folded into a
kinetic risk index per scale: ``SI + log10(max(1e-6, 1 - efficiency))``, the efficiency
capped at 0.999 and the worst index reported against the class's 0.5 threshold.

Spec: ``specs/models/eos/scale_control_assessor.toml``.
"""

from __future__ import annotations

import math
from typing import Any, cast

from azoth.core.range import apply_checks, checks_for
from azoth.core.record_outputs_gen import ScaleControlRiskRow, ScaleControlRow
from azoth.core.result import ScaleControlAssessorResult
from azoth.core.units import Q, from_si, input_records_to_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference.scale_inhibitor_performance import scale_inhibitor_performance
from azoth.eos.reference.scale_prediction import _state

CALC_ID = "eos.scale_control_assessor"

#: The scales, in the class's map order.
SCALES = ("CACO3", "BASO4", "SRSO4", "CASO4", "FECO3")


def scale_control_assessor(
    calcium: Q,
    barium: Q,
    strontium: Q,
    iron: Q,
    magnesium: Q,
    sodium: Q,
    bicarbonate: Q,
    sulphate: Q,
    total_dissolved_solids: Q,
    temperature: Q,
    pressure: Q,
    co2_partial_pressure: Q,
    ph: float,
    inhibitors: list[dict[str, Any]],
) -> ScaleControlAssessorResult:
    """The residual scaling risk after threshold inhibition.

    Args:
        calcium: calcium ion concentration.
        barium: barium ion concentration.
        strontium: strontium ion concentration.
        iron: ferrous iron concentration.
        magnesium: magnesium ion concentration.
        sodium: sodium ion concentration.
        bicarbonate: bicarbonate concentration.
        sulphate: sulphate ion concentration.
        total_dissolved_solids: TDS, which set the predictor's ionic strength.
        temperature: absolute temperature.
        pressure: absolute pressure.
        co2_partial_pressure: CO2 partial pressure; the predictor's auto-pH is off.
        ph: the water's pH.
        inhibitors: the configured inhibitor rows, one per scale at most.

    Returns:
        The uninhibited indices, the inhibited risk indices, the worst and the verdict.

    Raises:
        OutOfRangeError: if an input is outside its declared range.
        InvalidInputError: for an unknown scale type or chemistry in an inhibitor row.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> row = {
        ...     "scale_type": "caco3", "inhibitor_chemistry": "phosphonate",
        ...     "temperature": q(353.15, "K"), "saturation_ratio": 1.6851624454329173,
        ...     "total_dissolved_solids": q(80.0, "kg/m**3"),
        ...     "calcium": q(2.0, "kg/m**3"), "available_dose": q(0.02, "kg/m**3")}
        >>> r = scale_control_assessor(
        ...     q(2.0, "kg/m**3"), q(0.05, "kg/m**3"), q(0.0, "kg/m**3"),
        ...     q(0.0, "kg/m**3"), q(0.0, "kg/m**3"), q(20.0, "kg/m**3"),
        ...     q(0.8, "kg/m**3"), q(0.5, "kg/m**3"), q(80.0, "kg/m**3"),
        ...     q(353.15, "K"), q(1.0e5, "Pa"), q(0.0, "Pa"), 7.0, [row])
        >>> r.kinetically_controlled
        False
    """
    from azoth._models_gen import model as _model_for

    spec: dict[str, Any] = _model_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values: dict[str, float] = {
        name: input_to_si(spec, name, value)
        for name, value in {
            "calcium": calcium,
            "barium": barium,
            "strontium": strontium,
            "iron": iron,
            "magnesium": magnesium,
            "sodium": sodium,
            "bicarbonate": bicarbonate,
            "sulphate": sulphate,
            "total_dissolved_solids": total_dissolved_solids,
            "temperature": temperature,
            "pressure": pressure,
            "co2_partial_pressure": co2_partial_pressure,
        }.items()
    }
    values["ph"] = ph
    apply_checks(checks.on_input, values.get, warnings)

    rows = cast("list[dict[str, Any]]", input_records_to_si(spec, "inhibitors", inhibitors))

    point = _state(
        values,
        values["temperature"],
        values["pressure"],
        values["co2_partial_pressure"],
        ph,
        False,
        False,
    )
    indices = point["si"]

    thermodynamic_si: list[ScaleControlRow] = []
    kinetic_risk_index: list[ScaleControlRiskRow] = []
    worst: float | None = None

    for scale, si in zip(SCALES, indices, strict=True):
        thermodynamic = None if math.isnan(si) else si
        if math.isnan(si):
            risk = None
        else:
            match = next(
                (row for row in reversed(rows) if str(row["scale_type"]).lower() == scale.lower()),
                None,
            )
            if match is None:
                risk = si
            else:
                performance = scale_inhibitor_performance(
                    str(match["scale_type"]),
                    str(match["inhibitor_chemistry"]),
                    from_si(match["temperature"], "K"),
                    match["saturation_ratio"],
                    from_si(match["total_dissolved_solids"], "kg/m**3"),
                    from_si(match["calcium"], "kg/m**3"),
                    from_si(match["available_dose"], "kg/m**3"),
                )
                efficiency = min(0.999, max(0.0, performance.efficiency))
                risk = si + math.log10(max(1.0e-6, 1.0 - efficiency))

        if risk is not None:
            worst = risk if worst is None else max(worst, risk)
        thermodynamic_si.append(ScaleControlRow(scale=scale, saturation_index=thermodynamic))
        kinetic_risk_index.append(ScaleControlRiskRow(scale=scale, risk_index=risk))

    return ScaleControlAssessorResult(
        thermodynamic_si=tuple(thermodynamic_si),
        kinetic_risk_index=tuple(kinetic_risk_index),
        worst_kinetic_risk_index=worst,
        kinetically_controlled=worst is not None and worst < 0.5,
        warnings=tuple(warnings),
    )
