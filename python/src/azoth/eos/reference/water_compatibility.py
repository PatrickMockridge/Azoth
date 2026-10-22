"""``eos.water_compatibility`` - NeqSim's ``WaterCompatibilityScreener``.

Two waters mixed linearly at a list of injection-water percentages, each mixture evaluated by
the scale predictor with its auto-pH always on, and the worst saturation index reported.

Spec: ``specs/models/eos/water_compatibility.toml``.
"""

from __future__ import annotations

from typing import Any, cast

from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.record_outputs_gen import WaterMixingResult
from azoth.core.result import WaterCompatibilityResult
from azoth.core.units import Q, from_si, input_records_to_si
from azoth.core.warnings import Warning, WarningCode
from azoth.eos.reference.scale_prediction import scale_prediction as _predictor

CALC_ID = "eos.water_compatibility"

#: The class's calculator default pH; the two water records' own pH fields are never read.
DEFAULT_PH = 6.5

SCALES = ("CaCO3", "BaSO4", "SrSO4", "CaSO4", "FeCO3")


def water_compatibility(
    formation_water: list[dict[str, Q]],
    injection_water: list[dict[str, Q]],
    injection_water_percent: list[float],
) -> WaterCompatibilityResult:
    """The scale risk of mixing a formation water with an injection water.

    Args:
        formation_water: the formation water's analysis, exactly one row.
        injection_water: the injection water's analysis, exactly one row.
        injection_water_percent: the injection water's share at each mixture, in percent.

    Returns:
        One row per percentage and the worst case over them.

    Raises:
        InvalidInputError: for a water list that is not exactly one row.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> formation = [{
        ...     "calcium": q(0.4, "kg/m**3"), "barium": q(0.2, "kg/m**3"),
        ...     "strontium": q(0.005, "kg/m**3"), "iron": q(0.0, "kg/m**3"),
        ...     "magnesium": q(0.15, "kg/m**3"), "sodium": q(10.0, "kg/m**3"),
        ...     "bicarbonate": q(0.15, "kg/m**3"), "sulphate": q(0.01, "kg/m**3"),
        ...     "total_dissolved_solids": q(50.0, "kg/m**3"),
        ...     "temperature": q(353.15, "K"), "pressure": q(1.0e7, "Pa"),
        ...     "co2_partial_pressure": q(2.0e5, "Pa")}]
        >>> injection = [{
        ...     "calcium": q(0.02, "kg/m**3"), "barium": q(0.0, "kg/m**3"),
        ...     "strontium": q(0.0, "kg/m**3"), "iron": q(0.0, "kg/m**3"),
        ...     "magnesium": q(0.1, "kg/m**3"), "sodium": q(11.0, "kg/m**3"),
        ...     "bicarbonate": q(0.1, "kg/m**3"), "sulphate": q(2.7, "kg/m**3"),
        ...     "total_dissolved_solids": q(35.0, "kg/m**3"),
        ...     "temperature": q(293.15, "K"), "pressure": q(1.0e7, "Pa"),
        ...     "co2_partial_pressure": q(0.5e5, "Pa")}]
        >>> r = water_compatibility(formation, injection, [60.0])
        >>> r.worst_case_scale
        'BaSO4'
    """
    from azoth._models_gen import model as _model_for

    spec: dict[str, Any] = _model_for(CALC_ID)
    warnings: list[Warning] = []

    def single(rows: list[dict[str, Q]], role: str) -> dict[str, float]:
        if len(rows) != 1:
            raise InvalidInputError(role, f"`{role}` is one water and this is {len(rows)} row(s)")
        return cast("dict[str, float]", input_records_to_si(spec, role, rows)[0])

    formation = single(formation_water, "formation_water")
    injection = single(injection_water, "injection_water")

    for label, analysis in (
        ("formation_water", formation),
        ("injection_water", injection),
    ):
        for name, value, positive in (
            ("calcium", analysis["calcium"], False),
            ("barium", analysis["barium"], False),
            ("strontium", analysis["strontium"], False),
            ("iron", analysis["iron"], False),
            ("magnesium", analysis["magnesium"], False),
            ("sodium", analysis["sodium"], False),
            ("bicarbonate", analysis["bicarbonate"], False),
            ("sulphate", analysis["sulphate"], False),
            ("total_dissolved_solids", analysis["total_dissolved_solids"], False),
            ("temperature", analysis["temperature"], True),
            ("pressure", analysis["pressure"], True),
            ("co2_partial_pressure", analysis["co2_partial_pressure"], False),
        ):
            if positive and not value > 0.0:
                raise OutOfRangeError(f"{label}.{name}", value, "must be positive")
            if not positive and value < 0.0:
                raise OutOfRangeError(f"{label}.{name}", value, "must be non-negative")

    for percent in injection_water_percent:
        if not 0.0 <= percent <= 100.0:
            warnings.append(
                Warning(
                    code=WarningCode.OUT_OF_VALID_RANGE,
                    message=(
                        f"an injection-water percentage of {percent} is outside [0, 100]; "
                        "the class mixes with it as stated rather than clamping"
                    ),
                    field=None,
                )
            )

    def mix(name: str, fraction: float) -> float:
        return formation[name] * (1.0 - fraction) + injection[name] * fraction

    rows: list[WaterMixingResult] = []
    worst: tuple[float, str, float] | None = None
    for percent in injection_water_percent:
        fraction = percent / 100.0
        temperature_c = (
            formation["temperature"]
            - 273.15
            + fraction * ((injection["temperature"] - 273.15) - (formation["temperature"] - 273.15))
        )
        point = _predictor(
            from_si(mix("calcium", fraction), "kg/m**3"),
            from_si(mix("barium", fraction), "kg/m**3"),
            from_si(mix("strontium", fraction), "kg/m**3"),
            from_si(mix("iron", fraction), "kg/m**3"),
            from_si(mix("magnesium", fraction), "kg/m**3"),
            from_si(mix("sodium", fraction), "kg/m**3"),
            from_si(mix("bicarbonate", fraction), "kg/m**3"),
            from_si(mix("sulphate", fraction), "kg/m**3"),
            from_si(mix("total_dissolved_solids", fraction), "kg/m**3"),
            from_si(temperature_c + 273.15, "K"),
            from_si(mix("pressure", fraction), "Pa"),
            from_si(mix("co2_partial_pressure", fraction), "Pa"),
            DEFAULT_PH,
            True,
            False,
        )
        indices = (
            point.si_caco3,
            point.si_baso4,
            point.si_srso4,
            point.si_caso4,
            point.si_feco3,
        )
        for scale, index in zip(SCALES, indices, strict=True):
            if index is None:
                continue
            if worst is None or index > worst[2]:
                worst = (percent, scale, index)
        rows.append(
            WaterMixingResult(
                injection_water_percent=percent,
                si_caco3=indices[0],
                si_baso4=indices[1],
                si_srso4=indices[2],
                si_caso4=indices[3],
                si_feco3=indices[4],
            )
        )

    return WaterCompatibilityResult(
        mixtures=tuple(rows),
        worst_case_percent=None if worst is None else worst[0],
        worst_case_scale=None if worst is None else worst[1],
        worst_case_saturation_index=None if worst is None else worst[2],
        warnings=tuple(warnings),
    )
