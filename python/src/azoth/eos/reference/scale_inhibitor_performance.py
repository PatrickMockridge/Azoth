"""``eos.scale_inhibitor_performance`` - NeqSim's ``ScaleInhibitorPerformance``.

An empirical minimum-inhibitor-concentration correlation: a base MIC per mineral and
chemistry, four multiplicative factors over temperature, TDS, saturation and calcium,
and the efficiency at the dose actually available.

Spec: ``specs/calcs/eos/scale_inhibitor_performance.toml``.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ScaleInhibitorPerformanceResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.scale_inhibitor_performance"

_BASE = {"caco3": 5.0, "baso4": 10.0, "srso4": 8.0, "caso4": 6.0, "feco3": 12.0}
_CHEMISTRY_FACTOR = {
    "phosphonate": 1.0,
    "polymaleate": 1.2,
    "polyacrylate": 1.4,
    "phosphate_ester": 1.1,
    "vinyl_sulphonate": 1.3,
}
#: The temperature knee and slope, in degrees Celsius and per degree.
_KNEE = {
    "phosphonate": (100.0, 0.015),
    "polymaleate": (130.0, 0.012),
    "polyacrylate": (120.0, 0.013),
    "phosphate_ester": (150.0, 0.010),
    "vinyl_sulphonate": (100.0, 0.015),
}


def _base_mic(scale: str, chemistry: str) -> float:
    factor = _CHEMISTRY_FACTOR[chemistry]
    if scale == "baso4" and chemistry == "phosphonate":
        factor *= 0.8
    if scale == "caco3" and chemistry == "polymaleate":
        factor *= 0.9
    return _BASE[scale] * factor


def _temperature_factor(chemistry: str, temperature_c: float) -> float:
    knee, slope = _KNEE[chemistry]
    if temperature_c <= knee:
        return 1.0
    return math.exp(slope * (temperature_c - knee))


def _ionic_strength_factor(tds_mg_l: float) -> float:
    return max(0.5, 1.0 + 0.20 * math.log10(max(1.0, tds_mg_l) / 50000.0))


def _saturation_factor(saturation_ratio: float) -> float:
    return max(0.5, 1.0 + 0.6 * math.log10(max(1.0, saturation_ratio) / 5.0))


def _calcium_factor(chemistry: str, calcium_mg_l: float) -> float:
    if chemistry == "phosphonate" and calcium_mg_l > 1000.0:
        return 1.0 + 0.5 * math.log10(calcium_mg_l / 1000.0)
    return 1.0


def scale_inhibitor_performance(
    scale_type: str,
    inhibitor_chemistry: str,
    temperature: Q,
    saturation_ratio: float,
    total_dissolved_solids: Q,
    calcium: Q,
    available_dose: Q,
) -> ScaleInhibitorPerformanceResult:
    """The screening performance of an inhibitor at the stated dose.

    Args:
        scale_type: the mineral the inhibitor is dosed against.
        inhibitor_chemistry: the inhibitor family.
        temperature: the operating temperature.
        saturation_ratio: ``SR = IAP/Ksp``, a ratio and not its logarithm.
        total_dissolved_solids: the brine's TDS.
        calcium: calcium concentration, which carries the phosphonate penalty.
        available_dose: the active dose; zero reports zero efficiency.

    Returns:
        The MIC, the recommended dose, the efficiency and the class's four warnings.

    Raises:
        OutOfRangeError: if an input is outside its declared range.
        InvalidInputError: for an unknown scale type or chemistry.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = scale_inhibitor_performance(
        ...     "caco3", "phosphonate", q(333.15, "K"), 10.0,
        ...     q(50.0, "kg/m**3"), q(1.0, "kg/m**3"), q(0.0, "kg/m**3"))
        >>> abs(r.minimum_inhibitor_concentration.magnitude - 0.005903089986991944) < 1e-15
        True
    """
    from azoth._registry_gen import spec as _spec_for

    spec: dict[str, Any] = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    if scale_type not in _BASE:
        raise InvalidInputError(
            "scale_type",
            f"`{scale_type}` is not one of `caco3`, `baso4`, `srso4`, `caso4` or `feco3`",
        )
    if inhibitor_chemistry not in _CHEMISTRY_FACTOR:
        raise InvalidInputError(
            "inhibitor_chemistry",
            f"`{inhibitor_chemistry}` is not one of `phosphonate`, `polymaleate`, "
            "`polyacrylate`, `phosphate_ester` or `vinyl_sulphonate`",
        )

    values = {
        "temperature": input_to_si(spec, "temperature", temperature),
        "saturation_ratio": saturation_ratio,
        "total_dissolved_solids": input_to_si(
            spec, "total_dissolved_solids", total_dissolved_solids
        ),
        "calcium": input_to_si(spec, "calcium", calcium),
        "available_dose": input_to_si(spec, "available_dose", available_dose),
    }
    apply_checks(checks.on_input, values.get, warnings)

    temperature_c = values["temperature"] - 273.15
    tds_mg_l = values["total_dissolved_solids"] * 1000.0
    calcium_mg_l = values["calcium"] * 1000.0
    available_dose_mg_l = values["available_dose"] * 1000.0

    minimum = (
        _base_mic(scale_type, inhibitor_chemistry)
        * _temperature_factor(inhibitor_chemistry, temperature_c)
        * _ionic_strength_factor(tds_mg_l)
        * _saturation_factor(saturation_ratio)
        * _calcium_factor(inhibitor_chemistry, calcium_mg_l)
    )
    recommended = 1.5 * minimum

    if available_dose_mg_l > 0.0:
        ratio = available_dose_mg_l / minimum
        if ratio >= 1.0:
            efficiency = min(1.0, 0.85 + 0.15 * (1.0 - math.exp(-(ratio - 1.0))))
        else:
            efficiency = 0.85 * ratio
        adequate = ratio >= 1.0
    else:
        efficiency = 0.0
        adequate = False

    return ScaleInhibitorPerformanceResult(
        minimum_inhibitor_concentration=from_si(minimum * 1e-3, "kg/m**3"),
        recommended_dose=from_si(recommended * 1e-3, "kg/m**3"),
        efficiency=efficiency,
        adequate=adequate,
        thermal_degradation_risk=inhibitor_chemistry == "phosphonate" and temperature_c > 175.0,
        calcium_phosphonate_precipitation_risk=inhibitor_chemistry == "phosphonate"
        and calcium_mg_l > 2000.0,
        extreme_supersaturation_risk=saturation_ratio > 100.0,
        low_temperature_solubility_risk=inhibitor_chemistry == "phosphate_ester"
        and temperature_c < 40.0,
        warnings=tuple(warnings),
    )
