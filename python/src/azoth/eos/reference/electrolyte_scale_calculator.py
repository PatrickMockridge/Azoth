"""``eos.electrolyte_scale_calculator`` - NeqSim's ``ElectrolyteScaleCalculator``.

An ionic strength from the water's ions, Davies activity coefficients over a
temperature-adjusted Debye-Huckel A, and five saturation indices - calcite, barite,
gypsum, celestite and brucite - over the class's own temperature correlations.

Spec: ``specs/calcs/eos/electrolyte_scale_calculator.toml``.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ElectrolyteScaleCalculatorResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.electrolyte_scale_calculator"

#: The class's molar masses, in g/mol, as written upstream.
MM_CALCIUM = 40.08
MM_BARIUM = 137.33
MM_STRONTIUM = 87.62
MM_MAGNESIUM = 24.31
MM_SODIUM = 22.99
MM_POTASSIUM = 39.10
MM_IRON = 55.85
MM_CHLORIDE = 35.45
MM_SULPHATE = 96.06
MM_BICARBONATE = 61.02
MM_CARBONATE = 60.01


def _log10(value: float) -> float:
    """`log10`, NaN where the value is not positive - Rust's own answer, not a raise."""
    return math.log10(value) if value > 0.0 else math.nan


def _davies_log_gamma(charge: int, ionic_strength: float, debye_hueckel_a: float) -> float:
    if ionic_strength <= 0.0:
        return 0.0
    sqrt_i = math.sqrt(ionic_strength)
    return -debye_hueckel_a * charge * charge * (sqrt_i / (1.0 + sqrt_i) - 0.3 * ionic_strength)


def _saturation_index(activity_m: float, activity_x: float, p_ksp: float) -> float:
    if activity_m <= 0.0 or activity_x <= 0.0:
        return -10.0
    return _log10(activity_m * activity_x) + p_ksp


def _p_ksp(pressure_bara: float, temperature_k: float, ksp0: float, volume_delta: float) -> float:
    """`-log10(Ksp)` after the class's pressure correction."""
    if pressure_bara <= 1.013 or abs(volume_delta) < 1.0e-10:
        corrected = ksp0
    else:
        r_cm3_bar = 83.1446
        delta_p = pressure_bara - 1.01325
        corrected = ksp0 * math.exp(-volume_delta * delta_p / (r_cm3_bar * temperature_k))
    return -_log10(corrected)


def _caco3_ksp(temperature_k: float) -> float:
    log10_ksp = (
        -171.9065
        - 0.077993 * temperature_k
        + 2839.319 / temperature_k
        + 71.595 * _log10(temperature_k)
    )
    return math.pow(10.0, log10_ksp)


def _baso4_ksp(temperature_k: float) -> float:
    ln_ksp = (
        37588.0 / temperature_k
        - 747.610
        + math.log(temperature_k) * 119.280
        - temperature_k * 0.162830
        - 2880000.0 / (temperature_k * temperature_k)
    )
    return math.exp(ln_ksp)


def _caso4_ksp(temperature_k: float) -> float:
    ln_ksp = (
        -26309.9 / temperature_k
        + 815.978
        + math.log(temperature_k) * -138.361
        + temperature_k * 0.167863
        + 18.6143 / (temperature_k * temperature_k)
    )
    return math.exp(ln_ksp)


def _srso4_ksp(temperature_k: float) -> float:
    ln_ksp = -10452.9 / temperature_k + 19.790
    return math.exp(ln_ksp)


def _brucite_ksp(temperature_k: float) -> float:
    log10_ksp25 = -11.16
    delta_h_kj = -113.4
    gas_constant_kj = 0.0083145
    log10_ksp = log10_ksp25 - (delta_h_kj / (2.302585093 * gas_constant_kj)) * (
        1.0 / temperature_k - 1.0 / 298.15
    )
    return math.pow(10.0, log10_ksp)


def _p_kw(temperature_k: float) -> float:
    return 4470.99 / temperature_k - 6.0875 + 0.01706 * temperature_k


def _state(values: dict[str, float]) -> dict[str, float]:
    """`calculate()` in SI magnitudes, without the input range checks.

    The brine-mixing model calls this over its blends, whose values are not the caller's
    inputs; the class itself applies no range check anywhere.
    """
    m_ca = values["calcium"] * 1000.0 / MM_CALCIUM / 1000.0
    m_ba = values["barium"] * 1000.0 / MM_BARIUM / 1000.0
    m_sr = values["strontium"] * 1000.0 / MM_STRONTIUM / 1000.0
    m_mg = values["magnesium"] * 1000.0 / MM_MAGNESIUM / 1000.0
    m_na = values["sodium"] * 1000.0 / MM_SODIUM / 1000.0
    m_k = values["potassium"] * 1000.0 / MM_POTASSIUM / 1000.0
    m_fe = values["iron"] * 1000.0 / MM_IRON / 1000.0
    m_cl = values["chloride"] * 1000.0 / MM_CHLORIDE / 1000.0
    m_so4 = values["sulphate"] * 1000.0 / MM_SULPHATE / 1000.0
    m_hco3 = values["bicarbonate"] * 1000.0 / MM_BICARBONATE / 1000.0
    m_co3 = values["carbonate"] * 1000.0 / MM_CARBONATE / 1000.0

    ionic_strength = 0.5 * (
        m_ca * 4.0
        + m_ba * 4.0
        + m_sr * 4.0
        + m_mg * 4.0
        + m_na
        + m_k
        + m_fe * 4.0
        + m_cl
        + m_so4 * 4.0
        + m_hco3
        + m_co3 * 4.0
    )

    temperature_k = values["temperature"]
    pressure_bara = values["pressure"] / 1.0e5
    ph = values["ph"]

    debye_hueckel_a = 0.5092 + 8.5e-4 * (temperature_k - 273.15 - 25.0)
    if debye_hueckel_a < 0.5:
        debye_hueckel_a = 0.5

    gamma_divalent = math.pow(10.0, _davies_log_gamma(2, ionic_strength, debye_hueckel_a))
    gamma_monovalent = math.pow(10.0, _davies_log_gamma(1, ionic_strength, debye_hueckel_a))

    p_k2 = 10.329 + 0.013 * (temperature_k - 298.15) / 25.0
    m_co3_eff = m_hco3 * math.pow(10.0, ph - p_k2)

    p_ksp_caco3 = _p_ksp(pressure_bara, temperature_k, _caco3_ksp(temperature_k), -58.4)
    p_ksp_baso4 = _p_ksp(pressure_bara, temperature_k, _baso4_ksp(temperature_k), -46.4)
    p_ksp_caso4 = _p_ksp(pressure_bara, temperature_k, _caso4_ksp(temperature_k), -33.0)
    p_ksp_srso4 = _p_ksp(pressure_bara, temperature_k, _srso4_ksp(temperature_k), -47.0)

    si_caco3 = _saturation_index(
        gamma_divalent * max(m_ca, 1e-15),
        gamma_divalent * max(m_co3_eff + m_co3, 1e-15),
        p_ksp_caco3,
    )
    si_baso4 = _saturation_index(
        gamma_divalent * max(m_ba, 1e-15),
        gamma_divalent * max(m_so4, 1e-15),
        p_ksp_baso4,
    )
    si_caso4 = _saturation_index(
        gamma_divalent * max(m_ca, 1e-15),
        gamma_divalent * max(m_so4, 1e-15),
        p_ksp_caso4,
    )
    si_srso4 = _saturation_index(
        gamma_divalent * max(m_sr, 1e-15),
        gamma_divalent * max(m_so4, 1e-15),
        p_ksp_srso4,
    )

    hydroxide_activity = math.pow(10.0, -_p_kw(temperature_k) + ph)
    p_ksp_brucite = _p_ksp(pressure_bara, temperature_k, _brucite_ksp(temperature_k), -20.0)
    si_brucite = _saturation_index(
        gamma_divalent * max(m_mg, 1e-15),
        hydroxide_activity * hydroxide_activity,
        p_ksp_brucite,
    )

    return {
        "ionic_strength": ionic_strength,
        "debye_hueckel_a": debye_hueckel_a,
        "gamma_divalent": gamma_divalent,
        "gamma_monovalent": gamma_monovalent,
        "si_caco3": si_caco3,
        "si_baso4": si_baso4,
        "si_caso4": si_caso4,
        "si_srso4": si_srso4,
        "si_brucite": si_brucite,
        "hydroxide_activity": hydroxide_activity,
    }


def electrolyte_scale_calculator(
    calcium: Q,
    barium: Q,
    strontium: Q,
    magnesium: Q,
    sodium: Q,
    potassium: Q,
    iron: Q,
    chloride: Q,
    sulphate: Q,
    bicarbonate: Q,
    carbonate: Q,
    temperature: Q,
    pressure: Q,
    ph: float,
) -> ElectrolyteScaleCalculatorResult:
    """The activity-corrected saturation indices of a brine.

    Args:
        calcium: calcium ion concentration.
        barium: barium ion concentration.
        strontium: strontium ion concentration.
        magnesium: magnesium ion concentration.
        sodium: sodium ion concentration.
        potassium: potassium ion concentration.
        iron: ferrous iron concentration.
        chloride: chloride ion concentration.
        sulphate: sulphate ion concentration.
        bicarbonate: bicarbonate concentration, which sets carbonate through ``pK2``.
        carbonate: carbonate concentration, which adds to the derived carbonate.
        temperature: absolute temperature.
        pressure: absolute pressure, used by the Ksp pressure correction.
        ph: the water's pH, which sets carbonate and the hydroxide activity.

    Returns:
        The five indices and the class's ionic strength and coefficients.

    Raises:
        OutOfRangeError: if an input is outside its declared range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = electrolyte_scale_calculator(
        ...     q(0.412, "kg/m**3"), q(0.00002, "kg/m**3"), q(0.0079, "kg/m**3"),
        ...     q(1.29, "kg/m**3"), q(10.78, "kg/m**3"), q(0.399, "kg/m**3"),
        ...     q(0.000002, "kg/m**3"), q(19.35, "kg/m**3"), q(2.71, "kg/m**3"),
        ...     q(0.142, "kg/m**3"), q(0.0, "kg/m**3"), q(283.15, "K"),
        ...     q(1.013e5, "Pa"), 10.5)
        >>> abs(r.si_brucite - 0.27409434256769494) < 1e-12
        True
    """
    from azoth._registry_gen import spec as _spec_for

    spec: dict[str, Any] = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    concentrations = {
        "calcium": calcium,
        "barium": barium,
        "strontium": strontium,
        "magnesium": magnesium,
        "sodium": sodium,
        "potassium": potassium,
        "iron": iron,
        "chloride": chloride,
        "sulphate": sulphate,
        "bicarbonate": bicarbonate,
        "carbonate": carbonate,
    }
    values: dict[str, float] = {
        name: input_to_si(spec, name, value) for name, value in concentrations.items()
    }
    values["temperature"] = input_to_si(spec, "temperature", temperature)
    values["pressure"] = input_to_si(spec, "pressure", pressure)
    values["ph"] = ph
    apply_checks(checks.on_input, values.get, warnings)

    state = _state(values)
    return ElectrolyteScaleCalculatorResult(
        ionic_strength=from_si(state["ionic_strength"], "mol/kg"),
        debye_hueckel_a=state["debye_hueckel_a"],
        divalent_activity_coefficient=state["gamma_divalent"],
        monovalent_activity_coefficient=state["gamma_monovalent"],
        si_caco3=state["si_caco3"],
        si_baso4=state["si_baso4"],
        si_caso4=state["si_caso4"],
        si_srso4=state["si_srso4"],
        si_brucite=state["si_brucite"],
        hydroxide_activity=state["hydroxide_activity"],
        warnings=tuple(warnings),
    )
