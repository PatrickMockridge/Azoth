"""``eos.scale_prediction`` - NeqSim's ``ScalePredictionCalculator``.

The water-analysis route to mineral-scale saturation: an ionic strength from the TDS, a
Davies activity coefficient from a Debye-Huckel A built on water's density and dielectric,
six ion-pair corrections, and five saturation indices over the class's own Ksp set.

Spec: ``specs/calcs/eos/scale_prediction.toml``.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ScalePredictionResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.scale_prediction"

#: The class's molar masses, in mg/mol, and its two sodium constants.
MW_CALCIUM = 40078.0
MW_BARIUM = 137327.0
MW_STRONTIUM = 87620.0
MW_IRON = 55845.0
MW_MAGNESIUM = 24305.0
MW_SODIUM = 22990.0
MW_SODIUM_GETTER = 22989.77
MW_BICARBONATE = 61017.0
MW_SULPHATE = 96060.0


def _log10(value: float) -> float:
    """`log10`, NaN where the value is not positive - Rust's own answer, not a raise."""
    return math.log10(value) if value > 0.0 else math.nan


def _log(value: float) -> float:
    """`ln`, NaN where the value is not positive - Rust's own answer, not a raise."""
    return math.log(value) if value > 0.0 else math.nan


def _sqrt(value: float) -> float:
    """`sqrt`, NaN where the value is negative - Rust's own answer, not a raise."""
    return math.sqrt(value) if value >= 0.0 else math.nan


def _pow15(value: float) -> float:
    """`value**1.5`, NaN where the value is negative - Rust's own answer, not a raise."""
    return math.pow(value, 1.5) if value >= 0.0 else math.nan


def _debye_huckel_a(temperature_k: float) -> float:
    tc = temperature_k - 273.15
    rho = (
        0.99983
        + 5.0948e-5 * tc
        - 7.5722e-6 * tc * tc
        + 3.8907e-8 * tc * tc * tc
        - 1.2e-10 * tc * tc * tc * tc
    )
    rho = max(rho, 0.85)
    eps = 87.740 - 0.40008 * tc + 9.398e-4 * tc * tc - 1.410e-6 * tc * tc * tc
    eps = max(eps, 20.0)
    return 1.8246e6 * _sqrt(rho) / _pow15(eps * temperature_k)


def _activity_coefficient(charge: int, ionic_strength: float, deb: float) -> float:
    sqrt_i = _sqrt(ionic_strength)
    log_gamma = -deb * charge * charge * (sqrt_i / (1.0 + sqrt_i) - 0.3 * ionic_strength)
    return math.pow(10.0, log_gamma)


def _k1_carbonate(temperature_k: float) -> float:
    log_k1 = (
        -356.3094
        - 0.06091964 * temperature_k
        + 21834.37 / temperature_k
        + 126.8339 * _log10(temperature_k)
        - 1684915.0 / (temperature_k * temperature_k)
    )
    return math.pow(10.0, log_k1)


def _k2_carbonate(temperature_k: float) -> float:
    log_k2 = (
        -107.8871
        - 0.03252849 * temperature_k
        + 5151.79 / temperature_k
        + 38.92561 * _log10(temperature_k)
        - 563713.9 / (temperature_k * temperature_k)
    )
    return math.pow(10.0, log_k2)


def _henry_co2(temperature_k: float) -> float:
    log_kh = (
        -6.8346
        + 1684.88 / temperature_k
        + 21.6215 * _log10(temperature_k)
        - 0.012174 * temperature_k
    )
    return math.pow(10.0, log_kh)


def _pressure_correction(
    pressure_bara: float,
    second_order: bool,
    ksp0: float,
    temperature_k: float,
    volume_delta: float,
    kappa_delta: float,
) -> float:
    if pressure_bara <= 1.013 or abs(volume_delta) < 1.0e-10:
        return ksp0
    r_cm3_bar = 83.1446
    delta_p = pressure_bara - 1.01325
    ln_correction = -volume_delta * delta_p / (r_cm3_bar * temperature_k)
    if second_order and abs(kappa_delta) > 1.0e-12:
        ln_correction += 0.5 * kappa_delta * delta_p * delta_p / (r_cm3_bar * temperature_k)
    return ksp0 * math.exp(ln_correction)


def _ksp(temperature_k: float, pressure_bara: float, second_order: bool) -> list[float]:
    calcite = (
        -171.9065
        - 0.077993 * temperature_k
        + 2839.319 / temperature_k
        + 71.595 * _log10(temperature_k)
    )
    barite = 136.035 - 7680.41 / temperature_k - 48.595 * _log10(temperature_k)
    celestite = -10452.9 / temperature_k + 19.790
    anhydrite = -19966.8 / temperature_k + 454.860 - 69.840 * _log(temperature_k)
    siderite = (
        -59.3498
        - 0.041377 * temperature_k
        - 2.1963 / temperature_k
        + 24.5724 * _log10(temperature_k)
    )
    return [
        _pressure_correction(
            pressure_bara, second_order, math.pow(10.0, calcite), temperature_k, -58.4, -8.0e-3
        ),
        _pressure_correction(
            pressure_bara, second_order, math.pow(10.0, barite), temperature_k, -46.4, -6.0e-3
        ),
        _pressure_correction(
            pressure_bara, second_order, math.exp(celestite), temperature_k, -47.0, -6.5e-3
        ),
        _pressure_correction(
            pressure_bara, second_order, math.exp(anhydrite), temperature_k, -52.4, -7.0e-3
        ),
        _pressure_correction(
            pressure_bara, second_order, math.pow(10.0, siderite), temperature_k, -52.9, -7.5e-3
        ),
    ]


def _ion_pairing(
    c_ca: float,
    c_mg: float,
    c_na: float,
    c_so4: float,
    c_co3: float,
    c_hco3: float,
    temperature_k: float,
    gamma_divalent: float,
    gamma_monovalent: float,
) -> tuple[float, float, float, float]:
    def correction(log_k_25: float, enthalpy: float) -> float:
        return math.pow(
            10.0,
            log_k_25 + enthalpy / (2.303 * 8.314) * (1.0 / 298.15 - 1.0 / temperature_k),
        )

    k_caso4 = correction(2.31, 6100.0)
    k_mgso4 = correction(2.37, 5600.0)
    k_naso4 = correction(0.70, 3000.0)
    k_cahco3 = correction(1.11, 3600.0)
    k_mghco3 = correction(1.16, 2900.0)
    k_caco3 = correction(3.22, 14000.0)

    free_ca, free_mg, free_na = c_ca, c_mg, c_na
    free_so4, free_co3, free_hco3 = c_so4, c_co3, c_hco3
    g2 = gamma_divalent * gamma_divalent
    gm2 = gamma_monovalent * gamma_monovalent

    for _ in range(5):
        free_ca = c_ca / (
            1.0
            + k_caso4 * g2 * free_so4 * g2
            + k_cahco3 * g2 * free_hco3 * gamma_monovalent / gamma_monovalent
            + k_caco3 * g2 * free_co3 * g2
        )
        if free_ca < 0.0:
            free_ca = c_ca * 0.01
        free_mg = c_mg / (
            1.0 + k_mgso4 * g2 * free_so4 * g2 + k_mghco3 * g2 * free_hco3 * gm2 / gm2
        )
        if free_mg < 0.0:
            free_mg = c_mg * 0.01
        free_na = c_na / (1.0 + k_naso4 * gamma_monovalent * free_so4 * g2 / gamma_monovalent)
        if free_na < 0.0:
            free_na = c_na * 0.01
        free_so4 = c_so4 / (
            1.0
            + k_caso4 * free_ca * g2 * g2
            + k_mgso4 * free_mg * g2 * g2
            + k_naso4 * free_na * gamma_monovalent * g2 / gamma_monovalent
        )
        if free_so4 < 0.0:
            free_so4 = c_so4 * 0.01
        free_co3 = c_co3 / (1.0 + k_caco3 * free_ca * g2 * g2)
        if free_co3 < 0.0:
            free_co3 = c_co3 * 0.01
        free_hco3 = c_hco3 / (
            1.0
            + k_cahco3 * free_ca * g2 * gamma_monovalent / gamma_monovalent
            + k_mghco3 * free_mg * g2
        )
        if free_hco3 < 0.0:
            free_hco3 = c_hco3 * 0.01
    return free_ca, free_so4, free_co3, free_hco3


def _saturation_index(iap: float, ksp: float) -> float:
    return math.log10(iap / ksp) if iap > 0.0 and ksp > 0.0 else math.nan


def scale_prediction(
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
    auto_ph: bool,
    second_order_pressure_correction: bool,
) -> ScalePredictionResult:
    """The mineral-scale saturation indices of a water analysis.

    Args:
        calcium: calcium ion concentration.
        barium: barium ion concentration.
        strontium: strontium ion concentration.
        iron: ferrous iron concentration.
        magnesium: magnesium ion concentration.
        sodium: sodium ion concentration.
        bicarbonate: bicarbonate concentration.
        sulphate: sulphate ion concentration.
        total_dissolved_solids: TDS, which set the ionic strength.
        temperature: absolute temperature.
        pressure: absolute pressure.
        co2_partial_pressure: CO2 partial pressure, read by the auto-pH path.
        ph: the water's pH, ignored when ``auto_ph`` is on and its gate holds.
        auto_ph: whether to estimate the pH, **with the class's own Henry constant**.
        second_order_pressure_correction: whether to add the compressibility term.

    Returns:
        The five indices and every intermediate the class exposes.

    Raises:
        OutOfRangeError: if an input is outside its declared range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = scale_prediction(
        ...     q(0.4, "kg/m**3"), q(0.0, "kg/m**3"), q(0.0, "kg/m**3"), q(0.0, "kg/m**3"),
        ...     q(0.0, "kg/m**3"), q(0.0, "kg/m**3"), q(0.15, "kg/m**3"), q(0.01, "kg/m**3"),
        ...     q(35.0, "kg/m**3"), q(353.15, "K"), q(1.0e7, "Pa"), q(2.0e5, "Pa"),
        ...     6.5, False, False)
        >>> abs(r.si_caco3.magnitude - -0.48139366009468343) < 1e-12
        True
    """
    from azoth._registry_gen import spec as _spec_for

    spec: dict[str, Any] = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    mass = {
        "calcium": calcium,
        "barium": barium,
        "strontium": strontium,
        "iron": iron,
        "magnesium": magnesium,
        "sodium": sodium,
        "bicarbonate": bicarbonate,
        "sulphate": sulphate,
        "total_dissolved_solids": total_dissolved_solids,
    }
    values: dict[str, float] = {
        name: input_to_si(spec, name, value) for name, value in mass.items()
    }
    temperature_k = input_to_si(spec, "temperature", temperature)
    pressure_si = input_to_si(spec, "pressure", pressure)
    co2_si = input_to_si(spec, "co2_partial_pressure", co2_partial_pressure)
    values.update(
        {
            "temperature": temperature_k,
            "pressure": pressure_si,
            "co2_partial_pressure": co2_si,
            "ph": ph,
        }
    )
    apply_checks(checks.on_input, values.get, warnings)

    state = _state(
        values,
        temperature_k,
        pressure_si,
        co2_si,
        ph,
        auto_ph,
        second_order_pressure_correction,
    )
    ksp = state["ksp"]
    si = state["si"]
    risk_count = sum(1 for value in si if math.isfinite(value) and value > 0.0)

    def concentration(value_mol_l: float) -> Q:
        return from_si(value_mol_l * 1000.0, "mol/m**3")

    return ScalePredictionResult(
        ionic_strength=from_si(state["ionic_strength"] * 1000.0, "mol/m**3"),
        debye_huckel_a=state["debye_huckel_a"],
        divalent_activity_coefficient=state["gamma_divalent"],
        monovalent_activity_coefficient=state["gamma_monovalent"],
        ph_used=state["ph_used"],
        free_calcium=concentration(state["free_calcium"]),
        free_sulphate=concentration(state["free_sulphate"]),
        free_carbonate=concentration(state["free_carbonate"]),
        free_bicarbonate=concentration(state["free_bicarbonate"]),
        total_barium=concentration(state["total_barium"]),
        total_strontium=concentration(state["total_strontium"]),
        total_iron=concentration(state["total_iron"]),
        total_sodium=concentration(state["total_sodium"]),
        ksp_calcite=ksp[0],
        ksp_barite=ksp[1],
        ksp_celestite=ksp[2],
        ksp_anhydrite=ksp[3],
        ksp_siderite=ksp[4],
        log10_ksp_calcite=_log10(ksp[0]),
        log10_ksp_barite=_log10(ksp[1]),
        log10_ksp_celestite=_log10(ksp[2]),
        log10_ksp_anhydrite=_log10(ksp[3]),
        log10_ksp_siderite=_log10(ksp[4]),
        si_caco3=None if math.isnan(si[0]) else si[0],
        si_baso4=None if math.isnan(si[1]) else si[1],
        si_srso4=None if math.isnan(si[2]) else si[2],
        si_caso4=None if math.isnan(si[3]) else si[3],
        si_feco3=None if math.isnan(si[4]) else si[4],
        has_scaling_risk=risk_count > 0,
        scale_risk_count=risk_count,
        warnings=tuple(warnings),
    )


def _state(
    values: dict[str, float],
    temperature_k: float,
    pressure_si: float,
    co2_si: float,
    ph: float,
    auto_ph: bool,
    second_order_pressure_correction: bool,
) -> dict[str, Any]:
    """The class's `calculate()` in SI magnitudes, without the input range checks.

    The composers call this directly: they run it over mixtures and interpolated points whose
    values are not the caller's inputs, and the class itself applies no range check anywhere.
    """
    ionic_strength = values["total_dissolved_solids"] * 1000.0 / 40000.0
    deb = _debye_huckel_a(temperature_k)
    gamma_divalent = _activity_coefficient(2, ionic_strength, deb)
    gamma_monovalent = _activity_coefficient(1, ionic_strength, deb)

    pressure_bara = pressure_si / 1.0e5
    co2_bar = co2_si / 1.0e5
    ph_used = ph
    if auto_ph and co2_bar > 0.0 and values["bicarbonate"] > 0.0:
        co2_aq = _henry_co2(temperature_k) * co2_bar
        c_hco3 = values["bicarbonate"] * 1000.0 / MW_BICARBONATE
        if c_hco3 > 0.0 and co2_aq > 0.0:
            ph_used = -_log10(_k1_carbonate(temperature_k) * co2_aq / c_hco3)
        else:
            ph_used = 7.0

    def mol_per_l(name: str, molar_mass: float) -> float:
        return values[name] * 1000.0 / molar_mass

    c_ca = mol_per_l("calcium", MW_CALCIUM)
    c_ba = mol_per_l("barium", MW_BARIUM)
    c_sr = mol_per_l("strontium", MW_STRONTIUM)
    c_fe = mol_per_l("iron", MW_IRON)
    c_mg = mol_per_l("magnesium", MW_MAGNESIUM)
    c_na = mol_per_l("sodium", MW_SODIUM)
    c_hco3 = mol_per_l("bicarbonate", MW_BICARBONATE)
    c_so4 = mol_per_l("sulphate", MW_SULPHATE)
    h = math.pow(10.0, -ph_used)
    c_co3 = _k2_carbonate(temperature_k) * c_hco3 / h if h > 0.0 else 0.0

    free_ca, free_so4, free_co3, free_hco3 = _ion_pairing(
        c_ca, c_mg, c_na, c_so4, c_co3, c_hco3, temperature_k, gamma_divalent, gamma_monovalent
    )
    ksp = _ksp(temperature_k, pressure_bara, second_order_pressure_correction)

    si = [
        _saturation_index((free_ca * gamma_divalent) * (free_co3 * gamma_divalent), ksp[0]),
        _saturation_index((c_ba * gamma_divalent) * (free_so4 * gamma_divalent), ksp[1])
        if values["barium"] > 0.0 and values["sulphate"] > 0.0
        else math.nan,
        _saturation_index((c_sr * gamma_divalent) * (free_so4 * gamma_divalent), ksp[2])
        if values["strontium"] > 0.0 and values["sulphate"] > 0.0
        else math.nan,
        _saturation_index((free_ca * gamma_divalent) * (free_so4 * gamma_divalent), ksp[3])
        if values["sulphate"] > 0.0
        else math.nan,
        _saturation_index((c_fe * gamma_divalent) * (free_co3 * gamma_divalent), ksp[4])
        if values["iron"] > 0.0
        else math.nan,
    ]
    return {
        "ionic_strength": ionic_strength,
        "debye_huckel_a": deb,
        "gamma_divalent": gamma_divalent,
        "gamma_monovalent": gamma_monovalent,
        "ph_used": ph_used,
        "free_calcium": free_ca,
        "free_sulphate": free_so4,
        "free_carbonate": free_co3,
        "free_bicarbonate": free_hco3,
        "total_barium": c_ba,
        "total_strontium": c_sr,
        "total_iron": c_fe,
        "total_sodium": values["sodium"] * 1000.0 / MW_SODIUM_GETTER,
        "ksp": ksp,
        "si": si,
    }
