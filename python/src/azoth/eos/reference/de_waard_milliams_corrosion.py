"""``eos.de_waard_milliams_corrosion`` - NeqSim's ``DeWaardMilliamsCorrosion``.

The de Waard-Milliams 1991 CO2 rate with the de Waard-Lotz pH, scale, glycol and flow
corrections, plus the elemental-sulfur and H2S/FeS contributions and the NORSOK M-001
severity band.

Spec: ``specs/calcs/eos/de_waard_milliams_corrosion.toml``.
"""

from __future__ import annotations

import math

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import CorrosionSeverity, DeWaardMilliamsCorrosionResult, FeSMorphology
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.de_waard_milliams_corrosion"

SECONDS_PER_YEAR = 31_557_600.0
MM_PER_YEAR_TO_METRES_PER_SECOND = 1.0e-3 / SECONDS_PER_YEAR
BAR_PER_PASCAL = 1.0e-5


def de_waard_milliams_corrosion(
    temperature: Q,
    co2_partial_pressure: Q,
    ph: float,
    h2s_partial_pressure: Q,
    flow_velocity: Q,
    inhibitor_efficiency: float,
    glycol_fraction: float,
    sulfur_deposition_rate: Q,
    design_life: Q,
) -> DeWaardMilliamsCorrosionResult:
    """The de Waard-Milliams CO2/H2S/S8 corrosion rate and its corrections.

    Args:
        temperature: the operating temperature.
        co2_partial_pressure: the CO2 partial pressure, zero when no CO2 is present.
        ph: the aqueous phase's pH.
        h2s_partial_pressure: the H2S partial pressure, zero for a sweet service.
        flow_velocity: the bulk flow velocity.
        inhibitor_efficiency: the film-forming efficiency, clamped to 0..1.
        glycol_fraction: the glycol weight fraction, clamped to 0..1.
        sulfur_deposition_rate: the elemental-sulfur deposition rate.
        design_life: the exposure time the corrosion allowance is scaled by.

    Returns:
        The baseline and corrected rates, the four factors, the sulfur and FeS
        contributions, the FeCO3 index, the severity band, the morphology and the allowance.

    Raises:
        OutOfRangeError: if an input is outside the spec's hard range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = de_waard_milliams_corrosion(
        ...     q(333.15, "K"), q(2.0e5, "Pa"), 4.5, q(0.0, "Pa"),
        ...     q(2.0, "m/s"), 0.0, 0.0, q(0.0, "kg/(m**2*s)"),
        ...     q(788_940_000.0, "s"))
        >>> r.corrosion_severity
        <CorrosionSeverity.VERY_HIGH: 'very_high'>
    """
    from azoth._registry_gen import spec as _spec_for

    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values: dict[str, float] = {
        name: input_to_si(spec, name, value)
        for name, value in {
            "temperature": temperature,
            "co2_partial_pressure": co2_partial_pressure,
            "h2s_partial_pressure": h2s_partial_pressure,
            "flow_velocity": flow_velocity,
            "sulfur_deposition_rate": sulfur_deposition_rate,
            "design_life": design_life,
        }.items()
    }
    values["ph"] = ph
    values["inhibitor_efficiency"] = inhibitor_efficiency
    values["glycol_fraction"] = glycol_fraction
    apply_checks(checks.on_input, values.get, warnings)

    temperature_c = values["temperature"] - 273.15
    temperature_k = values["temperature"]
    pco2 = values["co2_partial_pressure"] * BAR_PER_PASCAL
    ph2s = values["h2s_partial_pressure"] * BAR_PER_PASCAL
    velocity = values["flow_velocity"]
    inhibitor = min(1.0, max(0.0, inhibitor_efficiency))
    glycol = min(1.0, max(0.0, glycol_fraction))
    sulfur_kg_m2_yr = max(0.0, values["sulfur_deposition_rate"]) * SECONDS_PER_YEAR

    baseline_mm_yr = (
        0.0 if pco2 <= 0.0 else 10.0 ** (5.8 - 1710.0 / temperature_k + 0.67 * math.log10(pco2))
    )

    log_pco2 = math.log10(pco2) if pco2 > 0.0 else -2.0
    ph_saturation = 3.71 - 0.5 * log_pco2 + 0.00133 * temperature_c
    ph_factor = 10.0 ** (0.32 * (ph_saturation - ph)) if ph >= ph_saturation else 1.0

    if temperature_c > 80.0 and ph > 5.0:
        reduction = min(0.7, (temperature_c - 80.0) * 0.01 * (ph - 4.0))
        scale_factor = max(0.1, 1.0 - reduction)
    else:
        scale_factor = 1.0

    if glycol <= 0.0:
        glycol_factor = 1.0
    elif glycol > 0.8:
        glycol_factor = 0.05
    else:
        glycol_factor = 1.0 - glycol

    flow_factor = 1.0 if velocity <= 1.0 else min(velocity**0.8, 5.0)
    inhibitor_factor = 1.0 - inhibitor

    co2_mm_yr = (
        baseline_mm_yr * ph_factor * scale_factor * glycol_factor * flow_factor * inhibitor_factor
    )
    sulfur_mm_yr = sulfur_kg_m2_yr * 1.74 / 7850.0 * 1000.0
    fes_mm_yr = (
        0.0
        if ph2s <= 0.0
        else min(math.exp(7.96 - 2320.0 / temperature_k + 0.67 * math.log(ph2s)), 100.0)
    )
    corrected_mm_yr = co2_mm_yr + sulfur_mm_yr
    total_mm_yr = corrected_mm_yr + fes_mm_yr

    log_ionic_product = 2.0 * (ph - 7.0) + math.log10(max(pco2, 0.001))
    log_ksp = -10.5 + 2000.0 / temperature_k
    feco3_saturation_index = 10.0 ** (log_ionic_product - log_ksp)

    if corrected_mm_yr < 0.1:
        severity = CorrosionSeverity.LOW
    elif corrected_mm_yr < 0.3:
        severity = CorrosionSeverity.MEDIUM
    elif corrected_mm_yr < 1.0:
        severity = CorrosionSeverity.HIGH
    else:
        severity = CorrosionSeverity.VERY_HIGH

    if ph2s <= 0.0 and sulfur_kg_m2_yr <= 0.0:
        morphology = FeSMorphology.NONE
    elif temperature_c < 60.0:
        morphology = FeSMorphology.MACKINAWITE
    elif temperature_c < 120.0:
        morphology = FeSMorphology.CUBIC_FES_TROILITE
    else:
        morphology = FeSMorphology.PYRRHOTITE

    def rate(mm_yr: float) -> Q:
        return from_si(mm_yr * MM_PER_YEAR_TO_METRES_PER_SECOND, "m/s")

    return DeWaardMilliamsCorrosionResult(
        baseline_rate=rate(baseline_mm_yr),
        ph_correction_factor=ph_factor,
        scale_correction_factor=scale_factor,
        glycol_correction_factor=glycol_factor,
        flow_correction_factor=flow_factor,
        inhibitor_factor=inhibitor_factor,
        corrected_rate=rate(corrected_mm_yr),
        sulfur_rate=rate(sulfur_mm_yr),
        fes_rate=rate(fes_mm_yr),
        co2_corrected_rate=rate(co2_mm_yr),
        total_corrosion_rate=rate(total_mm_yr),
        feco3_saturation_index=feco3_saturation_index,
        feco3_scale_forming=feco3_saturation_index > 1.0,
        is_sour_service=ph2s > 0.003,
        corrosion_severity=severity,
        fes_morphology=morphology,
        corrosion_allowance=from_si(
            corrected_mm_yr * 1.0e-3 * values["design_life"] / SECONDS_PER_YEAR, "m"
        ),
        warnings=tuple(warnings),
    )
