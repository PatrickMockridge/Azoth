//! `eos.de_waard_milliams_corrosion` - NeqSim's `DeWaardMilliamsCorrosion`.
//!
//! The de Waard-Milliams 1991 CO2 rate with the de Waard-Lotz pH, scale, glycol and flow
//! corrections, plus the elemental-sulfur and H2S/FeS contributions and the NORSOK M-001
//! severity band.

use azoth_core::units::{
    MassFlux, Pressure, ThermodynamicTemperature, Time, Velocity, meters, meters_per_second,
};
use azoth_core::{Result, apply_checks};

use crate::results::{CorrosionSeverity, DeWaardMilliamsCorrosionResult, FeSMorphology};
use crate::spec_gen;

/// Millimetres per year as metres per second, over the class's 365.25-day year.
const MM_PER_YEAR_TO_METRES_PER_SECOND: f64 = 1.0e-3 / 31_557_600.0;
/// Seconds in the class's year, for the sulfur input and the design life.
const SECONDS_PER_YEAR: f64 = 31_557_600.0;
/// Bars per pascal, for the correlations written in bar.
const BAR_PER_PASCAL: f64 = 1.0e-5;

/// The de Waard-Milliams CO2/H2S/S8 corrosion rate and its corrections.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's hard range.
///
/// # Example
/// ```
/// use azoth_core::units::{
///     kelvins, kilograms_per_square_meter_second, meters_per_second, pascals, seconds,
/// };
/// use azoth_eos::de_waard_milliams_corrosion::de_waard_milliams_corrosion;
///
/// let r = de_waard_milliams_corrosion(
///     kelvins(333.15),
///     pascals(2.0e5),
///     4.5,
///     pascals(0.0),
///     meters_per_second(2.0),
///     0.0,
///     0.0,
///     kilograms_per_square_meter_second(0.0),
///     seconds(788_940_000.0),
/// )?;
/// assert!(r.corrosion_severity.as_str() == "very_high");
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn de_waard_milliams_corrosion(
    temperature: ThermodynamicTemperature,
    co2_partial_pressure: Pressure,
    ph: f64,
    h2s_partial_pressure: Pressure,
    flow_velocity: Velocity,
    inhibitor_efficiency: f64,
    glycol_fraction: f64,
    sulfur_deposition_rate: MassFlux,
    design_life: Time,
) -> Result<DeWaardMilliamsCorrosionResult> {
    let spec = &spec_gen::DE_WAARD_MILLIAMS_CORROSION_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "temperature" => Some(temperature.value),
            "co2_partial_pressure" => Some(co2_partial_pressure.value),
            "h2s_partial_pressure" => Some(h2s_partial_pressure.value),
            "flow_velocity" => Some(flow_velocity.value),
            "design_life" => Some(design_life.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let temperature_c = temperature.value - 273.15;
    let temperature_k = temperature.value;
    let pco2 = co2_partial_pressure.value * BAR_PER_PASCAL;
    let ph2s = h2s_partial_pressure.value * BAR_PER_PASCAL;
    let velocity = flow_velocity.value;
    // numerics-ok: the class's setters clamp the efficiency and the glycol fraction to
    // [0, 1], so the clamp is the class's transcription rather than a default.
    let inhibitor = inhibitor_efficiency.clamp(0.0, 1.0);
    let glycol = glycol_fraction.clamp(0.0, 1.0);
    let sulfur_kg_m2_yr = sulfur_deposition_rate.value.max(0.0) * SECONDS_PER_YEAR;

    let baseline_mm_yr = if pco2 <= 0.0 {
        0.0
    } else {
        10.0_f64.powf(5.8 - 1710.0 / temperature_k + 0.67 * pco2.log10())
    };

    let log_pco2 = if pco2 > 0.0 { pco2.log10() } else { -2.0 };
    let ph_saturation = 3.71 - 0.5 * log_pco2 + 0.00133 * temperature_c;
    let ph_factor = if ph >= ph_saturation {
        10.0_f64.powf(0.32 * (ph_saturation - ph))
    } else {
        1.0
    };

    let scale_factor = if temperature_c > 80.0 && ph > 5.0 {
        let reduction = (0.7_f64).min((temperature_c - 80.0) * 0.01 * (ph - 4.0));
        (0.1_f64).max(1.0 - reduction)
    } else {
        1.0
    };

    let glycol_factor = if glycol <= 0.0 {
        1.0
    } else if glycol > 0.8 {
        0.05
    } else {
        1.0 - glycol
    };

    let flow_factor = if velocity <= 1.0 {
        1.0
    } else {
        velocity.powf(0.8).min(5.0)
    };

    let inhibitor_factor = 1.0 - inhibitor;
    let co2_mm_yr =
        baseline_mm_yr * ph_factor * scale_factor * glycol_factor * flow_factor * inhibitor_factor;
    let sulfur_mm_yr = sulfur_kg_m2_yr * 1.74 / 7850.0 * 1000.0;
    let fes_mm_yr = if ph2s <= 0.0 {
        0.0
    } else {
        (7.96 - 2320.0 / temperature_k + 0.67 * ph2s.ln())
            .exp()
            .min(100.0)
    };
    let corrected_mm_yr = co2_mm_yr + sulfur_mm_yr;
    let total_mm_yr = corrected_mm_yr + fes_mm_yr;

    let log_ionic_product = 2.0 * (ph - 7.0) + pco2.max(0.001).log10();
    let log_ksp = -10.5 + 2000.0 / temperature_k;
    let feco3_saturation_index = 10.0_f64.powf(log_ionic_product - log_ksp);

    let severity = if corrected_mm_yr < 0.1 {
        CorrosionSeverity::Low
    } else if corrected_mm_yr < 0.3 {
        CorrosionSeverity::Medium
    } else if corrected_mm_yr < 1.0 {
        CorrosionSeverity::High
    } else {
        CorrosionSeverity::VeryHigh
    };

    let morphology = if ph2s <= 0.0 && sulfur_kg_m2_yr <= 0.0 {
        FeSMorphology::None
    } else if temperature_c < 60.0 {
        FeSMorphology::Mackinawite
    } else if temperature_c < 120.0 {
        FeSMorphology::CubicFesTroilite
    } else {
        FeSMorphology::Pyrrhotite
    };

    let rate = |mm_yr: f64| meters_per_second(mm_yr * MM_PER_YEAR_TO_METRES_PER_SECOND);

    Ok(DeWaardMilliamsCorrosionResult {
        baseline_rate: rate(baseline_mm_yr),
        ph_correction_factor: ph_factor,
        scale_correction_factor: scale_factor,
        glycol_correction_factor: glycol_factor,
        flow_correction_factor: flow_factor,
        inhibitor_factor,
        corrected_rate: rate(corrected_mm_yr),
        sulfur_rate: rate(sulfur_mm_yr),
        fes_rate: rate(fes_mm_yr),
        co2_corrected_rate: rate(co2_mm_yr),
        total_corrosion_rate: rate(total_mm_yr),
        feco3_saturation_index,
        feco3_scale_forming: feco3_saturation_index > 1.0,
        is_sour_service: ph2s > 0.003,
        corrosion_severity: severity,
        fes_morphology: morphology,
        corrosion_allowance: meters(
            corrected_mm_yr * 1.0e-3 * design_life.value / SECONDS_PER_YEAR,
        ),
        warnings,
    })
}
