//! `eos.scale_inhibitor_performance` - NeqSim's `ScaleInhibitorPerformance`.
//!
//! An empirical minimum-inhibitor-concentration correlation: a base MIC per mineral and
//! chemistry, four multiplicative factors over temperature, TDS, saturation and calcium,
//! and the efficiency at the dose actually available.

use azoth_core::units::{MassDensity, ThermodynamicTemperature, kilograms_per_cubic_meter};
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::ScaleInhibitorPerformanceResult;
use crate::spec_gen;

#[derive(Clone, Copy, PartialEq)]
enum Mineral {
    Caco3,
    Baso4,
    Srso4,
    Caso4,
    Feco3,
}

#[derive(Clone, Copy, PartialEq)]
enum Chemistry {
    Phosphonate,
    Polymaleate,
    Polyacrylate,
    PhosphateEster,
    VinylSulphonate,
}

/// The screening performance of an inhibitor at the stated dose.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's range.
/// * [`azoth_core::AzothError::InvalidInput`] for an unknown scale type or chemistry.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, kilograms_per_cubic_meter};
/// use azoth_eos::scale_inhibitor_performance::scale_inhibitor_performance;
///
/// let r = scale_inhibitor_performance(
///     "caco3",
///     "phosphonate",
///     kelvins(333.15),
///     10.0,
///     kilograms_per_cubic_meter(50.0),
///     kilograms_per_cubic_meter(1.0),
///     kilograms_per_cubic_meter(0.0),
/// )?;
/// assert!((r.minimum_inhibitor_concentration.value - 0.005903089986991944).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn scale_inhibitor_performance(
    scale_type: &str,
    inhibitor_chemistry: &str,
    temperature: ThermodynamicTemperature,
    saturation_ratio: f64,
    total_dissolved_solids: MassDensity,
    calcium: MassDensity,
    available_dose: MassDensity,
) -> Result<ScaleInhibitorPerformanceResult> {
    let spec = &spec_gen::SCALE_INHIBITOR_PERFORMANCE_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "temperature" => Some(temperature.value),
            "saturation_ratio" => Some(saturation_ratio),
            "total_dissolved_solids" => Some(total_dissolved_solids.value),
            "calcium" => Some(calcium.value),
            "available_dose" => Some(available_dose.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let mineral = match scale_type {
        "caco3" => Mineral::Caco3,
        "baso4" => Mineral::Baso4,
        "srso4" => Mineral::Srso4,
        "caso4" => Mineral::Caso4,
        "feco3" => Mineral::Feco3,
        other => {
            return Err(AzothError::invalid_input(
                "scale_type",
                format!("`{other}` is not one of `caco3`, `baso4`, `srso4`, `caso4` or `feco3`"),
            ));
        }
    };
    let chemistry = match inhibitor_chemistry {
        "phosphonate" => Chemistry::Phosphonate,
        "polymaleate" => Chemistry::Polymaleate,
        "polyacrylate" => Chemistry::Polyacrylate,
        "phosphate_ester" => Chemistry::PhosphateEster,
        "vinyl_sulphonate" => Chemistry::VinylSulphonate,
        other => {
            return Err(AzothError::invalid_input(
                "inhibitor_chemistry",
                format!(
                    "`{other}` is not one of `phosphonate`, `polymaleate`, `polyacrylate`, \
                     `phosphate_ester` or `vinyl_sulphonate`"
                ),
            ));
        }
    };

    // The class's own units: mg/L and degrees Celsius.
    let temperature_c = temperature.value - 273.15;
    let tds_mg_l = total_dissolved_solids.value * 1000.0;
    let calcium_mg_l = calcium.value * 1000.0;
    let available_dose_mg_l = available_dose.value * 1000.0;

    let minimum = base_mic(mineral, chemistry)
        * temperature_factor(chemistry, temperature_c)
        * ionic_strength_factor(tds_mg_l)
        * saturation_factor(saturation_ratio)
        * calcium_factor(chemistry, calcium_mg_l);
    let recommended = 1.5 * minimum;

    let (efficiency, adequate) = if available_dose_mg_l > 0.0 {
        let ratio = available_dose_mg_l / minimum;
        let efficiency = if ratio >= 1.0 {
            (0.85 + 0.15 * (1.0 - (-(ratio - 1.0)).exp())).min(1.0)
        } else {
            0.85 * ratio
        };
        (efficiency, ratio >= 1.0)
    } else {
        (0.0, false)
    };

    Ok(ScaleInhibitorPerformanceResult {
        minimum_inhibitor_concentration: kilograms_per_cubic_meter(minimum * 1e-3),
        recommended_dose: kilograms_per_cubic_meter(recommended * 1e-3),
        efficiency,
        adequate,
        thermal_degradation_risk: chemistry == Chemistry::Phosphonate && temperature_c > 175.0,
        calcium_phosphonate_precipitation_risk: chemistry == Chemistry::Phosphonate
            && calcium_mg_l > 2000.0,
        extreme_supersaturation_risk: saturation_ratio > 100.0,
        low_temperature_solubility_risk: chemistry == Chemistry::PhosphateEster
            && temperature_c < 40.0,
        warnings,
    })
}

fn base_mic(scale: Mineral, chem: Chemistry) -> f64 {
    let base = match scale {
        Mineral::Caco3 => 5.0,
        Mineral::Baso4 => 10.0,
        Mineral::Srso4 => 8.0,
        Mineral::Caso4 => 6.0,
        Mineral::Feco3 => 12.0,
    };
    let mut chem_factor = match chem {
        Chemistry::Phosphonate => 1.0,
        Chemistry::Polymaleate => 1.2,
        Chemistry::Polyacrylate => 1.4,
        Chemistry::PhosphateEster => 1.1,
        Chemistry::VinylSulphonate => 1.3,
    };
    if scale == Mineral::Baso4 && chem == Chemistry::Phosphonate {
        chem_factor *= 0.8;
    }
    if scale == Mineral::Caco3 && chem == Chemistry::Polymaleate {
        chem_factor *= 0.9;
    }
    base * chem_factor
}

fn temperature_factor(chem: Chemistry, temperature_c: f64) -> f64 {
    let (knee, slope) = match chem {
        Chemistry::Phosphonate | Chemistry::VinylSulphonate => (100.0, 0.015),
        Chemistry::Polymaleate => (130.0, 0.012),
        Chemistry::Polyacrylate => (120.0, 0.013),
        Chemistry::PhosphateEster => (150.0, 0.010),
    };
    if temperature_c <= knee {
        return 1.0;
    }
    (slope * (temperature_c - knee)).exp()
}

fn ionic_strength_factor(tds_mg_l: f64) -> f64 {
    let tds_ref = 50000.0;
    let factor = 1.0 + 0.20 * (tds_mg_l.max(1.0) / tds_ref).log10();
    factor.max(0.5)
}

fn saturation_factor(saturation_ratio: f64) -> f64 {
    let factor = 1.0 + 0.6 * (saturation_ratio.max(1.0) / 5.0).log10();
    factor.max(0.5)
}

fn calcium_factor(chem: Chemistry, calcium_mg_l: f64) -> f64 {
    if chem == Chemistry::Phosphonate && calcium_mg_l > 1000.0 {
        return 1.0 + 0.5 * (calcium_mg_l / 1000.0).log10();
    }
    1.0
}
