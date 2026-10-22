//! `eos.brine_mixing_scale_evaluator` - NeqSim's `BrineMixingScaleEvaluator`.
//!
//! Two incompatible brines blended linearly at a set of mixing fractions, each blend evaluated
//! by the activity-corrected electrolyte calculator, and the worst saturation index reported.

use azoth_core::units::{MassDensity, Pressure, ThermodynamicTemperature};
use azoth_core::{Result, apply_checks};

use crate::electrolyte_scale_calculator::{ElectrolyteInputs, state};
use crate::model_gen;
use crate::results::{BrineMixingPointRow, BrineMixingScaleEvaluatorResult};

/// One brine's analysis, in the nine ions the class's composition carries.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrineAnalysis {
    /// Calcium.
    pub calcium: MassDensity,
    /// Barium.
    pub barium: MassDensity,
    /// Strontium.
    pub strontium: MassDensity,
    /// Magnesium.
    pub magnesium: MassDensity,
    /// Sodium.
    pub sodium: MassDensity,
    /// Potassium.
    pub potassium: MassDensity,
    /// Chloride.
    pub chloride: MassDensity,
    /// Sulphate.
    pub sulphate: MassDensity,
    /// Bicarbonate.
    pub bicarbonate: MassDensity,
}

/// The mineral saturation of two brines across their mixing range.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's range.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, kilograms_per_cubic_meter, pascals};
/// use azoth_eos::brine_mixing_scale_evaluator::{
///     BrineAnalysis, brine_mixing_scale_evaluator,
/// };
///
/// let formation = BrineAnalysis {
///     calcium: kilograms_per_cubic_meter(5.0),
///     barium: kilograms_per_cubic_meter(0.3),
///     strontium: kilograms_per_cubic_meter(0.8),
///     magnesium: kilograms_per_cubic_meter(0.5),
///     sodium: kilograms_per_cubic_meter(30.0),
///     potassium: kilograms_per_cubic_meter(0.4),
///     chloride: kilograms_per_cubic_meter(55.0),
///     sulphate: kilograms_per_cubic_meter(0.0),
///     bicarbonate: kilograms_per_cubic_meter(0.4),
/// };
/// let injection = BrineAnalysis {
///     calcium: kilograms_per_cubic_meter(0.4),
///     barium: kilograms_per_cubic_meter(0.0),
///     strontium: kilograms_per_cubic_meter(0.008),
///     magnesium: kilograms_per_cubic_meter(1.3),
///     sodium: kilograms_per_cubic_meter(11.0),
///     potassium: kilograms_per_cubic_meter(0.4),
///     chloride: kilograms_per_cubic_meter(19.0),
///     sulphate: kilograms_per_cubic_meter(2.7),
///     bicarbonate: kilograms_per_cubic_meter(0.14),
/// };
///
/// let r = brine_mixing_scale_evaluator(
///     &[formation], &[injection], kelvins(333.15), pascals(1.0e7), 6.0, 11,
/// )?;
/// assert_eq!(r.worst_mineral, "BaSO4");
/// assert_eq!(r.points.len(), 11);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn brine_mixing_scale_evaluator(
    brine_a: &[BrineAnalysis],
    brine_b: &[BrineAnalysis],
    temperature: ThermodynamicTemperature,
    pressure: Pressure,
    ph: f64,
    steps: usize,
) -> Result<BrineMixingScaleEvaluatorResult> {
    let spec = &model_gen::BRINE_MIXING_SCALE_EVALUATOR_SPEC;
    let mut warnings = Vec::new();

    let brine_a = single(brine_a, "brine_a")?;
    let brine_b = single(brine_b, "brine_b")?;

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "temperature" => Some(temperature.value),
            "pressure" => Some(pressure.value),
            "ph" => Some(ph),
            _ => None,
        },
        &mut warnings,
    )?;
    for (label, brine) in [("brine_a", brine_a), ("brine_b", brine_b)] {
        for (name, value) in [
            ("calcium", brine.calcium.value),
            ("barium", brine.barium.value),
            ("strontium", brine.strontium.value),
            ("magnesium", brine.magnesium.value),
            ("sodium", brine.sodium.value),
            ("potassium", brine.potassium.value),
            ("chloride", brine.chloride.value),
            ("sulphate", brine.sulphate.value),
            ("bicarbonate", brine.bicarbonate.value),
        ] {
            if value < 0.0 {
                return Err(azoth_core::AzothError::out_of_range(
                    format!("{label}.{name}"),
                    value,
                    "a concentration cannot be negative",
                ));
            }
        }
    }

    let steps = steps.max(2);
    let mut points = Vec::with_capacity(steps);
    let mut worst_si = f64::NEG_INFINITY;
    let mut worst_mineral = "NONE".to_owned();
    let mut worst_fraction_a = f64::NAN;

    for index in 0..steps {
        let fraction_a = index as f64 / (steps - 1) as f64;
        let fraction_b = 1.0 - fraction_a;
        let blend = |a: MassDensity, b: MassDensity| a.value * fraction_a + b.value * fraction_b;

        let inputs = ElectrolyteInputs {
            calcium_mg_l: blend(brine_a.calcium, brine_b.calcium) * 1000.0,
            barium_mg_l: blend(brine_a.barium, brine_b.barium) * 1000.0,
            strontium_mg_l: blend(brine_a.strontium, brine_b.strontium) * 1000.0,
            magnesium_mg_l: blend(brine_a.magnesium, brine_b.magnesium) * 1000.0,
            sodium_mg_l: blend(brine_a.sodium, brine_b.sodium) * 1000.0,
            potassium_mg_l: blend(brine_a.potassium, brine_b.potassium) * 1000.0,
            iron_mg_l: 0.0,
            chloride_mg_l: blend(brine_a.chloride, brine_b.chloride) * 1000.0,
            sulphate_mg_l: blend(brine_a.sulphate, brine_b.sulphate) * 1000.0,
            bicarbonate_mg_l: blend(brine_a.bicarbonate, brine_b.bicarbonate) * 1000.0,
            carbonate_mg_l: 0.0,
            temperature_c: temperature.value - 273.15,
            pressure_bara: pressure.value / 1.0e5,
            ph,
        };
        let state = state(&inputs);

        let row = BrineMixingPointRow {
            fraction_a,
            si_caco3: state.si_caco3,
            si_baso4: state.si_baso4,
            si_caso4: state.si_caso4,
            si_srso4: state.si_srso4,
        };
        for (mineral, si) in [
            ("CaCO3", row.si_caco3),
            ("BaSO4", row.si_baso4),
            ("CaSO4", row.si_caso4),
            ("SrSO4", row.si_srso4),
        ] {
            if si > worst_si {
                worst_si = si;
                worst_mineral = mineral.to_owned();
                worst_fraction_a = fraction_a;
            }
        }
        points.push(row);
    }

    Ok(BrineMixingScaleEvaluatorResult {
        points,
        worst_fraction_a,
        worst_mineral,
        worst_saturation_index: worst_si,
        warnings,
    })
}

fn single<'a>(rows: &'a [BrineAnalysis], label: &str) -> Result<&'a BrineAnalysis> {
    match rows {
        [row] => Ok(row),
        _ => Err(azoth_core::AzothError::invalid_input(
            label,
            format!("`{label}` is one brine and this is {} row(s)", rows.len()),
        )),
    }
}
