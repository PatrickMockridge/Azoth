//! `eos.scale_control_assessor` - NeqSim's `ScaleControlAssessor`.
//!
//! The scale predictor's five indices with a threshold inhibitor's efficiency folded into a
//! kinetic risk index per scale: `SI + log10(max(1e-6, 1 - efficiency))`, the efficiency
//! capped at 0.999 and the worst index reported against the class's 0.5 threshold.

use azoth_core::units::{MassDensity, Pressure, ThermodynamicTemperature};
use azoth_core::{Result, apply_checks};

use crate::model_gen;
use crate::results::{ScaleControlAssessorResult, ScaleControlRiskRow, ScaleControlRow};
use crate::scale_inhibitor_performance::scale_inhibitor_performance;
use crate::scale_prediction::{ScaleInputs, state};

/// The scales, in the class's map order.
const SCALES: [&str; 5] = ["CACO3", "BASO4", "SRSO4", "CASO4", "FECO3"];

/// One configured inhibitor, as the class's `addInhibitor` takes it.
#[derive(Debug, Clone, PartialEq)]
pub struct ScaleInhibitor {
    /// The scale the row is registered for.
    pub scale_type: String,
    /// The inhibitor family.
    pub inhibitor_chemistry: String,
    /// The temperature the inhibitor model was configured at.
    pub temperature: ThermodynamicTemperature,
    /// `SR = IAP/Ksp` the model was configured with.
    pub saturation_ratio: f64,
    /// The TDS the model was configured with.
    pub total_dissolved_solids: MassDensity,
    /// The calcium the model was configured with.
    pub calcium: MassDensity,
    /// The active dose available.
    pub available_dose: MassDensity,
}

/// The residual scaling risk after threshold inhibition.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's range.
/// * [`azoth_core::AzothError::InvalidInput`] for an unknown scale type or chemistry in an
///   inhibitor row.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, kilograms_per_cubic_meter, pascals};
/// use azoth_eos::scale_control_assessor::{ScaleInhibitor, scale_control_assessor};
///
/// let inhibitor = ScaleInhibitor {
///     scale_type: "caco3".to_owned(),
///     inhibitor_chemistry: "phosphonate".to_owned(),
///     temperature: kelvins(353.15),
///     saturation_ratio: 1.6851624454329173,
///     total_dissolved_solids: kilograms_per_cubic_meter(80.0),
///     calcium: kilograms_per_cubic_meter(2.0),
///     available_dose: kilograms_per_cubic_meter(0.02),
/// };
/// let r = scale_control_assessor(
///     kilograms_per_cubic_meter(2.0),
///     kilograms_per_cubic_meter(0.05),
///     kilograms_per_cubic_meter(0.0),
///     kilograms_per_cubic_meter(0.0),
///     kilograms_per_cubic_meter(0.0),
///     kilograms_per_cubic_meter(20.0),
///     kilograms_per_cubic_meter(0.8),
///     kilograms_per_cubic_meter(0.5),
///     kilograms_per_cubic_meter(80.0),
///     kelvins(353.15),
///     pascals(1.0e5),
///     pascals(0.0),
///     7.0,
///     &[inhibitor],
/// )?;
/// assert!(!r.kinetically_controlled);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn scale_control_assessor(
    calcium: MassDensity,
    barium: MassDensity,
    strontium: MassDensity,
    iron: MassDensity,
    magnesium: MassDensity,
    sodium: MassDensity,
    bicarbonate: MassDensity,
    sulphate: MassDensity,
    total_dissolved_solids: MassDensity,
    temperature: ThermodynamicTemperature,
    pressure: Pressure,
    co2_partial_pressure: Pressure,
    ph: f64,
    inhibitors: &[ScaleInhibitor],
) -> Result<ScaleControlAssessorResult> {
    let spec = &model_gen::SCALE_CONTROL_ASSESSOR_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "temperature" => Some(temperature.value),
            "pressure" => Some(pressure.value),
            "co2_partial_pressure" => Some(co2_partial_pressure.value),
            "calcium" => Some(calcium.value),
            "barium" => Some(barium.value),
            "strontium" => Some(strontium.value),
            "iron" => Some(iron.value),
            "magnesium" => Some(magnesium.value),
            "sodium" => Some(sodium.value),
            "bicarbonate" => Some(bicarbonate.value),
            "sulphate" => Some(sulphate.value),
            "total_dissolved_solids" => Some(total_dissolved_solids.value),
            "ph" => Some(ph),
            _ => None,
        },
        &mut warnings,
    )?;

    let state = state(&ScaleInputs {
        calcium_mg_l: calcium.value * 1000.0,
        barium_mg_l: barium.value * 1000.0,
        strontium_mg_l: strontium.value * 1000.0,
        iron_mg_l: iron.value * 1000.0,
        magnesium_mg_l: magnesium.value * 1000.0,
        sodium_mg_l: sodium.value * 1000.0,
        bicarbonate_mg_l: bicarbonate.value * 1000.0,
        sulphate_mg_l: sulphate.value * 1000.0,
        tds_mg_l: total_dissolved_solids.value * 1000.0,
        temperature_c: temperature.value - 273.15,
        pressure_bara: pressure.value / 1.0e5,
        co2_partial_pressure_bar: co2_partial_pressure.value / 1.0e5,
        ph,
        auto_ph: false,
        second_order_pressure_correction: false,
    });

    let mut thermodynamic_si = Vec::with_capacity(SCALES.len());
    let mut kinetic_risk_index = Vec::with_capacity(SCALES.len());
    let mut worst: Option<f64> = None;

    for (index, scale) in SCALES.iter().enumerate() {
        let si = state.si[index];
        let thermodynamic = (!si.is_nan()).then_some(si);
        let risk = if si.is_nan() {
            None
        } else if let Some(inhibitor) = inhibitors
            .iter()
            .rev()
            .find(|row| row.scale_type.eq_ignore_ascii_case(scale))
        {
            let performance = scale_inhibitor_performance(
                &inhibitor.scale_type,
                &inhibitor.inhibitor_chemistry,
                inhibitor.temperature,
                inhibitor.saturation_ratio,
                inhibitor.total_dissolved_solids,
                inhibitor.calcium,
                inhibitor.available_dose,
            )?;
            let efficiency = performance.efficiency.clamp(0.0, 0.999);
            Some(si + (1.0 - efficiency).max(1.0e-6).log10())
        } else {
            Some(si)
        };

        if let Some(value) = risk {
            worst = Some(worst.map_or(value, |current: f64| current.max(value)));
        }
        thermodynamic_si.push(ScaleControlRow {
            scale: (*scale).to_owned(),
            saturation_index: thermodynamic,
        });
        kinetic_risk_index.push(ScaleControlRiskRow {
            scale: (*scale).to_owned(),
            risk_index: risk,
        });
    }

    Ok(ScaleControlAssessorResult {
        thermodynamic_si,
        kinetic_risk_index,
        worst_kinetic_risk_index: worst,
        kinetically_controlled: worst.is_some_and(|value| value < 0.5),
        warnings,
    })
}
