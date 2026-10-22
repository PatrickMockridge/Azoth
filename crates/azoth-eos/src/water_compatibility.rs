//! `eos.water_compatibility` - NeqSim's `WaterCompatibilityScreener`.
//!
//! Two waters mixed linearly at a list of injection-water percentages, each mixture evaluated by
//! the scale predictor with its auto-pH always on, and the worst saturation index reported.

use azoth_core::Result;
use azoth_core::units::{MassDensity, Pressure, ThermodynamicTemperature};

use crate::results::{WaterCompatibilityResult, WaterMixingResult};
use crate::scale_prediction::{ScaleInputs, state};

/// The class's calculator default pH; the two water records' own pH fields are never read.
const DEFAULT_PH: f64 = 6.5;

/// One water's analysis, as the screener takes it.
///
/// The eleven quantities the two water setters carry except the pH, which is stored and never
/// passed to a calculator.
#[derive(Debug, Clone, PartialEq)]
pub struct WaterAnalysis {
    /// Calcium.
    pub calcium: MassDensity,
    /// Barium.
    pub barium: MassDensity,
    /// Strontium.
    pub strontium: MassDensity,
    /// Ferrous iron.
    pub iron: MassDensity,
    /// Magnesium.
    pub magnesium: MassDensity,
    /// Sodium.
    pub sodium: MassDensity,
    /// Bicarbonate.
    pub bicarbonate: MassDensity,
    /// Sulphate.
    pub sulphate: MassDensity,
    /// Total dissolved solids.
    pub total_dissolved_solids: MassDensity,
    /// Absolute temperature.
    pub temperature: ThermodynamicTemperature,
    /// Absolute pressure.
    pub pressure: Pressure,
    /// CO2 partial pressure.
    pub co2_partial_pressure: Pressure,
}

/// The scale risk of mixing two waters at a list of injection-water percentages.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's range.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, kilograms_per_cubic_meter, pascals};
/// use azoth_eos::water_compatibility::{WaterAnalysis, water_compatibility};
///
/// let formation = WaterAnalysis {
///     calcium: kilograms_per_cubic_meter(0.4),
///     barium: kilograms_per_cubic_meter(0.2),
///     strontium: kilograms_per_cubic_meter(0.005),
///     iron: kilograms_per_cubic_meter(0.0),
///     magnesium: kilograms_per_cubic_meter(0.15),
///     sodium: kilograms_per_cubic_meter(10.0),
///     bicarbonate: kilograms_per_cubic_meter(0.15),
///     sulphate: kilograms_per_cubic_meter(0.01),
///     total_dissolved_solids: kilograms_per_cubic_meter(50.0),
///     temperature: kelvins(353.15),
///     pressure: pascals(1.0e7),
///     co2_partial_pressure: pascals(2.0e5),
/// };
/// let mut injection = formation.clone();
/// injection.calcium = kilograms_per_cubic_meter(0.02);
/// injection.barium = kilograms_per_cubic_meter(0.0);
/// injection.strontium = kilograms_per_cubic_meter(0.0);
/// injection.bicarbonate = kilograms_per_cubic_meter(0.1);
/// injection.sulphate = kilograms_per_cubic_meter(2.7);
/// injection.total_dissolved_solids = kilograms_per_cubic_meter(35.0);
/// injection.magnesium = kilograms_per_cubic_meter(0.1);
/// injection.sodium = kilograms_per_cubic_meter(11.0);
/// injection.temperature = kelvins(293.15);
/// injection.co2_partial_pressure = pascals(0.5e5);
///
/// let r = water_compatibility(&[formation], &[injection], &[60.0])?;
/// assert_eq!(r.worst_case_scale.as_deref(), Some("BaSO4"));
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn water_compatibility(
    formation_water: &[WaterAnalysis],
    injection_water: &[WaterAnalysis],
    injection_water_percent: &[f64],
) -> Result<WaterCompatibilityResult> {
    let formation_water = single(formation_water, "formation_water")?;
    let injection_water = single(injection_water, "injection_water")?;
    let mut warnings = Vec::new();
    for (label, analysis) in [
        ("formation_water", formation_water),
        ("injection_water", injection_water),
    ] {
        for (name, value, positive) in [
            ("calcium", analysis.calcium.value, false),
            ("barium", analysis.barium.value, false),
            ("strontium", analysis.strontium.value, false),
            ("iron", analysis.iron.value, false),
            ("magnesium", analysis.magnesium.value, false),
            ("sodium", analysis.sodium.value, false),
            ("bicarbonate", analysis.bicarbonate.value, false),
            ("sulphate", analysis.sulphate.value, false),
            (
                "total_dissolved_solids",
                analysis.total_dissolved_solids.value,
                false,
            ),
            ("temperature", analysis.temperature.value, true),
            ("pressure", analysis.pressure.value, true),
            (
                "co2_partial_pressure",
                analysis.co2_partial_pressure.value,
                false,
            ),
        ] {
            if positive && (value.is_nan() || value <= 0.0) {
                return Err(azoth_core::AzothError::out_of_range(
                    format!("{label}.{name}"),
                    value,
                    "must be positive",
                ));
            }
            if !positive && value < 0.0 {
                return Err(azoth_core::AzothError::out_of_range(
                    format!("{label}.{name}"),
                    value,
                    "must be non-negative",
                ));
            }
        }
    }
    for percent in injection_water_percent {
        if !(0.0..=100.0).contains(percent) {
            warnings.push(azoth_core::Warning::new(
                azoth_core::WarningCode::OutOfValidRange,
                format!(
                    "an injection-water percentage of {percent} is outside [0, 100]; the class \
                     mixes with it as stated rather than clamping"
                ),
            ));
        }
    }

    let mix = |formation: f64, injection: f64, fraction: f64| {
        formation * (1.0 - fraction) + injection * fraction
    };
    let mut rows = Vec::with_capacity(injection_water_percent.len());
    let mut worst: Option<(f64, &'static str, f64)> = None;
    const SCALES: [&str; 5] = ["CaCO3", "BaSO4", "SrSO4", "CaSO4", "FeCO3"];
    for percent in injection_water_percent {
        let fraction = percent / 100.0;
        let temperature_c = mix(
            formation_water.temperature.value - 273.15,
            injection_water.temperature.value - 273.15,
            fraction,
        );
        let inputs = ScaleInputs {
            calcium_mg_l: mix(
                formation_water.calcium.value,
                injection_water.calcium.value,
                fraction,
            ) * 1000.0,
            barium_mg_l: mix(
                formation_water.barium.value,
                injection_water.barium.value,
                fraction,
            ) * 1000.0,
            strontium_mg_l: mix(
                formation_water.strontium.value,
                injection_water.strontium.value,
                fraction,
            ) * 1000.0,
            iron_mg_l: mix(
                formation_water.iron.value,
                injection_water.iron.value,
                fraction,
            ) * 1000.0,
            magnesium_mg_l: mix(
                formation_water.magnesium.value,
                injection_water.magnesium.value,
                fraction,
            ) * 1000.0,
            sodium_mg_l: mix(
                formation_water.sodium.value,
                injection_water.sodium.value,
                fraction,
            ) * 1000.0,
            bicarbonate_mg_l: mix(
                formation_water.bicarbonate.value,
                injection_water.bicarbonate.value,
                fraction,
            ) * 1000.0,
            sulphate_mg_l: mix(
                formation_water.sulphate.value,
                injection_water.sulphate.value,
                fraction,
            ) * 1000.0,
            tds_mg_l: mix(
                formation_water.total_dissolved_solids.value,
                injection_water.total_dissolved_solids.value,
                fraction,
            ) * 1000.0,
            temperature_c,
            pressure_bara: mix(
                formation_water.pressure.value,
                injection_water.pressure.value,
                fraction,
            ) / 1.0e5,
            co2_partial_pressure_bar: mix(
                formation_water.co2_partial_pressure.value,
                injection_water.co2_partial_pressure.value,
                fraction,
            ) / 1.0e5,
            ph: DEFAULT_PH,
            auto_ph: true,
            second_order_pressure_correction: false,
        };
        let point = state(&inputs);
        let si = point.si;
        let defined = |value: f64| (!value.is_nan()).then_some(value);
        for (name, value) in SCALES.iter().zip(si) {
            if value.is_nan() {
                continue;
            }
            if worst.is_none_or(|(_, _, current)| value > current) {
                worst = Some((*percent, name, value));
            }
        }
        rows.push(WaterMixingResult {
            injection_water_percent: *percent,
            si_caco3: defined(si[0]),
            si_baso4: defined(si[1]),
            si_srso4: defined(si[2]),
            si_caso4: defined(si[3]),
            si_feco3: defined(si[4]),
        });
    }

    Ok(WaterCompatibilityResult {
        mixtures: rows,
        worst_case_percent: worst.map(|(percent, _, _)| percent),
        worst_case_scale: worst.map(|(_, name, _)| name.to_string()),
        worst_case_saturation_index: worst.map(|(_, _, value)| value),
        warnings,
    })
}

/// The one water a role names, or a refusal when the list is not exactly one row.
fn single<'a>(rows: &'a [WaterAnalysis], role: &str) -> Result<&'a WaterAnalysis> {
    match rows {
        [row] => Ok(row),
        _ => Err(azoth_core::AzothError::invalid_input(
            role,
            format!("`{role}` is one water and this is {} row(s)", rows.len()),
        )),
    }
}
