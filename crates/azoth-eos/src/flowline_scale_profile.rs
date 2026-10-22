//! `eos.flowline_scale_profile` - NeqSim's `FlowlineScaleProfile`.
//!
//! A water analysis walked along a linear temperature and pressure profile: `N + 1` points from
//! the inlet fraction zero to the outlet fraction one, each evaluated by the scale predictor.

use azoth_core::units::{MassDensity, Pressure, ThermodynamicTemperature};
use azoth_core::{Result, apply_checks};

use crate::model_gen;
use crate::results::{FlowlineScaleProfileResult, FlowlineScaleSegment};
use crate::scale_prediction::{ScaleInputs, state};

/// The class's calculator default pH, which it never overwrites from `setWaterChemistry`.
const DEFAULT_PH: f64 = 6.5;

/// Mineral-scale saturation indices along a flowline's temperature and pressure profile.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's range.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, kilograms_per_cubic_meter, pascals};
/// use azoth_eos::flowline_scale_profile::flowline_scale_profile;
///
/// let r = flowline_scale_profile(
///     kilograms_per_cubic_meter(0.4), kilograms_per_cubic_meter(0.01),
///     kilograms_per_cubic_meter(0.005), kilograms_per_cubic_meter(0.0),
///     kilograms_per_cubic_meter(0.15), kilograms_per_cubic_meter(0.01),
///     kilograms_per_cubic_meter(0.15), kilograms_per_cubic_meter(0.01),
///     kilograms_per_cubic_meter(40.0), pascals(2.0e5), false,
///     kelvins(353.15), pascals(2.0e7), kelvins(293.15), pascals(5.0e6), 10,
/// )?;
/// assert_eq!(r.segments.len(), 11);
/// assert!((r.max_baso4.expect("barite has both ions") - 1.048081964483009).abs() < 1e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn flowline_scale_profile(
    calcium: MassDensity,
    barium: MassDensity,
    strontium: MassDensity,
    iron: MassDensity,
    magnesium: MassDensity,
    sodium: MassDensity,
    bicarbonate: MassDensity,
    sulphate: MassDensity,
    total_dissolved_solids: MassDensity,
    co2_partial_pressure: Pressure,
    auto_ph: bool,
    inlet_temperature: ThermodynamicTemperature,
    inlet_pressure: Pressure,
    outlet_temperature: ThermodynamicTemperature,
    outlet_pressure: Pressure,
    segments: usize,
) -> Result<FlowlineScaleProfileResult> {
    let spec = &model_gen::FLOWLINE_SCALE_PROFILE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "calcium" => Some(calcium.value),
            "barium" => Some(barium.value),
            "strontium" => Some(strontium.value),
            "iron" => Some(iron.value),
            "magnesium" => Some(magnesium.value),
            "sodium" => Some(sodium.value),
            "bicarbonate" => Some(bicarbonate.value),
            "sulphate" => Some(sulphate.value),
            "total_dissolved_solids" => Some(total_dissolved_solids.value),
            "co2_partial_pressure" => Some(co2_partial_pressure.value),
            "inlet_temperature" => Some(inlet_temperature.value),
            "inlet_pressure" => Some(inlet_pressure.value),
            "outlet_temperature" => Some(outlet_temperature.value),
            "outlet_pressure" => Some(outlet_pressure.value),
            "segments" => Some(segments as f64),
            _ => None,
        },
        &mut warnings,
    )?;

    let steps = segments.max(2);
    let mut rows = Vec::with_capacity(steps + 1);
    let mut maxima = [None::<f64>; 5];
    for index in 0..=steps {
        let fraction = index as f64 / steps as f64;
        let temperature_c = (inlet_temperature.value - 273.15)
            + fraction * ((outlet_temperature.value - 273.15) - (inlet_temperature.value - 273.15));
        let pressure_bara = inlet_pressure.value / 1.0e5
            + fraction * (outlet_pressure.value / 1.0e5 - inlet_pressure.value / 1.0e5);
        let inputs = ScaleInputs {
            calcium_mg_l: calcium.value * 1000.0,
            barium_mg_l: barium.value * 1000.0,
            strontium_mg_l: strontium.value * 1000.0,
            iron_mg_l: iron.value * 1000.0,
            magnesium_mg_l: magnesium.value * 1000.0,
            sodium_mg_l: sodium.value * 1000.0,
            bicarbonate_mg_l: bicarbonate.value * 1000.0,
            sulphate_mg_l: sulphate.value * 1000.0,
            tds_mg_l: total_dissolved_solids.value * 1000.0,
            temperature_c,
            pressure_bara,
            co2_partial_pressure_bar: co2_partial_pressure.value / 1.0e5,
            ph: DEFAULT_PH,
            auto_ph,
            second_order_pressure_correction: false,
        };
        let point = state(&inputs);
        let si = point.si;
        let defined = |value: f64| (!value.is_nan()).then_some(value);
        for (slot, value) in maxima.iter_mut().zip(si) {
            if value.is_nan() {
                continue;
            }
            if slot.is_none_or(|current| value > current) {
                *slot = Some(value);
            }
        }
        rows.push(FlowlineScaleSegment {
            distance_fraction: fraction,
            temperature: azoth_core::units::kelvins(temperature_c + 273.15),
            pressure: azoth_core::units::pascals(pressure_bara * 1.0e5),
            si_caco3: defined(si[0]),
            si_baso4: defined(si[1]),
            si_srso4: defined(si[2]),
            si_caso4: defined(si[3]),
            si_feco3: defined(si[4]),
        });
    }

    Ok(FlowlineScaleProfileResult {
        segments: rows,
        max_caco3: maxima[0],
        max_baso4: maxima[1],
        max_srso4: maxima[2],
        max_caso4: maxima[3],
        max_feco3: maxima[4],
        warnings,
    })
}
