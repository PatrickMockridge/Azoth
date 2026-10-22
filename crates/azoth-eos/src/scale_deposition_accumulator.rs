//! `eos.scale_deposition_accumulator` - NeqSim's `ScaleDepositionAccumulator`.
//!
//! A segmented pipeline's scale deposition from a temperature and pressure profile: the
//! CaCO3 index per segment, its deposit mass and wall thickness, the total, and the time to
//! half-block the worst segment.

use azoth_core::Result;
use azoth_core::apply_checks;
use azoth_core::units::{Length, MassDensity, MassFlux, Time, kilograms, meters, seconds};

use crate::model_gen;
use crate::results::{DepositionSegmentRow, ScaleDepositionAccumulatorResult};
use crate::scale_prediction::{ScaleInputs, state};

/// The segmented deposition of a pipeline profile.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's hard range.
///
/// # Example
/// ```
/// use azoth_core::units::{
///     kilograms_per_cubic_meter, kilograms_per_square_meter_second, meters, seconds,
/// };
/// use azoth_eos::scale_deposition_accumulator::scale_deposition_accumulator;
///
/// let r = scale_deposition_accumulator(
///     kilograms_per_cubic_meter(1.0),
///     kilograms_per_cubic_meter(0.5),
///     kilograms_per_cubic_meter(0.1),
///     kilograms_per_cubic_meter(0.01),
///     kilograms_per_cubic_meter(20.0),
///     kilograms_per_cubic_meter(50.0),
///     6.5,
///     meters(0.1),
///     meters(1000.0),
///     &[333.15],
///     &[1.0e7],
///     seconds(31_557_600.0),
///     kilograms_per_cubic_meter(2700.0),
///     kilograms_per_square_meter_second(1.5844043908606944e-8),
///     0.0,
/// )?;
/// assert!((r.total_deposition_mass.value - 5.112073003689538).abs() < 1e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn scale_deposition_accumulator(
    calcium: MassDensity,
    bicarbonate: MassDensity,
    sulphate: MassDensity,
    barium: MassDensity,
    sodium: MassDensity,
    total_dissolved_solids: MassDensity,
    ph: f64,
    diameter: Length,
    length: Length,
    segment_temperatures: &[f64],
    segment_pressures: &[f64],
    service_time: Time,
    scale_density: MassDensity,
    base_rate: MassFlux,
    inhibitor_efficiency: f64,
) -> Result<ScaleDepositionAccumulatorResult> {
    let spec = &model_gen::SCALE_DEPOSITION_ACCUMULATOR_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "ph" => Some(ph),
            "diameter" => Some(diameter.value),
            "length" => Some(length.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let service = service_time.value.max(0.0);
    let efficiency = inhibitor_efficiency.clamp(0.0, 1.0);
    let count = segment_temperatures.len().min(segment_pressures.len());

    if diameter.value <= 0.0 || length.value <= 0.0 || count == 0 {
        return Ok(ScaleDepositionAccumulatorResult {
            segments: Vec::new(),
            total_deposition_mass: kilograms(0.0),
            max_thickness: meters(0.0),
            time_to_blockage: seconds(f64::INFINITY),
            warnings,
        });
    }

    let segment_length = length.value / count as f64;
    let segment_area = std::f64::consts::PI * diameter.value * segment_length;

    let mut segments = Vec::with_capacity(count);
    let mut total_mass = 0.0;
    let mut max_thickness: f64 = 0.0;

    for index in 0..count {
        let segment = state(&ScaleInputs {
            calcium_mg_l: calcium.value * 1000.0,
            barium_mg_l: barium.value * 1000.0,
            strontium_mg_l: 0.0,
            iron_mg_l: 0.0,
            magnesium_mg_l: 0.0,
            sodium_mg_l: sodium.value * 1000.0,
            bicarbonate_mg_l: bicarbonate.value * 1000.0,
            sulphate_mg_l: sulphate.value * 1000.0,
            tds_mg_l: total_dissolved_solids.value * 1000.0,
            temperature_c: segment_temperatures[index] - 273.15,
            pressure_bara: segment_pressures[index] / 1.0e5,
            co2_partial_pressure_bar: 0.0,
            ph,
            auto_ph: false,
            second_order_pressure_correction: false,
        });
        let saturation_index = segment.si[0];
        let driving = (10.0f64.powf(saturation_index) - 1.0).max(0.0);
        let rate = base_rate.value * driving * (1.0 - efficiency);
        let mass = rate * segment_area * service;
        let thickness = mass / (scale_density.value * segment_area);

        total_mass += mass;
        max_thickness = max_thickness.max(thickness);
        segments.push(DepositionSegmentRow {
            saturation_index,
            deposition_mass: kilograms(mass),
            thickness: meters(thickness),
        });
    }

    let time_to_blockage = if max_thickness > 0.0 && service > 0.0 {
        // Half the inner radius, at the worst segment's linear growth.
        (diameter.value / 4.0) * service / max_thickness
    } else {
        f64::INFINITY
    };

    Ok(ScaleDepositionAccumulatorResult {
        segments,
        total_deposition_mass: kilograms(total_mass),
        max_thickness: meters(max_thickness),
        time_to_blockage: seconds(time_to_blockage),
        warnings,
    })
}
