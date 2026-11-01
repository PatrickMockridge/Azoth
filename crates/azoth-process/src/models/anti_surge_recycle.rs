//! `process.anti_surge_recycle` - the steady recycle as a registered id.
//!
//! Spec: `specs/models/process/anti_surge_recycle.toml`. The arithmetic is
//! [`crate::kernels::anti_surge_recycle`]; what is here is the boundary a case and a
//! cross-impl test address.
//!
//! **The calculator is a loop around `process.compressor`'s chart branch**, not a machine of
//! its own: the surge flow and the distance to surge are that branch's read-outs at each mixed
//! inlet state. The run-time recycle block is off throughout, because this loop owns it.

use azoth_core::units::{
    MassRate, Pressure, ThermodynamicTemperature, VolumeRate, cubic_meter_per_hour,
    cubic_meters_per_second, kilograms_per_second, meter,
};
use azoth_core::{Result, Warning, apply_checks};

use crate::compressor_chart::{CompressorChart, CurvePoint};
use crate::kernels::anti_surge_recycle::{
    DEFAULT_COOLER_TEMPERATURE, DEFAULT_DAMPING, DEFAULT_MARGIN, DEFAULT_MAX_ITERATIONS,
    DEFAULT_TOLERANCE, RecycleSettings, anti_surge_recycle as kernel,
};
use crate::model_gen;
use crate::models::compressor::{CompressorCurvePoint, CompressorSurgePoint};
use crate::stream::Stream;
use crate::surge_curve::SurgeCurve;

/// Result of `process.anti_surge_recycle`.
#[derive(Debug, Clone, PartialEq)]
pub struct AntiSurgeRecycleResult {
    /// Whether a recycle was opened.
    pub recycle_active: bool,
    /// Whether a pass met the gate; true on the inactive branch.
    pub converged: bool,
    /// Passes through the loop, zero when inactive.
    pub iterations: u32,
    /// The converged recycle mass.
    pub recycle_mass_flow: MassRate,
    /// The machine inlet's actual volumetric flow at the final state.
    pub inlet_volume_flow: VolumeRate,
    /// `getSurgeFlowRate()` at the final machine state.
    pub surge_flow: VolumeRate,
    /// The control line, `surge_flow * (1 + margin)`.
    pub control_flow: VolumeRate,
    /// `getDistanceToSurge()` at the final machine state.
    pub distance_to_surge: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

/// Converge the cooled recycle that holds a charted compressor on its control line.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] where a bound the spec declares is violated, the
/// shapes disagree, or the machine reports no surge flow; and whatever the chart, the mixer
/// or the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input
pub fn anti_surge_recycle(
    components: &[String],
    feed_n: f64,
    feed_z: &[f64],
    feed_p: Pressure,
    feed_t: ThermodynamicTemperature,
    curve: &[CompressorCurvePoint],
    surge_curve: &[CompressorSurgePoint],
    speed: Option<f64>,
    surge_control_margin: Option<f64>,
    recycle_cooler_temperature: Option<ThermodynamicTemperature>,
    damping: Option<f64>,
    tolerance: Option<f64>,
    max_iterations: Option<usize>,
) -> Result<AntiSurgeRecycleResult> {
    let spec = &model_gen::ANTI_SURGE_RECYCLE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "feed_n" => Some(feed_n),
            "feed_p" => Some(feed_p.value),
            "feed_t" => Some(feed_t.value),
            "speed" => speed,
            "surge_control_margin" => surge_control_margin,
            "recycle_cooler_temperature" => recycle_cooler_temperature.map(|value| value.value),
            "damping" => damping,
            "tolerance" => tolerance,
            "max_iterations" => max_iterations.map(|value| value as f64),
            _ => None,
        },
        &mut warnings,
    )?;

    let feed = Stream::from_pt(components.to_vec(), feed_z.to_vec(), feed_n, feed_p, feed_t)?;

    let points: Vec<CurvePoint> = curve
        .iter()
        .map(|point| CurvePoint {
            speed: point.speed,
            flow: point.flow.get::<cubic_meter_per_hour>(),
            head: point.head.get::<meter>(),
            efficiency: point.efficiency,
        })
        .collect();
    let (surge_flows, surge_heads): (Vec<f64>, Vec<f64>) = surge_curve
        .iter()
        .map(|point| {
            (
                point.flow.get::<cubic_meter_per_hour>(),
                point.head.get::<meter>(),
            )
        })
        .unzip();
    let mut chart = CompressorChart::new(&points)?;
    chart.set_surge_curve(SurgeCurve::new(&surge_flows, &surge_heads)?);

    let settings = RecycleSettings {
        control_margin: surge_control_margin.unwrap_or(DEFAULT_MARGIN),
        cooler_temperature: recycle_cooler_temperature
            .map_or(DEFAULT_COOLER_TEMPERATURE, |value| value.value),
        damping: damping.unwrap_or(DEFAULT_DAMPING),
        tolerance: tolerance.unwrap_or(DEFAULT_TOLERANCE),
        max_iterations: max_iterations.unwrap_or(DEFAULT_MAX_ITERATIONS),
    };
    let outcome = kernel(&feed, &chart, speed.unwrap_or(3000.0), &settings)?;

    Ok(AntiSurgeRecycleResult {
        recycle_active: outcome.recycle_active,
        converged: outcome.converged,
        iterations: outcome.iterations,
        recycle_mass_flow: kilograms_per_second(outcome.recycle_mass_flow),
        inlet_volume_flow: cubic_meters_per_second(outcome.inlet_volume_flow / 3600.0),
        surge_flow: cubic_meters_per_second(outcome.surge_flow / 3600.0),
        control_flow: cubic_meters_per_second(outcome.control_flow / 3600.0),
        distance_to_surge: outcome.distance_to_surge,
        warnings,
    })
}
