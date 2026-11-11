//! `process.anti_surge_recycle` - the steady recycle that holds a machine on its control line.

use azoth_core::units::kelvins;
use azoth_core::{AzothError, Result};

use crate::compressor_chart::CompressorChart;
use crate::kernels::compressor::{SurgeSettings, compressor_with_chart, feed_volumetric_flow};
use crate::kernels::mixer::mixer;
use crate::stream::Stream;

/// `AntiSurgeRecycleCalculator.surgeControlMargin`'s default.
pub const DEFAULT_MARGIN: f64 = 0.05;
/// `recycleCoolerTemperature`'s default: 40 C.
pub const DEFAULT_COOLER_TEMPERATURE: f64 = 313.15;
/// `AntiSurgeRecycleCalculator.damping`'s default.
pub const DEFAULT_DAMPING: f64 = 0.8;
/// `AntiSurgeRecycleCalculator.tolerance`'s default, as a fraction of the control flow.
pub const DEFAULT_TOLERANCE: f64 = 0.005;
/// `AntiSurgeRecycleCalculator.maxIterations`'s default.
pub const DEFAULT_MAX_ITERATIONS: usize = 20;

/// The calculator's settable state, with the defaults its fields carry.
#[derive(Debug, Clone, Copy)]
pub struct RecycleSettings {
    /// `setSurgeControlMargin`: the control line as a fraction above the surge flow.
    pub control_margin: f64,
    /// `setRecycleCoolerTemperature`: the recycle clone's temperature, K.
    pub cooler_temperature: f64,
    /// `setDamping`: the fixed-point relaxation.
    pub damping: f64,
    /// `setTolerance`: the convergence gate, as a fraction of the control flow.
    pub tolerance: f64,
    /// `setMaxIterations`: the fixed-point cap.
    pub max_iterations: usize,
}

impl Default for RecycleSettings {
    fn default() -> Self {
        Self {
            control_margin: DEFAULT_MARGIN,
            cooler_temperature: DEFAULT_COOLER_TEMPERATURE,
            damping: DEFAULT_DAMPING,
            tolerance: DEFAULT_TOLERANCE,
            max_iterations: DEFAULT_MAX_ITERATIONS,
        }
    }
}

/// `AntiSurgeRecycleCalculator.Result`, in the units the kernel works in.
#[derive(Debug, Clone, PartialEq)]
pub struct RecycleOutcome {
    /// `isRecycleActive`: the natural inlet volume flow was below the control line.
    pub recycle_active: bool,
    /// `isConverged`: a pass met the gate. True on the inactive branch, which performs none.
    pub converged: bool,
    /// `getIterations`: passes through the loop, zero when inactive.
    pub iterations: u32,
    /// `getRecycleMassFlow`, kg/s - the class reports kg/hr.
    pub recycle_mass_flow: f64,
    /// `getInletVolumeFlow`, m3/hr.
    pub inlet_volume_flow: f64,
    /// `getSurgeFlow`, m3/hr.
    pub surge_flow: f64,
    /// `getControlFlow`, m3/hr: `surge_flow * (1 + margin)`.
    pub control_flow: f64,
    /// `getDistanceToSurge` at the final machine state.
    pub distance_to_surge: f64,
}

/// Converge the cooled recycle that holds a charted machine on its control line.
///
/// # Errors
/// `InvalidInput` where the machine reports no surge flow, and whatever the chart, the mixer
/// or the flash refuses.
pub fn anti_surge_recycle(
    feed: &Stream,
    chart: &CompressorChart,
    speed: f64,
    settings: &RecycleSettings,
) -> Result<RecycleOutcome> {
    let natural = compressor_with_chart(feed, chart, speed, &SurgeSettings::default())?;
    let surge = natural.surge_flow_rate.ok_or_else(no_surge_flow)?;
    let control = surge * (1.0 + settings.control_margin);
    let distance = natural.distance_to_surge.ok_or_else(no_surge_flow)?;
    let feed_volume = feed_volumetric_flow(feed)? * 3600.0;

    // The class's own test is `inletVol < controlFlow`, and the inactive branch is its
    // negation: a NaN inlet is inactive rather than silently on the control line.
    let short_of_control = feed_volume < control;
    if !short_of_control {
        return Ok(RecycleOutcome {
            recycle_active: false,
            converged: true,
            iterations: 0,
            recycle_mass_flow: 0.0,
            inlet_volume_flow: feed_volume,
            surge_flow: surge,
            control_flow: control,
            distance_to_surge: distance,
        });
    }

    // `feedMass` in the class: the natural suction's mass flow, kg/hr.
    let molar_mass = feed.molar_mass()?.value;
    let feed_mass = feed.mass_flow()? * 3600.0;
    let mut recycle_mass = ((control - feed_volume) / feed_volume.max(1.0e-6) * feed_mass).max(0.0);

    let mut inlet_volume = feed_volume;
    let mut surge_flow = surge;
    let mut control_flow = control;
    let mut distance_to_surge = distance;
    let mut converged = false;
    let mut iterations = 0u32;

    for index in 0..settings.max_iterations {
        iterations = u32::try_from(index + 1).unwrap_or(u32::MAX);
        let recycle_n = recycle_mass / 3600.0 / molar_mass;
        let recycle = Stream::from_pt(
            feed.components.clone(),
            feed.z.clone(),
            recycle_n,
            feed.p,
            kelvins(settings.cooler_temperature),
        )?;
        let mixed = mixer(&[feed.clone(), recycle], None)?;
        let machine = compressor_with_chart(&mixed, chart, speed, &SurgeSettings::default())?;

        inlet_volume = feed_volumetric_flow(&mixed)? * 3600.0;
        surge_flow = machine.surge_flow_rate.ok_or_else(no_surge_flow)?;
        control_flow = surge_flow * (1.0 + settings.control_margin);
        distance_to_surge = machine.distance_to_surge.ok_or_else(no_surge_flow)?;

        let error = control_flow - inlet_volume;
        if error.abs() < settings.tolerance * control_flow {
            converged = true;
            break;
        }
        let total_mass = mixed.mass_flow()? * 3600.0;
        let volume_per_mass = inlet_volume / total_mass.max(1.0e-6);
        recycle_mass =
            (recycle_mass + settings.damping * error / volume_per_mass.max(1.0e-9)).max(0.0);
    }

    Ok(RecycleOutcome {
        recycle_active: true,
        converged,
        iterations,
        recycle_mass_flow: recycle_mass / 3600.0,
        inlet_volume_flow: inlet_volume,
        surge_flow,
        control_flow,
        distance_to_surge,
    })
}

/// The refusal a machine without a usable surge curve has to make.
fn no_surge_flow() -> AzothError {
    AzothError::invalid_input(
        "surge_curve",
        "the machine reports no surge flow, and a control line needs one",
    )
}
