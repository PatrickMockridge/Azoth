//! `MapTurboExpanderCompressor` - an expander and a map-driven compressor on one shaft.
//!
//! **Two modes, one identity.** The expander's power, less the mechanical loss and the bearing
//! loss, is what the compressor may consume. In `SPECIFIED_PRESSURES` the compressor runs
//! independently at its own setting and the residual is reported; in `BALANCED_SPEED` the
//! shaft speed is **bisected** over the compressor map's own measured speed range until the
//! compressor's absorbed power meets the available power - pegged at the lower bound with
//! `UNDER_POWER_SURGE` when even the slowest map line draws more than is available, and at the
//! upper with `OVER_POWER_MAX_SPEED` when even the fastest draws less.
//!
//! **This is the second shaft.** `process.turbo_expander_compressor` is NeqSim's older
//! `TurboExpanderCompressor` (a U/C efficiency map and Q/N correction curves, solved by a
//! damped Newton over a fixed 1000-9000 rpm band); this one drives a real brake compressor on
//! a `CompressorChart` - the Khader map in the oracle - and its speed range is the map's.

use azoth_core::units::Pressure;
use azoth_core::{AzothError, Result};

use crate::compressor_chart::CompressorChart;
use crate::kernels::compressor::{SurgeSettings, compressor, compressor_with_chart};
use crate::kernels::expander::expander;
use crate::stream::Stream;

/// `MapTurboExpanderCompressor.mechanicalEfficiency`'s default.
pub const DEFAULT_MECHANICAL_EFFICIENCY: f64 = 0.98;
/// Its `bearingLossPower`'s default.
pub const DEFAULT_BEARING_LOSS: f64 = 0.0;
/// `Compressor.speed`'s initialiser, rpm - the specified mode reports it when no chart is
/// present.
pub const DEFAULT_COMPRESSOR_SPEED: f64 = 3000.0;
/// The bisection's cap and its interval stop, rpm.
const BISECTION_ITERATIONS: usize = 80;
const SPEED_TOLERANCE: f64 = 1.0e-3;

/// Which solve the shaft runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShaftMode {
    /// The compressor runs at its own setting and the residual is reported.
    SpecifiedPressures,
    /// The speed is bisected until the compressor's power meets the available power.
    BalancedSpeed,
}

/// What the shaft found.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OperatingStatus {
    /// The balance closed, or a specified-pressure run that did not need one.
    Balanced,
    /// Even the slowest map line draws more than the expander produces.
    UnderPowerSurge,
    /// Even the fastest map line draws less.
    OverPowerMaxSpeed,
}

impl OperatingStatus {
    /// The spec's `values` spelling, which the transport and the result carry.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Balanced => "balanced",
            Self::UnderPowerSurge => "under_power_surge",
            Self::OverPowerMaxSpeed => "over_power_max_speed",
        }
    }
}

/// The machine's configuration.
pub struct MapMachineSettings<'a> {
    /// The expander's stated outlet pressure.
    pub expander_out_pressure: Pressure,
    /// The expander's isentropic efficiency.
    pub expander_efficiency: f64,
    /// Which solve runs.
    pub mode: ShaftMode,
    /// The expander's power that reaches the shaft.
    pub mechanical_efficiency: f64,
    /// The bearing loss, W.
    pub bearing_loss: f64,
    /// The shaft-speed lower bound, rpm, where the chart's is not wanted.
    pub min_speed: Option<f64>,
    /// The shaft-speed upper bound, rpm.
    pub max_speed: Option<f64>,
    /// The compressor's map; the balanced mode requires one.
    pub chart: Option<&'a CompressorChart>,
    /// The speed the specified mode reports, rpm.
    pub compressor_speed: f64,
    /// The specified mode's outlet pressure, for a compressor without a map.
    pub compressor_outlet_pressure: Option<Pressure>,
    /// The specified mode's isentropic efficiency, for a compressor without a map.
    pub compressor_isentropic_efficiency: f64,
}

/// What the machine solved.
pub struct MapMachineOutcome {
    /// The expander's outlet.
    pub expander_outlet: Stream,
    /// The compressor's outlet.
    pub compressor_outlet: Stream,
    /// The shaft speed, rpm.
    pub shaft_speed: f64,
    /// The expander's power that reaches the shaft, W.
    pub available_shaft_power: f64,
    /// The compressor's absorbed power, W.
    pub consumed_compressor_power: f64,
    /// `available - consumed`, W.
    pub power_balance_residual: f64,
    /// Whether the speed hit a bound.
    pub speed_limited: bool,
    /// Which mode outcome this is.
    pub status: OperatingStatus,
}

/// Run the map-based single shaft.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a mechanical efficiency outside `(0, 1]`, a negative
/// bearing loss, a balanced-speed run with no map or no outlet pressure in a map-less
/// specified run; and whatever the expander, the compressor or the flash refuses.
pub fn map_turbo_expander(
    expander_feed: &Stream,
    compressor_feed: &Stream,
    settings: &MapMachineSettings<'_>,
) -> Result<MapMachineOutcome> {
    if !(settings.mechanical_efficiency > 0.0 && settings.mechanical_efficiency <= 1.0) {
        return Err(AzothError::invalid_input(
            "mechanical_efficiency",
            format!(
                "the mechanical efficiency must be in (0, 1], and {} is not",
                settings.mechanical_efficiency
            ),
        ));
    }
    if settings.bearing_loss.is_nan() || settings.bearing_loss < 0.0 {
        return Err(AzothError::invalid_input(
            "bearing_loss",
            format!(
                "the bearing loss is {} W, and it cannot be negative",
                settings.bearing_loss
            ),
        ));
    }

    let expander_outlet = expander(
        expander_feed,
        settings.expander_out_pressure,
        settings.expander_efficiency,
    )?;
    // The expander's power is negative - it produces - and the shaft sees its magnitude less
    // the losses.
    let expander_power = expander_feed.n * (expander_outlet.h.value - expander_feed.h.value);
    let available =
        (expander_power.abs() * settings.mechanical_efficiency - settings.bearing_loss).max(0.0);

    let power_at = |speed: f64| -> Result<(f64, Stream)> {
        if let Some(chart) = settings.chart {
            let outcome =
                compressor_with_chart(compressor_feed, chart, speed, &SurgeSettings::default())?;
            return Ok((outcome.power.value, outcome.outlet));
        }
        let pressure = settings.compressor_outlet_pressure.ok_or_else(|| {
            AzothError::invalid_input(
                "compressor_outlet_pressure",
                "a machine without a compressor map reaches only the pressure it is given, and \
                 this one states none",
            )
        })?;
        let outlet = compressor(
            compressor_feed,
            pressure,
            settings.compressor_isentropic_efficiency,
        )?;
        let power = compressor_feed.n * (outlet.h.value - compressor_feed.h.value);
        Ok((power, outlet))
    };

    let (shaft_speed, speed_limited, status, consumed, compressor_outlet) = match settings.mode {
        ShaftMode::SpecifiedPressures => {
            let speed = settings.compressor_speed;
            let (power, outlet) = power_at(speed)?;
            (speed, false, OperatingStatus::Balanced, power, outlet)
        }
        ShaftMode::BalancedSpeed => {
            let chart = settings.chart.ok_or_else(|| {
                AzothError::invalid_input(
                    "curve",
                    "the balanced-speed mode bisects the compressor map, and this machine \
                     states none",
                )
            })?;
            let low = settings
                .min_speed
                .unwrap_or_else(|| chart.min_curve_speed());
            let high = settings
                .max_speed
                .unwrap_or_else(|| chart.max_curve_speed());
            if !low.is_finite() || !high.is_finite() || high <= low {
                return Err(AzothError::invalid_input(
                    "curve",
                    format!(
                        "the map's speed range is {low} to {high} rpm, and a bisection needs one"
                    ),
                ));
            }
            let (low_power, low_outlet) = power_at(low)?;
            if low_power - available > 0.0 {
                (
                    low,
                    true,
                    OperatingStatus::UnderPowerSurge,
                    low_power,
                    low_outlet,
                )
            } else {
                let (high_power, high_outlet) = power_at(high)?;
                if high_power - available < 0.0 {
                    (
                        high,
                        true,
                        OperatingStatus::OverPowerMaxSpeed,
                        high_power,
                        high_outlet,
                    )
                } else {
                    let (mut a, mut b) = (low, high);
                    for _ in 0..BISECTION_ITERATIONS {
                        if (b - a).abs() < SPEED_TOLERANCE {
                            break;
                        }
                        let mid = 0.5 * (a + b);
                        let (power, _) = power_at(mid)?;
                        if power - available > 0.0 {
                            b = mid;
                        } else {
                            a = mid;
                        }
                    }
                    let speed = 0.5 * (a + b);
                    let (power, outlet) = power_at(speed)?;
                    (speed, false, OperatingStatus::Balanced, power, outlet)
                }
            }
        }
    };

    Ok(MapMachineOutcome {
        expander_outlet,
        compressor_outlet,
        shaft_speed,
        available_shaft_power: available,
        consumed_compressor_power: consumed,
        power_balance_residual: available - consumed,
        speed_limited,
        status,
    })
}
