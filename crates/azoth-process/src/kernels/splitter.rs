//! `unit_ops.splitter` - one stream, split into several with the same state.

use azoth_core::units::{Length, Pressure, Velocity, meters, meters_per_second, pascals};
use azoth_core::{AzothError, Result};

use crate::segment::phase::system_volumetric_flow_from_feed;
use crate::stream::Stream;

/// `CapacityConstraint.MAX_UTILIZATION`: the ceiling on any utilization.
pub const MAX_UTILIZATION: f64 = 9.99;
/// `Splitter.designPressureDrop`, the class's `0.05` bar, in pascals.
pub const DEFAULT_DESIGN_PRESSURE_DROP_PA: f64 = 5000.0;
/// `Splitter.maxDesignVelocity`, m/s.
pub const DEFAULT_MAX_DESIGN_VELOCITY_M_S: f64 = 30.0;

/// Split a feed into several outlets with the same composition, pressure,
/// temperature and molar enthalpy, scaled by the given fractions.
///
/// The fractions are normalised, so `[0.3, 0.7]` and `[3.0, 7.0]` mean the same
/// split. Nothing else changes - a splitter does not flash, mix or lose pressure.
pub fn splitter(feed: &Stream, fractions: &[f64]) -> Result<Vec<Stream>> {
    if fractions.is_empty() {
        return Err(AzothError::invalid_input(
            "fractions",
            "a splitter needs at least one outlet",
        ));
    }
    let total: f64 = fractions.iter().sum();
    if total <= 0.0 {
        return Err(AzothError::invalid_input(
            "fractions",
            "the split fractions must sum to a positive value",
        ));
    }

    Ok(fractions
        .iter()
        .map(|&f| Stream {
            components: feed.components.clone(),
            z: feed.z.clone(),
            n: feed.n * f / total,
            p: feed.p,
            t: feed.t,
            h: feed.h,
            // A split moves no composition, so it moves no phase: each product is on the side its
            // feed was on, and a product of a feed nothing flashed is one nothing knows about.
            vapour_fraction: feed.vapour_fraction,
        })
        .collect())
}

/// `Splitter`'s capacity map and its three summaries, in SI.
///
/// The two constraints are the class's own: a `pressureDrop` whose value is always zero and a
/// `velocity` read from the inlet's system volume flow over the header area. The six per-constraint
/// fields are `None` when that constraint is not in the map, which is the class's lazy
/// initialization, not an absent answer.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SplitterCapacity {
    /// `pressureDrop`'s current value, Pa: always zero.
    pub pressure_drop_value: Option<Pressure>,
    /// `pressureDrop`'s design value, Pa.
    pub pressure_drop_design: Option<Pressure>,
    /// `pressureDrop`'s utilization.
    pub pressure_drop_utilization: Option<f64>,
    /// `velocity`'s current value, m/s.
    pub velocity_value: Option<Velocity>,
    /// `velocity`'s design value, m/s.
    pub velocity_design: Option<Velocity>,
    /// `velocity`'s utilization.
    pub velocity_utilization: Option<f64>,
    /// `getMaxUtilization`, 0.0 on an empty map.
    pub max_utilization: f64,
    /// `isCapacityExceeded`.
    pub capacity_exceeded: bool,
    /// `getBottleneckConstraint`'s class-side name.
    pub bottleneck_constraint: Option<&'static str>,
}

/// The capacity map `Splitter.splitterCapacityConstraints` holds, evaluated on `feed`.
///
/// **The map is lazy.** The class's constructor creates the mechanical design but no constraints;
/// `setDesignPressureDrop` and `setMaxDesignVelocity` each initialize it, and
/// `setCapacityAnalysisEnabled(true)` initializes it only when empty. An input that is `Some` is a
/// setter that was called, so it initializes the map; a non-positive limit suppresses its own
/// constraint but does not suppress the initialization.
///
/// # Errors
/// Whatever the feed's system volume flow raises - the capacity read flashes the feed for the
/// density `getFlowRate("m3/hr")` uses. The flash is skipped when the velocity constraint is
/// absent or the diameter is not positive, exactly the class's supplier guard.
pub fn capacity(
    feed: &Stream,
    design_pressure_drop: Option<Pressure>,
    max_design_velocity: Option<Velocity>,
    header_diameter: Option<Length>,
    analysis_enabled: bool,
) -> Result<SplitterCapacity> {
    if design_pressure_drop.is_none() && max_design_velocity.is_none() && !analysis_enabled {
        return Ok(SplitterCapacity::default());
    }

    let design_pressure_drop =
        design_pressure_drop.unwrap_or_else(|| pascals(DEFAULT_DESIGN_PRESSURE_DROP_PA));
    let max_design_velocity =
        max_design_velocity.unwrap_or_else(|| meters_per_second(DEFAULT_MAX_DESIGN_VELOCITY_M_S));

    let pressure_drop = (design_pressure_drop.value > 0.0).then(|| {
        let value = pascals(0.0);
        (
            value,
            design_pressure_drop,
            utilization(value.value, design_pressure_drop.value),
        )
    });

    let velocity = if max_design_velocity.value > 0.0 {
        let diameter = header_diameter.unwrap_or_else(|| meters(0.0));
        let value = if diameter.value > 0.0 {
            let volume_flow = system_volumetric_flow_from_feed(feed)?;
            let area = std::f64::consts::PI * diameter.value * diameter.value / 4.0;
            meters_per_second(volume_flow / area)
        } else {
            meters_per_second(0.0)
        };
        Some((
            value,
            max_design_velocity,
            utilization(value.value, max_design_velocity.value),
        ))
    } else {
        None
    };

    // `getBottleneckConstraint` walks the insertion-ordered map - pressureDrop first, then
    // velocity - and keeps a constraint only when its utilization is strictly above the running
    // maximum, which starts at zero. The always-zero pressure drop therefore never wins, and a
    // zero-utilization map has no bottleneck at all.
    let mut max_utilization = 0.0;
    let mut bottleneck_constraint = None;
    if let Some((_, _, value)) = pressure_drop {
        if value > max_utilization {
            max_utilization = value;
            bottleneck_constraint = Some("pressureDrop");
        }
    }
    if let Some((_, _, value)) = velocity {
        if value > max_utilization {
            max_utilization = value;
            bottleneck_constraint = Some("velocity");
        }
    }

    Ok(SplitterCapacity {
        pressure_drop_value: pressure_drop.map(|(value, _, _)| value),
        pressure_drop_design: pressure_drop.map(|(_, design, _)| design),
        pressure_drop_utilization: pressure_drop.map(|(_, _, value)| value),
        velocity_value: velocity.map(|(value, _, _)| value),
        velocity_design: velocity.map(|(_, design, _)| design),
        velocity_utilization: velocity.map(|(_, _, value)| value),
        max_utilization,
        capacity_exceeded: pressure_drop.is_some_and(|(_, _, value)| value > 1.0)
            || velocity.is_some_and(|(_, _, value)| value > 1.0),
        bottleneck_constraint,
    })
}

/// `CapacityConstraint.getUtilization` for one of the splitter's maximum constraints.
///
/// A design value that is not positive answers zero, as the class's does; the cap is the class's
/// own `MAX_UTILIZATION`.
fn utilization(value: f64, design: f64) -> f64 {
    if design <= 0.0 {
        0.0
    } else {
        (value / design).min(MAX_UTILIZATION)
    }
}
