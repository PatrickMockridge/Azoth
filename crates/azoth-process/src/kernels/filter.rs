//! `unit_ops.filter` - a fixed pressure drop, flashed at the feed's temperature.

use azoth_core::Result;
use azoth_core::units::{Pressure, pascals};

use crate::stream::Stream;

/// What a filter did: the outlet record, and the drop it actually applied.
pub struct FilterOutcome {
    /// The outlet record.
    pub outlet: Stream,
    /// The drop applied, Pa - which is the requested one unless it was clamped.
    pub applied_drop: Pressure,
}

/// Drop a stream's pressure by a fixed amount at a constant temperature.
///
/// **It holds the temperature, which is what separates this from `throttling_valve`.**
/// `Filter.run` sets the reduced pressure and runs a `TPflash`, so the outlet is at the
/// *feed's* temperature and its enthalpy moves with the pressure; a valve is a `PHflash` and
/// holds the enthalpy instead. On the captured fluid the two differ by about `25` J/mol over
/// one bar, so a port that took the ordinary isenthalpic reading of "a pressure drop across a
/// restriction" would be wrong by a visible amount.
///
/// **A drop larger than the inlet pressure is clamped, not refused.** `run` applies
/// `min(max(0, dP), max(0, P_in - 1e-6 bar))` and logs a warning, so the outlet lands at
/// `1e-6` bar above vacuum rather than at a negative pressure. That is the class's own
/// behaviour and this reproduces it, returning the applied drop so the caller can say so.
pub fn filter(feed: &Stream, pressure_drop: Pressure) -> Result<FilterOutcome> {
    // `Math.max(0.0, inletPressure - 1.0e-6)` in the class's bara: a millionth of a bar.
    let ceiling = (feed.p.value - 0.1).max(0.0);
    let applied = pressure_drop.value.clamp(0.0, ceiling);
    let p_out = pascals(feed.p.value - applied);

    // The class flashes only when the drop is above its own `1e-10` bar threshold; below it
    // the outlet is the cloned inlet untouched. A TP flash at the same state is that state,
    // so the two agree and this needs no second branch.
    Ok(FilterOutcome {
        outlet: Stream::from_pt(
            feed.components.clone(),
            feed.z.clone(),
            feed.n,
            p_out,
            feed.t,
        )?,
        applied_drop: pascals(applied),
    })
}
