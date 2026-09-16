//! `unit_ops.throttling_valve` - an isenthalpic pressure drop.

use azoth_core::units::Pressure;

use crate::stream::Stream;

/// Drop a stream to a lower pressure without heat or work, so its molar enthalpy
/// is unchanged and its temperature is the one the mixture reaches at that state.
pub fn throttling_valve(feed: &Stream, outlet_pressure: Pressure) -> azoth_core::Result<Stream> {
    Stream::from_ph(
        feed.components.clone(),
        feed.z.clone(),
        feed.n,
        outlet_pressure,
        feed.h,
    )
}
