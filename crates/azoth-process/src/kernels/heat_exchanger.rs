//! `unit_ops.heat_exchanger` - a duty moved from a hot stream to a cold one.

use azoth_core::units::{Power, joules_per_mole};

use crate::stream::Stream;

/// Transfer a duty from the hot stream to the cold stream, at constant pressure
/// on each side. Each outlet's temperature is the one its mixture reaches at the
/// shifted enthalpy.
///
/// The two streams may be different fluids; each resolves its own mixture. The
/// duty is positive when heat flows from hot to cold.
pub fn heat_exchanger(
    hot: &Stream,
    cold: &Stream,
    duty: Power,
) -> azoth_core::Result<(Stream, Stream)> {
    let hot_out = Stream::from_ph(
        hot.components.clone(),
        hot.z.clone(),
        hot.n,
        hot.p,
        joules_per_mole(hot.h.value - duty.value / hot.n),
    )?;
    let cold_out = Stream::from_ph(
        cold.components.clone(),
        cold.z.clone(),
        cold.n,
        cold.p,
        joules_per_mole(cold.h.value + duty.value / cold.n),
    )?;
    Ok((hot_out, cold_out))
}
