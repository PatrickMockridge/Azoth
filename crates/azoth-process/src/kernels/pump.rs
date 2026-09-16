//! `unit_ops.pump` - a liquid raised to a higher pressure.

use azoth_core::AzothError;
use azoth_core::units::{Pressure, joules_per_mole};
use azoth_eos::{pr_molar_volume, pt_flash};

use crate::stream::Stream;

/// Raise a liquid stream's pressure, adding the pump's work as enthalpy.
///
/// A pump moves a liquid, whose molar volume is nearly constant, so the reversible
/// shaft work is `v·(P_out - P_in)` and the actual work divides it by the isentropic
/// efficiency. The outlet temperature is the one the mixture reaches at the shifted
/// enthalpy.
pub fn pump(
    feed: &Stream,
    outlet_pressure: Pressure,
    efficiency: f64,
) -> azoth_core::Result<Stream> {
    if !(efficiency > 0.0 && efficiency <= 1.0) {
        return Err(AzothError::invalid_input(
            "efficiency",
            "the pump's isentropic efficiency must be in (0, 1]",
        ));
    }

    let (mixture, _) = feed.mixture()?;
    let flash = pt_flash(&mixture, feed.t, feed.p, &feed.z)?;
    let molar_volume = pr_molar_volume(flash.z_liquid, feed.t, feed.p)?.v.value;

    let dh = molar_volume * (outlet_pressure.value - feed.p.value) / efficiency;
    Stream::from_ph(
        feed.components.clone(),
        feed.z.clone(),
        feed.n,
        outlet_pressure,
        joules_per_mole(feed.h.value + dh),
    )
}
