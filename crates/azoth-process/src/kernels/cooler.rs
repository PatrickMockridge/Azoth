//! `unit_ops.cooler` - `unit_ops.heater`'s arithmetic, with the sign a cooler implies.

use azoth_core::units::{Power, Pressure, ThermodynamicTemperature};
use azoth_core::Result;

use crate::kernels::heater::{HeaterOutcome, heater};
use crate::stream::Stream;

/// Cool a stream to a stated temperature, or by a stated duty.
///
/// **This is a delegation and the delegation is the port.** `Cooler extends Heater` and
/// overrides `runTransient`, the mechanical design and a handful of getters - **not `run`**,
/// which is the whole of the steady-state arithmetic. Its own 258 lines are a dynamic
/// valve/NTU temperature-control model that does not touch the outlet record, and the
/// caption to that claim is a measurement rather than a reading: `ProcessProbe cooler` runs
/// the same six rows as `ProcessProbe heater` through this class and the two captures are
/// **byte-identical**.
///
/// Writing the arithmetic out again here would be a second implementation of one machine's
/// physics and the two could then disagree - which is the same reason `GasScrubber`'s
/// stream-side kernel is `Separator`'s and not a copy of it.
pub fn cooler(
    feed: &Stream,
    outlet_temperature: Option<ThermodynamicTemperature>,
    duty: Option<Power>,
    pressure_drop: Option<Pressure>,
) -> Result<HeaterOutcome> {
    heater(feed, outlet_temperature, duty, pressure_drop)
}
