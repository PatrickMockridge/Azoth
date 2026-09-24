//! `unit_ops.gas_scrubber` - `unit_ops.separator`'s arithmetic, under the other entry.

use azoth_core::Result;
use azoth_core::units::{Power, Pressure};

use crate::kernels::separator::separator;
use crate::stream::Stream;

/// Flash a feed into vapour and liquid outlets, as a scrubber does.
///
/// **This is a delegation and the delegation is the port.** `GasScrubber extends Separator`
/// and **does not override `run`**: its own 112 lines are constructors, a mechanical design
/// and the capacity metric, and the steady state is the separator's flash exactly. So the
/// arithmetic is [`separator`]'s, and writing it out again would be a second implementation
/// of one vessel.
///
/// **What the class adds is the Souders-Brown capacity metric, and it is not ported.**
/// `getCapacityUtilization` is a function of the vapour's volumetric flow, the liquid's
/// density and two *mechanical* parameters - the internal diameter (`setInternalDiameter`)
/// and a design gas load factor (`setDesignGasLoadFactor`, `0.04`-`0.10` m/s for a vertical
/// scrubber) - which the palette entry declares neither of. It is a check on whether the
/// vessel is big enough for the stream, not a statement about the stream, and it belongs to
/// the mechanical design the tranche declares out.
///
/// # Errors
/// The separator's own: a pressure drop that leaves a non-positive pressure, an entrainment
/// fraction outside `[0, 1]`, and whatever the flash refuses.
pub fn gas_scrubber(
    feed: &Stream,
    pressure_drop: Pressure,
    gas_in_liquid: f64,
    heat_input: Option<Power>,
) -> Result<(Stream, Stream)> {
    separator(feed, pressure_drop, gas_in_liquid, heat_input)
}
