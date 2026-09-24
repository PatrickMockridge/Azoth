//! `unit_ops.tank` - the two-phase flash a tank's steady state is.

use azoth_core::Result;
use azoth_core::units::pascals;

use crate::kernels::mixer::mixer;
use crate::kernels::separator::separator;
use crate::stream::Stream;

/// Split a tank's contents into a gas and a liquid outlet.
///
/// **The steady state is a `VUflash` at the fluid's own volume and internal energy, and a
/// fluid at its stable state re-imposing its own `V` and `U` returns that state.**
/// `Tank.run` calls `ops.VUflash(thermoSystem2.getVolume(), thermoSystem2.getInternalEnergy())`
/// on the *fluid's* volume rather than the vessel's, and the capture's two-phase row is
/// `process.separator`'s first row to the last digit, at every field of both outlets.
///
/// **The design volume never reaches it.** `setVolume` is read by `getVolume`,
/// `validateMechanicalDesign`, `toJson` and `displayResult`, and by nothing in `run`. Two
/// capture rows that differ only in `setVolume` are identical, which is the measurement
/// that says so - and is why the palette entry declares no `volume`.
///
/// **`run` does not enable `multiPhaseCheck`, and a tank is not a three-phase machine.**
/// Its two tests name "gas" and "oil", and its flash is the one place in this family that
/// does not ask for a third phase. On the capture's three-phase-capable feed the VU flash
/// finds two: the "oil" outlet carries the water, and the gas leaves at a water mole
/// fraction of `0.234` where a `ThreePhaseSeparator` at the same state leaves at `0.0016`.
///
/// **The single-phase branch is not reproduced, and the reason is the class's own
/// construction.** `run`'s absent-oil branch writes `gasOutStream` where `liquidOutStream`
/// was meant, and `setInletStream` builds the liquid outlet from `getPhases()[1]` whatever
/// the phase count is - so on one phase the two outlets come back *swapped* against
/// `Separator`'s on the same feed. A kernel that is a pure function of its inlets has no
/// construction-time state to swap, so the port answers with the flash, and the capture's
/// two single-phase rows are the declared divergence.
///
/// # Errors
/// Whatever [`mixer`] and [`separator`] raise: feeds carrying no flow in total, and
/// whatever the databank or the flash refuses.
pub fn tank(feeds: &[Stream]) -> Result<(Stream, Stream)> {
    let joined = mixer(feeds, None)?;
    // Zero drop, no entrainment, no heat input: the whole of what a tank's steady state
    // adds to a flash is that it holds the feed's pressure and temperature while doing it.
    separator(&joined, pascals(0.0), 0.0, None)
}
