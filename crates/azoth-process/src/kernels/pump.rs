//! `unit_ops.pump` - a liquid raised to a higher pressure.

use azoth_core::AzothError;
use azoth_core::units::{Pressure, joules_per_mole, joules_per_mole_kelvin};
use azoth_eos::{ph_flash, ps_flash};

use crate::stream::Stream;

/// Raise a liquid stream's pressure, adding the pump's work as enthalpy.
///
/// **The work is the isentropic head, taken from a flash.** NeqSim's `Pump.run` clones the
/// fluid to the outlet pressure and runs a `PSflash` at the inlet entropy, so the work is
/// `(H(P_out, s_in) - H(P_in)) / eta` — the exact enthalpy a real fluid needs to reach that
/// pressure along an isentrope, not the incompressible approximation `v (P_out - P_in)`.
/// The two agree to `1.2e-3` on pure n-butane across 5 to 20 bar, and this is the one
/// NeqSim computes: `calculateAsCompressor` defaults to `true`, which selects exactly that
/// branch. The incompressible form is the linearisation `H(P_out, s_in) - H(P_in)` is not,
/// and the difference is what an earlier version of this kernel carried.
///
/// The outlet temperature is then the one at which the mixture carries the shifted
/// enthalpy, which is a second flash and not an assumption of isothermal compression.
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

    let (mixture, ideal_gas) = feed.mixture()?;
    let (entropy, _) = ps_flash::entropy_at(&mixture, &ideal_gas, feed.t, feed.p, &feed.z)?;

    // The isentropic outlet: the temperature the fluid reaches at the raised pressure
    // carrying the inlet's entropy.
    let isentropic = ps_flash::ps_flash(
        &mixture,
        &ideal_gas,
        outlet_pressure,
        joules_per_mole_kelvin(entropy),
        &feed.z,
    )?;
    let (h_isentropic, _) = ph_flash::enthalpy_at(
        &mixture,
        &ideal_gas,
        isentropic.temperature,
        outlet_pressure,
        &feed.z,
    )?;

    let dh = (h_isentropic - feed.h.value) / efficiency;
    Stream::from_ph(
        feed.components.clone(),
        feed.z.clone(),
        feed.n,
        outlet_pressure,
        joules_per_mole(feed.h.value + dh),
    )
}
