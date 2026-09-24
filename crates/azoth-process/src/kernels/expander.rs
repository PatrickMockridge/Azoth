//! `unit_ops.expander` - the isentropic route with the efficiency on the other side.

use azoth_core::units::{Pressure, joules_per_mole, joules_per_mole_kelvin};
use azoth_core::{AzothError, Result};
use azoth_eos::{ph_flash, ps_flash};

use crate::stream::Stream;

/// Drop a stream's pressure along an isentrope, multiplying the step by the efficiency.
///
/// **One operator away from `unit_ops.compressor`, and the operator is the point.**
/// `Expander.run`'s steady-state branch is `Compressor.run`'s with `* isentropicEfficiency`
/// where the compressor has `/ isentropicEfficiency` - and the two are the same physical rule,
/// because an expansion's isentropic enthalpy difference is *negative*: dividing it would
/// make the machine recover more work than the reversible one, with an efficiency below one.
/// At `eta = 1` both reduce to the isentropic step, which is the row each capture pins with an
/// entropy production of zero.
///
/// The route itself is `Pump.run`'s: the inlet's own entropy - derived from `(T, P, z)` and
/// not carried on the record - flashed to the outlet pressure for the reversible enthalpy, and
/// a second flash at the enthalpy the efficiency leaves.
///
/// Named out: the polytropic branch's stepped integration (`numbersteps` increments of `dp`,
/// each one a `PSflash` and a `PHflash`), the `powerSet` mode that takes a shaft power instead
/// of an efficiency, the `useOutTemperature` mode, and the same chart, anti-surge and property
/// overrides the compressor names.
pub fn expander(feed: &Stream, outlet_pressure: Pressure, efficiency: f64) -> Result<Stream> {
    if !(efficiency > 0.0 && efficiency <= 1.0) {
        return Err(AzothError::invalid_input(
            "isentropic_efficiency",
            "the expander's isentropic efficiency must be in (0, 1]",
        ));
    }

    let (mixture, ideal_gas) = feed.mixture()?;
    let (entropy, _) = ps_flash::entropy_at(&mixture, &ideal_gas, feed.t, feed.p, &feed.z)?;

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

    let dh = (h_isentropic - feed.h.value) * efficiency;
    Stream::from_ph(
        feed.components.clone(),
        feed.z.clone(),
        feed.n,
        outlet_pressure,
        joules_per_mole(feed.h.value + dh),
    )
}
