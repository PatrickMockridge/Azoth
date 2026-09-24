//! `unit_ops.heater` - a stream moved to a stated temperature or duty.

use azoth_core::units::{
    Power, Pressure, ThermodynamicTemperature, joules_per_mole, pascals, watts,
};
use azoth_core::{AzothError, Result};

use crate::stream::Stream;

/// What a heater did: the outlet record, and the duty it moved.
pub struct HeaterOutcome {
    /// The outlet record: the inlet's flow and composition at the new state.
    pub outlet: Stream,
    /// The duty it moved, W. Negative is heat removed, which is what a cooler is.
    pub duty: Power,
}

/// Heat or cool a stream to a stated temperature, or by a stated duty.
///
/// **The duty is in W and the enthalpy is per mole, and they meet because a NeqSim process
/// fluid's moles are one second of flow.** `Heater.run` adds `getDuty()` to
/// `system.getEnthalpy()` directly - W to J - which balances only while the fluid's total
/// moles equal its mol/s rate. Here the conversion is explicit, `duty / n` being the shift
/// in J/mol, so the convention is stated rather than inherited.
///
/// **The reported duty is the enthalpy the outlet reached, not the one the caller stated.**
/// `run` overwrites its own `energyInput` with `newH - oldH` after the flash, whatever branch
/// it took, so a temperature-specified heater still reports a duty - the state's, computed
/// from what the flash did.
///
/// **A temperature and a duty together are refused, and that is a measurement rather than a
/// preference.** `setOutletTemperature` sets `setTemperature` and clears `setEnergyInput`;
/// `setEnergyInput` (and so `setDuty`) sets `setEnergyInput` and clears `setTemperature`. The
/// pair is therefore *order-dependent* - the last setter called decides - which is a state
/// the declaration cannot express. Refusing it keeps this model a function of its declared
/// inputs; picking one would be inventing a rule the class does not have.
///
/// **With neither, the drop is isothermal.** `run`'s else branch is `T_in + dT` with `dT`
/// defaulting to zero, so a heater that specifies nothing still moves the pressure.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if a temperature and a duty are both given, if the
/// pressure drop would leave a non-positive pressure, or if a duty is asked of a stream
/// carrying no flow, and whatever the flash refuses.
pub fn heater(
    feed: &Stream,
    outlet_temperature: Option<ThermodynamicTemperature>,
    duty: Option<Power>,
    pressure_drop: Option<Pressure>,
) -> Result<HeaterOutcome> {
    if outlet_temperature.is_some() && duty.is_some() {
        return Err(AzothError::invalid_input(
            "duty",
            "a heater holds a temperature specification or a duty specification and not \
             both: `Heater.setOutletTemperature` clears the duty flag and `setDuty` clears \
             the temperature flag, so which one the class runs is the order they were set in \
             and not a function of these inputs",
        ));
    }
    let p_out = match pressure_drop {
        Some(drop) => pascals(feed.p.value - drop.value),
        None => feed.p,
    };
    if p_out.value <= 0.0 {
        return Err(AzothError::invalid_input(
            "pressure_drop",
            format!(
                "a pressure drop of {} Pa takes a {} Pa inlet to {} Pa, which is not a pressure",
                pressure_drop.map_or(0.0, |drop| drop.value),
                feed.p.value,
                p_out.value
            ),
        ));
    }

    let outlet = match (outlet_temperature, duty) {
        (Some(t), _) => Stream::from_pt(feed.components.clone(), feed.z.clone(), feed.n, p_out, t)?,
        (None, Some(power)) => {
            if feed.n == 0.0 {
                return Err(AzothError::invalid_input(
                    "duty",
                    "a duty is moved per mole of flow, and this stream carries none",
                ));
            }
            Stream::from_ph(
                feed.components.clone(),
                feed.z.clone(),
                feed.n,
                p_out,
                joules_per_mole(feed.h.value + power.value / feed.n),
            )?
        }
        (None, None) => Stream::from_pt(
            feed.components.clone(),
            feed.z.clone(),
            feed.n,
            p_out,
            feed.t,
        )?,
    };

    Ok(HeaterOutcome {
        duty: watts(feed.n * (outlet.h.value - feed.h.value)),
        outlet,
    })
}
