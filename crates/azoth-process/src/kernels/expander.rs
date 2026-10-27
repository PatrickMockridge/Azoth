//! `unit_ops.expander` - the isentropic route with the efficiency on the other side.

use azoth_core::units::{Pressure, joules_per_mole, joules_per_mole_kelvin, pascals};
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
/// Named out: the `powerSet` mode that takes a shaft power instead of an efficiency, the
/// `useOutTemperature` mode, and the same chart, anti-surge and property overrides the
/// compressor names. The polytropic branch is [`expander_polytropic`].
pub fn expander(feed: &Stream, outlet_pressure: Pressure, efficiency: f64) -> Result<Stream> {
    if !(efficiency > 0.0 && efficiency <= 1.0) {
        return Err(AzothError::invalid_input(
            "isentropic_efficiency",
            "the expander's isentropic efficiency must be in (0, 1]",
        ));
    }
    expander_replay(feed, outlet_pressure, efficiency)
}

/// The isentropic route without the upper bound, for the turbo-expander's replay.
///
/// **The coupled machine assigns its solved efficiency to the field directly**, bypassing
/// `setIsentropicEfficiency`'s `(0, 1]` clamp because the outlet-temperature inverse needs a
/// required efficiency that can exceed one on the way to a back-calculated design value. The
/// public [`expander`] keeps the clamp; this is the path the machine's replay takes.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] for a non-positive efficiency, and whatever the
/// flash refuses.
pub fn expander_replay(
    feed: &Stream,
    outlet_pressure: Pressure,
    efficiency: f64,
) -> Result<Stream> {
    if efficiency.is_nan() || efficiency <= 0.0 {
        return Err(AzothError::invalid_input(
            "isentropic_efficiency",
            format!("the replay needs a positive efficiency, and {efficiency} is not"),
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

/// Drop a stream's pressure in steps, multiplying each isentropic step by the efficiency.
///
/// **`Expander.run`'s polytropic branch, `Expander.java:620-633`.** The pressure falls in
/// `numbersteps` equal increments; each step takes the *current* entropy and enthalpy, flashes
/// isentropically to the next pressure, and sets the real enthalpy to
/// `h + (h_is - h) eta` before a `PHflash` lands the state. The expander **multiplies** where
/// the compressor divides, because an expansion's isentropic difference is negative.
///
/// **The step count is floored at one**, which is the class's `Math.max(1, ...)`; a count of
/// zero is one step rather than an error.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the efficiency is outside `(0, 1]`, and
/// whatever the flash refuses.
pub fn expander_polytropic(
    feed: &Stream,
    outlet_pressure: Pressure,
    polytropic_efficiency: f64,
    number_of_steps: usize,
) -> Result<Stream> {
    if !(polytropic_efficiency > 0.0 && polytropic_efficiency <= 1.0) {
        return Err(AzothError::invalid_input(
            "polytropic_efficiency",
            "the expander's polytropic efficiency must be in (0, 1]",
        ));
    }

    let steps = number_of_steps.max(1);
    let (mixture, ideal_gas) = feed.mixture()?;
    let increment = (outlet_pressure.value - feed.p.value) / steps as f64;
    let mut state = feed.clone();
    for _ in 0..steps {
        let (entropy, _) = ps_flash::entropy_at(&mixture, &ideal_gas, state.t, state.p, &state.z)?;
        let next_pressure = pascals(state.p.value + increment);
        let isentropic = ps_flash::ps_flash(
            &mixture,
            &ideal_gas,
            next_pressure,
            joules_per_mole_kelvin(entropy),
            &state.z,
        )?;
        let (h_isentropic, _) = ph_flash::enthalpy_at(
            &mixture,
            &ideal_gas,
            isentropic.temperature,
            next_pressure,
            &state.z,
        )?;
        let enthalpy = state.h.value + (h_isentropic - state.h.value) * polytropic_efficiency;
        state = Stream::from_ph(
            state.components.clone(),
            state.z.clone(),
            state.n,
            next_pressure,
            joules_per_mole(enthalpy),
        )?;
    }
    Ok(state)
}
