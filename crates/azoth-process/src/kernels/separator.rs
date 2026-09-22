//! `unit_ops.separator` - a flash, split into vapour and liquid outlets.

use azoth_core::units::{Power, Pressure, joules_per_mole, pascals};
use azoth_core::{AzothError, Result};
use azoth_eos::{Phase, ph_flash, pt_flash};

use crate::stream::Stream;

/// Flash a feed into vapour and liquid outlets.
///
/// The flash is at `P_in - pressure_drop` and at the **feed's own temperature**, which is
/// the vessel's: `Separator.run` reduces the pressure and runs a `TPflash`, so a pressure
/// drop is not a throttling and the outlets do not carry the inlet's enthalpy. A
/// `heat_input`, when given, moves the flash to the enthalpy it implies instead and the
/// outlet temperature is solved for.
///
/// `gas_in_liquid` carries that fraction of the vapour's moles into the liquid outlet,
/// which is NeqSim's `gasInLiquid` entrainment on its `feed`/`mole` basis. It moves
/// material between the phases and not energy: both outlets are then flashed at the same
/// temperature, so the balance across the transfer is not closed - which is what the class
/// does and is stated in the spec rather than corrected here.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the pressure drop would take the outlet to
/// a non-positive pressure or the entrainment fraction is outside `[0, 1]`, and whatever
/// the flash refuses.
pub fn separator(
    feed: &Stream,
    pressure_drop: Pressure,
    gas_in_liquid: f64,
    heat_input: Option<Power>,
) -> Result<(Stream, Stream)> {
    if !(0.0..=1.0).contains(&gas_in_liquid) {
        return Err(AzothError::invalid_input(
            "gas_in_liquid",
            format!("an entrainment fraction is a fraction, and {gas_in_liquid} is not in [0, 1]"),
        ));
    }
    let p_out = pascals(feed.p.value - pressure_drop.value);
    if p_out.value <= 0.0 {
        return Err(AzothError::invalid_input(
            "pressure_drop",
            format!(
                "a pressure drop of {} Pa takes a {} Pa inlet to {} Pa, which is not a pressure",
                pressure_drop.value, feed.p.value, p_out.value
            ),
        ));
    }

    let (mixture, ideal_gas) = feed.mixture()?;

    let (temperature, x, y, beta) = match heat_input.filter(|power| power.value != 0.0) {
        Some(power) => {
            // One second of flow, for the reason `Heater.run` gives: the duty is W and the
            // enthalpy it moves is per mole, so it divides by the same flow a case states
            // in mol/s.
            let target = feed.h.value + power.value / feed.n;
            let flash = ph_flash::ph_flash(
                &mixture,
                &ideal_gas,
                p_out,
                joules_per_mole(target),
                &feed.z,
            )?;
            (
                flash.temperature,
                flash.x,
                flash.y,
                vapour_fraction(flash.phase, flash.beta),
            )
        }
        None => {
            let flash = pt_flash(&mixture, feed.t, p_out, &feed.z)?;
            (
                feed.t,
                flash.x,
                flash.y,
                vapour_fraction(flash.phase, flash.beta),
            )
        }
    };

    let n_vapour = feed.n * beta;
    let n_liquid = feed.n * (1.0 - beta);
    // NeqSim moves material between the phases only where both exist: `addPhaseFractionToPhase`
    // returns unchanged unless the from-phase and the to-phase are both present, so a
    // single-phase feed is untouched whatever the fraction says.
    let moved = if n_vapour > 0.0 && n_liquid > 0.0 {
        n_vapour * gas_in_liquid
    } else {
        0.0
    };

    let liquid_n = n_liquid + moved;
    let liquid_z: Vec<f64> = if moved > 0.0 {
        (0..x.len())
            .map(|i| (n_liquid * x[i] + moved * y[i]) / liquid_n)
            .collect()
    } else {
        x
    };

    // Each outlet is a *stream*, so its enthalpy is the state's at its own composition,
    // pressure and temperature - not the phase root the flash happened to report.
    Ok((
        Stream::from_pt(
            feed.components.clone(),
            y,
            n_vapour - moved,
            p_out,
            temperature,
        )?,
        Stream::from_pt(
            feed.components.clone(),
            liquid_z,
            liquid_n,
            p_out,
            temperature,
        )?,
    ))
}

/// The vapour fraction a flash's phase and `beta` imply.
///
/// `beta` is absent for a single-phase answer - the flash reports `None` rather than a
/// number outside `[0, 1]`, because there is no split to report - so the phase is what
/// says which side of the interval the outlet is on.
fn vapour_fraction(phase: Phase, beta: Option<f64>) -> f64 {
    match phase {
        Phase::TwoPhase => beta.unwrap_or(0.0),
        Phase::AllVapour => 1.0,
        Phase::AllLiquid | Phase::Trivial => 0.0,
    }
}
