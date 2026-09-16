//! `unit_ops.separator` - a flash, split into vapour and liquid outlets.

use azoth_core::units::ThermodynamicTemperature;
use azoth_eos::{Phase, molar_enthalpy_entropy, pt_flash};

use crate::stream::Stream;

/// Flash a feed into vapour and liquid outlets at the given temperature.
///
/// The split is the equilibrium one from [`azoth_eos::pt_flash`] at the feed's
/// pressure and the separator's temperature. The vapour and liquid outlets each
/// carry their phase's composition and molar enthalpy; a single-phase feed empties
/// one outlet rather than inventing a split that does not exist.
pub fn separator(
    feed: &Stream,
    temperature: ThermodynamicTemperature,
) -> azoth_core::Result<(Stream, Stream)> {
    let (mixture, ideal_gas) = feed.mixture()?;
    let flash = pt_flash(&mixture, temperature, feed.p, &feed.z)?;

    let beta = match flash.phase {
        Phase::TwoPhase => flash.beta.unwrap_or(0.0),
        Phase::AllVapour => 1.0,
        Phase::AllLiquid | Phase::Trivial => 0.0,
    };

    let vapour = Stream {
        components: feed.components.clone(),
        z: flash.y.clone(),
        n: feed.n * beta,
        p: feed.p,
        t: temperature,
        h: molar_enthalpy_entropy(
            &mixture,
            &ideal_gas,
            temperature,
            feed.p,
            &flash.y,
            flash.z_vapour,
        )?
        .h,
    };
    let liquid = Stream {
        components: feed.components.clone(),
        z: flash.x.clone(),
        n: feed.n * (1.0 - beta),
        p: feed.p,
        t: temperature,
        h: molar_enthalpy_entropy(
            &mixture,
            &ideal_gas,
            temperature,
            feed.p,
            &flash.x,
            flash.z_liquid,
        )?
        .h,
    };

    Ok((vapour, liquid))
}
