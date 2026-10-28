//! `eos.sound_speed` - the cubic's isentropic speed of sound.
//!
//! ```text
//! c^2 = -gamma (v^2 / M) (dP/dv)_T        gamma = cp/cv
//! cp - cv = -T (dP/dT)^2_v / (dP/dv)_T
//! ```
//!
//! **NeqSim's formula, in SI.** `PhaseEos.getSoundSpeed` is
//! `bs = -(1/V)(Cv/Cp)/(dPdVTn)`, `sqrt(V/Mw/bs)` in the class's internal volume units;
//! the 1e5 volume and the 1e-10 inside `dPdVTn` cancel, so the answer is the expression
//! above with the phase's molar volume and molar mass. The two derivatives are closed forms
//! of the same residual the departure functions already evaluate:
//!
//! ```text
//! P = R T/(v - b) - a/((v + d1 b)(v + d2 b))
//! (dP/dv)_T = -R T/(v - b)^2 + a (2v + (d1 + d2) b) / f^2,   f = (v + d1 b)(v + d2 b)
//! (dP/dT)_v = R/(v - b) - (da/dT)/f
//! da/dT = a psi_bar/T
//! ```
//!
//! `psi_bar` is `dln a/dln T` for the mixture - the composition-weighted average of the
//! components' `dln alpha/dln T`, which `PhaseState` already carries for the departure heat
//! capacity. The reduced `A = aP/(RT)^2` carries the `T^-2` of `Tr` and cancels the `(RT)^2`
//! of `a`, so the dimensional attraction's logarithmic slope is `psi_bar` itself.
//!
//! **The system value is the class's beta-weighted mean**, not a mixture speed - see the
//! spec's assumptions.

use azoth_core::units::{Pressure, ThermodynamicTemperature, meters_per_second};
use azoth_core::{AzothError, Result, apply_checks};

use crate::cubic::Cubic;
use crate::databank;
use crate::mixture::{Mixture, ReducedParameters, RootSide};
use crate::molar_enthalpy_entropy;
use crate::pr_molar_volume::MOLAR_GAS_CONSTANT;
use crate::pt_flash;

use crate::model_gen;
use crate::results::{Phase, SoundSpeedResult};

/// The speed of sound in a mixture at a state.
///
/// # Errors
/// [`AzothError::InvalidInput`] if the composition's length disagrees with the components,
/// or the state gives no admissible derivative; [`AzothError::OutOfRange`] for a
/// non-positive temperature or pressure; and whatever the flash or the databank refuses.
pub fn sound_speed(
    components: &[String],
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<SoundSpeedResult> {
    let spec = &model_gen::SOUND_SPEED_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "P" => Some(p.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, ideal_gas) = databank::mixture_of(&names, Cubic::Pr, None)?;
    if z.len() != mixture.len() {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "a mixture of {} components needs {} mole fractions, and {z:?} has {}",
                mixture.len(),
                mixture.len(),
                z.len()
            ),
        ));
    }

    let reduced = mixture.reduced_parameters(t, p)?;
    let flash = pt_flash::pt_flash(&mixture, t, p, z)?;
    let (phase_speed, system_speed) = match flash.phase {
        Phase::TwoPhase => {
            let beta = flash.vapour_fraction.ok_or_else(|| {
                AzothError::invalid_input(
                    "state",
                    "the flash split the feed and reported no vapour fraction, so the system \
                     speed has no weight for either phase",
                )
            })?;
            let vapour = phase_speed_of(
                &mixture,
                &ideal_gas,
                &reduced,
                t,
                p,
                &flash.y,
                flash.z_vapour,
            )?;
            let liquid = phase_speed_of(
                &mixture,
                &ideal_gas,
                &reduced,
                t,
                p,
                &flash.x,
                flash.z_liquid,
            )?;
            (vapour, beta * vapour + (1.0 - beta) * liquid)
        }
        Phase::AllVapour | Phase::AllLiquid | Phase::Trivial => {
            // **The feed's own root, at the feed's composition.** `pt_flash` on a
            // superheated stream reports the *incipient* phase's root, whose composition can
            // sit far from the feed's; the class reads the one phase it holds.
            let side = if matches!(flash.phase, Phase::AllVapour) {
                RootSide::Vapour
            } else {
                RootSide::Liquid
            };
            let root = mixture.phase_state(&reduced, z, side)?.z;
            let speed = phase_speed_of(&mixture, &ideal_gas, &reduced, t, p, z, root)?;
            (speed, speed)
        }
    };

    Ok(SoundSpeedResult {
        sound_speed: meters_per_second(phase_speed),
        system_sound_speed: meters_per_second(system_speed),
        phase: flash.phase,
        warnings,
    })
}

/// One phase's isentropic speed, m/s.
fn phase_speed_of(
    mixture: &Mixture,
    ideal_gas: &crate::molar_enthalpy_entropy::IdealGasModel,
    reduced: &ReducedParameters,
    t: ThermodynamicTemperature,
    p: Pressure,
    composition: &[f64],
    root: f64,
) -> Result<f64> {
    let state = mixture.phase_state_at(reduced, composition, root)?;
    let molar_mass = molar_mass_of(mixture, composition)?;
    let constant = MOLAR_GAS_CONSTANT;
    let v = root * constant * t.value / p.value;
    let (delta1, delta2) = (mixture.cubic().delta1(), mixture.cubic().delta2());

    // `A` is dimensionless (`aP/(RT)^2`), so the dimensional `a` and `b` are these.
    let a = state.a_mix * (constant * t.value) * (constant * t.value) / p.value;
    let b = state.b_mix * constant * t.value / p.value;
    let d_a_d_t = a * state.psi_bar / t.value;
    let f = (v + delta1 * b) * (v + delta2 * b);
    let d_f = 2.0 * v + (delta1 + delta2) * b;
    let d_p_d_v = -constant * t.value / ((v - b) * (v - b)) + a * d_f / (f * f);
    let d_p_d_t = constant / (v - b) - d_a_d_t / f;
    if d_p_d_v.is_nan() || d_p_d_v >= 0.0 || !d_p_d_v.is_finite() || !d_p_d_t.is_finite() {
        return Err(AzothError::invalid_input(
            "state",
            format!(
                "the cubic's pressure-volume derivative is {d_p_d_v} Pa m3/mol at {:.3} K and \
                 {:.3} Pa, and a sound speed needs a stable phase",
                t.value, p.value
            ),
        ));
    }

    let cp = molar_enthalpy_entropy::molar_enthalpy_entropy(
        mixture,
        ideal_gas,
        t,
        p,
        composition,
        root,
    )?
    .cp
    .value;
    let cv = cp + t.value * d_p_d_t * d_p_d_t / d_p_d_v;
    if cv.is_nan() || cv <= 0.0 || !cp.is_finite() {
        return Err(AzothError::invalid_input(
            "state",
            format!("the phase's heat capacities are cp={cp}, cv={cv}, and both must be positive"),
        ));
    }
    let gamma = cp / cv;
    Ok(v * (-gamma * d_p_d_v / molar_mass).sqrt())
}

/// The composition's molar mass, kg/mol.
fn molar_mass_of(mixture: &Mixture, composition: &[f64]) -> Result<f64> {
    let mut total = 0.0;
    for (component, fraction) in mixture.components().iter().zip(composition) {
        let mass = component.molar_mass.ok_or_else(|| {
            AzothError::property_unavailable(
                "components",
                "molar_mass",
                "a component of this mixture carries no molar mass, so the phase has none",
            )
        })?;
        total += fraction * mass;
    }
    Ok(total)
}
