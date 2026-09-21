//! `eos.tp_multiflash_wax` - how much of a feed is wax at a state.
//!
//! ```text
//! seed   the two-phase flash, plus a wax phase
//! solve  Q(beta) = sum_k beta_k - sum_i z_i ln E_i      E_i = sum_k beta_k / phi_ik
//! keep   the wax only where the solve converges and it is above the floor
//! ```
//!
//! Spec: `specs/models/eos/tp_multiflash_wax.toml`, which carries the seeding, the floor's
//! rule and the four states it is checked at.
//!
//! NeqSim's `TPmultiflashWAX`. **It adds one phase to [`crate::tp_multiflash`] and reuses the
//! solve** rather than restating it: the fraction Newton is Michelsen's on the whole vector,
//! and a phase whose coefficients are known is a phase it can carry whatever their source.
//!
//! # What is different about a solid
//!
//! A cubic phase's coefficients are a function of *its own composition* on a chosen root, so
//! they are rebuilt at every trial. The wax solid's are
//! [`crate::wax_solid_fugacity`]'s - `phi_liq exp(...)`, with the mole fraction cancelled out
//! of `SolidFug/(P x)` - so they depend on the state and the component alone and are one
//! vector throughout. **And a substance that is not a wax former is excluded by a number**:
//! NeqSim's `ComponentWax.fugcoef` returns `1e50` for one, and `x_i = z_i/(E_i phi_i)` is then
//! zero to any precision the answer is read at.
//!
//! # The two-phase fallback
//!
//! Where no wax forms, the three-phase set does not converge - the Newton has a phase in it
//! that has nowhere to go - and the answer is the two-phase flash's, which is the same state.
//! `converged` is false in that case and says which of the two was solved.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result, apply_checks};

use crate::algorithm_of;
use crate::mixture::{Mixture, RootSide};
use crate::model_gen;
use crate::multiphase::{MultiphasePhase, PhaseKind, solve_phase_fractions};
use crate::pt_flash::pt_flash;
use crate::results::TpMultiflashWaxResult;

/// Below this a phase is not there.
///
/// The fraction solve drives a phase that should not exist to its floor rather than to a
/// negative amount, so this is the line between "the wax is a phase" and "the wax is what the
/// solve could not remove": the same `1.1e-12` `eos.tp_multiflash`'s merge uses.
const WAX_FLOOR: f64 = 1.1e-12;

/// The fraction of a feed that is wax at a temperature and pressure.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `z` is the wrong length, negative, or does not sum to
///   one; or if a wax former carries no melt data.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
#[allow(non_snake_case)] // `T`, `P` and `z` are the symbols in the chemistry
pub fn tp_multiflash_wax(
    mixture: &Mixture,
    T: ThermodynamicTemperature,
    P: Pressure,
    z: &[f64],
) -> Result<TpMultiflashWaxResult> {
    let spec = &model_gen::TP_MULTIFLASH_WAX_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T.value),
            "P" => Some(P.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let n = mixture.len();
    if z.len() != n {
        return Err(AzothError::invalid_input(
            "z",
            format!("a feed for {n} components has {} entries", z.len()),
        ));
    }
    if let Some(bad) = z.iter().position(|&value| value < 0.0) {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "z[{bad}] is {} but a mole fraction cannot be negative",
                z[bad]
            ),
        ));
    }
    let sum: f64 = z.iter().sum();
    if (sum - 1.0).abs() > 1.0e-09 {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "the mole fractions sum to {sum}, not to one. Renormalising them here would \
                 make a composition error invisible in every number downstream, so it is \
                 refused instead"
            ),
        ));
    }

    let algorithm = algorithm_of(spec)?;
    let reduced = mixture.reduced_parameters(T, P)?;

    // **The two-phase flash seeds it**, and is also the answer if the wax does not survive.
    let flash = pt_flash(mixture, T, P, z)?;
    warnings.extend(flash.warnings.iter().cloned());
    let split_of_flash = flash.beta.unwrap_or(0.0);
    let two_phase = |iterations: u32, residual: f64, converged: bool| TpMultiflashWaxResult {
        wax_fraction: 0.0,
        phase_count: 2,
        beta: vec![1.0 - split_of_flash, split_of_flash],
        x: vec![flash.x.clone(), flash.y.clone()],
        iterations,
        residual,
        converged,
        warnings: warnings.clone(),
    };

    // A feed with no wax former in it has no solid to find, and the two-phase answer is the
    // whole of it. NeqSim reaches the same place through a phase whose every coefficient is
    // its `1e50` marker, which the solve can only take to zero.
    if !mixture.components().iter().any(|c| c.wax_former) {
        return Ok(two_phase(0, 0.0, true));
    }

    let mut phases = vec![
        MultiphasePhase {
            fraction: 1.0 - split_of_flash,
            composition: flash.x.clone(),
            kind: PhaseKind::Cubic(RootSide::Liquid),
        },
        MultiphasePhase {
            fraction: split_of_flash,
            composition: flash.y.clone(),
            kind: PhaseKind::Cubic(RootSide::Vapour),
        },
        MultiphasePhase {
            // The wax's starting share is the feed's own composition: there is no trial
            // composition to make, because the phase's coefficients do not depend on one.
            fraction: WAX_FLOOR,
            composition: z.to_vec(),
            kind: PhaseKind::Wax,
        },
    ];

    let split = solve_phase_fractions(mixture, &reduced, z, &mut phases, algorithm)?;
    let wax = split.fractions[2];
    if !split.converged || wax <= WAX_FLOOR {
        return Ok(two_phase(split.iterations, split.residual, split.converged));
    }

    Ok(TpMultiflashWaxResult {
        wax_fraction: wax,
        phase_count: split.fractions.len() as u32,
        beta: split.fractions.clone(),
        x: split.compositions.clone(),
        iterations: split.iterations,
        residual: split.residual,
        converged: split.converged,
        warnings,
    })
}
