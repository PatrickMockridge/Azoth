//! `process.splitter` - one feed divided into branches.
//!
//! Spec: `specs/models/process/splitter.yaml`
//!
//! # The port
//!
//! NeqSim's `Splitter.run` is lines 376-438 of an 827-line file, and its physics:
//!
//! ```text
//! n_k = f_k * n_in                   (Splitter.java:420-425, addComponent per component)
//! T_k = T_in, P_k = P_in             (copied by the clone at :407, never set)
//! state_k = TPflash(T_k, P_k)        (:427)
//! ```
//!
//! The rest of the file is the split-factor bookkeeping (`:385-404`, which converts a
//! fixed-flow specification into fractions and normalises them) and the transient mixer
//! path.
//!
//! # One flash, not `S` flashes
//!
//! NeqSim flashes every branch (`:427`), and it has to: its streams are mutable objects
//! that a caller may have written to between the split and the flash. Here every branch
//! is at the **same** temperature, pressure and composition by construction, and an
//! isothermal flash is a function of exactly those three - so `S` flashes would return
//! `S` copies of one answer.
//!
//! This model runs the flash **once** and reports one `phase` and one `beta`, with the
//! branch flows as a vector. That is not a shortcut: it is the same answer, stated in a
//! form that cannot disagree with itself. A splitter whose branches were at different
//! states would be `S` separators, not one splitter.
//!
//! # No ideal-gas model
//!
//! This is the one unit operation whose spec carries no `cp_a`..`P_ref` block, and the
//! reason is that it needs none: a splitter does no energy balance. Its only
//! thermodynamic call is a `TPflash`, which is a function of the temperature, the
//! pressure and the composition. Carrying an ideal-gas datum nothing reads would be a
//! field a caller has to supply and no test could justify.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result, apply_checks};
use azoth_eos::mixture::Mixture;
use azoth_eos::pt_flash::pt_flash;

use crate::model_gen;
use crate::results::SplitterResult;

/// One feed divided into branches at the feed's own temperature and pressure.
///
/// `fractions` is the share of the feed each branch takes. They must **sum to one**,
/// and that is checked rather than corrected - the rule everywhere else in this
/// library, and the reason is the same here: silently rescaling a caller's fractions
/// would make their error invisible while changing every number downstream.
///
/// NeqSim normalises instead (`Splitter.java:388-404`), and its own guard is worth
/// reading: it zeroes negatives, and if the total is not positive it sets every factor
/// to zero and the first to one. Both are defensible for a transient solver that must
/// keep running; neither is defensible for a model whose whole output is those numbers.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if `fractions` is empty, holds a negative
///   entry, or does not sum to one to within `1e-9`.
/// * [`azoth_core::AzothError::OutOfRange`] if the feed's state is non-positive, or the
///   flow is not positive.
/// * Whatever the flash raises.
pub fn splitter(
    mixture: &Mixture,
    t_in: ThermodynamicTemperature,
    p_in: Pressure,
    n_in: f64,
    z: &[f64],
    fractions: &[f64],
) -> Result<SplitterResult> {
    let spec = &model_gen::SPLITTER_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t_in.value),
            "P" => Some(p_in.value),
            "n" => Some(n_in),
            _ => None,
        },
        &mut warnings,
    )?;

    if fractions.is_empty() {
        return Err(AzothError::InvalidInput {
            field: "fractions".to_string(),
            reason: "a splitter with no branches has nothing to split into".to_string(),
        });
    }
    if let Some(negative) = fractions.iter().find(|f| **f < 0.0) {
        return Err(AzothError::InvalidInput {
            field: "fractions".to_string(),
            reason: format!(
                "a negative share ({negative}) asks a branch to give back more than it \
                 takes, which is a mixer written backwards"
            ),
        });
    }
    let total: f64 = fractions.iter().sum();
    if (total - 1.0).abs() > 1.0e-9 {
        return Err(AzothError::InvalidInput {
            field: "fractions".to_string(),
            reason: format!(
                "the shares must account for the whole feed and sum to {total}, not 1. \
                 A splitter that loses or invents material is not a splitter, and \
                 normalising here would hide the arithmetic error that produced this"
            ),
        });
    }

    // One flash: every branch is at the same state by construction.
    let flash = pt_flash(mixture, t_in, p_in, z)?;
    warnings.extend(flash.warnings);

    Ok(SplitterResult {
        temperature: t_in,
        pressure: p_in,
        phase: flash.phase,
        beta: flash.beta,
        flows: fractions.iter().map(|f| f * n_in).collect(),
        iterations: flash.iterations,
        warnings,
    })
}
