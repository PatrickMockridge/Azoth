//! `process.splitter` - one feed divided into branches.
//!
//! Spec: `specs/models/process/splitter.yaml`, which carries the provenance, the reason
//! there is no ideal-gas block, and the divergence over normalising the fractions.
//!
//! Every branch leaves at the feed's own temperature, pressure and composition, so one
//! isothermal flash describes all of them. The model runs it **once** and reports one
//! `phase` and one `beta`, with the branch flows as a vector - not a shortcut, but the
//! same answer stated in a form that cannot disagree with itself.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result, apply_checks};
use azoth_eos::mixture::Mixture;
use azoth_eos::pt_flash::pt_flash;

use crate::model_gen;
use crate::results::SplitterResult;

/// One feed divided into branches at the feed's own temperature and pressure.
///
/// `fractions` is the share of the feed each branch takes. They must **sum to one**, and
/// that is checked rather than corrected: silently rescaling a caller's fractions would
/// make their error invisible while changing every number downstream.
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
