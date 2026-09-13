//! `hydraulics.crane_k_factors` - fitting losses by the equivalent-length
//! method.
//!
//! ```text
//! K = f_t * sum(n_ld for each fitting)
//! ```
//!
//! Spec: `specs/calcs/hydraulics/crane_k_factors.yaml`
//!
//! # This calc cannot validate its own inputs, and says so
//!
//! The *method* is standard: a fitting's resistance coefficient is its
//! equivalent length ratio times the friction factor. The **coefficients** come
//! from `data/fittings/crane_k_factors.csv`, where every row is currently an
//! estimated dummy value - a placeholder of plausible magnitude, not from Crane
//! TP-410 or any other standard. See that file's header.
//!
//! The consequence is unusual and worth stating plainly: because the
//! coefficients are placeholders, there is no correct value for this calc to be
//! checked against, and **no test in this repository can detect a wrong
//! coefficient**. The tests validate the arithmetic and the registry lookup.
//! That is why every result built from estimated rows carries an
//! [`chemeng_core::WarningCode::EstimatedData`] warning: a pressure drop
//! computed from this data can be wrong by a factor of two and still look
//! entirely reasonable.

use crate::fittings::{self, Fitting};
use crate::results::{KComponent, KFactorsResult};
use crate::spec_gen;
use chemeng_core::{Result, Warning, WarningCode, apply_checks};

/// Total resistance coefficient for a list of fittings.
///
/// `f_t` is the Darcy friction factor used as the basis. Crane specifies `f_T`,
/// the fully turbulent friction factor at the nominal fitting size; using the
/// actual friction factor at the flow Reynolds number instead is a common and
/// slightly more accurate variant, and is what the `chemeng pipe` CLI does.
/// Both are defensible and they give different answers, which is why the value
/// is an explicit input rather than something this function guesses.
///
/// # Errors
/// * [`chemeng_core::ChemEngError::UnknownFitting`] for an id not in the registry. An error
///   rather than a skip: treating an unknown fitting as zero loss would
///   under-report pressure drop.
/// * [`chemeng_core::ChemEngError::OutOfRange`] if `f_t` is not positive, or if the list is
///   empty - an empty list would silently return zero loss, which reads as "no
///   fittings" when the caller may have meant to supply some.
///
/// # Example
/// ```
/// use chemeng_hydraulics::crane_k_factors;
/// use chemeng_core::CalcResult; // for has_warning
///
/// let r = crane_k_factors(&["90_elbow", "gate_valve_open"], 0.018)?;
/// assert!((r.k_total - 0.684).abs() < 1e-12);
/// assert_eq!(r.components.len(), 2);
/// // The registry holds placeholder data, so the result says so.
/// assert!(r.has_warning(chemeng_core::WarningCode::EstimatedData));
/// # Ok::<(), chemeng_core::ChemEngError>(())
/// ```
pub fn crane_k_factors(fittings: &[&str], f_t: f64) -> Result<KFactorsResult> {
    let spec = &spec_gen::CRANE_K_FACTORS_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |q| (q == "f_t").then_some(f_t),
        &mut warnings,
    )?;

    let count = fittings.len() as f64;
    apply_checks(
        spec.derived_checks(),
        |q| (q == "n_fittings").then_some(count),
        &mut warnings,
    )?;

    let mut components = Vec::with_capacity(fittings.len());
    let mut estimated = Vec::new();
    for id in fittings {
        let fitting: &'static Fitting = fittings::find(id)?;
        if fitting.is_estimated() {
            estimated.push(fitting.id.as_str());
        }
        components.push(KComponent {
            fitting_id: fitting.id.clone(),
            n_ld: fitting.n_ld,
            k: f_t * fitting.n_ld,
        });
    }

    // Provenance warning, driven by the data rather than hardcoded: promote a
    // row to `verified` and this warning stops firing for that row with no code
    // change.
    if !estimated.is_empty() {
        warnings.push(Warning::new(
            WarningCode::EstimatedData,
            format!(
                "{} of {} fitting(s) use ESTIMATED DUMMY coefficients that are not \
                 engineering data ({}). This resistance coefficient is a placeholder \
                 and must not be used to size equipment. Populate \
                 data/fittings/crane_k_factors.csv from the primary standard and set \
                 verify_status=verified.",
                estimated.len(),
                fittings.len(),
                estimated.join(", ")
            ),
        ));
    }

    let k_total = components.iter().map(|c| c.k).sum();

    Ok(KFactorsResult {
        k_total,
        f_t,
        components,
        warnings,
    })
}

/// Every fitting id in the registry, for CLI completion and diagnostics.
///
/// # Errors
/// Propagates a malformed embedded registry.
pub fn known_fittings() -> Result<Vec<&'static str>> {
    Ok(fittings::registry()?
        .iter()
        .map(|f| f.id.as_str())
        .collect())
}
