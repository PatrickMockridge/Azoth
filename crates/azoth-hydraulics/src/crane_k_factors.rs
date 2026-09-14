//! `hydraulics.crane_k_factors` - fitting losses by the equivalent-length method.
//!
//! ```text
//! K = f_t * sum(n_ld for each fitting)
//! ```
//!
//! Spec: `specs/calcs/hydraulics/crane_k_factors.yaml`, which records that the
//! coefficients this sums are placeholders and what that means for testing.
//!
//! The method is standard; the coefficients come from
//! `data/fittings/crane_k_factors.csv`, and every row there is an estimated dummy value
//! rather than a value from Crane TP-410 or any other standard.

use crate::fittings::{self, Fitting};
use crate::results::{KComponent, KFactorsResult};
use crate::spec_gen;
use azoth_core::{Result, apply_checks};

/// Total resistance coefficient for a list of fittings.
///
/// `f_t` is the Darcy friction factor used as the basis. Crane specifies `f_T`,
/// the fully turbulent friction factor at the nominal fitting size; using the
/// actual friction factor at the flow Reynolds number instead is a common and
/// slightly more accurate variant, and is what the `azoth pipe` CLI does.
/// Both are defensible and they give different answers, which is why the value
/// is an explicit input rather than something this function guesses.
///
/// # Errors
/// * [`azoth_core::AzothError::UnknownFitting`] for an id not in the registry. An error
///   rather than a skip: treating an unknown fitting as zero loss would
///   under-report pressure drop.
/// * [`azoth_core::AzothError::OutOfRange`] if `f_t` is not positive, or if the list is
///   empty - an empty list would silently return zero loss, which reads as "no
///   fittings" when the caller may have meant to supply some.
///
/// # Example
/// ```
/// use azoth_hydraulics::crane_k_factors;
///
/// let r = crane_k_factors(&["90_elbow", "gate_valve_open"], 0.018)?;
/// assert!((r.k_total - 0.684).abs() < 1e-12);
/// assert_eq!(r.components.len(), 2);
/// # Ok::<(), azoth_core::AzothError>(())
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
    for id in fittings {
        let fitting: &'static Fitting = fittings::find(id)?;
        components.push(KComponent {
            fitting_id: fitting.id.clone(),
            n_ld: fitting.n_ld,
            k: f_t * fitting.n_ld,
        });
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
