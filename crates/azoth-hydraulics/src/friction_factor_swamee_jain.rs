//! `hydraulics.friction_factor_swamee_jain` - the explicit Swamee-Jain friction
//! factor.
//!
//! ```text
//! f = 0.25 / (log10(relative_roughness / 3.7 + 5.74 / Re**0.9))**2
//! ```
//!
//! Spec: `specs/calcs/hydraulics/friction_factor_swamee_jain.yaml`, which carries the
//! citation.
//!
//! This is an *approximation to* Colebrook, not a more correct alternative to
//! it. It is offered so a caller can trade about 1% accuracy for the absence of
//! an iteration, and the wording matters: presenting the two as equivalent
//! methods would be wrong. A test compares them against the paper's own stated
//! accuracy, so a regression fails against the published claim rather than
//! against a number chosen here.

use crate::results::SwameeJainResult;
use crate::spec_gen;
use azoth_core::{Result, apply_checks};

/// Swamee-Jain explicit friction factor. `relative_roughness` is `epsilon / D`.
///
/// Returns immediately: no iteration, and so no convergence report. Where the
/// Colebrook result carries `iterations` and `converged`, this carries only `f`.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Re <= 0` (`Re**0.9` is zero
///   at the origin, making the `5.74/Re**0.9` term singular) or if
///   `relative_roughness < 0`.
///
/// Outside the paper's stated range the value is still returned, carrying an
/// `OutOfValidRange` warning: the equation is well defined there, it simply is
/// not covered by the accuracy claim.
///
/// # Example
/// ```
/// use azoth_hydraulics::friction_factor_swamee_jain;
///
/// let r = friction_factor_swamee_jain(100_000.0, 4.6e-4)?;
/// assert!((r.f - 0.02024003).abs() < 1e-6);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn friction_factor_swamee_jain(re: f64, relative_roughness: f64) -> Result<SwameeJainResult> {
    let spec = &spec_gen::FRICTION_FACTOR_SWAMEE_JAIN_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "re" => Some(re),
            "relative_roughness" => Some(relative_roughness),
            _ => None,
        },
        &mut warnings,
    )?;

    let inner = (relative_roughness / 3.7 + 5.74 / re.powf(0.9)).log10();
    let f = 0.25 / (inner * inner);

    apply_checks(
        spec.derived_checks(),
        |q| (q == "f").then_some(f),
        &mut warnings,
    )?;

    Ok(SwameeJainResult { f, warnings })
}
