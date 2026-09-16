//! `eos.pr_danesh_alpha` - the Danesh alpha function.
//!
//! ```text
//! m = 0.37464 + 1.54226*omega - 0.26992*omega**2
//! m_mod = 1.21*m above Tr = 1 else m
//! alpha = (1 + m_mod*(1 - sqrt(Tr)))**2
//! ```
//!
//! Spec: `specs/calcs/eos/pr_danesh_alpha.toml`, which carries the provenance and why the
//! 1.21 factor damps the attraction past the critical temperature.

use azoth_core::{Result, apply_checks};

use crate::results::PrDaneshAlphaResult;
use crate::spec_gen;

/// The Danesh alpha function for a pure component.
///
/// `omega` is the acentric factor; `Tr` the reduced temperature. Above `Tr = 1` the
/// Peng-Robinson Soave coefficient is multiplied by 1.21, which damps the attraction
/// into the supercritical region.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tr <= 0` - the formula takes a square
///   root of `Tr`.
///
/// # Example
/// ```
/// use azoth_eos::pr_danesh_alpha;
///
/// let r = pr_danesh_alpha(0.152, 0.7)?;
/// assert!((r.alpha - 1.2066270990034567).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn pr_danesh_alpha(omega: f64, Tr: f64) -> Result<PrDaneshAlphaResult> {
    let spec = &spec_gen::PR_DANESH_ALPHA_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "omega" => Some(omega),
            "Tr" => Some(Tr),
            _ => None,
        },
        &mut warnings,
    )?;

    let m = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
    let m_mod = if Tr > 1.0 { 1.21 * m } else { m };
    let alpha = (1.0 + m_mod * (1.0 - Tr.sqrt())) * (1.0 + m_mod * (1.0 - Tr.sqrt()));

    apply_checks(
        spec.derived_checks(),
        |name| (name == "alpha").then_some(alpha),
        &mut warnings,
    )?;

    Ok(PrDaneshAlphaResult { alpha, warnings })
}
