//! `eos.pr_lee_kesler_alpha` - the Peng-Robinson alpha function with a Soave-form
//! m-factor.
//!
//! ```text
//! alpha = (1 + m*(1 - sqrt(Tr)))**2
//! m     = 0.480 + 1.574*omega - 0.176*omega**2
//! ```
//!
//! Spec: `specs/calcs/eos/pr_lee_kesler_alpha.toml`, which carries the provenance and why
//! the m-factor is Soave's rather than the Lee-Kesler method the class name suggests.

use azoth_core::{Result, apply_checks};

use crate::results::PrLeeKeslerAlphaResult;
use crate::spec_gen;

/// The Peng-Robinson alpha function with a Soave-form m-factor.
///
/// `omega` is the acentric factor; `Tr` the reduced temperature. The alpha is
/// cubic-agnostic - the `a` and `b` it scales are the cubic's own.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tr <= 0` - the formula takes a square
///   root of `Tr`.
///
/// # Example
/// ```
/// use azoth_eos::pr_lee_kesler_alpha;
///
/// let r = pr_lee_kesler_alpha(0.1, 0.7)?;
/// assert!((r.alpha - 1.2184305594583278).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn pr_lee_kesler_alpha(omega: f64, Tr: f64) -> Result<PrLeeKeslerAlphaResult> {
    let spec = &spec_gen::PR_LEE_KESLER_ALPHA_SPEC;
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

    let m = 0.480 + 1.574 * omega - 0.176 * omega * omega;
    let t = 1.0 + m * (1.0 - Tr.sqrt());
    let alpha = t * t;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "alpha").then_some(alpha),
        &mut warnings,
    )?;

    Ok(PrLeeKeslerAlphaResult { alpha, warnings })
}
