//! `eos.matcop_alpha` - the Mathias-Copeman alpha function.
//!
//! ```text
//! alpha = (1 + mc1*(1 - sqrt(Tr)) + mc2*(1 - sqrt(Tr))**2 + mc3*(1 - sqrt(Tr))**3)**2
//! ```
//!
//! Spec: `specs/calcs/eos/matcop_alpha.toml`, which carries the provenance and why the
//! three coefficients are the caller's rather than derived.

use azoth_core::{Result, apply_checks};

use crate::results::MatcopAlphaResult;
use crate::spec_gen;

/// The Mathias-Copeman alpha function for a pure component.
///
/// `mc1`, `mc2` and `mc3` are the fitted Mathias-Copeman parameters; `Tr` the reduced
/// temperature. The alpha is cubic-agnostic - the `a` and `b` it scales are the cubic's
/// own, so one correlation serves Peng-Robinson and Soave-Redlich-Kwong alike.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tr <= 0` - the formula takes a square
///   root of `Tr`.
///
/// # Example
/// ```
/// use azoth_eos::matcop_alpha;
///
/// let r = matcop_alpha(0.1, 0.05, 0.02, 0.7)?;
/// assert!((r.alpha - 1.0358255509077803).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn matcop_alpha(mc1: f64, mc2: f64, mc3: f64, Tr: f64) -> Result<MatcopAlphaResult> {
    let spec = &spec_gen::MATCOP_ALPHA_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "mc1" => Some(mc1),
            "mc2" => Some(mc2),
            "mc3" => Some(mc3),
            "Tr" => Some(Tr),
            _ => None,
        },
        &mut warnings,
    )?;

    let one_minus = 1.0 - Tr.sqrt();
    let inner = 1.0
        + mc1 * one_minus
        + mc2 * one_minus * one_minus
        + mc3 * one_minus * one_minus * one_minus;
    let alpha = inner * inner;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "alpha").then_some(alpha),
        &mut warnings,
    )?;

    Ok(MatcopAlphaResult { alpha, warnings })
}
