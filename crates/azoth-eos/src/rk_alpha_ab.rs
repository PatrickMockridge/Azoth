//! `eos.rk_alpha_ab` - the Redlich-Kwong alpha function and reduced parameters.
//!
//! ```text
//! alpha     = 1/Tr**0.5
//! a_reduced = Omega_a * alpha * Pr / Tr**2
//! b_reduced = Omega_b * Pr / Tr
//! ```
//!
//! Spec: `specs/calcs/eos/rk_alpha_ab.toml`. The original Redlich-Kwong temperature
//! dependence, kappa-free - which is the whole difference from
//! [`crate::srk_alpha_ab`], whose Soave correlation adds a coefficient to this shape.

use azoth_core::{Result, apply_checks};

use crate::cubic::Cubic;
use crate::results::RkAlphaAbResult;
use crate::spec_gen;

/// The Redlich-Kwong alpha function and the reduced attraction parameters.
///
/// `Tr` and `Pr` are reduced against the caller's critical point. There is no
/// `kappa`: RK's alpha is `1/sqrt(Tr)`, with no fitted coefficient.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tr <= 0` or `Pr <= 0`.
///
/// # Example
/// ```
/// use azoth_eos::rk_alpha_ab;
///
/// let r = rk_alpha_ab(0.8, 0.25)?;
/// assert!((r.alpha - 1.118033988749895).abs() < 1e-12);
/// assert!((r.a_reduced - 0.1866943088347048).abs() < 1e-12);
/// assert!((r.b_reduced - 0.02707510936404929).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` and `Pr` are the symbols in the published equation
pub fn rk_alpha_ab(Tr: f64, Pr: f64) -> Result<RkAlphaAbResult> {
    let spec = &spec_gen::RK_ALPHA_AB_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "Tr" => Some(Tr),
            "Pr" => Some(Pr),
            _ => None,
        },
        &mut warnings,
    )?;

    let alpha = 1.0 / Tr.sqrt();
    // Guarded by the checks above, so Tr is positive here.
    let a_reduced = Cubic::Rk.omega_a() * alpha * Pr / (Tr * Tr);
    let b_reduced = Cubic::Rk.omega_b() * Pr / Tr;

    apply_checks(
        spec.derived_checks(),
        |name| match name {
            "alpha" => Some(alpha),
            "a_reduced" => Some(a_reduced),
            "b_reduced" => Some(b_reduced),
            _ => None,
        },
        &mut warnings,
    )?;

    Ok(RkAlphaAbResult {
        alpha,
        a_reduced,
        b_reduced,
        warnings,
    })
}
