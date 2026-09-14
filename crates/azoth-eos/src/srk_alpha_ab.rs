//! `eos.srk_alpha_ab` - the Soave-Redlich-Kwong alpha function and reduced parameters.
//!
//! ```text
//! alpha     = (1 + kappa*(1 - Tr**0.5))**2
//! a_reduced = Omega_a * alpha * Pr / Tr**2
//! b_reduced = Omega_b * Pr / Tr
//! ```
//!
//! Spec: `specs/calcs/eos/srk_alpha_ab.toml`, which carries the provenance and the two
//! Omega constants' derivation from the cubic's triple-root condition.
//!
//! The alpha function is Soave's, the same shape as Peng-Robinson's with a different
//! coefficient; the Omega pair is Redlich-Kwong's original. Only the constants differ
//! from [`crate::pr_alpha_ab`].

use azoth_core::{Result, apply_checks};

use crate::cubic::Cubic;
use crate::results::SrkAlphaAbResult;
use crate::spec_gen;

/// The Soave-Redlich-Kwong alpha function and the reduced attraction parameters.
///
/// `kappa` comes from [`crate::srk_kappa`]; `Tr` and `Pr` are reduced against the
/// caller's critical point. All inputs and outputs are dimensionless.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tr <= 0` or `Pr <= 0`.
///
/// # Example
/// ```
/// use azoth_eos::{srk_alpha_ab, srk_kappa};
///
/// let kappa = srk_kappa(0.152)?.kappa;
/// let r = srk_alpha_ab(kappa, 0.8, 0.25)?;
/// assert!((r.alpha - 1.1567082960277375).abs() < 1e-12);
/// assert!((r.a_reduced - 0.19315231739218255).abs() < 1e-12);
/// assert!((r.b_reduced - 0.02707510936404929).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` and `Pr` are the symbols in the published equation
pub fn srk_alpha_ab(kappa: f64, Tr: f64, Pr: f64) -> Result<SrkAlphaAbResult> {
    let spec = &spec_gen::SRK_ALPHA_AB_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "kappa" => Some(kappa),
            "Tr" => Some(Tr),
            "Pr" => Some(Pr),
            _ => None,
        },
        &mut warnings,
    )?;

    // A square, so never negative; exactly 1 at Tr = 1 whatever kappa is, which
    // is what makes the critical point special.
    let attraction = 1.0 + kappa * (1.0 - Tr.sqrt());
    let alpha = attraction * attraction;
    // Guarded by the checks above, so Tr is positive here.
    let a_reduced = Cubic::Srk.omega_a() * alpha * Pr / (Tr * Tr);
    let b_reduced = Cubic::Srk.omega_b() * Pr / Tr;

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

    Ok(SrkAlphaAbResult {
        alpha,
        a_reduced,
        b_reduced,
        warnings,
    })
}
