//! `eos.pr_alpha_ab` - the Peng-Robinson alpha function and reduced parameters.
//!
//! ```text
//! alpha     = (1 + kappa*(1 - Tr**0.5))**2
//! a_reduced = Omega_a * alpha * Pr / Tr**2
//! b_reduced = Omega_b * Pr / Tr
//! ```
//!
//! Peng, D. Y.; Robinson, D. B. (1976). "A New Two-Constant Equation of State."
//! Ind. Eng. Chem. Fundam. 15(1), 59-64. DOI 10.1021/i160057a011
//!
//! Spec: `specs/calcs/eos/pr_alpha_ab.yaml`
//!
//! # Why the Omega constants are not the printed ones
//!
//! The paper prints `Omega_a = 0.45724` and `Omega_b = 0.07780`. Those are
//! roundings, and using them makes the cubic fail its own critical point: at
//! `Tr = Pr = 1` the equation's real root comes out 0.321379025174 instead of
//! Peng-Robinson's critical compressibility 0.307401308699 - an error of 4.55%,
//! exactly where the equation is supposed to be anchored.
//!
//! The cause is conditioning. The constants exist to place a *triple* root at the
//! critical point, and a triple root is cubically ill-conditioned: perturbing the
//! coefficients by `epsilon` moves the roots by about `epsilon**(1/3)`. A
//! rounding error of 5e-6 in `Omega_b` therefore moves the critical root by
//! around 1.7e-2.
//!
//! So this module carries the full-precision pair. They are not transcribed from
//! anywhere - they are the unique solution of the triple-root condition, which is
//! the check `every_omega_is_the_triple_root_solution` in this crate's tests
//! performs. See the spec's `verification` notes for what remains unconfirmed.

use azoth_core::{Result, apply_checks};

use crate::results::PrAlphaAbResult;
use crate::spec_gen;

/// `Omega_a`, the attraction constant of the Peng-Robinson cubic.
///
/// Full precision, and load-bearing - see the module documentation. The paper
/// prints `0.45724`, which is this rounded, and which puts the cubic's critical
/// point 4.55% out.
pub const OMEGA_A: f64 = 0.4572355289213822;

/// `Omega_b`, the repulsion constant of the Peng-Robinson cubic.
///
/// Full precision. The paper prints `0.07780`.
pub const OMEGA_B: f64 = 0.07779607390388846;

/// The Peng-Robinson alpha function and the reduced attraction parameters.
///
/// `kappa` comes from [`crate::pr_kappa`]; `Tr` and `Pr` are reduced against the
/// caller's critical point, because this library ships no component databank.
///
/// All three inputs and all three outputs are dimensionless, so nothing here
/// touches units - which is the point of writing the cubic in `A` and `B` rather
/// than in `a` and `b`.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tr <= 0` (a square root and a
///   squared divisor) or if `Pr <= 0` (`B` appears as a divisor in the fugacity
///   expression this feeds).
///
/// # Example
/// ```
/// use azoth_eos::{pr_alpha_ab, pr_kappa};
///
/// let kappa = pr_kappa(0.152)?.kappa;
/// let r = pr_alpha_ab(kappa, 0.8, 0.25)?;
/// assert!((r.alpha - 1.1313346661636197).abs() < 1e-15);
/// assert!((r.a_reduced - 0.20206500174625697).abs() < 1e-15);
/// assert!((r.b_reduced - 0.02431127309496514).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` and `Pr` are the symbols in the published equation
pub fn pr_alpha_ab(kappa: f64, Tr: f64, Pr: f64) -> Result<PrAlphaAbResult> {
    let spec = &spec_gen::PR_ALPHA_AB_SPEC;
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
    //
    // Written as `t * t` rather than `.powi(2)` to mirror the Python side's
    // `t ** 2`, following the same convention `darcy_weisbach` uses for
    // `v**2`/`v * v`.
    //
    // Measured: the two backends produce bit-identical results for all three
    // outputs on the worked example, as they do for `darcy_weisbach`. That is not
    // *guaranteed* - IEEE-754 pins `+ - * /` and `sqrt`, and does not pin what a
    // language's `x ** 2` lowering does - so this calc claims agreement within the
    // spec's tolerance and not bit-equality. The measurement is recorded because a
    // claim that is true but unstated is worth less than one that is checked.
    let attraction = 1.0 + kappa * (1.0 - Tr.sqrt());
    let alpha = attraction * attraction;
    // Guarded by the checks above, so Tr is positive here.
    let a_reduced = OMEGA_A * alpha * Pr / (Tr * Tr);
    let b_reduced = OMEGA_B * Pr / Tr;

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

    Ok(PrAlphaAbResult {
        alpha,
        a_reduced,
        b_reduced,
        warnings,
    })
}
