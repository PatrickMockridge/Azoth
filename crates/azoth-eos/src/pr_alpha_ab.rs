//! `eos.pr_alpha_ab` - the Peng-Robinson alpha function and reduced parameters.
//!
//! ```text
//! alpha     = (1 + kappa*(1 - Tr**0.5))**2
//! a_reduced = Omega_a * alpha * Pr / Tr**2
//! b_reduced = Omega_b * Pr / Tr
//! ```
//!
//! Spec: `specs/calcs/eos/pr_alpha_ab.toml`, which carries the provenance, the
//! derivation of the two Omega constants from the triple-root condition, and what is
//! not claimed about the source.

use azoth_core::{Result, apply_checks};

use crate::cubic::Cubic;
use crate::results::PrAlphaAbResult;
use crate::spec_gen;

/// `Omega_a`, the attraction constant of the Peng-Robinson cubic.
///
/// **NeqSim 3.20.0's value, not the paper's, and carrying it is the point.** This
/// library is a port, and NeqSim's `ComponentPR` constructor sets
/// `a = .45724333333 * R**2 * Tc**2 / Pc`. The Peng-Robinson paper prints `0.45724`
/// and the cubic's triple-root condition gives `0.4572355289213822` exactly;
/// NeqSim's is neither. It is the paper's printed value plus 3.3333e-6 - and so is
/// its `Omega_b`, by the same offset - which is what makes the pair read as a
/// transcription artefact carried forward rather than as a refit.
///
/// Substituting NeqSim's pair for the exact one reproduces NeqSim's own `TPflash` to
/// twelve significant figures, where the exact pair leaves a 1.4e-4 residue - see
/// `validation/eos/methane_butane_flash_against_neqsim.json`, which records both the
/// measurement and how to repeat it. A port that corrected its upstream would
/// disagree with it by 1.4e-4 forever, and could never be validated against it at
/// all. The spec's `assumptions` carries the argument in full.
pub const OMEGA_A: f64 = Cubic::Pr.omega_a();

/// `Omega_b`, the repulsion constant of the Peng-Robinson cubic.
///
/// NeqSim's value, for the reason [`OMEGA_A`] gives: `ComponentPR` sets
/// `b = .077803333 * R * Tc / Pc`, against the cubic's exact `0.07779607390388846`
/// and the paper's printed `0.07780`. It is 7.26e-6 above the exact value and 3.3333e-6
/// above the printed one - the same offset `Omega_a` carries, and the larger of the
/// two departures from the exact pair.
pub const OMEGA_B: f64 = Cubic::Pr.omega_b();

/// The Peng-Robinson alpha function and the reduced attraction parameters.
///
/// `kappa` comes from [`crate::pr_kappa`]; `Tr` and `Pr` are reduced against the
/// caller's critical point. The library does ship a component databank
/// ([`crate::databank`]) - this kernel takes reduced variables, so a caller holding
/// a name reduces there first.
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
/// assert!((r.a_reduced - 0.20206845072985785).abs() < 1e-15);
/// assert!((r.b_reduced - 0.0243135415625).abs() < 1e-15);
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
