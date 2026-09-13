//! `eos.pr_departure` - the Peng-Robinson fugacity coefficient and departures.
//!
//! ```text
//! psi      = -kappa*sqrt(Tr) / (1 + kappa*(1 - sqrt(Tr)))
//! I        = ln((z + (1 + sqrt(2))*B) / (z + (1 - sqrt(2))*B))
//! ln_phi   = z - 1 - ln(z - B) - C*I
//! h_dep_rt = (z - 1) + C*(psi - 1)*I
//! s_dep_r  = ln(z - B) + C*psi*I          where C = A / (2*sqrt(2)*B)
//! ```
//!
//! Peng, D. Y.; Robinson, D. B. (1976). "A New Two-Constant Equation of State."
//! Ind. Eng. Chem. Fundam. 15(1), 59-64. DOI 10.1021/i160057a011
//!
//! Spec: `specs/calcs/eos/pr_departure.yaml`
//!
//! # The Gibbs identity
//!
//! `h_dep_rt - s_dep_r` equals `ln_phi`, exactly, in real arithmetic - the two
//! `psi` terms cancel:
//!
//! ```text
//! (z - 1) + C*(psi - 1)*I - ln(z - B) - C*psi*I
//!   = (z - 1) - ln(z - B) - C*I
//!   = ln_phi
//! ```
//!
//! For a pure component that is not a coincidence: the departure Gibbs energy
//! divided by `RT` **is** the logarithm of the fugacity coefficient, because
//! `G = H - TS`. So the identity is a consequence of what the three functions mean,
//! and any disagreement between them is rounding rather than a residual to be
//! tolerated.
//!
//! That is what makes it the strongest test here. A sign error in either departure
//! function leaves `ln_phi` untouched and moves the other to a value that is still
//! entirely plausible as an enthalpy or entropy departure - the point-value cases
//! would not necessarily notice. The identity does.
//!
//! # What this calc cannot check
//!
//! `z` must be a root of the cubic that `A` and `B` define, and `kappa` and `Tr`
//! must describe the same state. Neither is checkable here - the cubic is solved by
//! [`crate::pr_z_factor`], and the state is the caller's - so a mismatched input set
//! produces departure functions that are internally consistent and describe a state
//! that does not exist. The spec records both as assumptions.
//!
//! The failure is quieter for `kappa` and `Tr` than for `z`: they enter only through
//! `psi`, so a stale coefficient shifts the enthalpy and entropy departures without
//! touching `ln_phi` at all.

use azoth_core::{Result, apply_checks};

use crate::results::PrDepartureResult;
use crate::spec_gen;

/// The Peng-Robinson fugacity coefficient and departure functions, for one state.
///
/// All five arguments are the dimensionless quantities this namespace works in, and
/// all three outputs are dimensionless for the same reason: `h_dep_rt` is the
/// departure enthalpy over `R*T` and `s_dep_r` the departure entropy over `R`, so
/// the multiplication by `R` and `T` happens where those live - the model layer -
/// rather than here.
///
/// `kappa` may come from [`crate::pr_kappa`] or [`crate::prsv_kappa`]: this calc uses
/// it only through `psi`, the logarithmic derivative of the alpha function, and both
/// correlations supply a coefficient for the same alpha function.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `b_reduced <= 0` (it is a divisor) or
///   if `z <= b_reduced` (which makes `ln(z - B)` the logarithm of a negative
///   number).
///
/// # Example
/// ```
/// use azoth_eos::{pr_alpha_ab, pr_departure, pr_kappa, pr_z_factor};
///
/// let kappa = pr_kappa(0.152)?.kappa;
/// let ab = pr_alpha_ab(kappa, 0.8, 0.25)?;
/// let z = pr_z_factor(ab.a_reduced, ab.b_reduced)?;
/// let d = pr_departure(ab.a_reduced, ab.b_reduced, z.z_max, kappa, 0.8)?;
/// assert!((d.ln_phi + 0.19131055684257678).abs() < 1e-12);
/// // The identity, which is what ties the three together.
/// assert!((d.h_dep_rt - d.s_dep_r - d.ln_phi).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn pr_departure(
    a_reduced: f64,
    b_reduced: f64,
    z: f64,
    kappa: f64,
    Tr: f64,
) -> Result<PrDepartureResult> {
    let spec = &spec_gen::PR_DEPARTURE_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "a_reduced" => Some(a_reduced),
            "b_reduced" => Some(b_reduced),
            "z" => Some(z),
            "kappa" => Some(kappa),
            "Tr" => Some(Tr),
            _ => None,
        },
        &mut warnings,
    )?;

    // The logarithmic derivative of the alpha function. Hoisted rather than
    // recomputed for each use so the two languages cannot evaluate it twice in
    // different orders.
    let sqrt_tr = Tr.sqrt();
    let psi = -kappa * sqrt_tr / (1.0 + kappa * (1.0 - sqrt_tr));

    let i_term = ((z + (1.0 + std::f64::consts::SQRT_2) * b_reduced)
        / (z + (1.0 - std::f64::consts::SQRT_2) * b_reduced))
        .ln();
    let coefficient = a_reduced / (2.0 * std::f64::consts::SQRT_2 * b_reduced);
    let ln_z_minus_b = (z - b_reduced).ln();

    let ln_phi = z - 1.0 - ln_z_minus_b - coefficient * i_term;
    let h_dep_rt = (z - 1.0) + coefficient * (psi - 1.0) * i_term;
    let s_dep_r = ln_z_minus_b + coefficient * psi * i_term;

    apply_checks(
        spec.derived_checks(),
        |quantity| match quantity {
            // The bound that matters is on the difference, not on `z`: `z = B` is
            // the zero-volume limit and nothing about `z` alone says where it is.
            "z_minus_b_reduced" => Some(z - b_reduced),
            _ => None,
        },
        &mut warnings,
    )?;

    Ok(PrDepartureResult {
        ln_phi,
        h_dep_rt,
        s_dep_r,
        warnings,
    })
}
