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
//! Spec: `specs/calcs/eos/pr_departure.toml`, which carries the provenance, the
//! derivation of the Gibbs identity `h_dep_rt - s_dep_r = ln_phi`, and the assumptions
//! about `z`, `kappa` and `Tr` this calc cannot check.

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
/// assert!((d.ln_phi + 0.19131116744197554).abs() < 1e-12);
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

    // The heat-capacity departure, from `Cp^R/R = y + T*(dy/dT)_P` with `y = h_dep_rt`.
    //
    // Every derivative below is already multiplied by `T` - that is what the `t_`
    // prefix means - so the `1/T` that each of them carries cancels and no term here
    // needs the absolute temperature, only `Tr`. That is what lets this be an output
    // of the same call as the three above, and why it declares no input they do not.
    let t_da = a_reduced * (psi - 2.0);
    let t_db = -b_reduced;
    let t_dpsi =
        -kappa * (1.0 + kappa) * Tr / (2.0 * sqrt_tr * (1.0 + kappa * (1.0 - sqrt_tr)).powi(2));
    let t_dc = coefficient * (psi - 1.0);

    // `z` is a root of `F(z, T) = 0`, so `dz/dT = -(dF/dT)/(dF/dz)` and the chain
    // rule carries the two reduced-parameter derivatives through both brackets.
    let d_f_dz = 3.0 * z * z
        + 2.0 * (b_reduced - 1.0) * z
        + (a_reduced - 3.0 * b_reduced * b_reduced - 2.0 * b_reduced);
    let t_dfdt = t_db * z * z
        + (t_da - 6.0 * b_reduced * t_db - 2.0 * t_db) * z
        + (3.0 * b_reduced * b_reduced * t_db + 2.0 * b_reduced * t_db
            - t_da * b_reduced
            - a_reduced * t_db);
    let t_dz = -t_dfdt / d_f_dz;

    let n_plus = z + (1.0 + std::f64::consts::SQRT_2) * b_reduced;
    let n_minus = z + (1.0 - std::f64::consts::SQRT_2) * b_reduced;
    let t_di = (t_dz + (1.0 + std::f64::consts::SQRT_2) * t_db) / n_plus
        - (t_dz + (1.0 - std::f64::consts::SQRT_2) * t_db) / n_minus;

    let cp_dep_r = h_dep_rt
        + t_dz
        + t_dc * (psi - 1.0) * i_term
        + coefficient * t_dpsi * i_term
        + coefficient * (psi - 1.0) * t_di;

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
        cp_dep_r,
        warnings,
    })
}
