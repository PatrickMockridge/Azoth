//! `eos.uniquac_activity_coefficients` - the activity coefficients from the UNIQUAC
//! model.
//!
//! Spec: `specs/models/eos/uniquac_activity_coefficients.toml`. A *direct* model: no
//! iteration, so no `algorithm` block.

use azoth_core::{AzothError, Result, apply_checks};

use crate::model_gen;
use crate::results::UniquacActivityCoefficientsResult;

/// The activity coefficients of a mixture, from UNIQUAC (Abrams-Prausnitz).
///
/// `phi_i = r_i x_i / sum r_j x_j` and `theta_i = q_i x_i / sum q_j x_j` are the
/// volume and surface fractions, `l_i = 5 (r_i - q_i) - (r_i - 1)`, and
/// `tau_ij = exp(-aij[i][j] / T)`. The combinatorial term is
/// `ln gamma^C_i = ln(phi_i/x_i) + 5 q_i ln(theta_i/phi_i) + l_i - (phi_i/x_i) sum
/// x_j l_j`, and the residual is `q_i (1 - ln(sum_j theta_j tau_ji) - sum_j theta_j
/// tau_ij / sum_k theta_k tau_kj)`. `aij` is directional with a zero diagonal; `x` is
/// checked rather than renormalised.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the vectors disagree in length, `aij` is not
///   `N x N` with a zero diagonal, or `x` is not a composition.
/// * [`AzothError::OutOfRange`] if `T` is not positive.
///
/// # Example
/// ```
/// use azoth_eos::uniquac_activity_coefficients;
///
/// let r = uniquac_activity_coefficients(
///     298.15,
///     &[0.5, 0.5],
///     &[1.4311, 0.92],
///     &[1.432, 1.4],
///     &[vec![0.0, -71.0], vec![209.0, 0.0]],
/// )?;
/// assert!((r.gamma[0] - 1.2185441848728196).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T` and `x` are the symbols in the chemistry
pub fn uniquac_activity_coefficients(
    T: f64,
    x: &[f64],
    r: &[f64],
    q: &[f64],
    aij: &[Vec<f64>],
) -> Result<UniquacActivityCoefficientsResult> {
    let spec = &model_gen::UNIQUAC_ACTIVITY_COEFFICIENTS_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T),
            _ => None,
        },
        &mut warnings,
    )?;

    let n = x.len();
    if n == 0 {
        return Err(AzothError::invalid_input(
            "x",
            "a mixture of zero components has no activity coefficient",
        ));
    }
    if r.len() != n || q.len() != n {
        return Err(AzothError::invalid_input(
            "r",
            format!(
                "the per-component vectors disagree in length: `x` has {n} entries, `r` \
                 {} and `q` {}",
                r.len(),
                q.len()
            ),
        ));
    }
    if aij.len() != n {
        return Err(AzothError::invalid_input(
            "aij",
            format!("`aij` has {} rows but must be {n} x {n}", aij.len()),
        ));
    }
    for (i, row) in aij.iter().enumerate() {
        if row.len() != n {
            return Err(AzothError::invalid_input(
                "aij",
                format!("row {i} of `aij` has {} entries but must be {n}", row.len()),
            ));
        }
        if row[i] != 0.0 {
            return Err(AzothError::invalid_input(
                "aij",
                format!("the diagonal must be zero, but aij[{i}][{i}] = {}", row[i]),
            ));
        }
    }
    if let Some(bad) = x.iter().position(|&value| value < 0.0) {
        return Err(AzothError::invalid_input(
            "x",
            format!(
                "x[{bad}] is {} but a mole fraction cannot be negative",
                x[bad]
            ),
        ));
    }
    let sum: f64 = x.iter().sum();
    if (sum - 1.0).abs() > 1.0e-09 {
        return Err(AzothError::invalid_input(
            "x",
            format!(
                "the mole fractions sum to {sum}, not to one. Renormalising them here \
                 would make a composition error invisible in every number downstream, \
                 so it is refused instead"
            ),
        ));
    }

    let sum_r: f64 = x.iter().zip(r).map(|(&xi, &ri)| xi * ri).sum();
    let sum_q: f64 = x.iter().zip(q).map(|(&xi, &qi)| xi * qi).sum();
    let l: Vec<f64> = r
        .iter()
        .zip(q)
        .map(|(&ri, &qi)| 5.0 * (ri - qi) - (ri - 1.0))
        .collect();
    let sum_l: f64 = x.iter().zip(&l).map(|(&xi, &li)| xi * li).sum();
    let theta: Vec<f64> = q.iter().zip(x).map(|(&qi, &xi)| qi * xi / sum_q).collect();
    let tau: Vec<Vec<f64>> = (0..n)
        .map(|m| (0..n).map(|k| (-aij[m][k] / T).exp()).collect())
        .collect();

    let mut ln_gamma = vec![0.0; n];
    let mut gamma = vec![0.0; n];
    for i in 0..n {
        let phi = r[i] * x[i] / sum_r;
        let lng_c =
            (phi / x[i]).ln() + 5.0 * q[i] * (theta[i] / phi).ln() + l[i] - (phi / x[i]) * sum_l;

        let s1: f64 = (0..n).map(|j| theta[j] * tau[j][i]).sum();
        let mut s3 = 0.0;
        for j in 0..n {
            let denom: f64 = (0..n).map(|k| theta[k] * tau[k][j]).sum();
            s3 += theta[j] * tau[i][j] / denom;
        }
        let lng_r = q[i] * (1.0 - s1.ln() - s3);

        let lng = lng_c + lng_r;
        ln_gamma[i] = lng;
        gamma[i] = lng.exp();
    }

    Ok(UniquacActivityCoefficientsResult {
        ln_gamma,
        gamma,
        warnings,
    })
}
