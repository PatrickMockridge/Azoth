//! `eos.nrtl_activity_coefficients` - the activity coefficients from the NRTL
//! local-composition model.
//!
//! Spec: `specs/models/eos/nrtl_activity_coefficients.toml`. A *direct* model: no
//! iteration, so no `algorithm` block.

use azoth_core::{AzothError, Result, apply_checks};

use crate::model_gen;
use crate::results::NrtlActivityCoefficientsResult;

/// The activity coefficients of a mixture, from NRTL (Renon-Prausnitz).
///
/// `Dij[i][j] = g_ij` is in Kelvin so that `tau_ij = Dij[i][j] / T` is dimensionless,
/// `G_ij = exp(-alpha[i][j] * tau_ij)`, and
/// `ln gamma_i = (sum_j tau_ji G_ji x_j) / (sum_j G_ji x_j)
/// + sum_j (x_j G_ij / C_j) (tau_ij - D_j / C_j)` where `C_j = sum_l G_lj x_l` and
/// `D_j = sum_l tau_lj G_lj x_l`. `alpha` is symmetric and both matrices have a zero
/// diagonal; `x` is checked rather than renormalised.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the matrices are not `N x N` for the `N` of `x`,
///   `alpha` is not symmetric, a diagonal is not zero, or `x` is not a composition.
/// * [`AzothError::OutOfRange`] if `T` is not positive.
///
/// # Example
/// ```
/// use azoth_eos::nrtl_activity_coefficients;
///
/// let r = nrtl_activity_coefficients(
///     350.0,
///     &[0.5, 0.5],
///     &[vec![0.0, -48.68], vec![610.6, 0.0]],
///     &[vec![0.0, 0.303], vec![0.303, 0.0]],
/// )?;
/// assert!((r.gamma[0] - 1.2277269290049744).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Dij`, `T` and `x` are the symbols in the chemistry
pub fn nrtl_activity_coefficients(
    T: f64,
    x: &[f64],
    Dij: &[Vec<f64>],
    alpha: &[Vec<f64>],
) -> Result<NrtlActivityCoefficientsResult> {
    let spec = &model_gen::NRTL_ACTIVITY_COEFFICIENTS_SPEC;
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
    if Dij.len() != n || alpha.len() != n {
        return Err(AzothError::invalid_input(
            "Dij",
            format!(
                "the matrices are {} x {} (Dij) and {} x {} (alpha) but `x` has {n} \
                 entries; a matrix input must be N x N for the same N as `x`",
                Dij.len(),
                Dij.first().map_or(0, Vec::len),
                alpha.len(),
                alpha.first().map_or(0, Vec::len)
            ),
        ));
    }
    for i in 0..n {
        if Dij[i].len() != n || alpha[i].len() != n {
            return Err(AzothError::invalid_input(
                "Dij",
                format!(
                    "row {i} of a matrix is {} (Dij) / {} (alpha) entries but must be {n}",
                    Dij[i].len(),
                    alpha[i].len()
                ),
            ));
        }
    }
    for i in 0..n {
        for j in 0..n {
            if i == j {
                if Dij[i][j] != 0.0 || alpha[i][j] != 0.0 {
                    return Err(AzothError::invalid_input(
                        "Dij",
                        format!(
                            "the diagonal must be zero, but Dij[{i}][{i}] = {} and \
                             alpha[{i}][{i}] = {}",
                            Dij[i][j], alpha[i][j]
                        ),
                    ));
                }
            } else if alpha[i][j] != alpha[j][i] {
                return Err(AzothError::invalid_input(
                    "alpha",
                    format!(
                        "`alpha` must be symmetric, but alpha[{i}][{j}] = {} and \
                         alpha[{j}][{i}] = {}",
                        alpha[i][j], alpha[j][i]
                    ),
                ));
            }
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

    let mut tau = vec![vec![0.0; n]; n];
    let mut g = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            tau[i][j] = Dij[i][j] / T;
            g[i][j] = (-alpha[i][j] * tau[i][j]).exp();
        }
    }

    let mut ln_gamma = vec![0.0; n];
    let mut gamma = vec![0.0; n];
    for i in 0..n {
        let mut a = 0.0;
        let mut b = 0.0;
        for j in 0..n {
            a += tau[j][i] * g[j][i] * x[j];
            b += g[j][i] * x[j];
        }
        let mut f = 0.0;
        for j in 0..n {
            let mut c = 0.0;
            let mut d = 0.0;
            for l in 0..n {
                c += g[l][j] * x[l];
                d += tau[l][j] * g[l][j] * x[l];
            }
            f += (x[j] * g[i][j] / c) * (tau[i][j] - d / c);
        }
        let lng = a / b + f;
        ln_gamma[i] = lng;
        gamma[i] = lng.exp();
    }

    Ok(NrtlActivityCoefficientsResult {
        ln_gamma,
        gamma,
        warnings,
    })
}
