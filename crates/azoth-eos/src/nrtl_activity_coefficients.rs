//! `eos.nrtl_activity_coefficients` - the activity coefficients from the NRTL
//! local-composition model.
//!
//! Spec: `specs/models/eos/nrtl_activity_coefficients.toml`. A *direct* model: no
//! iteration, so no `algorithm` block.

use azoth_core::{AzothError, Result, apply_checks};

use crate::databank::NrtlParameters;
use crate::model_gen;
use crate::results::NrtlActivityCoefficientsResult;

/// The activity coefficients of a mixture, from NRTL (Renon-Prausnitz).
///
/// `params.dij[i][j] = g_ij` is in Kelvin so that `tau_ij = params.dij[i][j] / T` is
/// dimensionless, `G_ij = exp(-params.alpha[i][j] * tau_ij)`, and
/// `ln gamma_i = (sum_j tau_ji G_ji x_j) / (sum_j G_ji x_j)
/// + sum_j (x_j G_ij / C_j) (tau_ij - D_j / C_j)` where `C_j = sum_l G_lj x_l` and
/// `D_j = sum_l tau_lj G_lj x_l`. `params` is resolved by name through
/// [`crate::databank::nrtl_parameters`], which is symmetric with a zero diagonal in
/// `alpha` by construction; `x` is checked rather than renormalised.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if either matrix is not `N x N` for the `N` of `x`,
///   or `x` is not a composition.
/// * [`AzothError::OutOfRange`] if `T` is not positive.
///
/// # Example
/// ```
/// use azoth_eos::databank::NrtlParameters;
/// use azoth_eos::nrtl_activity_coefficients;
///
/// let params = NrtlParameters {
///     alpha: vec![0.0, 0.303, 0.303, 0.0],
///     dij: vec![0.0, -48.68, 610.6, 0.0],
/// };
/// let r = nrtl_activity_coefficients(&params, 350.0, &[0.5, 0.5])?;
/// assert!((r.gamma[0] - 1.2277269290049744).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T` and `x` are the symbols in the chemistry
pub fn nrtl_activity_coefficients(
    params: &NrtlParameters,
    T: f64,
    x: &[f64],
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
    for (field, matrix) in [("alpha", &params.alpha), ("dij", &params.dij)] {
        if matrix.len() != n * n {
            return Err(AzothError::invalid_input(
                "components",
                format!(
                    "the {field} matrix has {} entries but `x` has {n} components, \
                     which needs {} - an N x N matrix row-major",
                    matrix.len(),
                    n * n
                ),
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

    let mut tau = vec![vec![0.0; n]; n];
    let mut g = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            tau[i][j] = params.dij[i * n + j] / T;
            g[i][j] = (-params.alpha[i * n + j] * tau[i][j]).exp();
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
