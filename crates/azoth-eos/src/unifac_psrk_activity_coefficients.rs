//! `eos.unifac_psrk_activity_coefficients` - the activity coefficients from UNIFAC
//! with the PSRK temperature-dependent interaction parameters.
//!
//! Spec: `specs/models/eos/unifac_psrk_activity_coefficients.toml`. A *direct* model:
//! no iteration, so no `algorithm` block.

use azoth_core::{AzothError, Result, apply_checks};

use crate::databank::UnifacPsrkParameters;
use crate::model_gen;
use crate::results::UnifacPsrkActivityCoefficientsResult;
use crate::unifac_activity_coefficients::{Combinatorial, unifac_ln_gamma};

/// The activity coefficients of a mixture, from UNIFAC with PSRK's interaction
/// parameters.
///
/// Identical to `eos.unifac_activity_coefficients` except that the main-group
/// interaction is a function of temperature:
/// `a_mn(T) = a_mn + b_mn T + c_mn T^2`, NeqSim's
/// `ComponentGEUnifacPSRK.calcaij`. The group basis, the combinatorial term and the
/// residual are the same, so this evaluates the interaction matrix at `T` and calls the
/// shared body. `params` is resolved by name through
/// [`crate::databank::unifac_psrk_parameters`]; `x` is checked rather than renormalised.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the shapes disagree (`groups` is `N x G`, `aij` is
///   `G x G`), a group count is negative, or `x` is not a composition.
/// * [`AzothError::OutOfRange`] if `T` is not positive.
///
/// # Example
/// ```
/// use azoth_eos::databank::unifac_psrk_parameters;
/// use azoth_eos::unifac_psrk_activity_coefficients;
///
/// // Water/methane is a pair whose b and c are non-zero, so the temperature matters.
/// let params = unifac_psrk_parameters(&["water", "methane"])?;
/// let cold = unifac_psrk_activity_coefficients(&params, 250.0, &[0.5, 0.5])?;
/// let hot = unifac_psrk_activity_coefficients(&params, 400.0, &[0.5, 0.5])?;
/// assert_ne!(cold.gamma[0], hot.gamma[0]);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T` and `x` are the symbols in the chemistry
pub fn unifac_psrk_activity_coefficients(
    params: &UnifacPsrkParameters,
    T: f64,
    x: &[f64],
) -> Result<UnifacPsrkActivityCoefficientsResult> {
    let spec = &model_gen::UNIFAC_PSRK_ACTIVITY_COEFFICIENTS_SPEC;
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
    let g = params.group_r.len();
    if g == 0 || params.group_q.len() != g {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "`group_r` has {g} entries and `group_q` has {}; both must be the same \
                 non-empty length",
                params.group_q.len()
            ),
        ));
    }
    if params.groups.len() != n * g {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "`groups` has {} entries but `x` has {n} components and there are {g} \
                 groups, which needs {}",
                params.groups.len(),
                n * g
            ),
        ));
    }
    if let Some(bad) = params.groups.iter().position(|&count| count < 0.0) {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "group count {bad} is {} but a count cannot be negative",
                bad
            ),
        ));
    }
    for (field, matrix) in [
        ("aij", &params.aij),
        ("bij", &params.bij),
        ("cij", &params.cij),
    ] {
        if matrix.len() != g * g {
            return Err(AzothError::invalid_input(
                "components",
                format!(
                    "`{field}` has {} entries but must be {g} x {g}",
                    matrix.len()
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

    let aij: Vec<f64> = params
        .aij
        .iter()
        .zip(&params.bij)
        .zip(&params.cij)
        .map(|((&a, &b), &c)| a + b * T + c * T * T)
        .collect();

    let (ln_gamma, gamma) = unifac_ln_gamma(
        &params.groups,
        &params.group_r,
        &params.group_q,
        &aij,
        T,
        x,
        Combinatorial::StavermanGuggenheim,
    );

    Ok(UnifacPsrkActivityCoefficientsResult {
        ln_gamma,
        gamma,
        warnings,
    })
}
