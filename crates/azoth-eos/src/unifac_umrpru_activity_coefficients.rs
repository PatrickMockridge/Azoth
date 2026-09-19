//! `eos.unifac_umrpru_activity_coefficients` - the activity coefficients from UNIFAC
//! with UMR-PRU's group-interaction parameters.
//!
//! Spec: `specs/models/eos/unifac_umrpru_activity_coefficients.toml`. A *direct*
//! model: no iteration, so no `algorithm` block.

use azoth_core::{AzothError, Result, apply_checks};

use crate::databank::UnifacUmrpruParameters;
use crate::model_gen;
use crate::results::UnifacUmrpruActivityCoefficientsResult;
use crate::unifac_activity_coefficients::{Combinatorial, unifac_ln_gamma};

/// The temperature the UMR-PRU interaction is fitted about: the group parameters are
/// `a + b (T - 298.15) + c (T - 298.15)^2`, not `a + b T + c T^2`.
///
/// A constant of the model rather than of the state.
const REFERENCE_TEMPERATURE: f64 = 298.15;

/// The interaction matrix at a temperature, `a_mn(T) = a_mn + b_mn (T - 298.15) + c_mn
/// (T - 298.15)^2`, which is NeqSim's `ComponentGEUnifacUMRPRU.calcaij`.
///
/// Shared with the UMR mixing rule, which reads the same matrix for its `alpha_mix`:
/// the formula is one statement, so a change to it moves both.
#[must_use]
pub fn umrpru_aij(params: &UnifacUmrpruParameters, t: f64) -> Vec<f64> {
    let dt = t - REFERENCE_TEMPERATURE;
    params
        .aij
        .iter()
        .zip(&params.bij)
        .zip(&params.cij)
        .map(|((&a, &b), &c)| a + b * dt + c * dt * dt)
        .collect()
}

/// The activity coefficients of a mixture, from UNIFAC with UMR-PRU's parameters.
///
/// The group decomposition is NeqSim's `UNIFACcompUMRPRU` - 139 subgroups rather than
/// the 133 of the plain UNIFAC table - and the interaction is
/// `a_mn(T) = a_mn + b_mn (T - 298.15) + c_mn (T - 298.15)^2`, NeqSim's
/// `ComponentGEUnifacUMRPRU.calcaij`. **The combinatorial term is not the one
/// `eos.unifac_activity_coefficients` uses**: `ComponentGEUnifacUMRPRU` takes the
/// Flory-Huggins term alone where `ComponentGEUnifac` takes Staverman-Guggenheim. See
/// [`Combinatorial`]. `params` is resolved by name through
/// [`crate::databank::unifac_umrpru_parameters`]; `x` is checked rather than
/// renormalised.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the shapes disagree (`groups` is `N x G`, `aij` is
///   `G x G`), a group count is negative, or `x` is not a composition.
/// * [`AzothError::OutOfRange`] if `T` is not positive.
///
/// # Example
/// ```
/// use azoth_eos::databank::{UmrpruSet, unifac_umrpru_parameters};
/// use azoth_eos::unifac_umrpru_activity_coefficients;
///
/// let params = unifac_umrpru_parameters(&["methanol", "water"], UmrpruSet::Umr)?;
/// let r = unifac_umrpru_activity_coefficients(&params, 298.15, &[0.5, 0.5])?;
/// assert!(r.gamma[0] > 0.0);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T` and `x` are the symbols in the chemistry
pub fn unifac_umrpru_activity_coefficients(
    params: &UnifacUmrpruParameters,
    T: f64,
    x: &[f64],
) -> Result<UnifacUmrpruActivityCoefficientsResult> {
    let spec = &model_gen::UNIFAC_UMRPRU_ACTIVITY_COEFFICIENTS_SPEC;
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

    let aij = umrpru_aij(params, T);

    let (ln_gamma, gamma) = unifac_ln_gamma(
        &params.groups,
        &params.group_r,
        &params.group_q,
        &aij,
        T,
        x,
        Combinatorial::FloryHuggins,
    );

    Ok(UnifacUmrpruActivityCoefficientsResult {
        ln_gamma,
        gamma,
        warnings,
    })
}
