//! `eos.ge_uniquac_phase` - the fugacity coefficients of a UNIQUAC
//! activity-coefficient liquid.
//!
//! Spec: `specs/models/eos/ge_uniquac_phase.toml`. A *direct* model: no iteration, so no
//! `algorithm` block.
//!
//! NeqSim's `PhaseGEUniquac` is the port source. Its liquid is `ComponentGE.fugcoef`'s
//! `phi_i = gamma_i P0_i / P`, the same composition the other GE phases make - every GE
//! phase inherits that method and none overrides it - with `gamma_i` from UNIQUAC.
//! Both halves are already ported and this is their composition rather than a third
//! piece of physics.
//!
//! NeqSim's own UNIQUAC `gamma` does not exist: `ComponentGEUniquac.getGamma` returns
//! `0.0` and its formula is commented out, so unlike the other phases this one has no
//! oracle for its composition. The `r`/`q` half is checked instead, against
//! `ComponentGEUnifac.getR`/`getQ`.

use azoth_core::units::pascals;
use azoth_core::{AzothError, Result, apply_checks};

use crate::databank::{GeUniquacPhaseParameters, UniquacParameters};
use crate::ge_phase::ge_fugacities;
use crate::model_gen;
use crate::results::GeUniquacPhaseResult;

/// The fugacity coefficients of a liquid whose non-ideality is UNIQUAC's.
///
/// `phi_i = gamma_i P0_i / P`, which is NeqSim's `ComponentGE.fugcoef` for this phase.
/// `gamma_i` comes from UNIQUAC on `params`' `r` and `q` with the caller's `aij`;
/// `P0_i` is Antoine's correlation on the component's own coefficients at `T`. `params`
/// is resolved by name through [`crate::databank::ge_uniquac_phase_parameters`]; `x` is
/// checked rather than renormalised.
///
/// `aij` is `N x N`, in kelvin, directional with a zero diagonal - the same argument
/// `eos.uniquac_activity_coefficients` takes, and for the same reason: no upstream table
/// carries it.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `params` and `x` disagree in length, if `x` is not a
///   composition, or if `aij` is not an `N x N` matrix with a zero diagonal.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
#[allow(non_snake_case)] // `T`, `P` and `x` are the symbols in the chemistry
pub fn ge_uniquac_phase(
    params: &GeUniquacPhaseParameters,
    T: f64,
    P: f64,
    x: &[f64],
    aij: &[Vec<f64>],
) -> Result<GeUniquacPhaseResult> {
    let spec = &model_gen::GE_UNIQUAC_PHASE_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T),
            "P" => Some(P),
            _ => None,
        },
        &mut warnings,
    )?;

    let n = x.len();
    if n == 0 {
        return Err(AzothError::invalid_input(
            "x",
            "a mixture of zero components has no fugacity coefficient",
        ));
    }
    if params.antoine.len() != n {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "the phase has {} component(s) but `x` has {n} entries",
                params.antoine.len()
            ),
        ));
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
    if aij.len() != n || aij.iter().any(|row| row.len() != n) {
        return Err(AzothError::invalid_input(
            "aij",
            format!(
                "the interaction matrix is {} x {} but the mixture has {n} components",
                aij.len(),
                aij.first().map_or(0, Vec::len)
            ),
        ));
    }

    // The activity coefficients, from the same r and q and the same interaction matrix
    // `eos.uniquac_activity_coefficients` uses - resolved through the same function, so the
    // two cannot disagree.
    let uniquac = UniquacParameters {
        r: params.r.clone(),
        q: params.q.clone(),
    };
    let activity =
        crate::uniquac_activity_coefficients::uniquac_activity_coefficients(&uniquac, T, x, aij)?;
    warnings.extend(activity.warnings);
    let gamma = activity.gamma;

    // The saturation pressures and the composition, in the one place every GE phase
    // shares them.
    let saturated = ge_fugacities(&gamma, &params.antoine, T, P)?;
    warnings.extend(saturated.warnings);

    let ln_gamma: Vec<f64> = gamma.iter().map(|g| g.ln()).collect();
    Ok(GeUniquacPhaseResult {
        gamma,
        ln_gamma,
        ln_phi: saturated.ln_phi,
        p_sat: saturated.p_sat.into_iter().map(pascals).collect(),
        warnings,
    })
}
