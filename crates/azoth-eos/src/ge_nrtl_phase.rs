//! `eos.ge_nrtl_phase` - the fugacity coefficients of an NRTL activity-coefficient
//! liquid.
//!
//! Spec: `specs/models/eos/ge_nrtl_phase.toml`. A *direct* model: no iteration, so no
//! `algorithm` block.
//!
//! NeqSim's `PhaseGENRTL` is the port source. Its liquid does not use a cubic at all:
//! `ComponentGE.fugcoef` sets `phi_i = gamma_i P0_i / P` for a `solvent` component and
//! `phi_i = (gamma_i / gamma_i^inf) H_i / P` for every other one. `gamma_i` is the NRTL
//! activity coefficient; `P0_i` is the pure-component saturation pressure and `H_i` the
//! Henry coefficient. The activity model, Antoine and the Henry reference state are all
//! ported, and this is their composition rather than a third piece of physics.

use azoth_core::units::pascals;
use azoth_core::{AzothError, Result, apply_checks};

use crate::databank::GeNrtlPhaseParameters;
use crate::model_gen;
use crate::results::GeNrtlPhaseResult;

/// The fugacity coefficients of a liquid whose non-ideality is NRTL's.
///
/// `ComponentGE.fugcoef` for this phase: `gamma_i P0_i / P` on the solvent arm and
/// `(gamma_i / gamma_i^inf) H_i / P` on the Henry arm. `gamma_i` comes from NRTL on
/// `params` at `T` and `x`; the reference `gamma_i^inf` is the same model on the pair the
/// Henry arm constructs. `params` is resolved by name through
/// [`crate::databank::ge_nrtl_phase_parameters`]; `x` is checked rather than
/// renormalised.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `params` and `x` disagree in length, if `x` is not
///   a composition, or if a component's Antoine label names no form this build evaluates.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
///
/// # Example
/// ```
/// use azoth_eos::databank::ge_nrtl_phase_parameters;
/// use azoth_eos::ge_nrtl_phase;
///
/// let params = ge_nrtl_phase_parameters(&["methanol", "water"], None)?;
/// let r = ge_nrtl_phase(&params, 298.15, 100_000.0, &[0.5, 0.5])?;
/// assert!((r.gamma[0] - 1.2331788561677834).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T`, `P` and `x` are the symbols in the chemistry
pub fn ge_nrtl_phase(
    params: &GeNrtlPhaseParameters,
    T: f64,
    P: f64,
    x: &[f64],
) -> Result<GeNrtlPhaseResult> {
    let spec = &model_gen::GE_NRTL_PHASE_SPEC;
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
    if params.references.len() != n {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "the phase has {} component(s) but `x` has {n} entries",
                params.references.len()
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

    // The activity coefficients, from the same matrices `eos.nrtl_activity_coefficients`
    // uses - resolved through the same function, so the two cannot disagree.
    let nrtl = crate::databank::NrtlParameters {
        alpha: params.alpha.clone(),
        dij: params.dij.clone(),
    };
    let gamma = crate::nrtl_activity_coefficients::nrtl_activity_coefficients(&nrtl, T, x)?.gamma;

    // The reference states and the composition, in the one place every GE phase shares
    // them. The reference `gamma^inf` is the phase's own NRTL over the 2 x 2 submatrix of
    // the pair - or, without water, over the diluted full composition.
    let water = params
        .names
        .iter()
        .position(|name| name.eq_ignore_ascii_case("water"));
    let fugacities = crate::ge_phase::ge_fugacities_with_references(
        &gamma,
        &params.references,
        x,
        water,
        T,
        P,
        |_, reference| {
            let m = reference.indices.len();
            let mut alpha = Vec::with_capacity(m * m);
            let mut dij = Vec::with_capacity(m * m);
            for &i in &reference.indices {
                for &j in &reference.indices {
                    alpha.push(params.alpha[i * n + j]);
                    dij.push(params.dij[i * n + j]);
                }
            }
            let subset = crate::databank::NrtlParameters { alpha, dij };
            let evaluated = crate::nrtl_activity_coefficients::nrtl_activity_coefficients(
                &subset,
                T,
                &reference.x,
            )?;
            Ok(evaluated.gamma[0])
        },
    )?;
    warnings.extend(fugacities.warnings);

    let ln_gamma: Vec<f64> = gamma.iter().map(|g| g.ln()).collect();
    Ok(GeNrtlPhaseResult {
        gamma,
        ln_gamma,
        ln_phi: fugacities.ln_phi,
        p_sat: fugacities.p_sat.into_iter().map(pascals).collect(),
        henry: fugacities.henry.into_iter().map(pascals).collect(),
        gamma_inf: fugacities.gamma_inf,
        warnings,
    })
}
