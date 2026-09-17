//! `eos.ge_nrtl_phase` - the fugacity coefficients of an NRTL activity-coefficient
//! liquid.
//!
//! Spec: `specs/models/eos/ge_nrtl_phase.toml`. A *direct* model: no iteration, so no
//! `algorithm` block.
//!
//! NeqSim's `PhaseGENRTL` is the port source. Its liquid does not use a cubic at all:
//! `ComponentGE.fugcoef` sets `phi_i = gamma_i P0_i / P`, where `gamma_i` is the NRTL
//! activity coefficient and `P0_i` the pure-component saturation pressure. Both halves
//! are already ported - [`crate::nrtl_activity_coefficients`] and
//! [`crate::antoine_vapor_pressure`] - and this is their composition rather than a third
//! piece of physics.

use azoth_core::units::pascals;
use azoth_core::{AzothError, Result, apply_checks};

use crate::databank::GeNrtlPhaseParameters;
use crate::model_gen;
use crate::results::GeNrtlPhaseResult;

/// The fugacity coefficients of a liquid whose non-ideality is NRTL's.
///
/// `phi_i = gamma_i P0_i / P`, which is NeqSim's `ComponentGE.fugcoef` for this phase.
/// `gamma_i` comes from NRTL on `params` at `T` and `x`; `P0_i` is Antoine's correlation
/// on the component's own coefficients at `T`. `params` is resolved by name through
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

    // The activity coefficients, from the same matrices `eos.nrtl_activity_coefficients`
    // uses - resolved through the same function, so the two cannot disagree.
    let nrtl = crate::databank::NrtlParameters {
        alpha: params.alpha.clone(),
        dij: params.dij.clone(),
    };
    let gamma = crate::nrtl_activity_coefficients::nrtl_activity_coefficients(&nrtl, T, x)?.gamma;

    // The saturation pressures and the composition, in the one place every GE phase
    // shares them.
    let saturated = crate::ge_phase::ge_fugacities(&gamma, &params.antoine, T, P)?;
    warnings.extend(saturated.warnings);

    let ln_gamma: Vec<f64> = gamma.iter().map(|g| g.ln()).collect();
    Ok(GeNrtlPhaseResult {
        gamma,
        ln_gamma,
        ln_phi: saturated.ln_phi,
        p_sat: saturated.p_sat.into_iter().map(pascals).collect(),
        warnings,
    })
}
