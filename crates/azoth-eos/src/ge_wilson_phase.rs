//! `eos.ge_wilson_phase` - the fugacity coefficients of a Wilson activity-coefficient
//! liquid.
//!
//! Spec: `specs/models/eos/ge_wilson_phase.toml`. A *direct* model: no iteration, so no
//! `algorithm` block.
//!
//! NeqSim's `PhaseGEWilson` is the port source. Its liquid is `ComponentGE.fugcoef`'s
//! `phi_i = gamma_i P0_i / P`, the same composition `eos.ge_nrtl_phase` and
//! `eos.ge_unifac_phase` make, with `gamma_i` from the paraffin-wax Wilson correlation.
//! Every GE phase inherits that method and none overrides it, so this is a composition of
//! two already-ported pieces rather than a third piece of physics.

use azoth_core::units::pascals;
use azoth_core::{AzothError, Result, apply_checks};

use crate::databank::GeWilsonPhaseParameters;
use crate::ge_phase::ge_fugacities;
use crate::mixture::Mixture;
use crate::model_gen;
use crate::results::GeWilsonPhaseResult;

/// The fugacity coefficients of a liquid whose non-ideality is Wilson's.
///
/// `phi_i = gamma_i P0_i / P`, which is NeqSim's `ComponentGE.fugcoef` for this phase.
/// `gamma_i` comes from the Wilson correlation on `mixture` at `T` and `x`; `P0_i` is
/// Antoine's correlation on the component's own coefficients at `T`. `params` is resolved
/// by name through [`crate::databank::ge_wilson_phase_parameters`] and `mixture` through
/// [`crate::databank::mixture_of`], from the same names; `x` is checked rather than
/// renormalised.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `params`, `mixture` and `x` disagree in length, or if
///   `x` is not a composition.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
///
/// # Example
/// ```
/// use azoth_eos::Cubic;
/// use azoth_eos::databank::{ge_wilson_phase_parameters, mixture_of};
/// use azoth_eos::ge_wilson_phase;
///
/// let names = ["n-octane", "nc10"];
/// let params = ge_wilson_phase_parameters(&names, None)?;
/// let (mixture, _) = mixture_of(&names, Cubic::Pr, None)?;
/// let r = ge_wilson_phase(&params, &mixture, 350.0, 100_000.0, &[0.4, 0.6])?;
/// assert!((r.gamma[0] - 1.3244982894719637).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T`, `P` and `x` are the symbols in the chemistry
pub fn ge_wilson_phase(
    params: &GeWilsonPhaseParameters,
    mixture: &Mixture,
    T: f64,
    P: f64,
    x: &[f64],
) -> Result<GeWilsonPhaseResult> {
    let spec = &model_gen::GE_WILSON_PHASE_SPEC;
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

    let n = mixture.len();
    if params.antoine.len() != n {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "the phase has {} component(s) but the mixture has {n}",
                params.antoine.len()
            ),
        ));
    }
    if x.len() != n {
        return Err(AzothError::invalid_input(
            "x",
            format!("a mixture of {n} components has {} mole fractions", x.len()),
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

    // The activity coefficients, from the same correlation
    // `eos.wilson_activity_coefficients` uses - resolved through the same function, so the
    // two cannot disagree.
    let activity =
        crate::wilson_activity_coefficients::wilson_activity_coefficients(mixture, T, x)?;
    warnings.extend(activity.warnings);
    let gamma = activity.gamma;

    // The saturation pressures and the composition, in the one place every GE phase
    // shares them.
    let saturated = ge_fugacities(&gamma, &params.antoine, T, P)?;
    warnings.extend(saturated.warnings);

    let ln_gamma: Vec<f64> = gamma.iter().map(|g| g.ln()).collect();
    Ok(GeWilsonPhaseResult {
        gamma,
        ln_gamma,
        ln_phi: saturated.ln_phi,
        p_sat: saturated.p_sat.into_iter().map(pascals).collect(),
        warnings,
    })
}
