//! `eos.ge_wilson_phase` - the fugacity coefficients of a Wilson activity-coefficient
//! liquid.
//!
//! Spec: `specs/models/eos/ge_wilson_phase.toml`. A *direct* model: no iteration, so no
//! `algorithm` block.
//!
//! NeqSim's `PhaseGEWilson` is the port source. Its liquid is `ComponentGE.fugcoef`'s
//! two reference arms, the same composition `eos.ge_nrtl_phase` and `eos.ge_unifac_phase`
//! make, with `gamma_i` from the paraffin-wax Wilson correlation. Every GE phase inherits
//! that method and none overrides it.

use azoth_core::units::pascals;
use azoth_core::{AzothError, Result, apply_checks};

use crate::databank::GeWilsonPhaseParameters;
use crate::mixture::Mixture;
use crate::model_gen;
use crate::results::GeWilsonPhaseResult;

/// The fugacity coefficients of a liquid whose non-ideality is Wilson's.
///
/// `ComponentGE.fugcoef` for this phase: `gamma_i P0_i / P` on the solvent arm and
/// `(gamma_i / gamma_i^inf) H_i / P` on the Henry arm. `gamma_i` comes from the Wilson
/// correlation on `mixture` at `T` and `x`; the reference `gamma_i^inf` is the same
/// correlation on the pair or subset the arm names. `params` is resolved by name through
/// [`crate::databank::ge_wilson_phase_parameters`] and `mixture` through
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
    if params.references.len() != n {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "the phase has {} component(s) but the mixture has {n}",
                params.references.len()
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

    // The reference states and the composition. The Wilson correlation reads only each
    // component's molar mass and critical temperature, so the reference pair is a
    // `Mixture` over the two cloned components.
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
            let components: Vec<crate::mixture::Component> = reference
                .indices
                .iter()
                .map(|&i| mixture.components()[i].clone())
                .collect();
            let m = components.len();
            let subset = Mixture::new(components, vec![0.0; m * m])?;
            let evaluated = crate::wilson_activity_coefficients::wilson_activity_coefficients(
                &subset,
                T,
                &reference.x,
            )?;
            Ok(evaluated.gamma[0])
        },
    )?;
    warnings.extend(fugacities.warnings);

    let ln_gamma: Vec<f64> = gamma.iter().map(|g| g.ln()).collect();
    Ok(GeWilsonPhaseResult {
        gamma,
        ln_gamma,
        ln_phi: fugacities.ln_phi,
        p_sat: fugacities.p_sat.into_iter().map(pascals).collect(),
        henry: fugacities.henry.into_iter().map(pascals).collect(),
        gamma_inf: fugacities.gamma_inf,
        warnings,
    })
}
