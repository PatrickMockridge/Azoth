//! `eos.ge_unifac_phase` - the fugacity coefficients of a UNIFAC activity-coefficient
//! liquid.
//!
//! Spec: `specs/models/eos/ge_unifac_phase.toml`. A *direct* model: no iteration, so no
//! `algorithm` block.
//!
//! NeqSim's `PhaseGEUnifac` is the port source. Its liquid is `ComponentGE.fugcoef`'s
//! two reference arms, the composition `eos.ge_nrtl_phase` makes - every GE phase
//! inherits that method and none overrides it - with `gamma_i` from UNIFAC's group tables
//! instead of NRTL's matrices.

use azoth_core::units::pascals;
use azoth_core::{AzothError, Result, apply_checks};

use crate::databank::{GeUnifacPhaseParameters, UnifacParameters};
use crate::model_gen;
use crate::results::GeUnifacPhaseResult;

/// The fugacity coefficients of a liquid whose non-ideality is UNIFAC's.
///
/// `ComponentGE.fugcoef` for this phase: `gamma_i P0_i / P` on the solvent arm and
/// `(gamma_i / gamma_i^inf) H_i / P` on the Henry arm. `gamma_i` comes from UNIFAC on
/// `params`' group tables at `T` and `x`; the reference `gamma_i^inf` is the same model on
/// the pair or subset the arm names. `params` is resolved by name through
/// [`crate::databank::ge_unifac_phase_parameters`]; `x` is checked rather than
/// renormalised.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `params` and `x` disagree in length, or if `x` is not
///   a composition.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
///
/// # Example
/// ```
/// use azoth_eos::databank::ge_unifac_phase_parameters;
/// use azoth_eos::ge_unifac_phase;
///
/// let params = ge_unifac_phase_parameters(&["methanol", "water"], None)?;
/// let r = ge_unifac_phase(&params, 298.15, 100_000.0, &[0.5, 0.5])?;
/// assert!((r.gamma[0] - 1.1156815062468024).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T`, `P` and `x` are the symbols in the chemistry
pub fn ge_unifac_phase(
    params: &GeUnifacPhaseParameters,
    T: f64,
    P: f64,
    x: &[f64],
) -> Result<GeUnifacPhaseResult> {
    let spec = &model_gen::GE_UNIFAC_PHASE_SPEC;
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

    // The activity coefficients, from the same group tables
    // `eos.unifac_activity_coefficients` uses - resolved through the same function, so the
    // two cannot disagree.
    let unifac = UnifacParameters {
        groups: params.groups.clone(),
        group_r: params.group_r.clone(),
        group_q: params.group_q.clone(),
        aij: params.aij.clone(),
    };
    let activity =
        crate::unifac_activity_coefficients::unifac_activity_coefficients(&unifac, T, x)?;
    warnings.extend(activity.warnings);
    let gamma = activity.gamma;

    // The reference states and the composition. UNIFAC divides by `x_i` in its
    // combinatorial term, so the reference is evaluated on the subset the arm names,
    // never on the full list with the other components zeroed.
    let water = params
        .names
        .iter()
        .position(|name| name.eq_ignore_ascii_case("water"));
    let groups_per_component = params.group_r.len();
    let fugacities = crate::ge_phase::ge_fugacities_with_references(
        &gamma,
        &params.references,
        x,
        water,
        T,
        P,
        |_, reference| {
            let mut groups = Vec::with_capacity(reference.indices.len() * groups_per_component);
            for &i in &reference.indices {
                groups.extend_from_slice(
                    &params.groups[i * groups_per_component..(i + 1) * groups_per_component],
                );
            }
            let subset = UnifacParameters {
                groups,
                group_r: params.group_r.clone(),
                group_q: params.group_q.clone(),
                aij: params.aij.clone(),
            };
            let evaluated = crate::unifac_activity_coefficients::unifac_activity_coefficients(
                &subset,
                T,
                &reference.x,
            )?;
            Ok(evaluated.gamma[0])
        },
    )?;
    warnings.extend(fugacities.warnings);

    let ln_gamma: Vec<f64> = gamma.iter().map(|g| g.ln()).collect();
    Ok(GeUnifacPhaseResult {
        gamma,
        ln_gamma,
        ln_phi: fugacities.ln_phi,
        p_sat: fugacities.p_sat.into_iter().map(pascals).collect(),
        henry: fugacities.henry.into_iter().map(pascals).collect(),
        gamma_inf: fugacities.gamma_inf,
        warnings,
    })
}
