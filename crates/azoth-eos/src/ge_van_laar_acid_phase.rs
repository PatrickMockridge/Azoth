//! `eos.ge_van_laar_acid_phase` - the fugacity coefficients of the water-nitric-sulfuric
//! acid liquid.
//!
//! Spec: `specs/models/eos/ge_van_laar_acid_phase.toml`. A *direct* model: no iteration,
//! so no `algorithm` block.
//!
//! NeqSim's `ComponentGEVanLaarAcid.fugcoef` is the port source, and it is short:
//!
//! ```text
//! computeGamma(phase);
//! if (acidIndex == 0) return NON_MODELED_COMPONENT_PENALTY;   // 1.0e12
//! fugacityCoefficient = gamma * pureVaporPressureBar(T) / P;
//! ```
//!
//! It is the one GE phase whose component overrides `fugcoef`, and it does so to **ignore**
//! `referenceStateType`: both acids are tagged `solute` in the component database, so the
//! inherited method would give them a Henry's-law coefficient, and the model's whole point
//! is the symmetric Raoult identity `f_i = gamma_i x_i P0_i`.
//!
//! Its `P0` is not Antoine's either. For the three modelled acids it comes from
//! [`crate::nitric_sulfuric_acid_vapor_pressure`]; only a species the model does not cover
//! falls back to the database Antoine. That is why [`crate::ge_phase::ge_fugacities`] is
//! split: the composition is shared, the `P0` is not.

use azoth_core::units::pascals;
use azoth_core::{AzothError, Result, apply_checks};

use crate::databank::{AntoineRecord, GeVanLaarAcidPhaseParameters, VanLaarAcidParameters};
use crate::ge_phase::{combine, saturation};
use crate::model_gen;
use crate::results::GeVanLaarAcidPhaseResult;

/// The fugacity coefficient NeqSim gives a species the model does not cover.
///
/// A *penalty* rather than physics, and a large finite number rather than an infinity:
/// driving the component out of the liquid is the intent, and a finite value leaves a
/// flash over the mixture able to iterate while it does so.
const NON_MODELED_COMPONENT_PENALTY: f64 = 1.0e12;

/// The fugacity coefficients of the water-nitric-sulfuric acid liquid.
///
/// `phi_i = gamma_i P0_i / P`, where `gamma_i` is
/// [`crate::van_laar_acid_activity_coefficients`]' and `P0_i` is the acid correlation for
/// the three modelled species and the database Antoine for anything else. `params` is
/// resolved by name through [`crate::databank::ge_van_laar_acid_phase_parameters`]; `x` is
/// checked rather than renormalised.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `params` and `x` disagree in length, or if `x` is not
///   a composition.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive, or outside the acid
///   correlation's own bounds.
#[allow(non_snake_case)] // `T`, `P` and `x` are the symbols in the chemistry
pub fn ge_van_laar_acid_phase(
    params: &GeVanLaarAcidPhaseParameters,
    T: f64,
    P: f64,
    x: &[f64],
) -> Result<GeVanLaarAcidPhaseResult> {
    let spec = &model_gen::GE_VAN_LAAR_ACID_PHASE_SPEC;
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
    if params.antoine.len() != n || params.acid_index.len() != n {
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

    // The activity coefficients, from the same acid identities and the same renormalised
    // acid basis `eos.van_laar_acid_activity_coefficients` uses.
    let acid = VanLaarAcidParameters {
        acid_index: params.acid_index.clone(),
    };
    let activity = crate::van_laar_acid_activity_coefficients::van_laar_acid_activity_coefficients(
        &acid, T, x,
    )?;
    warnings.extend(activity.warnings);
    let gamma = activity.gamma;

    // `P0` per component: the acid correlation for the three modelled species, the database
    // Antoine for anything else. **`saturation` is asked only about the components that
    // actually take the fallback**, because a component the acid correlation covers must not
    // be refused for having no Antoine row of its own - the three of them carry upstream's
    // `none` marker, and asking would refuse this model for its own components.
    let fallback: Vec<AntoineRecord> = params
        .antoine
        .iter()
        .zip(&params.acid_index)
        .filter(|(_, index)| **index == 0)
        .map(|(record, _)| record.clone())
        .collect();
    let (fallback_p_sat, antoine_warnings) = saturation(&fallback, T)?;
    warnings.extend(antoine_warnings);
    let acids = crate::nitric_sulfuric_acid_vapor_pressure::nitric_sulfuric_acid_vapor_pressure(
        azoth_core::units::kelvins(T),
    )?;
    warnings.extend(acids.warnings);

    let mut next = fallback_p_sat.into_iter();
    let p_sat: Vec<f64> = params
        .acid_index
        .iter()
        .map(|&index| match index {
            1 => acids.p_water.value,
            2 => acids.p_nitric_acid.value,
            3 => acids.p_sulfuric_acid.value,
            _ => next
                .next()
                .expect("one fallback pressure per component the acid correlation does not cover"),
        })
        .collect();

    // The penalty replaces the whole coefficient, not just `P0`: NeqSim returns it before
    // computing anything else, so an uncovered component's `ln_phi` is `ln(1e12)` whatever
    // its composition.
    let ln_phi: Vec<f64> = params
        .acid_index
        .iter()
        .zip(combine(&gamma, &p_sat, P))
        .map(|(&index, composed)| {
            if index == 0 {
                NON_MODELED_COMPONENT_PENALTY.ln()
            } else {
                composed
            }
        })
        .collect();

    let ln_gamma: Vec<f64> = gamma.iter().map(|g| g.ln()).collect();
    Ok(GeVanLaarAcidPhaseResult {
        gamma,
        ln_gamma,
        ln_phi,
        p_sat: p_sat.into_iter().map(pascals).collect(),
        warnings,
    })
}
