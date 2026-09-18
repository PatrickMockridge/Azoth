//! The Soave-Redlich-Kwong CPA phase state.
//!
//! `PhaseSrkCPA` as a *model*: the cubic's attraction and covolume replaced by each
//! component's fitted `aCPA`/`bCPA`, mixed with the `cpakij_SRK` column, and the Wertheim
//! association contribution added to the residual Helmholtz energy.
//!
//! **The fluid is resolved from names here, and again in Python.** `eos.eos_cg_phase` sets
//! the precedent: the names cross the boundary unresolved and each implementation looks
//! them up in its own databank. That is what makes the two-kernel comparison cover the
//! *resolution* as well as the arithmetic, which matters more for this model than for any
//! other in the library - an associating mixture mixes with `cpakij_SRK` and a classical
//! one with `KIJPR`, and on water/methanol those differ by a factor of two.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result};

use crate::association::R;
use crate::databank;
use crate::mixture::{Mixture, RootSide};
use crate::model_gen;
use crate::results::SrkCpaPhaseResult;

/// One CPA phase's state at a temperature, pressure and composition.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `z` is not one entry per component, is not a
///   composition, or a name is not in the databank.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive, or no volume root exists
///   above the mixture's covolume.
pub fn srk_cpa_phase(
    components: &[String],
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    compressed_phase: &str,
) -> Result<SrkCpaPhaseResult> {
    let spec = &model_gen::SRK_CPA_PHASE_SPEC;
    let mut warnings = Vec::new();

    azoth_core::range::apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "P" => Some(p.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let n = components.len();
    if z.len() != n {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "a mixture of {n} components needs {n} mole fractions, but z has {}",
                z.len()
            ),
        ));
    }
    let sum: f64 = z.iter().sum();
    if (sum - 1.0).abs() > 1.0e-9 {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "the mole fractions sum to {sum}, not to one. Renormalising them here would \
                 make a composition error invisible in every number downstream, so it is \
                 refused instead"
            ),
        ));
    }

    let side = match compressed_phase {
        "liquid" => RootSide::Liquid,
        "vapour" => RootSide::Vapour,
        other => {
            return Err(AzothError::invalid_input(
                "compressed_phase",
                format!("is {other:?}; the roots this model has are \"liquid\" and \"vapour\""),
            ));
        }
    };

    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, _) = databank::associating_mixture_of(&names, crate::Cubic::Srk, None)?;
    let reduced = mixture.reduced_parameters(t, p)?;
    warnings.extend(reduced.warnings.iter().cloned());

    let state = mixture.phase_state(&reduced, z, side)?;
    let r_t = R * t.value;
    Ok(SrkCpaPhaseResult {
        z_factor: state.z,
        ln_phi: state.ln_phi,
        h_res: azoth_core::units::joules_per_mole(state.h_dep_rt * r_t),
        s_res: azoth_core::units::joules_per_mole_kelvin(state.s_dep_r * R),
        warnings,
    })
}

/// A mixture's phase state without the databank, for a caller that has one already.
///
/// The same arithmetic as [`srk_cpa_phase`] with the resolution lifted out, which is what
/// the transport calls when a caller has built the fluid themselves and what a test uses
/// to pin the two apart.
///
/// # Errors
/// As [`srk_cpa_phase`], without the names.
pub fn srk_cpa_phase_of(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    side: RootSide,
) -> Result<SrkCpaPhaseResult> {
    let reduced = mixture.reduced_parameters(t, p)?;
    let state = mixture.phase_state(&reduced, z, side)?;
    let r_t = R * t.value;
    Ok(SrkCpaPhaseResult {
        z_factor: state.z,
        ln_phi: state.ln_phi,
        h_res: azoth_core::units::joules_per_mole(state.h_dep_rt * r_t),
        s_res: azoth_core::units::joules_per_mole_kelvin(state.s_dep_r * R),
        warnings: reduced.warnings,
    })
}
