//! The CPA phase state, cubic family aside.
//!
//! `eos.srk_cpa_phase` and `eos.pr_cpa_phase` are the same model with a different fitted
//! parameter set: the cubic's attraction and covolume replaced by each component's fitted
//! `aCPA`/`bCPA`, mixed with that family's `cpakij` column, and the Wertheim association
//! contribution added to the residual Helmholtz energy. Everything else - the validation,
//! the root choice, the departure functions - is one thing, so it lives here once and the
//! two models are the two lines that name the family and the result type.
//!
//! **The fluid is resolved from names here, and again in Python.** `eos.eos_cg_phase` sets
//! the precedent: the names cross the boundary unresolved and each implementation looks
//! them up in its own databank. That is what makes the two-kernel comparison cover the
//! *resolution* as well as the arithmetic, which matters more for this model than for any
//! other in the library - an associating mixture mixes with the `cpa` column and a
//! classical one with `KIJPR`, and on water/methanol those differ by a factor of two.
//!
//! **The family is not a mixing rule.** `AssociationCubic::of` reads it off the cubic, and
//! the two fitted sets are separate fits rather than one converted into the other - water's
//! `kappa_AB` is 0.0692 for SRK against 0.046473789 for PR - so choosing is a selection and
//! reading the wrong family is a different fluid.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result};

use crate::association::R;
use crate::databank;
use crate::mixture::{Mixture, RootSide};
use azoth_core::spec::ModelSpec;

/// One CPA phase's state, before it is wrapped in the family's own result type.
pub(crate) struct CpaPhaseState {
    /// The compressibility factor at the chosen root.
    pub(crate) z_factor: f64,
    /// `ln phi_i` at the chosen root.
    pub(crate) ln_phi: Vec<f64>,
    /// The residual enthalpy, in J/mol.
    pub(crate) h_res: f64,
    /// The residual entropy, in J/(mol*K).
    pub(crate) s_res: f64,
    /// Whatever the range checks and the reduction raised.
    pub(crate) warnings: Vec<azoth_core::Warning>,
}

/// One CPA phase's state at a temperature, pressure and composition.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `z` is not one entry per component, is not a
///   composition, or a name is not in the databank.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive, or no volume root exists
///   above the mixture's covolume.
#[allow(clippy::too_many_arguments)] // The spec's declared inputs plus the family.
pub(crate) fn cpa_phase(
    spec: &ModelSpec,
    components: &[String],
    cubic: crate::Cubic,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    compressed_phase: &str,
) -> Result<CpaPhaseState> {
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

    let side = side_of(compressed_phase)?;
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, _) = databank::associating_mixture_of(&names, cubic, None)?;
    let state = phase_state_of(&mixture, t, p, z, side)?;
    warnings.extend(state.warnings);
    Ok(CpaPhaseState { warnings, ..state })
}

/// A mixture's phase state without the databank, for a caller that has one already.
///
/// The same arithmetic as [`cpa_phase`] with the resolution lifted out, which is what the
/// transport calls when a caller has built the fluid themselves and what a test uses to
/// pin the two apart.
///
/// # Errors
/// As [`cpa_phase`], without the names.
pub(crate) fn phase_state_of(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    side: RootSide,
) -> Result<CpaPhaseState> {
    let reduced = mixture.reduced_parameters(t, p)?;
    let state = mixture.phase_state(&reduced, z, side)?;
    let r_t = R * t.value;
    Ok(CpaPhaseState {
        z_factor: state.z,
        ln_phi: state.ln_phi,
        h_res: state.h_dep_rt * r_t,
        s_res: state.s_dep_r * R,
        warnings: reduced.warnings,
    })
}

/// The root a `compressed_phase` string names.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the spelling is neither `"liquid"` nor `"vapour"`.
pub(crate) fn side_of(compressed_phase: &str) -> Result<RootSide> {
    match compressed_phase {
        "liquid" => Ok(RootSide::Liquid),
        "vapour" | "vapor" => Ok(RootSide::Vapour),
        other => Err(AzothError::invalid_input(
            "compressed_phase",
            format!("is {other:?}; the roots this model has are \"liquid\" and \"vapour\""),
        )),
    }
}

/// The same, for a model whose mixture is not a cubic family plus a `kij` column.
///
/// `eos.umr_cpa_phase` resolves its fluid through the UMR mixing rule rather than through
/// `associating_mixture_of`, and every other line - the range checks, the composition
/// check, the root choice, the departure - is the same. Taking the resolver as an argument
/// keeps one expression of those rather than a second copy that would have to be kept in
/// step.
///
/// # Errors
/// As [`cpa_phase`], plus whatever the resolver refuses.
pub(crate) fn cpa_phase_with(
    spec: &ModelSpec,
    components: &[String],
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    compressed_phase: &str,
    resolve: impl FnOnce(&[&str]) -> Result<Mixture>,
) -> Result<CpaPhaseState> {
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

    let side = side_of(compressed_phase)?;
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let mixture = resolve(&names)?;
    let state = phase_state_of(&mixture, t, p, z, side)?;
    warnings.extend(state.warnings);
    Ok(CpaPhaseState { warnings, ..state })
}
