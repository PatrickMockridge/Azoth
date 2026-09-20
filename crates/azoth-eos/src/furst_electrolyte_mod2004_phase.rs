//! `eos.furst_electrolyte_mod2004_phase` - the 2004 revision of the Fürst phase state.
//!
//! Spec: `specs/models/eos/furst_electrolyte_mod2004_phase.toml`. A *direct* model: the root is
//! solved by the shared cubic-plus-pressure iteration and there is no outer loop, so no
//! `algorithm` block.
//!
//! NeqSim's `SystemFurstElectrolyteEos`, which is an SRK cubic with Schwartzentruber's
//! attractive term for every component, carrying **three additive Helmholtz terms** whose
//! pressure moves the root:
//!
//! ```text
//! F = F_SRK + FSR2 sr2On + FLR lrOn + FBorn bornOn
//! ```
//!
//! the short-range `W` term, the MSA long-range term, and the Born solvation term. The
//! arithmetic is [`crate::furst_electrolyte`]'s; this module is the databank resolution, the
//! mixing rule the system installs, and the root.
//!
//! # What this model does not have
//!
//! **No flash.** `eos.pt_flash`'s Jacobian is `d ln phi / d n`, which the Huron-Vidal mixing
//! rule this system installs refuses - that derivative is the excess Gibbs energy's second
//! derivative. The boundary is the one `eos.umr_cpa_phase` and `eos.soreide_whitson_phase`
//! sit behind, so the cases record phase states at compositions NeqSim reports.
//!
//! **And no `h_res` or `s_res`.** The departure is not assembled, so a model reporting one
//! would be reporting the cubic's.

use azoth_core::Result;
use azoth_core::units::{Pressure, ThermodynamicTemperature};

use crate::model_gen;
use crate::results::FurstElectrolyteMod2004PhaseResult;

/// One Fürst electrolyte phase's state, as the **2004 revision** computes it.
///
/// The same kernels as `eos.furst_electrolyte_phase` with five quantities zeroed and one
/// term added - see the spec. The zeroing is the source's own: `volInit` multiplies the
/// solvent dielectric constant's temperature derivatives by zero rather than omitting them,
/// and `calcSolventdiElectricdn` returns `0.0` with its body commented out.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if `x` is not one entry per component or is
///   not a composition, a name is not in the databank, or `compressed_phase` is neither
///   root's spelling.
/// * [`azoth_core::AzothError::PropertyUnavailable`] if a component carries no heat-capacity
///   coefficients.
/// * [`azoth_core::AzothError::OutOfRange`] if `T` or `P` is not positive, or no volume root
///   exists above the mixture's covolume.
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn furst_electrolyte_mod2004_phase(
    components: &[String],
    T: ThermodynamicTemperature,
    P: Pressure,
    x: &[f64],
    compressed_phase: &str,
) -> Result<FurstElectrolyteMod2004PhaseResult> {
    let spec = &model_gen::FURST_ELECTROLYTE_MOD2004_PHASE_SPEC;
    let mut warnings = Vec::new();
    azoth_core::apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T.value),
            "P" => Some(P.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let n = components.len();
    if x.len() != n {
        return Err(azoth_core::AzothError::invalid_input(
            "x",
            format!("a mixture of {n} components has {} mole fractions", x.len()),
        ));
    }
    let sum: f64 = x.iter().sum();
    if (sum - 1.0).abs() > 1.0e-09 {
        return Err(azoth_core::AzothError::invalid_input(
            "x",
            format!(
                "the mole fractions sum to {sum}, not to one. Renormalising them here \
                 would make a composition error invisible in every number downstream, so \
                 it is refused instead"
            ),
        ));
    }

    let side = crate::cpa_phase::side_of(compressed_phase)?;
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, _) = crate::furst_electrolyte::furst_mod2004_mixture_of(
        &names,
        crate::furst_dielectric::MixingRule::default_for_the_model(),
        None,
    )?;
    let reduced = mixture.reduced_parameters(T, P)?;
    warnings.extend(reduced.warnings.iter().cloned());
    let state = mixture.phase_state(&reduced, x, side)?;

    Ok(FurstElectrolyteMod2004PhaseResult {
        z_factor: state.z,
        ln_phi: state.ln_phi,
        warnings,
    })
}
