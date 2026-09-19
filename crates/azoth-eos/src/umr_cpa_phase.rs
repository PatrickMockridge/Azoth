//! The UMR-CPA phase state.
//!
//! `PhaseUMRCPA` as a *model*: a Peng-Robinson cubic whose attraction is mixed by the UMR
//! universal rule over the UNIFAC-UMR-PRU activity coefficients, with the Wertheim
//! association contribution added to the residual Helmholtz energy.
//!
//! **The only model in this library whose attraction is not mixed by an interaction
//! matrix.** `PhaseUMRCPA extends PhasePrEos`, and its pressure is `P_PR(UMR mixing) +
//! P_assoc`: the cubic's `A` is `n B R T alpha_mix` with `alpha_mix = sum_i x_i
//! (a_i^T/(b_i R T) + hwfc ln gamma_i)`, so a `kij` column is not read at all. What the
//! model needs instead is each component's UNIFAC group decomposition, which
//! [`crate::databank::umr_cpa_mixture_of`] resolves.
//!
//! **Three parameter sets, not one.** A component's `a` and `b` are the `UMRCPA_a0` and
//! `UMRCPA_b` columns where it carries them - a third fitted set beside the SRK and PR
//! families - its bounding volume and energy are `UMRCPA_assocVolume` and
//! `UMRCPA_assocEnergy`, and its attractive term is the five-parameter Mathias-Copeman form
//! seeded with `UMRCPA_MC1..5`. Water is the check: `kappa_AB` 0.125 and `eps` 14177 J/mol
//! here against the PR family's, so a model reading the wrong set would be a different
//! fluid that still converged.
//!
//! **No flash.** `eos.pt_flash`'s Jacobian is `d ln phi / d n`, which
//! [`crate::mixture::Mixture::phase_derivatives`] refuses for every excess-Gibbs rule
//! because that derivative is the excess Gibbs energy's own Hessian. The boundary is the
//! one `MixingRule::HuronVidal` and `MixingRule::WongSandler` already sit behind.

use azoth_core::Result;
use azoth_core::units::{Pressure, ThermodynamicTemperature};

use crate::model_gen;
use crate::results::UmrCpaPhaseResult;

/// One UMR-CPA phase's state at a temperature, pressure and composition.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if `z` is not one entry per component, is
///   not a composition, a name is not in the databank, or `compressed_phase` is neither
///   root's spelling.
/// * [`azoth_core::AzothError::PropertyUnavailable`] if a component carries no
///   `UMRCPA_MC1..5` set or no `UNIFACcompUMRPRU` group decomposition.
/// * [`azoth_core::AzothError::OutOfRange`] if `T` or `P` is not positive, or no volume
///   root exists above the mixture's covolume.
pub fn umr_cpa_phase(
    components: &[String],
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    compressed_phase: &str,
) -> Result<UmrCpaPhaseResult> {
    let state = crate::cpa_phase::cpa_phase_with(
        &model_gen::UMR_CPA_PHASE_SPEC,
        components,
        t,
        p,
        z,
        compressed_phase,
        |names| crate::databank::umr_cpa_mixture_of(names).map(|(mixture, _)| mixture),
    )?;
    Ok(UmrCpaPhaseResult {
        z_factor: state.z_factor,
        ln_phi: state.ln_phi,
        h_res: azoth_core::units::joules_per_mole(state.h_res),
        s_res: azoth_core::units::joules_per_mole_kelvin(state.s_res),
        warnings: state.warnings,
    })
}
