//! The Peng-Robinson CPA phase state.
//!
//! `PhasePrCPA` as a *model*: the cubic's attraction and covolume replaced by each
//! component's fitted `aCPA_PR`/`bCPA_PR`, mixed with the `cpakij_PR` column, and the
//! Wertheim association contribution added to the residual Helmholtz energy.
//!
//! **The first model to read the PR family's CPA columns**, and the reason the two families
//! are separate selections rather than one converted into the other: water's `kappa_AB` is
//! 0.0692 for SRK against 0.046473789 for PR, and its fitted covolume is 1.4515 against
//! 1.456360879 - so a model that read the wrong set would be a different fluid.
//!
//! **There is no NeqSim oracle for this model, and that is a finding rather than an
//! omission.** Against the pinned 3.20.0 jar, `SystemPrCPA` builds `ComponentSrkCPA`
//! components that carry their sites, and its phase never sums them: `PhaseSrkCPA` does that
//! in its init path and `PhasePrCPA` has the field and the setter and no block that sets it,
//! so water reports four sites and the phase's own total is zero. Every association term is
//! then computed over nothing and the flash is a Peng-Robinson run wearing the name.
//! `validation/neqsim/PrCpaFlash.java` prints both counts for exactly that reason. The
//! finding is recorded with the tranche's other two; this model's case carries its own
//! numbers and says so.

use azoth_core::Result;
use azoth_core::units::{Pressure, ThermodynamicTemperature};

use crate::mixture::{Mixture, RootSide};
use crate::model_gen;
use crate::results::PrCpaPhaseResult;

/// One CPA phase's state at a temperature, pressure and composition.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if `z` is not one entry per component, is not
///   a composition, or a name is not in the databank.
/// * [`azoth_core::AzothError::OutOfRange`] if `T` or `P` is not positive, or no volume root
///   exists above the mixture's covolume.
pub fn pr_cpa_phase(
    components: &[String],
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    compressed_phase: &str,
) -> Result<PrCpaPhaseResult> {
    let state = crate::cpa_phase::cpa_phase(
        &model_gen::PR_CPA_PHASE_SPEC,
        components,
        crate::Cubic::Pr,
        t,
        p,
        z,
        compressed_phase,
    )?;
    Ok(PrCpaPhaseResult {
        z_factor: state.z_factor,
        ln_phi: state.ln_phi,
        h_res: azoth_core::units::joules_per_mole(state.h_res),
        s_res: azoth_core::units::joules_per_mole_kelvin(state.s_res),
        warnings: state.warnings,
    })
}

/// A mixture's phase state without the databank, for a caller that has one already.
///
/// # Errors
/// As [`pr_cpa_phase`], without the names.
pub fn pr_cpa_phase_of(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    side: RootSide,
) -> Result<PrCpaPhaseResult> {
    let state = crate::cpa_phase::phase_state_of(mixture, t, p, z, side)?;
    Ok(PrCpaPhaseResult {
        z_factor: state.z_factor,
        ln_phi: state.ln_phi,
        h_res: azoth_core::units::joules_per_mole(state.h_res),
        s_res: azoth_core::units::joules_per_mole_kelvin(state.s_res),
        warnings: state.warnings,
    })
}
