//! `eos.electrolyte_cpa_statoil_phase` - the phase state of an electrolyte-CPA fluid.
//!
//! Spec: `specs/models/eos/electrolyte_cpa_statoil_phase.toml`. A *direct* model: the root
//! is solved by the shared cubic-plus-pressure iteration and there is no outer loop, so no
//! `algorithm` block.
//!
//! NeqSim's `SystemElectrolyteCPAstatoil`, whose phase is `PhaseElectrolyteCPAstatoil` and
//! whose components are `ComponentElectrolyteCPAstatoil`. It is `SystemFurstElectrolyteEos`
//! with the Wertheim association added - the fluid is [`crate::furst_electrolyte`]'s, and
//! [`crate::furst_electrolyte::electrolyte_cpa_statoil_mixture_of`] is the resolver.
//!
//! # What this model does not have
//!
//! **No flash.** `eos.pt_flash`'s Jacobian is `d ln phi / d n`, which the Huron-Vidal mixing
//! rule this system installs refuses - the same boundary
//! [`crate::furst_electrolyte_phase`] sits behind. The cases record phase states at
//! compositions NeqSim reports.
//!
//! **And no `h_res` or `s_res`.** The association's departure is `A/(RT) - T d(A/RT)/dT` and
//! the electrolyte terms' is not assembled, so a model reporting one would be reporting part
//! of it.

use azoth_core::Result;
use azoth_core::units::{Pressure, ThermodynamicTemperature};

use crate::model_gen;
use crate::results::ElectrolyteCpaStatoilPhaseResult;

/// One electrolyte-CPA phase's state at a temperature, pressure and composition.
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
pub fn electrolyte_cpa_statoil_phase(
    components: &[String],
    T: ThermodynamicTemperature,
    P: Pressure,
    x: &[f64],
    compressed_phase: &str,
) -> Result<ElectrolyteCpaStatoilPhaseResult> {
    let spec = &model_gen::ELECTROLYTE_CPA_STATOIL_PHASE_SPEC;
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
    let (mixture, _) = crate::furst_electrolyte::electrolyte_cpa_statoil_mixture_of(&names, None)?;
    let reduced = mixture.reduced_parameters(T, P)?;
    warnings.extend(reduced.warnings.iter().cloned());
    let state = mixture.phase_state(&reduced, x, side)?;

    Ok(ElectrolyteCpaStatoilPhaseResult {
        z_factor: state.z,
        ln_phi: state.ln_phi,
        warnings,
    })
}
