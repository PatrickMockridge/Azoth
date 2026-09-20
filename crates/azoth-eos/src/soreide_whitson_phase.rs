//! `eos.soreide_whitson_phase` - the phase state of a Soreide-Whitson fluid.
//!
//! Spec: `specs/models/eos/soreide_whitson_phase.toml`. A *direct* model: the root is the
//! cubic's, chosen by ordering, so no iteration and no `algorithm` block.
//!
//! NeqSim's `SystemSoreideWhitson`, whose distance from Peng-Robinson 1978 is exactly two
//! things:
//!
//! ```text
//! kij(i, water) = the salinity correlation, in a water-rich phase
//! alpha(water)  = AttractiveTermSoreideWhitson, salinity-dependent
//! ```
//!
//! Everything else - every other pair, every other component, every phase whose water is
//! below a mole fraction of `0.8` - is PR78. The two are independent paths and this model
//! is where they meet: [`crate::soreide_whitson_alpha`] and
//! [`crate::mixing_rule::MixingRule::phase_kij`] are the pieces, and this is the assembly.
//!
//! # What this model does not have
//!
//! **No flash.** `eos.pt_flash`'s Jacobian is `d ln phi / d n`, which
//! [`crate::mixture::Mixture::phase_derivatives`] refuses for this rule because its
//! interaction matrix moves with the composition being differentiated against. The
//! boundary is the one `MixingRule::HuronVidal` and `MixingRule::WongSandler` sit behind.
//!
//! **And the salinity is the caller's.** NeqSim's `calcSalinity` recomputes the phase's
//! molality from the flashed aqueous phase and pushes it into the alpha through
//! `setSalinityFromPhase` - a call its own `getClass().getName()` test never reaches,
//! which is the defect reported as the alpha issue. A caller here states the brine, and
//! the model states what it read.

use azoth_core::Result;
use azoth_core::units::{Pressure, ThermodynamicTemperature};

use crate::model_gen;
use crate::results::SoreideWhitsonPhaseResult;

/// One Soreide-Whitson phase's state at a temperature, pressure and composition.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if `x` is not one entry per component or is
///   not a composition, a name is not in the databank, an ion is named, the salinity is
///   negative, or `compressed_phase` is neither root's spelling.
/// * [`azoth_core::AzothError::PropertyUnavailable`] if a component carries no
///   heat-capacity coefficients.
/// * [`azoth_core::AzothError::OutOfRange`] if `T` or `P` is not positive, or no volume
///   root exists above the mixture's covolume.
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn soreide_whitson_phase(
    components: &[String],
    T: ThermodynamicTemperature,
    P: Pressure,
    x: &[f64],
    salinity: f64,
    compressed_phase: &str,
) -> Result<SoreideWhitsonPhaseResult> {
    let spec = &model_gen::SOREIDE_WHITSON_PHASE_SPEC;
    let mut warnings = Vec::new();
    azoth_core::apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T.value),
            "P" => Some(P.value),
            "salinity" => Some(salinity),
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
    let (mixture, _) = crate::databank::soreide_whitson_mixture_of(&names, salinity, None)?;
    let reduced = mixture.reduced_parameters(T, P)?;
    warnings.extend(reduced.warnings.iter().cloned());
    let state = mixture.phase_state(&reduced, x, side)?;

    Ok(SoreideWhitsonPhaseResult {
        z_factor: state.z,
        ln_phi: state.ln_phi,
        warnings,
    })
}
