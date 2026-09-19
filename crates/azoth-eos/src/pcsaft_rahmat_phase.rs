//! `eos.pcsaft_rahmat_phase` - the PC-SAFT phase state.
//!
//! ```text
//! solve P_calc(v) = P   then   Z = P v/(R T)   and   ln phi_i = d(nF)/dn_i - ln Z
//! ```
//!
//! Gross and Sadowski's perturbed-chain statistical associating fluid theory without its
//! association term, as NeqSim's `PhasePCSAFTRahmat` carries it. The layers are
//! [`crate::pcsaft`]'s and the solve is [`crate::pcsaft_phase`]'s, checked one at a time
//! against `validation/neqsim/PcsaftProbe.java`; what is here is the assembly and the
//! names.
//!
//! # Why the class name is in the model's
//!
//! **`PhasePCSAFTRahmat` is the PC-SAFT that NeqSim 3.20.0 runs.** `SystemPCSAFT` builds
//! it, and the sibling it extends - `PhasePCSAFT` - is never instantiated anywhere in the
//! checkout: there is no `new PhasePCSAFT()` in it, and the only two subclasses are this
//! one and the associating `PhasePCSAFTa`. The two non-associating classes differ in
//! `dF_HC_SAFTdVdV`, whose chain term is short a factor of `dnSAFT/dV` in the base class,
//! and in the volume solve's seed - a denominator and a seed, neither of which moves the
//! fixed point a converged solve reaches. So they answer the same volumes and the same
//! fugacity coefficients, and one model names the class a caller can actually reach
//! rather than the shell around it. `validation/neqsim/PcsaftDVdVProbe.java` measures the
//! difference in the derivative.

use azoth_core::units::{Pressure, ThermodynamicTemperature, cubic_meters_per_mole};
use azoth_core::{AzothError, Result};

use crate::mixture::RootSide;
use crate::pcsaft::PcsaftComponent;
use crate::pcsaft_phase::{molar_volume, parameters_of, side_of};
use crate::results::PcsaftRahmatPhaseResult;

/// The PC-SAFT phase state at a temperature, a pressure and a composition.
///
/// The names cross unresolved and are looked up here, so the two-kernel comparison covers
/// the resolution as well as the arithmetic - which for this model means the
/// `KIJPCSAFT` column and the zero-means-absent convention on the parameters.
///
/// # Errors
/// * [`AzothError::InvalidInput`] as [`parameters_of`], or if `z` is not one entry per
///   component or does not sum to one.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive, or the wanted branch has no
///   root at this state.
pub fn pcsaft_rahmat_phase(
    components: &[String],
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    compressed_phase: &str,
) -> Result<PcsaftRahmatPhaseResult> {
    let spec = &crate::model_gen::PCSAFT_RAHMAT_PHASE_SPEC;
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

    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (parameters, kij) = parameters_of(&names)?;
    let side = side_of(compressed_phase)?;
    let mut result = phase_state_of(&parameters, &kij, t, p, z, side)?;
    warnings.extend(result.warnings);
    result.warnings = warnings;
    Ok(result)
}

/// The same state, for a caller that resolved the fluid itself.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `z` is not one entry per component or does not sum to
///   one.
/// * [`AzothError::OutOfRange`] as [`molar_volume`] and
///   [`crate::pcsaft::ln_fugacity_coefficients`].
pub fn phase_state_of(
    components: &[PcsaftComponent],
    kij: &[f64],
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    side: RootSide,
) -> Result<PcsaftRahmatPhaseResult> {
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

    let solved = molar_volume(components, kij, z, t.value, p.value, side)?;
    Ok(PcsaftRahmatPhaseResult {
        z_factor: solved.z,
        ln_phi: crate::pcsaft::ln_fugacity_coefficients(components, kij, z, t.value, solved.v)?,
        v: cubic_meters_per_mole(solved.v),
        warnings: Vec::new(),
    })
}
