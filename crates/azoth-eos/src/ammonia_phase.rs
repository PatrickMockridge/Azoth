//! `eos.ammonia_phase` - the ammonia reference phase state at a temperature and pressure.
//!
//! ```text
//! solve P(rho) = P for the molar density (gas root), then the Helmholtz derivatives
//! give Z, u, h, s, cv, cp and g.
//! ```
//!
//! Spec: `specs/models/eos/ammonia_phase.toml`, which carries the assembly and what the
//! caller is responsible for. Pure ammonia, so there is no composition input; the
//! density solve is Newton from the ideal-gas guess, which finds the gas-like root.

use azoth_core::units::{
    Pressure, ThermodynamicTemperature, joules_per_mole, joules_per_mole_kelvin,
};
use azoth_core::{Result, apply_checks};

use crate::ammonia;
use crate::model_gen;
use crate::results::AmmoniaPhaseResult;

/// The ammonia phase state at a temperature and pressure.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `T` or `P` is not positive.
pub fn ammonia_phase(t: ThermodynamicTemperature, p: Pressure) -> Result<AmmoniaPhaseResult> {
    let spec = &model_gen::AMMONIA_PHASE_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "P" => Some(p.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let tk = t.value;
    let p_pa = p.value;
    let rho = ammonia::solve_density(tk, p_pa, false);
    let props = ammonia::properties(tk, rho);

    Ok(AmmoniaPhaseResult {
        z_factor: props.z,
        u: joules_per_mole(props.u),
        h: joules_per_mole(props.h),
        s: joules_per_mole_kelvin(props.s),
        cv: joules_per_mole_kelvin(props.cv),
        cp: joules_per_mole_kelvin(props.cp),
        g: joules_per_mole(props.g),
        warnings,
    })
}
