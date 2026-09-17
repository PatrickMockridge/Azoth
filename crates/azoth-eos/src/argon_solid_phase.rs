//! `eos.argon_solid_phase` - the solid argon phase state at a temperature and pressure.
//!
//! ```text
//! solve -a_V = P for the molar volume, then the second-order Helmholtz derivatives
//! give Z, u, h, s, cv, cp and g.
//! ```
//!
//! Spec: `specs/models/eos/argon_solid_phase.toml`. Pure solid argon, so there is no
//! composition input; the volume solve is bracketing plus bisection/Newton in logarithmic
//! volume, which finds the solid root.

use azoth_core::units::{
    Pressure, ThermodynamicTemperature, joules_per_mole, joules_per_mole_kelvin,
};
use azoth_core::{Result, apply_checks};

use crate::argon_solid;
use crate::model_gen;
use crate::results::ArgonSolidPhaseResult;

/// The solid argon phase state at a temperature and pressure.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `T` or `P` is outside the published range.
pub fn argon_solid_phase(
    t: ThermodynamicTemperature,
    p: Pressure,
) -> Result<ArgonSolidPhaseResult> {
    let spec = &model_gen::ARGON_SOLID_PHASE_SPEC;
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

    let state = argon_solid::properties(t.value, p.value);

    Ok(ArgonSolidPhaseResult {
        z_factor: p.value * state.v / (argon_solid::R * t.value),
        u: joules_per_mole(state.u),
        h: joules_per_mole(state.h),
        s: joules_per_mole_kelvin(state.s),
        cv: joules_per_mole_kelvin(state.cv),
        cp: joules_per_mole_kelvin(state.cp),
        g: joules_per_mole(state.g),
        warnings,
    })
}
