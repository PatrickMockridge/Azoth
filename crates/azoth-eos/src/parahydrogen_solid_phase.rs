//! `eos.parahydrogen_solid_phase` - the solid para-hydrogen phase state at a temperature
//! and pressure.
//!
//! ```text
//! solve -a_V = P for the molar volume, then the second-order Helmholtz derivatives
//! give Z, u, h, s, cv, cp and g.
//! ```
//!
//! Spec: `specs/models/eos/parahydrogen_solid_phase.toml`. Pure solid para-hydrogen, so
//! there is no composition input; the volume solve is bracketing plus bisection/Newton in
//! logarithmic volume, which finds the solid root.

use azoth_core::units::{
    Pressure, ThermodynamicTemperature, joules_per_mole, joules_per_mole_kelvin,
};
use azoth_core::{Result, apply_checks};

use crate::model_gen;
use crate::parahydrogen_solid;
use crate::results::ParahydrogenSolidPhaseResult;

/// The solid para-hydrogen phase state at a temperature and pressure.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `T` or `P` is outside the published range.
pub fn parahydrogen_solid_phase(
    t: ThermodynamicTemperature,
    p: Pressure,
) -> Result<ParahydrogenSolidPhaseResult> {
    let spec = &model_gen::PARAHYDROGEN_SOLID_PHASE_SPEC;
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

    let state = parahydrogen_solid::properties(t.value, p.value);

    Ok(ParahydrogenSolidPhaseResult {
        z_factor: p.value * state.v / (parahydrogen_solid::R * t.value),
        u: joules_per_mole(state.u),
        h: joules_per_mole(state.h),
        s: joules_per_mole_kelvin(state.s),
        cv: joules_per_mole_kelvin(state.cv),
        cp: joules_per_mole_kelvin(state.cp),
        g: joules_per_mole(state.g),
        warnings,
    })
}
