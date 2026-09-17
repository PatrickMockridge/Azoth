//! `eos.water_phase` - the IAPWS-IF97 water phase state at a temperature and pressure.
//!
//! ```text
//! select the region by T vs Tsat(p), evaluate the region's Gibbs function, and
//! rearrange its derivatives into Z, u, h, s, cv, cp and g, mass-to-molar.
//! ```
//!
//! Spec: `specs/models/eos/water_phase.toml`. Pure water, so there is no composition
//! input, and no density solve: IF97 is a Gibbs formulation in `(p, T)`, so the region
//! is chosen directly rather than by a root of `P(rho)`.

use azoth_core::units::{
    Pressure, ThermodynamicTemperature, joules_per_mole, joules_per_mole_kelvin,
};
use azoth_core::{Result, apply_checks};

use crate::iapws_if97;
use crate::model_gen;
use crate::results::WaterPhaseResult;

/// The IAPWS-IF97 water phase state at a temperature and pressure.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `T` or `P` is not positive.
pub fn water_phase(t: ThermodynamicTemperature, p: Pressure) -> Result<WaterPhaseResult> {
    let spec = &model_gen::WATER_PHASE_SPEC;
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
    let p_mpa = p.value / 1.0e6;
    let props = iapws_if97::properties(p_mpa, tk);

    // IF97 is mass-based (kJ/kg); the model is molar SI (J/mol), so a specific
    // property becomes molar on multiplying by the molar mass and by the kJ -> J
    // factor, the single place the two conventions meet.
    let molar = iapws_if97::MOLAR_MASS * 1.0e3;

    Ok(WaterPhaseResult {
        z_factor: props.z,
        u: joules_per_mole(props.u * molar),
        h: joules_per_mole(props.h * molar),
        s: joules_per_mole_kelvin(props.s * molar),
        cv: joules_per_mole_kelvin(props.cv * molar),
        cp: joules_per_mole_kelvin(props.cp * molar),
        g: joules_per_mole(props.g * molar),
        warnings,
    })
}
