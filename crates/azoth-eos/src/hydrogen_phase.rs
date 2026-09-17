//! `eos.hydrogen_phase` - the Leachman hydrogen phase state at a temperature and pressure.
//!
//! ```text
//! solve P(rho) = P for the molar density (gas root), then the Helmholtz derivatives
//! give Z, u, h, s, cv, cp and g.
//! ```
//!
//! Spec: `specs/models/eos/hydrogen_phase.toml`. Pure hydrogen, so there is no composition
//! input; the spin-isomer is a boundary-only choice (normal, para or ortho) that does not
//! change the declared inputs.

use azoth_core::units::{
    Pressure, ThermodynamicTemperature, joules_per_mole, joules_per_mole_kelvin,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::leachman::{self, HydrogenType};
use crate::model_gen;
use crate::results::HydrogenPhaseResult;

/// The Leachman hydrogen phase state at a temperature and pressure.
///
/// `hydrogen_type` names the spin-isomer the equation is parameterised for. It is
/// boundary-only: the spec declares `T` and `P`, and the isomer is a caller's choice that
/// the bridge carries as a string.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
/// * [`AzothError::InvalidInput`] if `hydrogen_type` is not `normal`, `para` or `ortho`.
pub fn hydrogen_phase(
    t: ThermodynamicTemperature,
    p: Pressure,
    hydrogen_type: &str,
) -> Result<HydrogenPhaseResult> {
    let spec = &model_gen::HYDROGEN_PHASE_SPEC;
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

    let ht: HydrogenType = hydrogen_type
        .parse()
        .map_err(|e| AzothError::invalid_input("hydrogen_type", e))?;

    let tk = t.value;
    let p_pa = p.value;
    let rho = leachman::solve_density(tk, p_pa, ht);
    let props = leachman::properties(tk, rho, ht);

    Ok(HydrogenPhaseResult {
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
