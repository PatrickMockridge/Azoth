//! `eos.eos_cg_phase` - the EOS-CG phase state at a temperature, pressure and composition.
//!
//! ```text
//! map the component names to the EOS-CG order, solve the density in logarithmic
//! volume, then the Helmholtz derivatives give Z, u, h, s, cv, cp and g.
//! ```
//!
//! Spec: `specs/models/eos/eos_cg_phase.toml`. The composition is the EOS-CG component
//! set, mapped by name to the model's fixed 28-component order.

use azoth_core::units::{
    Pressure, ThermodynamicTemperature, joules_per_mole, joules_per_mole_kelvin,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::eos_cg;
use crate::model_gen;
use crate::results::EosCgPhaseResult;

/// The EOS-CG phase state of a combustion-gas mixture.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
/// * [`AzothError::InvalidInput`] if a component name is not an EOS-CG component, the
///   fractions do not match the names in length, or the fractions do not sum to one.
pub fn eos_cg_phase(
    components: &[String],
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<EosCgPhaseResult> {
    let spec = &model_gen::EOS_CG_PHASE_SPEC;
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

    if components.len() != z.len() {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "{} components but {} fractions; the two must match",
                components.len(),
                z.len()
            ),
        ));
    }

    let mut x = [0.0; eos_cg::NCOMP + 1];
    for (name, fraction) in components.iter().zip(z) {
        let index = eos_cg::component_index(name).ok_or_else(|| {
            AzothError::invalid_input("components", format!("unknown EOS-CG component {name:?}"))
        })?;
        x[index] += fraction;
    }
    let sum: f64 = x.iter().skip(1).sum();
    if (sum - 1.0).abs() > 1.0e-9 {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "the composition sums to {sum}, not to one; renormalising here would hide a caller's error"
            ),
        ));
    }

    let p_kpa = p.value / 1000.0;
    let density = eos_cg::solve_density(t.value, p_kpa, &x);
    let props = eos_cg::properties(t.value, density, &x);

    Ok(EosCgPhaseResult {
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
