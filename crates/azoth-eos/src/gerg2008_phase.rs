//! `eos.gerg2008_phase` - the GERG-2008 phase state at a temperature, pressure and composition.
//!
//! ```text
//! map the component names to the GERG-2008 order, solve the density in logarithmic
//! volume, then the Helmholtz derivatives give Z, u, h, s, cv, cp and g.
//! ```
//!
//! Spec: `specs/models/eos/gerg2008_phase.toml`. The composition is the GERG-2008
//! component set, mapped by name to the model's fixed 21-component order.

use azoth_core::units::{
    Pressure, ThermodynamicTemperature, joules_per_mole, joules_per_mole_kelvin,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::gerg2008;
use crate::model_gen;
use crate::results::Gerg2008PhaseResult;

/// The GERG-2008 phase state of a natural-gas mixture.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
/// * [`AzothError::InvalidInput`] if a component name is not a GERG-2008 component, the
///   fractions do not match the names in length, or the fractions do not sum to one.
pub fn gerg2008_phase(
    components: &[String],
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<Gerg2008PhaseResult> {
    let spec = &model_gen::GERG2008_PHASE_SPEC;
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

    let mut x = [0.0; gerg2008::NCOMP + 1];
    for (name, fraction) in components.iter().zip(z) {
        let index = gerg2008::component_index(name).ok_or_else(|| {
            AzothError::invalid_input(
                "components",
                format!("unknown GERG-2008 component {name:?}"),
            )
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
    let density = gerg2008::solve_density(t.value, p_kpa, &x);
    let props = gerg2008::properties(t.value, density, &x);

    Ok(Gerg2008PhaseResult {
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
