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
/// `compressed_phase` selects the **root**, and below the critical temperature the two are
/// different states at the same temperature and pressure - the spec declares it, because a
/// caller cannot supply a compressibility factor instead without solving this model first.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive, or if the dense root was
///   asked for at a state whose dense root is outside the range the equation is fitted to.
/// * [`AzothError::InvalidInput`] if `hydrogen_type` is not `normal`, `para` or `ortho`, or
///   `compressed_phase` is neither `liquid` nor `vapour`.
pub fn hydrogen_phase(
    t: ThermodynamicTemperature,
    p: Pressure,
    hydrogen_type: &str,
    compressed_phase: &str,
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
    let rho = match compressed_phase {
        "vapour" => leachman::solve_density(tk, p_pa, ht),
        "liquid" => leachman::solve_density_dense(tk, p_pa, ht).ok_or_else(|| {
            AzothError::out_of_range(
                "compressed_phase",
                p_pa,
                "the dense root was asked for, and this state has none within the density \
                 range the equation is fitted to: the isotherm does not cross the pressure \
                 there below the ceiling the solve brackets from. `vapour` is the root such a \
                 state is on.",
            )
        })?,
        other => {
            return Err(AzothError::invalid_input(
                "compressed_phase",
                format!("`{other}`; expected `liquid` or `vapour`"),
            ));
        }
    };
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
