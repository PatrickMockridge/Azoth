//! `eos.solid_fugacity` - a pure solid's fugacity coefficient, from tabulated properties.
//!
//! ```text
//! phi_solid = phi_liq(T, P) exp( -dH_fus/(R T) (1 - T/T_tp)
//!                              + dCp_SL/(R T) (T_tp - T)
//!                              - dCp_SL/R ln(T_tp/T)
//!                              - dV_sl (P_bar - 1)/(R T) )
//! ```
//!
//! Spec: `specs/calcs/eos/solid_fugacity.toml`, which carries the three input routes for the
//! heat-capacity difference and the measurement that pins the whole expression.
//!
//! NeqSim's `ComponentSolid.fugcoef2`. **The reference is a *liquid*** - `fugcoef2`
//! initialises it `PhaseType.LIQUID`, where the class's other entry (`fugcoef`) uses a gas at
//! the component's solid vapour pressure and is a different model - and it is the host's own
//! class, so `eos` is an input and not a label.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result, apply_checks};

use crate::Alpha;
use crate::Cubic;
use crate::mixture::{Component, Mixture, RootSide};
use crate::results::SolidFugacityResult;
use crate::spec_gen;

/// NeqSim's `R`, which `Component` states as a literal.
pub const R: f64 = 8.314_462_1;

/// One atmosphere in bar: the pressure the volume term is referred to in NeqSim's own code.
const REFERENCE_PRESSURE_BAR: f64 = 1.0;

/// A pure solid's fugacity coefficient at a state.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `eos` is not a cubic this reaches.
/// * [`AzothError::OutOfRange`] if any of the constants or the state is not positive.
#[allow(clippy::too_many_arguments, non_snake_case)]
pub fn solid_fugacity(
    heat_of_fusion: f64,
    triple_point_temperature: f64,
    delta_cp_sl: f64,
    delta_solid_volume: f64,
    tc: ThermodynamicTemperature,
    pc: Pressure,
    omega: f64,
    T: ThermodynamicTemperature,
    P: Pressure,
    eos: &str,
) -> Result<SolidFugacityResult> {
    let spec = &spec_gen::SOLID_FUGACITY_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "heat_of_fusion" => Some(heat_of_fusion),
            "triple_point_temperature" => Some(triple_point_temperature),
            "T" => Some(T.value),
            "P" => Some(P.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let cubic = match eos {
        "pr" => Cubic::Pr,
        "srk" => Cubic::Srk,
        other => {
            return Err(AzothError::invalid_input(
                "eos",
                format!(
                    "`{other}` is not a cubic this reaches: NeqSim builds the reference liquid \
                     from the *host phase's* class, and only `pr` and `srk` are ported"
                ),
            ));
        }
    };

    // **The reference liquid**: this component alone, on the fluid's own cubic, on the liquid
    // root.
    let alpha = match cubic {
        Cubic::Pr => Alpha::Pr,
        _ => Alpha::Srk,
    };
    let reference = Mixture::new(vec![Component::new(tc, pc, omega)?], vec![0.0])?
        .with_cubic(cubic)
        .with_alpha(alpha);
    let reduced = reference.reduced_parameters(T, P)?;
    let state = reference.phase_state(&reduced, &[1.0], RootSide::Liquid)?;
    let phi_liquid = state.ln_phi[0].exp();

    let pressure_bar = P.value / 1.0e5;
    let fusion = -heat_of_fusion / (R * T.value) * (1.0 - T.value / triple_point_temperature);
    let heat_capacity = delta_cp_sl / (R * T.value) * (triple_point_temperature - T.value)
        - delta_cp_sl / R * (triple_point_temperature / T.value).ln();
    let volume = -delta_solid_volume * (pressure_bar - REFERENCE_PRESSURE_BAR) / (R * T.value);

    Ok(SolidFugacityResult {
        fugacity_coefficient: phi_liquid * (fusion + heat_capacity + volume).exp(),
        warnings,
    })
}
