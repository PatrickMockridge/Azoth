//! `eos.wax_solid_fugacity` - a wax cut's solid fugacity coefficient.
//!
//! ```text
//! phi_wax = phi_liq(T, P) exp( -dH_fus/(R T) (1 - T/T_tp)
//!                            + dCp_SL/R (T_tp/T - 1 - ln(T_tp/T))
//!                            - (v_liq - v_sol)(P - P_ref)/(R T) )
//! ```
//!
//! Spec: `specs/calcs/eos/wax_solid_fugacity.toml`, which carries the fits, the one-bar
//! reference pressure, and the measurement that pins all of it.
//!
//! NeqSim's `ComponentWax.fugcoef2`, the default `PhaseWax` component model. **The coefficient
//! is a pure-component quantity**: `SolidFug = x f_liq exp(...)` and the reported coefficient
//! is `SolidFug/(P x)`, so the mole fraction cancels and no composition survives.
//!
//! # The units, which are where this goes wrong
//!
//! `v_liq` is the reference liquid's molar volume in **m³/mol**, and the pressure term is
//! evaluated against **one bar** - `P_ref = 1.0` in NeqSim's bar-valued code, so `1.0e5` Pa
//! here. Getting either wrong is not a rounding: with the pressure in bar against an
//! m³/mol volume the term comes out `1e5` too small, and the coefficient is out by a per cent
//! that grows with the cut's molar mass - which is what this port did first, and what the
//! reference probe's separate capture of the two halves caught.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result, apply_checks};

use crate::Alpha;
use crate::Cubic;
use crate::mixture::{Component, Mixture, RootSide};
use crate::results::WaxSolidFugacityResult;
use crate::spec_gen;

/// NeqSim's `R`, which the class's compiled form states as a literal.
pub const R: f64 = 8.314_462_1;

/// One bar, the pressure the volume term is referred to.
const REFERENCE_PRESSURE: f64 = 1.0e5;

/// The solid's molar volume over the liquid's, which is NeqSim's own shortcut.
const SOLID_VOLUME_RATIO: f64 = 0.9;

/// A wax cut's solid fugacity coefficient at a state.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if any of the cut's constants or the state is not positive.
#[allow(clippy::too_many_arguments)]
pub fn wax_solid_fugacity(
    molar_mass: f64,
    tc: ThermodynamicTemperature,
    pc: Pressure,
    omega: f64,
    heat_of_fusion: f64,
    triple_point_temperature: ThermodynamicTemperature,
    t: ThermodynamicTemperature,
    p: Pressure,
    eos: &str,
) -> Result<WaxSolidFugacityResult> {
    let spec = &spec_gen::WAX_SOLID_FUGACITY_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "molar_mass" => Some(molar_mass),
            "tc" => Some(tc.value),
            "pc" => Some(pc.value),
            "heat_of_fusion" => Some(heat_of_fusion),
            "triple_point_temperature" => Some(triple_point_temperature.value),
            "T" => Some(t.value),
            "P" => Some(p.value),
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
                     from the *host phase's* class, and only `pr` and `srk` have a wax route"
                ),
            ));
        }
    };

    // **The reference liquid: this component alone, on the fluid's own cubic, on the liquid
    // root.** NeqSim's `setSolidRefFluidPhase` clones the host phase's class and adds one
    // component, so the alpha is that cubic's own Soave form at the cut's acentric factor.
    let alpha = match cubic {
        Cubic::Pr => Alpha::Pr,
        _ => Alpha::Srk,
    };
    let reference = Mixture::new(
        vec![
            Component::new(tc, pc, omega)
                .map_err(|_| AzothError::out_of_range("tc/pc", tc.value, "the cut's constants"))?
                .with_molar_mass(Some(molar_mass)),
        ],
        vec![0.0],
    )?
    .with_cubic(cubic)
    .with_alpha(alpha);
    let reduced = reference.reduced_parameters(t, p)?;
    let state = reference.phase_state(&reduced, &[1.0], RootSide::Liquid)?;
    let phi_liquid = state.ln_phi[0].exp();
    let v_liquid = state.z * R * t.value / p.value;

    let v_solid = SOLID_VOLUME_RATIO * v_liquid;
    let pressure_term = -(v_liquid - v_solid) * (p.value - REFERENCE_PRESSURE) / R / t.value;
    let molar_mass_grams = molar_mass * 1000.0;
    let delta_cp_sl = (0.3033 * molar_mass_grams - 4.635e-4 * molar_mass_grams * t.value) * 4.184;
    let triple_point_ratio = triple_point_temperature.value / t.value;
    let heat_capacity_term = delta_cp_sl / R * (triple_point_ratio - 1.0 - triple_point_ratio.ln());
    let fusion_term =
        -heat_of_fusion / (R * t.value) * (1.0 - t.value / triple_point_temperature.value);

    let coefficient = phi_liquid * (fusion_term + heat_capacity_term + pressure_term).exp();

    Ok(WaxSolidFugacityResult {
        fugacity_coefficient: coefficient,
        warnings,
    })
}
