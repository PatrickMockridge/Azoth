//! The bottom end: `Reboiler`, which is a stage with a boilup ratio.

use azoth_core::units::{Power, Pressure, ThermodynamicTemperature, watts};
use azoth_core::{AzothError, Result};
use azoth_eos::pv_reflux_flash::{RefluxPhase, pv_reflux_flash};

use crate::kernels::mixer;
use crate::stream::Stream;

use super::phase_fractions;
use super::tray::{self, SideDraws, TrayOutcome};

/// How the reboiler is specified.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ReboilerMode {
    /// No ratio: the `super.run` branch, which is the stage's own flash. An *equilibrium*
    /// reboiler, and the mode `ReboilerMode.EQUILIBRIUM` names.
    Equilibrium,
    /// `PVrefluxflash(ratio, 1)`: the boilup `V/B`, at the temperature that gives it.
    /// `ReboilerMode.VAPOR_BOILUP_RATIO`.
    VaporBoilupRatio(f64),
}

/// What one end hands the column.
///
/// **The two products are named for the column and not for the class.** NeqSim's reboiler
/// publishes a gas and a liquid and its condenser publishes a *distillate* through the gas
/// getter and a *reflux* through the liquid one; the names here are the column's, because
/// that is what the caller is assembling.
#[derive(Debug, Clone)]
pub struct ReboilerOutcome {
    /// The end's temperature, which is its flash's.
    pub temperature: ThermodynamicTemperature,
    /// The end's pressure.
    pub pressure: Pressure,
    /// The vapour rising into the column, or `None` where the flash found none.
    pub vapour: Option<Stream>,
    /// The liquid leaving the column, or `None` where the flash found none.
    pub liquid: Option<Stream>,
    /// The heat the end needs, W. Positive for a reboiler: it is the enthalpy the boilup
    /// carries out less the enthalpy the downcomer brought in.
    pub duty: Power,
}

/// The column's bottom end.
///
/// `Reboiler.run` has two branches and no more: with no ratio set it is `super.run`, the
/// stage's own flash, and with one it is `PVrefluxflash(ratio, 1)`. That flash is a
/// **temperature** search for the state whose *liquid* fraction satisfies `ratio = 1/beta_L -
/// 1`, which makes the ratio the boilup `V/B`; the condenser's is the same search on the
/// other phase, which is why the two ratios are reciprocals and why the choice is not
/// cosmetic.
///
/// **The duty is the outlet enthalpy less the inlets'**, in W, which is the class's own
/// `getMaterialOutletEnthalpy() - calcMixStreamEnthalpy0()`. On the equilibrium branch it is
/// zero by construction - the stage's flash is at the mixed enthalpy - and the capture
/// measures `1.4e-11` W there, which is the residual of the flash rather than a duty.
///
/// **`heat_input` reaches the equilibrium branch only.** The ratio branch solves for the
/// temperature from the ratio, so a duty and a boilup ratio are two specifications of the
/// same end and the class reads the ratio.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] for a negative or non-finite ratio, and whatever
/// the stage's or the flash's own refusals are.
pub fn reboiler(
    inlets: &[Stream],
    tray_pressure: Option<Pressure>,
    out_temperature: Option<ThermodynamicTemperature>,
    heat_input: Power,
    mode: ReboilerMode,
) -> Result<ReboilerOutcome> {
    let inlets_enthalpy: f64 = inlets.iter().map(|s| s.n * s.h.value).sum();

    match mode {
        ReboilerMode::Equilibrium => {
            let out = tray::tray(
                inlets,
                tray_pressure,
                out_temperature,
                heat_input,
                SideDraws::NONE,
                false,
            )?;
            let duty = outlet_enthalpy(&out) - inlets_enthalpy;
            Ok(ReboilerOutcome {
                temperature: out.temperature,
                pressure: out.pressure,
                vapour: out.gas,
                liquid: out.liquid,
                duty: watts(duty),
            })
        }
        ReboilerMode::VaporBoilupRatio(ratio) => {
            if !ratio.is_finite() || ratio < 0.0 {
                return Err(AzothError::invalid_input(
                    "boilup_ratio",
                    format!("a boilup ratio is a ratio of two flows, and {ratio} is not one"),
                ));
            }
            let mixed = mixer(inlets, tray_pressure)?;
            let (mixture, _ideal_gas) = mixed.mixture()?;
            let flash = pv_reflux_flash(
                &mixture,
                mixed.p,
                ratio,
                RefluxPhase::Liquid,
                mixed.t,
                &mixed.z,
            )?;
            let (vapour_fraction, liquid_fraction) = phase_fractions(flash.phase, flash.beta)?;

            let vapour = if vapour_fraction > 0.0 {
                Some(Stream::from_pt(
                    mixed.components.clone(),
                    flash.y,
                    mixed.n * vapour_fraction,
                    mixed.p,
                    flash.t,
                )?)
            } else {
                None
            };
            let liquid = if liquid_fraction > 0.0 {
                Some(Stream::from_pt(
                    mixed.components.clone(),
                    flash.x,
                    mixed.n * liquid_fraction,
                    mixed.p,
                    flash.t,
                )?)
            } else {
                None
            };
            let duty =
                enthalpy_of(vapour.as_ref()) + enthalpy_of(liquid.as_ref()) - inlets_enthalpy;

            Ok(ReboilerOutcome {
                temperature: flash.t,
                pressure: mixed.p,
                vapour,
                liquid,
                duty: watts(duty),
            })
        }
    }
}

/// The total enthalpy of a stage's two outlets, W.
fn outlet_enthalpy(out: &TrayOutcome) -> f64 {
    enthalpy_of(out.gas.as_ref()) + enthalpy_of(out.liquid.as_ref())
}

/// The total enthalpy of one outlet, W: zero where the phase is absent, which is NeqSim's
/// own `getMaterialStreamEnthalpy` on a zero-flow stream.
fn enthalpy_of(stream: Option<&Stream>) -> f64 {
    stream.map_or(0.0, |s| s.n * s.h.value)
}
