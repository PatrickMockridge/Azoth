//! The stage: `SimpleTray`, which is a mixer, a duty and a flash.

use azoth_core::units::{Power, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{AzothError, Result};
use azoth_eos::{Phase, ph_flash, pt_flash};

use crate::kernels::mixer;
use crate::stream::Stream;

/// What one stage hands up and down.
///
/// **An absent outlet is `None`, not a zero-flow stream.** `SimpleTray` reports a zero-flow
/// stream carrying the *tray's* composition, and the capture shows what that costs: a
/// subcooled feed's vapour outlet has `n = 0` and `h = -Infinity`, because NeqSim divides an
/// empty system's enthalpy by its zero moles. That is not a state this library can express -
/// a `Stream`'s molar enthalpy is finite by construction - and inventing one would be a
/// fabricated state that looks like a real one. The absence of a phase is a fact, so it is
/// spelled as one.
#[derive(Debug, Clone)]
pub struct TrayOutcome {
    /// The tray's temperature, which is the flash's.
    pub temperature: ThermodynamicTemperature,
    /// The tray's pressure.
    pub pressure: Pressure,
    /// The vapour leaving the stage, or `None` where the flash found no vapour.
    pub gas: Option<Stream>,
    /// The liquid leaving the stage, or `None` where the flash found no liquid.
    pub liquid: Option<Stream>,
}

/// One equilibrium stage: mix the inlets, add the duty, flash.
///
/// The arithmetic is `SimpleTray.run`'s. Its inlets are mixed by [`crate::kernels::mixer`],
/// which is the same arithmetic the mixer's own kernel runs - **and the tray's
/// `calcMixStreamEnthalpy` overrides the mixer's**, starting from the tray's `heatInput` so
/// that the duty enters the flash as a raised enthalpy rather than as a separate term.
///
/// **The duty is divided by the flow, and that was measured rather than assumed.** NeqSim's
/// `PHflash(enthalpy, 0)` takes a *total* enthalpy: 5000 W raises one mol/s by 15.32 K and
/// two mol/s by 8.85 K, in the capture's last two rows. A molar reading would have raised
/// both by the same amount, and no unit-flow case can tell the two apart - which is why the
/// capture carries a two-mol row.
///
/// With `out_temperature` given the flash is a `TPflash` at that temperature, which is
/// `SimpleTray`'s other branch; without it the flash is at the mixed enthalpy.
///
/// **The temperature is not what the capture is oracled on, and that is a NeqSim finding.**
/// Its `PHflash(h, 0)` returns a state whose own enthalpy differs from the one it was asked
/// for, measured on a fresh fluid with no tray in the picture: the capture's
/// `neqsim_phflash_honours_its_enthalpy` block asks for `+2500`, `+5000` and `-2000` J/mol and
/// reads back `+179.5`, `-27.3` and `+0.00002` J/mol of error. That is the solver's
/// stopping rule and not its equations - a 179 J/mol error on a 2500 J/mol step is not a
/// property of any equation of state. So a case here asserts the *outlets* against the
/// capture, where the two libraries agree to 0.06 J/mol, and asserts the temperature against
/// the stage's own enthalpy balance, where azoth's flash is consistent to `5.6e-7` J/mol.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if there are no inlets, if they carry nothing,
/// or if the flash converges trivially - where NeqSim reads a phase type this library's
/// model does not prove; and whatever the flash refuses.
pub fn tray(
    inlets: &[Stream],
    tray_pressure: Option<Pressure>,
    out_temperature: Option<ThermodynamicTemperature>,
    heat_input: Power,
) -> Result<TrayOutcome> {
    let mixed = mixer(inlets, tray_pressure)?;
    let (mixture, ideal_gas) = mixed.mixture()?;

    let (temperature, phase, x, y, gas_fraction, liquid_fraction) = match out_temperature {
        Some(t) => {
            let flash = pt_flash(&mixture, t, mixed.p, &mixed.z)?;
            let (gas, liquid) = phase_fractions(flash.phase, flash.beta)?;
            (t, flash.phase, flash.x, flash.y, gas, liquid)
        }
        None => {
            // `calcMixStreamEnthalpy`: the inlets' total enthalpy plus the duty, on the
            // same flow its own flash divides by.
            let target = mixed.h.value + heat_input.value / mixed.n;
            let flash = ph_flash(
                &mixture,
                &ideal_gas,
                mixed.p,
                joules_per_mole(target),
                &mixed.z,
            )?;
            let (gas, liquid) = phase_fractions(flash.phase, flash.beta)?;
            (
                flash.temperature,
                flash.phase,
                flash.x,
                flash.y,
                gas,
                liquid,
            )
        }
    };

    // **A single phase has the tray's composition, and the flash's `x` or `y` is not it.**
    // Measured: on the capture's subcooled row the flash's `x` comes back
    // `[0.2310, 0.2649, 0.3369, 0.1672]` where the tray's own composition is
    // `[0.1, 0.3, 0.4, 0.2]` - a trial phase rather than the phase that is there - and NeqSim
    // reports the feed's, because the one phase *is* the feed.
    let single = matches!(phase, Phase::AllVapour | Phase::AllLiquid);
    let gas_z = if single { mixed.z.clone() } else { y };
    let liquid_z = if single { mixed.z.clone() } else { x };

    let gas = if gas_fraction > 0.0 {
        Some(Stream::from_pt(
            mixed.components.clone(),
            gas_z,
            mixed.n * gas_fraction,
            mixed.p,
            temperature,
        )?)
    } else {
        None
    };
    let liquid = if liquid_fraction > 0.0 {
        Some(Stream::from_pt(
            mixed.components.clone(),
            liquid_z,
            mixed.n * liquid_fraction,
            mixed.p,
            temperature,
        )?)
    } else {
        None
    };

    Ok(TrayOutcome {
        temperature,
        pressure: mixed.p,
        gas,
        liquid,
    })
}

/// What share of the tray's flow leaves as vapour and as liquid.
///
/// **The phase decides, and `beta` does not.** NeqSim's two outlet getters look for a phase
/// *type*: a subcooled feed has no gas phase and reports a vapour of zero flow while its
/// `getBeta()` still answers `1.0`. Reading the split off `beta` would give that feed all of
/// its flow as vapour, which is the opposite of the measurement.
///
/// A trivial solution is refused for the reason `unit_ops.shortcut_distillation_column`
/// gives: every K-value straddled one, so nothing here proves which single phase the tray
/// holds, and NeqSim reads a phase type this model does not derive.
fn phase_fractions(phase: Phase, beta: Option<f64>) -> Result<(f64, f64)> {
    match phase {
        Phase::TwoPhase => {
            let vapour = beta.ok_or_else(|| {
                AzothError::invalid_input(
                    "flash",
                    "a two-phase answer with no vapour fraction, so the split is unknown",
                )
            })?;
            Ok((vapour, 1.0 - vapour))
        }
        Phase::AllVapour => Ok((1.0, 0.0)),
        Phase::AllLiquid => Ok((0.0, 1.0)),
        Phase::Trivial => Err(AzothError::invalid_input(
            "flash",
            "the tray's flash is a trivial solution - every K-value straddled one, so nothing \
             here proves which single phase the tray holds",
        )),
    }
}
