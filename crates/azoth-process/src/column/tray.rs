//! The stage: `SimpleTray`, which is a mixer, a duty and a flash.

use azoth_core::units::{Power, Pressure, ThermodynamicTemperature, joules_per_mole, kelvins};
use azoth_core::warning::Warning;
use azoth_core::{AzothError, Result};
use azoth_eos::{Cubic, Phase, ph_flash, pt_flash};
use azoth_reactions::databank::formation_properties;
use azoth_reactions::reactive_ph_flash::reactive_ph_flash;
use azoth_reactions::reactive_tp_flash::{ReactiveTpFlashResult, reactive_tp_flash};

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
    /// **What the stage had to fall back on.** NeqSim's reactive route falls back from a
    /// reactive PH flash to a reactive TP one and then to a plain TP flash, logging each step;
    /// a stage says which it took rather than doing it silently.
    pub warnings: Vec<Warning>,
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
    reactive: bool,
) -> Result<TrayOutcome> {
    let mixed = mixer(inlets, tray_pressure)?;
    if reactive {
        return reactive_stage(&mixed, out_temperature, heat_input);
    }
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
        warnings: Vec::new(),
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

/// **A reactive stage: `SimpleTray.run` with `useReactiveFlash` set.**
///
/// The class's own route, in its own order. A tray whose outlet temperature is stated takes a
/// **reactive TP flash** there; a tray that flashes at its own enthalpy takes a **reactive PH
/// flash**, and - because that answer is a temperature and not a phase list - a reactive TP
/// flash at it, which is what NeqSim's own system holds afterwards. A reactive PH flash that
/// fails falls back the way the class falls back: a reactive TP flash, then a plain one, and
/// the stage says which it took.
///
/// **The fluid is the process layer's cubic**, which is the whole reason
/// `reactions.reactive_tp_flash` takes one: this column is PR, so a reactive tray that flashed
/// it with SRK would be a different machine from the trays beside it.
fn reactive_stage(
    mixed: &Stream,
    out_temperature: Option<ThermodynamicTemperature>,
    heat_input: Power,
) -> Result<TrayOutcome> {
    let names = mixed.components.clone();
    let moles: Vec<f64> = mixed.z.iter().map(|z| z * mixed.n).collect();
    let pressure = mixed.p.value;
    let mut warnings = Vec::new();

    let (temperature, tp) = match out_temperature {
        Some(pin) => {
            let t = pin.value;
            (t, reactive_tp(&names, t, pressure, &moles)?)
        }
        None => {
            // `calcMixStreamEnthalpy`: the inlets' total enthalpy plus the duty, and then the
            // *formation inventory* the reactive specification is stated against - the class's
            // constructor forms that sum and a sensible enthalpy without it is a different
            // number, which finds a different temperature.
            let target =
                mixed.h.value * mixed.n + heat_input.value + formation_inventory(&names, &moles)?;
            match reactive_ph_flash(
                &names,
                Cubic::Pr,
                mixed.t.value,
                pressure,
                &moles,
                target,
                2,
            ) {
                Ok(out) => {
                    let t = out.temperature;
                    // The PH answer is a temperature; the split is the TP flash at it, which is
                    // the state NeqSim's system is left holding.
                    (t, reactive_tp(&names, t, pressure, &moles)?)
                }
                Err(error) => {
                    let t = out_temperature.map_or(mixed.t.value, |pin| pin.value);
                    warnings.push(reactive_fallback(format!(
                        "the reactive PH flash did not answer ({error}), so the tray flashed \
                         reactively at {t} K instead"
                    )));
                    (t, reactive_tp(&names, t, pressure, &moles)?)
                }
            }
        }
    };

    let (gas, liquid) = reactive_outlets(&names, mixed.p, temperature, &tp, &mut warnings)?;
    Ok(TrayOutcome {
        temperature: kelvins(temperature),
        pressure: mixed.p,
        gas,
        liquid,
        warnings,
    })
}

/// `testOps.reactiveTPflash()`: the reactive split at a stated state.
fn reactive_tp(
    names: &[String],
    temperature: f64,
    pressure: f64,
    moles: &[f64],
) -> Result<ReactiveTpFlashResult> {
    reactive_tp_flash(names, Cubic::Pr, temperature, pressure, moles, 2)
}

/// `sum_i n_i dHf_i`: the formation inventory the reactive enthalpy specification carries.
fn formation_inventory(names: &[String], moles: &[f64]) -> Result<f64> {
    let mut total = 0.0;
    for (name, amount) in names.iter().zip(moles) {
        let properties = formation_properties(name)?.ok_or_else(|| {
            AzothError::invalid_input(
                "components",
                format!("the component databank carries no formation row for `{name}`"),
            )
        })?;
        total += amount * properties.enthalpy_of_formation;
    }
    Ok(total)
}

/// The two outlets a reactive answer holds, or a refusal where it holds no way to tell them
/// apart.
///
/// **The label decides where there is one**, and the flash emits it exactly where it built a
/// vapour/liquid pair. Where it did not - a reacting solve's two phases converge to the *same*
/// composition on every captured state - a tray has no basis for sending one row up and the
/// other down, so it refuses rather than picking a row.
///
/// **A single phase needs no pair and is still typed**: the port's own flash at the answer's
/// state says whether that one phase is a vapour or a liquid, which is what NeqSim's
/// `getGasOutStream` reads off its system.
fn reactive_outlets(
    names: &[String],
    pressure: Pressure,
    temperature: f64,
    tp: &ReactiveTpFlashResult,
    warnings: &mut Vec<Warning>,
) -> Result<(Option<Stream>, Option<Stream>)> {
    // **A phase the class calls negligible is no outlet.** `removeNegligiblePhases` drops one
    // below `MIN_PHASE_FRACTION`, and the delegation's answer on a fluid that turns out to be
    // single-phase carries exactly such a row: measured on methane/ethane at 300 K and 15 bar,
    // the pairs come back `1e-15` and `1.0` mol with the *feed's* composition in both, because
    // the class's constructor leaves two phase objects each holding the whole feed. Reporting a
    // stage that hands 1e-15 mol/s down a column would be reporting a stream the port refuses to
    // fabricate anyway.
    let total: f64 = tp
        .phase_moles
        .iter()
        .map(|row| row.iter().sum::<f64>())
        .sum();
    let held: Vec<(usize, f64)> = tp
        .phase_moles
        .iter()
        .enumerate()
        .map(|(index, row)| (index, row.iter().sum::<f64>()))
        .filter(|(_, moles)| {
            total > 0.0 && moles / total >= azoth_reactions::reactive_flash::MIN_PHASE_FRACTION
        })
        .collect();

    let stream_of = |index: usize, amount: f64| -> Result<Stream> {
        let row = &tp.phase_moles[index];
        let z: Vec<f64> = row.iter().map(|n| n / amount).collect();
        Stream::from_pt(names.to_vec(), z, amount, pressure, kelvins(temperature))
    };

    match held.as_slice() {
        // **One phase, and the port's own flash says which it is.** The delegation's indices are
        // its bookkeeping rather than a type - a single-phase fluid comes back as a full second
        // row whatever it is - so the character is measured at the answer's own state instead,
        // which is what `getGasOutStream` reads off NeqSim's system.
        [(index, amount)] => {
            let row = &tp.phase_moles[*index];
            let z: Vec<f64> = row.iter().map(|n| n / amount).collect();
            let names_ref: Vec<&str> = names.iter().map(String::as_str).collect();
            let (mixture, _) = azoth_eos::databank::mixture_of(&names_ref, Cubic::Pr, None)?;
            let flash = pt_flash(&mixture, kelvins(temperature), pressure, &z)?;
            let stream = stream_of(*index, *amount)?;
            match flash.phase {
                Phase::AllVapour => Ok((Some(stream), None)),
                Phase::AllLiquid => Ok((None, Some(stream))),
                _ => Err(AzothError::invalid_input(
                    "flash",
                    "a reactive answer of one phase that this flash does not hold as one, so \
                     nothing here says which outlet it leaves by",
                )),
            }
        }
        // **Two phases, and the label decides** - which is only asked for where both are real.
        [(first, a), (second, b)] => {
            let (gas_index, liquid_index, gas_moles, liquid_moles) =
                match (tp.phase_type.get(*first), tp.phase_type.get(*second)) {
                    (Some(f), Some(s)) if f == "vapour" && s == "liquid" => {
                        (*first, *second, *a, *b)
                    }
                    (Some(f), Some(s)) if f == "liquid" && s == "vapour" => {
                        (*second, *first, *b, *a)
                    }
                    _ => {
                        warnings.push(reactive_fallback(
                            "the reactive answer's two phases carry no vapour/liquid label, so \
                             this stage cannot form two outlets from it"
                                .to_string(),
                        ));
                        return Err(AzothError::invalid_input(
                            "flash",
                            "the reactive flash returned two phases without a vapour/liquid \
                             label: where its two phases converge to the same composition, naming \
                             one of them the vapour would be picking a row rather than measuring \
                             one",
                        ));
                    }
                };
            Ok((
                Some(stream_of(gas_index, gas_moles)?),
                Some(stream_of(liquid_index, liquid_moles)?),
            ))
        }
        _ => Err(AzothError::invalid_input(
            "flash",
            "the reactive answer holds no phase above the class's own negligible fraction",
        )),
    }
}

/// A stage's caveat. **The code is `SolverNotConverged`**, which is what happened: the flash
/// the class asked for did not answer, and the stage took the class's own next step. There is no
/// separate code for a fallback, and the vocabulary is a cross-language contract.
fn reactive_fallback(note: String) -> Warning {
    Warning {
        code: azoth_core::warning::WarningCode::SolverNotConverged,
        message: note,
        field: None,
    }
}
