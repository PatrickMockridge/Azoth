//! `process.separator` - one feed split into a gas and a liquid at one state.
//!
//! Spec: `specs/models/process/separator.yaml`
//!
//! # The port
//!
//! NeqSim's `Separator.run` is lines 674-782 of a 4,511-line
//! `process/equipment/separator/Separator.java`. Read out of it, the thermodynamics is:
//!
//! ```text
//! P_out   = P_in - pressure_drop                     (Separator.java:694)
//! state   = PHflash(P_out, H_in + heat_input)        (:697-706)   if a duty is set
//!         | TPflash(T_in, P_out)                     (:707-711)   otherwise
//! n_gas   = beta       * n_in                        (:740-749)
//! n_liq   = (1 - beta) * n_in                        (:750-759)
//! outlet T = outlet P = the flash's                  (Stream.setThermoSystemFromPhase)
//! ```
//!
//! That is the whole of it. Everything else in those 109 lines is a low-flow bypass, a
//! memoization test on `|dH/H| < 1e-6`, an entrainment carry-over model
//! (`addPhaseFractionToPhase`) and a re-flash of the phase outlets; and everything else
//! in the file is geometry, level control and mechanical design. None of it is ported,
//! and `docs/src/roadmap.md` records that boundary.
//!
//! # What is *not* a transcription, and why
//!
//! **NeqSim runs its own internal `Mixer` over its inlet streams** (`:675`) before any
//! of the above. A separator with one feed does not need it, and a separator with
//! several is a mixer upstream of a separator - which is what `process.mixer` is. This
//! takes one feed, mixed by the caller.
//!
//! **The duty branch measures its enthalpy from the feed's own state.** NeqSim reads
//! `getEnthalpy()` off a system whose pressure has been changed but which has not been
//! re-flashed, so what it adds the duty to is whatever the previous flash left on the
//! object. This evaluates [`enthalpy_at`] at `(T_in, P_in)` explicitly, so the answer
//! depends on the inputs and not on what ran before - the same decision
//! `eos.stability_test` records about NeqSim's `Component.getK()`.
//!
//! # The split is decided on `phase`, never on `beta`
//!
//! `beta` is absent for a single-phase feed, and present-but-extrapolated for a
//! two-phase *verdict* reached by a negative-flash root - the value can lie outside
//! `[0, 1]`. Multiplying the feed by 1.888 would produce a wrong answer shaped exactly
//! like a right one. So the branch is on [`Phase`], and `beta` is only read where the
//! flash has said there is a genuine split.
//!
//! [`Phase::Trivial`] is the case that has no good answer. The flash converged to
//! `x = y = z`, which says the feed is single phase and **not which one**. The whole
//! feed is routed to the gas outlet, which is a convention rather than a finding, and
//! the result carries `phase` so a caller can see that it happened. A separator that
//! needs to know which phase it is holding has to ask `eos.stability_test`, or a
//! phase-boundary model, before it asks this one.

use azoth_core::units::{Power, Pressure, ThermodynamicTemperature, joules_per_mole, pascals};
use azoth_core::{Result, Warning, apply_checks};
use azoth_eos::Phase;
use azoth_eos::mixture::Mixture;
use azoth_eos::molar_enthalpy_entropy::IdealGasModel;
use azoth_eos::ph_flash::{enthalpy_at, ph_flash};
use azoth_eos::pt_flash::pt_flash;

use crate::beta_of;
use crate::model_gen;
use crate::results::SeparatorResult;

/// A flash's answer, normalised across the two models the two branches call.
///
/// `eos.pt_flash` reports no temperature - it was given one - and `eos.ph_flash` solves
/// for it. Everything below the split reads the same six fields from either, so the two
/// branches differ in one line each and the split itself is written once.
struct Split {
    temperature: ThermodynamicTemperature,
    beta: Option<f64>,
    x: Vec<f64>,
    y: Vec<f64>,
    phase: Phase,
    iterations: u32,
    warnings: Vec<Warning>,
}

/// A feed split into a gas and a liquid at one temperature and pressure.
///
/// `n_in` is a molar flow in mol/s, and it is the one argument here carried as a bare
/// `f64` rather than a `uom` quantity - see the crate note in `lib.rs` for why.
///
/// `heat_duty` is the one specification this model applies: zero means the separator is
/// isothermal at its inlet temperature, which is NeqSim's default and the only branch
/// that needs no enthalpy. A non-zero duty makes it a `PHflash` instead, and the
/// temperature becomes an answer rather than an input.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's declared range -
///   a zero flow, a negative pressure drop, a non-positive absolute state.
/// * [`azoth_core::AzothError::InvalidInput`] if the flash reports a two-phase split with no
///   vapour fraction, which is a contradiction rather than a state.
/// * Whatever the flash raises: an outlet pressure at or below zero, or a duty that no
///   temperature on `eos.ph_flash`'s bracket reaches.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, pascals, watts};
/// use azoth_eos::mixture::{Component, Mixture};
/// use azoth_eos::molar_enthalpy_entropy::IdealGasModel;
/// use azoth_process::separator;
///
/// let mixture = Mixture::new(
///     vec![
///         Component::new(kelvins(190.56), pascals(4_599_000.0), 0.0115)?,
///         Component::new(kelvins(425.12), pascals(3_796_000.0), 0.2002)?,
///     ],
///     vec![0.0, 0.01289789, 0.01289789, 0.0],
/// )?;
/// let ideal_gas = IdealGasModel {
///     cp_a: vec![3.0, 5.0],
///     cp_b: vec![0.0, 0.0],
///     cp_c: vec![0.0, 0.0],
///     cp_d: vec![0.0, 0.0],
///     h_ref: vec![0.0, 0.0],
///     s_ref: vec![0.0, 0.0],
///     t_ref: kelvins(300.0),
///     p_ref: pascals(100_000.0),
/// };
/// // A two-phase feed, so the split is a real one and the flows add back up.
/// let r = separator(
///     &mixture,
///     &ideal_gas,
///     kelvins(300.0),
///     pascals(2_000_000.0),
///     10.0,
///     &[0.6, 0.4],
///     pascals(0.0),
///     watts(0.0),
/// )?;
/// assert!((r.gas_flow + r.liquid_flow - 10.0).abs() < 1.0e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn separator(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t_in: ThermodynamicTemperature,
    p_in: Pressure,
    n_in: f64,
    z: &[f64],
    pressure_drop: Pressure,
    heat_duty: Power,
) -> Result<SeparatorResult> {
    let spec = &model_gen::SEPARATOR_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t_in.value),
            "P" => Some(p_in.value),
            "n" => Some(n_in),
            "pressure_drop" => Some(pressure_drop.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let p_out = pascals(p_in.value - pressure_drop.value);

    let split = if heat_duty.value == 0.0 {
        let flash = pt_flash(mixture, t_in, p_out, z)?;
        Split {
            temperature: t_in,
            beta: flash.beta,
            x: flash.x,
            y: flash.y,
            phase: flash.phase,
            iterations: flash.iterations,
            warnings: flash.warnings,
        }
    } else {
        let (h_in, _) = enthalpy_at(mixture, ideal_gas, t_in, p_in, z)?;
        let h_target = joules_per_mole(h_in + heat_duty.value / n_in);
        let flash = ph_flash(mixture, ideal_gas, p_out, h_target, z)?;
        Split {
            temperature: flash.temperature,
            beta: flash.beta,
            x: flash.x,
            y: flash.y,
            phase: flash.phase,
            iterations: flash.iterations,
            warnings: flash.warnings,
        }
    };
    warnings.extend(split.warnings);

    let (gas_flow, liquid_flow, gas_z, liquid_z) = match split.phase {
        Phase::TwoPhase => {
            let beta = beta_of(spec, split.beta)?;
            (
                beta * n_in,
                (1.0 - beta) * n_in,
                split.y.clone(),
                split.x.clone(),
            )
        }
        // An empty outlet carries the feed's composition and a flow of zero. A stream
        // with no flow has no phase composition of its own, and a row of zeros is not a
        // composition - it would not sum to one, and the first downstream unit to read
        // it would produce a plausible wrong answer.
        Phase::AllLiquid => (0.0, n_in, z.to_vec(), split.x.clone()),
        Phase::AllVapour | Phase::Trivial => (n_in, 0.0, split.y.clone(), z.to_vec()),
    };

    Ok(SeparatorResult {
        temperature: split.temperature,
        pressure: p_out,
        beta: split.beta,
        gas_flow,
        liquid_flow,
        gas_z,
        liquid_z,
        phase: split.phase,
        iterations: split.iterations,
        warnings,
    })
}
