//! `process.separator` - one feed split into a gas and a liquid at one state.
//!
//! Spec: `specs/models/process/separator.yaml`

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
/// `n_in` is a molar flow in mol/s. A `heat_duty` of zero flashes isothermally at the
/// inlet temperature; a non-zero duty flashes at constant enthalpy, so the temperature
/// becomes an answer.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's declared range.
/// * [`azoth_core::AzothError::InvalidInput`] if the flash reports a two-phase split with no
///   vapour fraction.
/// * Whatever the flash raises.
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
/// // Two-phase at this state, so both flows are non-zero and they sum to the feed.
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
        // A zero-flow outlet carries the feed's composition: a row of zeros is not a
        // composition, and a unit reading one would produce a plausible wrong answer.
        Phase::AllLiquid => (0.0, n_in, z.to_vec(), split.x.clone()),
        // `Trivial` means the flash established that the feed is single phase but not
        // which one, so the routing is a convention and `phase` is what says so.
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
