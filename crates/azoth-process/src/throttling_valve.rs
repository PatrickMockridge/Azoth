//! `process.throttling_valve` - pressure dropped at constant enthalpy.
//!
//! Spec: `specs/models/process/throttling_valve.yaml`
//!
//! # The port
//!
//! NeqSim's `ThrottlingValve.run` is lines 363-453 of a 1,894-line file, and the
//! thermodynamics in it is two statements:
//!
//! ```text
//! P_out = P_in - deltaP                     (ThrottlingValve.java:407-412, clamped at :424-432)
//! T_out = PHflash(P_out, H_in)              (runPHflashWithNaNRetry, :467-495)
//! ```
//!
//! The rest of the file is valve sizing - the IEC 60534 capacity calculation from `Kv`
//! and opening, the `findOutletPressureForFixedKv` solve, the opening actuator - and the
//! transient response. None of it is here, and none of it needs to be: **`Kv` sizing is
//! already in this library** as `hydraulics.control_valve_cv`, which is the same
//! standard applied to a valve rather than to a stream.
//!
//! # Why there is no flow in this signature
//!
//! A valve does not change how much is flowing or what it is, and the outlet state does
//! not depend on the flow rate at all - an isenthalpic flash is a *molar* property, and
//! dividing by the flow changes nothing. So the flow is genuinely absent from this
//! model's inputs rather than merely unused, and the same is true of the composition
//! beyond the fact that it fixes the mixture.
//!
//! NeqSim computes a `molarFlow` in its `run` too (`:371`) and applies it only on the
//! transient path. A steady-state caller has no use for it.
//!
//! # The temperature really does change
//!
//! A throttling valve is the classic Joule-Thomson device, and a model that returned
//! the inlet temperature would be wrong for every real gas. The `PHflash` is what makes
//! it right, and it is the reason this unit operation could not be written before
//! `eos.ph_flash` existed.

use azoth_core::units::{Pressure, ThermodynamicTemperature, joules_per_mole, pascals};
use azoth_core::{Result, apply_checks};
use azoth_eos::mixture::Mixture;
use azoth_eos::molar_enthalpy_entropy::IdealGasModel;
use azoth_eos::ph_flash::{enthalpy_at, ph_flash};

use crate::model_gen;
use crate::results::ThrottlingValveResult;

/// A stream's pressure dropped at constant enthalpy.
///
/// `pressure_drop` is positive and is subtracted from `p_in`. There is no
/// `accept_negative_dp` flag as NeqSim has (`:424-432`): the spec's range check refuses
/// a negative drop instead, because a valve that raises the pressure is not a valve and
/// admitting one would let a caller model a compressor by accident.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `p_in` is not positive or the drop would take
///   the outlet to or below zero - the flash catches that last one, and the check it
///   already applies is the one that matters.
/// * Whatever the flashes raise.
pub fn throttling_valve(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t_in: ThermodynamicTemperature,
    p_in: Pressure,
    z: &[f64],
    pressure_drop: Pressure,
) -> Result<ThrottlingValveResult> {
    let spec = &model_gen::THROTTLING_VALVE_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t_in.value),
            "P" => Some(p_in.value),
            "pressure_drop" => Some(pressure_drop.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let p_out = pascals(p_in.value - pressure_drop.value);

    // The enthalpy the fluid arrived with. This is the whole of the model: what leaves
    // has the same enthalpy and less pressure, and everything else follows from that.
    let (h_in, _) = enthalpy_at(mixture, ideal_gas, t_in, p_in, z)?;
    let flash = ph_flash(mixture, ideal_gas, p_out, joules_per_mole(h_in), z)?;
    warnings.extend(flash.warnings);

    Ok(ThrottlingValveResult {
        temperature: flash.temperature,
        pressure: p_out,
        phase: flash.phase,
        beta: flash.beta,
        iterations: flash.iterations,
        warnings,
    })
}
