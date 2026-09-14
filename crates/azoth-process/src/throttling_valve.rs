//! `process.throttling_valve` - pressure dropped at constant enthalpy.
//!
//! Spec: `specs/models/process/throttling_valve.yaml`, which carries the provenance, the
//! absent flow rate and the divergence over a negative drop.
//!
//! There is no flow in this signature and none is needed: an isenthalpic flash is a
//! *molar* property, so the outlet state does not depend on the flow rate, and a valve
//! changes neither the flow nor the composition. The `PHflash` is the whole of the
//! thermodynamics - a valve is the classic Joule-Thomson device, and a model returning the
//! inlet temperature would be wrong for every real gas.

use azoth_core::units::{Pressure, ThermodynamicTemperature, joules_per_mole, pascals};
use azoth_core::{Result, apply_checks};
use azoth_eos::mixture::Mixture;
use azoth_eos::molar_enthalpy_entropy::IdealGasModel;
use azoth_eos::ph_flash::{enthalpy_at, ph_flash};

use crate::model_gen;
use crate::results::ThrottlingValveResult;

/// A stream's pressure dropped at constant enthalpy.
///
/// `pressure_drop` is positive and is subtracted from `p_in`. The spec's range check
/// refuses a negative drop, because a valve that raises the pressure is not a valve.
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
