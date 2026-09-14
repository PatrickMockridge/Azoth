//! `process.heater` - a duty applied at a fixed pressure.
//!
//! Spec: `specs/models/process/heater.yaml`, which carries the provenance and the reason
//! one of the four specifications is ported and the others are not.
//!
//! A cooler is this model with a negative `heat_duty`; the sign is the whole difference.

use azoth_core::units::{Power, Pressure, ThermodynamicTemperature, joules_per_mole, pascals};
use azoth_core::{Result, apply_checks};
use azoth_eos::mixture::Mixture;
use azoth_eos::molar_enthalpy_entropy::IdealGasModel;
use azoth_eos::ph_flash::{enthalpy_at, ph_flash};

use crate::model_gen;
use crate::results::HeaterResult;

/// A duty applied to a stream at a fixed pressure.
///
/// `heat_duty` is in watts and is **signed**: positive adds energy and raises the
/// temperature, negative removes it and is a cooler. `n` is the molar flow in mol/s, and
/// it is needed rather than decorative - the duty is an *extensive* quantity and the
/// flash inverts a *molar* enthalpy, so the conversion `Q / n` is the model.
///
/// That division is safe because the spec's range check puts `n` strictly above zero,
/// and it is the reason the check is there rather than a nicety. A duty too large for any
/// temperature on `eos.ph_flash`'s bracket is a solver failure and is reported as one,
/// not clamped.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's declared range -
///   a zero flow, a non-positive absolute state - or if the flash cannot find the
///   requested enthalpy on its bracket.
/// * [`azoth_core::AzothError::InvalidInput`] if the composition is not a composition.
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn heater(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t_in: ThermodynamicTemperature,
    p_in: Pressure,
    n: f64,
    z: &[f64],
    pressure_drop: Pressure,
    heat_duty: Power,
) -> Result<HeaterResult> {
    let spec = &model_gen::HEATER_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t_in.value),
            "P" => Some(p_in.value),
            "n" => Some(n),
            "pressure_drop" => Some(pressure_drop.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let p_out = pascals(p_in.value - pressure_drop.value);

    let (h_in, _) = enthalpy_at(mixture, ideal_gas, t_in, p_in, z)?;
    let h_out = joules_per_mole(h_in + heat_duty.value / n);
    let flash = ph_flash(mixture, ideal_gas, p_out, h_out, z)?;
    warnings.extend(flash.warnings);

    Ok(HeaterResult {
        temperature: flash.temperature,
        pressure: p_out,
        phase: flash.phase,
        beta: flash.beta,
        iterations: flash.iterations,
        warnings,
    })
}
