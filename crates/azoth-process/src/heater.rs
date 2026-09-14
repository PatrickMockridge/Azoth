//! `process.heater` - a duty applied at a fixed pressure.
//!
//! Spec: `specs/models/process/heater.yaml`
//!
//! # The port, and the unit it replaces
//!
//! NeqSim has `Heater` and `Cooler` as two classes, and **`Cooler` has no `run()` of its
//! own** - `Cooler.java` is 258 lines and inherits `Heater.run` for its steady state.
//! The two differ only in the sign of the duty and in what their dynamics do. So a
//! cooler here is this model with a negative `heat_duty`, which is what NeqSim's own
//! code says and is one fewer unit operation to wire, spec and test.
//!
//! `Heater.run` is lines 402-471 of a 1,113-line file, and its physics:
//!
//! ```text
//! P_out   = P_in - pressureDrop                  (Heater.java:432-435)
//! H_out   = H_in + energyInput                   (:431)
//! T_out   = PHflash(P_out, H_out)                (:445-446)
//! Q       = H_out - H_in                         (:457)
//! ```
//!
//! Everything else in the file is the specification switch (`out stream`,
//! `setTemperature`, `setEnergyInput`, `deltaT` - `:437-450`), the energy port and the
//! caches.
//!
//! # Why only the duty specification is ported
//!
//! NeqSim's other three specifications are not here, and the reason is that **none of
//! them is a unit operation**. A heater with a specified outlet temperature is a
//! `TPflash`, which is `eos.pt_flash` - the temperature is already an input there and
//! there is no energy balance to solve. A specified `deltaT` is the same call with the
//! temperature added first. This model exists to answer the question those cannot: given
//! the duty, what temperature results.
//!
//! So the Pareto set is one specification, and the other two are compositions of models
//! this library already has.

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
