//! `process.pump` - a pressure rise in a liquid, at a stated isentropic efficiency.
//!
//! Spec: `specs/models/process/pump.yaml`
//!
//! Physics, and everything not ported, is in [`crate::isentropic`].
//!
//! # Which of NeqSim's three paths this is
//!
//! `Pump.run` is lines 515-677 of a 1,700-line file and it has four ways through it.
//! This is the **default** one - `:558-572`, taken when `calculateAsCompressor` is true,
//! which it is by default (`Pump.java:114`) - and that path is literally the compressor's:
//! a `PSflash` at the inlet entropy, then a `PHflash` at the actual enthalpy.
//!
//! The other three are not ported. The pump-curve path (`:573-644`) needs a
//! manufacturer's head-flow curve, which is data this library does not ship and cannot
//! check. The simple pressure-rise path (`:645-664`) computes a hydraulic power from
//! `dP * volumetric flow` and divides by the efficiency, which is the incompressible
//! approximation of the same thing - it is a *different model* rather than a different
//! answer, and a caller who wants it can compose `hydraulics.pump_power` with
//! `process.heater`. The fixed-outlet-temperature path (`:549-556`) is a `TPflash` and so
//! is `eos.pt_flash`.
//!
//! # Why a liquid pump needs an isentropic flash at all
//!
//! Because the temperature rises, slightly, and the model should say by how much rather
//! than assume zero. For water at ordinary conditions it is a few hundredths of a kelvin
//! per hundred bar - small, real, and free here, since the flash this runs is the same
//! one a compressor runs.

use azoth_core::units::{Pressure, ThermodynamicTemperature, pascals};
use azoth_core::{Result, apply_checks};
use azoth_eos::mixture::Mixture;
use azoth_eos::molar_enthalpy_entropy::IdealGasModel;

use crate::isentropic::{self, Direction};
use crate::model_gen;
use crate::results::PumpResult;

/// A stream pumped to a stated outlet pressure at a stated isentropic efficiency.
///
/// `efficiency` is required and range-checked for the reason
/// [`crate::compressor::compressor`] gives: NeqSim defaults it to one and clamps, and a
/// silently ideal pump is a pump that understates its power.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if `outlet_pressure` is not above the inlet
///   pressure.
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's declared range.
/// * Whatever the flashes raise.
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn pump(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t_in: ThermodynamicTemperature,
    p_in: Pressure,
    n: f64,
    z: &[f64],
    outlet_pressure: Pressure,
    efficiency: f64,
) -> Result<PumpResult> {
    let spec = &model_gen::PUMP_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t_in.value),
            "P" => Some(p_in.value),
            "n" => Some(n),
            "outlet_pressure" => Some(outlet_pressure.value),
            "efficiency" => Some(efficiency),
            _ => None,
        },
        &mut warnings,
    )?;

    let solved = isentropic::run(
        mixture,
        ideal_gas,
        t_in,
        p_in,
        pascals(outlet_pressure.value),
        n,
        z,
        efficiency,
        Direction::Consuming,
    )?;
    warnings.extend(solved.warnings);

    Ok(PumpResult {
        temperature: solved.temperature,
        pressure: pascals(outlet_pressure.value),
        power: solved.power,
        beta: solved.beta,
        phase: solved.phase,
        isentropic_temperature: solved.isentropic_temperature,
        iterations: solved.iterations,
        warnings,
    })
}
