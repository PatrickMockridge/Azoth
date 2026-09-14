//! `process.pump` - a pressure rise in a liquid, at a stated isentropic efficiency.
//!
//! Spec: `specs/models/process/pump.yaml`, which carries the provenance, which of the four
//! paths through NeqSim's `Pump.run` this is, and why a liquid pump runs an isentropic
//! flash at all. The procedure is [`crate::isentropic`]'s.

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
/// [`crate::compressor::compressor`] gives: a silently ideal pump is a pump that
/// understates its power.
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
