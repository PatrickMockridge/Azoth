//! `process.compressor` - a pressure rise at a stated isentropic efficiency.
//!
//! Spec: `specs/models/process/compressor.yaml`, which carries the provenance and what is
//! not ported. The procedure is [`crate::isentropic`]'s - the three machines share it and
//! differ only in which way the efficiency scales the ideal enthalpy change.

use azoth_core::units::{Pressure, ThermodynamicTemperature, pascals};
use azoth_core::{Result, apply_checks};
use azoth_eos::mixture::Mixture;
use azoth_eos::molar_enthalpy_entropy::IdealGasModel;

use crate::isentropic::{self, Direction};
use crate::model_gen;
use crate::results::CompressorResult;

/// A stream compressed to a stated outlet pressure at a stated isentropic efficiency.
///
/// `efficiency` is the **isentropic** efficiency, in `(0, 1]`, and is a required input
/// with a range check: an efficiency of one is a machine with no losses, and a model that
/// assumes it silently understates every duty it is asked for.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if `outlet_pressure` is not above the inlet
///   pressure - a compressor that lowers the pressure is an expander, and there is one.
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's declared range,
///   including an efficiency above one or below zero.
/// * Whatever the flashes raise.
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn compressor(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t_in: ThermodynamicTemperature,
    p_in: Pressure,
    n: f64,
    z: &[f64],
    outlet_pressure: Pressure,
    efficiency: f64,
) -> Result<CompressorResult> {
    let spec = &model_gen::COMPRESSOR_SPEC;
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

    Ok(CompressorResult {
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
