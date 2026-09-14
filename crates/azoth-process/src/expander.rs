//! `process.expander` - a pressure drop that produces work.
//!
//! Spec: `specs/models/process/expander.yaml`, which carries the provenance and what is
//! not ported. The procedure is [`crate::isentropic`]'s.
//!
//! What this model adds is the one expression that makes it an expander rather than a
//! compressor: [`Direction::Producing`] **multiplies** the ideal enthalpy change where the
//! compressor divides it. An expansion's enthalpy change is negative, so multiplying makes
//! it less negative - the real outlet is warmer than the ideal one, which is the
//! irreversibility the efficiency stands for.

use azoth_core::units::{Pressure, ThermodynamicTemperature, pascals};
use azoth_core::{Result, apply_checks};
use azoth_eos::mixture::Mixture;
use azoth_eos::molar_enthalpy_entropy::IdealGasModel;

use crate::isentropic::{self, Direction};
use crate::model_gen;
use crate::results::ExpanderResult;

/// A stream expanded to a stated outlet pressure, producing work.
///
/// `outlet_pressure` must be **below** the inlet pressure - the mirror of the check the
/// compressor makes, and the reason there are two models rather than one with a sign.
///
/// `power` on the result is **negative**, because the fluid is doing the work: the
/// convention across all three machines is that the shaft power is positive into the
/// fluid.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if `outlet_pressure` is not below the inlet
///   pressure.
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's declared range.
/// * Whatever the flashes raise.
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn expander(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t_in: ThermodynamicTemperature,
    p_in: Pressure,
    n: f64,
    z: &[f64],
    outlet_pressure: Pressure,
    efficiency: f64,
) -> Result<ExpanderResult> {
    let spec = &model_gen::EXPANDER_SPEC;
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
        Direction::Producing,
    )?;
    warnings.extend(solved.warnings);

    Ok(ExpanderResult {
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
