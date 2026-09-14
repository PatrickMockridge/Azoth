//! `process.expander` - a pressure drop that produces work.
//!
//! Spec: `specs/models/process/expander.yaml`
//!
//! Physics, and everything not ported, is in [`crate::isentropic`].
//!
//! # The one line that makes it an expander
//!
//! NeqSim's `Expander` **extends `Compressor`** and overrides `run` - `Expander.java` is
//! 669 lines, of which `run` is 61 (`:608-668`). Its isentropic path differs from the
//! compressor's in exactly one expression:
//!
//! ```text
//! Compressor.java:1738   dH = (H(P_out, s_in) - H_in) / isentropicEfficiency
//! Expander.java:653      dH = (H(P_out, s_in) - H_in) * isentropicEfficiency
//! ```
//!
//! Dividing gives an outlet further from the inlet than the ideal one, which is a machine
//! that consumes work; multiplying gives one closer to it, which is a machine that
//! produces it. The enthalpy change is negative for an expansion, so the multiplication
//! makes it *less* negative - the real outlet is warmer than the ideal one, which is the
//! irreversibility the efficiency stands for.
//!
//! That difference is the whole of this file. [`crate::isentropic::Direction`] names it.
//!
//! # The polytropic path is not ported
//!
//! `Expander.java:620-633` steps the pressure in five stages and applies the polytropic
//! efficiency to each. It is a different model of the same machine rather than a more
//! accurate one, and like the compressor's polytropic paths it needs `usePolytropicCalc`
//! set. The default is the isentropic path (`usePolytropicCalc = false`, inherited from
//! `Compressor.java:111`) and that is what is here.

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
/// fluid. NeqSim's expander reports the opposite sign on its energy port
/// (`Expander.java:661`, `setDuty(-dH)`), so a duty read from here and a duty read from
/// there differ by a sign - which is exactly the kind of difference worth stating rather
/// than discovering.
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
