//! `eos.ts_flash` - the pressure (or temperature) a mixture reaches at a given molar entropy.
//!
//! Spec: `specs/models/eos/ts_flash.toml`
//!
//! A thin procedure model over the flash-property solver, mirroring `eos.tv_flash` /
//! `eos.pv_flash` with the property exchanged.

use azoth_core::units::{MolarHeatCapacity, ThermodynamicTemperature, pascals};
use azoth_core::{Result, apply_checks};

use crate::algorithm_of;
use crate::flash_property::{Property, solve_pressure};
use crate::mixture::Mixture;
use crate::model_gen;
use crate::molar_enthalpy_entropy::IdealGasModel;
use crate::results::TsFlashResult;

/// The pressure at which a mixture has a given molar entropy.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if a state input is not positive.
/// * [`azoth_core::AzothError::SolverNotConverged`] if the iteration reaches its cap.
pub fn ts_flash(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t: ThermodynamicTemperature,
    s: MolarHeatCapacity,
    z: &[f64],
) -> Result<TsFlashResult> {
    let spec = &model_gen::TS_FLASH_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "S" => Some(s.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let algorithm = algorithm_of(spec)?;
    let start = 1.0e5;
    let solved = solve_pressure(
        mixture,
        ideal_gas,
        t,
        s.value,
        z,
        Property::Entropy,
        algorithm,
        start,
    )?;

    apply_checks(
        spec.derived_checks(),
        |quantity| (quantity == "P").then_some(solved.pressure),
        &mut warnings,
    )?;

    Ok(TsFlashResult {
        pressure: pascals(solved.pressure),
        beta: solved.flash.beta,
        x: solved.flash.x,
        y: solved.flash.y,
        k: solved.flash.k,
        phase: solved.flash.phase,
        z_liquid: solved.flash.z_liquid,
        z_vapour: solved.flash.z_vapour,
        iterations: solved.iterations,
        residual: solved.residual,
        warnings: solved.warnings,
    })
}
