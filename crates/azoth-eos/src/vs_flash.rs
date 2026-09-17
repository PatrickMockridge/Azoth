//! `eos.vs_flash` - the pressure and temperature a mixture reaches at a given volume and
//! entropy.
//!
//! Spec: `specs/models/eos/vs_flash.toml`
//!
//! A vessel at fixed volume and entropy: both state variables are solved for, by the same
//! decoupled 2x2 Newton [`crate::vu_flash`] uses, with an entropy asked for instead of an
//! energy - NeqSim's `VSflash`, which is its `OptimizedVUflash` loop with a Q function
//! in the entropy rather than in the enthalpy.

use azoth_core::units::{MolarHeatCapacity, MolarVolume, kelvins, pascals};
use azoth_core::{Result, apply_checks};

use crate::algorithm_of;
use crate::flash_property::{EnergyTarget, solve_pressure_temperature};
use crate::mixture::Mixture;
use crate::model_gen;
use crate::molar_enthalpy_entropy::IdealGasModel;
use crate::pr_molar_volume::MOLAR_GAS_CONSTANT;
use crate::results::VsFlashResult;

/// The pressure and temperature at which a mixture has a given molar volume and entropy.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `v` is not positive.
/// * [`azoth_core::AzothError::SolverNotConverged`] if the iteration reaches its cap.
pub fn vs_flash(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    v: MolarVolume,
    s: MolarHeatCapacity,
    z: &[f64],
) -> Result<VsFlashResult> {
    let spec = &model_gen::VS_FLASH_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "V" => Some(v.value),
            "S" => Some(s.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let algorithm = algorithm_of(spec)?;
    let start_t = algorithm.initial_temperature.unwrap_or(300.0);
    let start_p = MOLAR_GAS_CONSTANT * start_t / v.value;
    let solved = solve_pressure_temperature(
        mixture,
        ideal_gas,
        v.value,
        EnergyTarget::Entropy(s.value),
        z,
        algorithm,
        start_p,
        start_t,
    )?;

    Ok(VsFlashResult {
        pressure: pascals(solved.pressure),
        temperature: kelvins(solved.temperature),
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
