//! `eos.vh_flash` - the pressure and temperature a mixture reaches at a given volume and
//! enthalpy.
//!
//! Spec: `specs/models/eos/vh_flash.toml`
//!
//! A vessel at fixed volume and enthalpy: both state variables are solved for, by the same
//! decoupled 2x2 Newton [`crate::vu_flash`] uses, with the energy target supplied whole -
//! NeqSim's `VHflashQfunc`, which is its `OptimizedVUflash` with `Hspec` in place of
//! `Uspec + P Vspec`.

use azoth_core::units::{MolarEnergy, MolarVolume, kelvins, pascals};
use azoth_core::{Result, apply_checks};

use crate::algorithm_of;
use crate::flash_property::{EnergyTarget, solve_pressure_temperature};
use crate::mixture::Mixture;
use crate::model_gen;
use crate::molar_enthalpy_entropy::IdealGasModel;
use crate::pr_molar_volume::MOLAR_GAS_CONSTANT;
use crate::results::VhFlashResult;

/// The pressure and temperature at which a mixture has a given molar volume and enthalpy.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `v` is not positive.
/// * [`azoth_core::AzothError::SolverNotConverged`] if the iteration reaches its cap.
pub fn vh_flash(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    v: MolarVolume,
    h: MolarEnergy,
    z: &[f64],
) -> Result<VhFlashResult> {
    let spec = &model_gen::VH_FLASH_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "V" => Some(v.value),
            "H" => Some(h.value),
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
        EnergyTarget::Enthalpy(h.value),
        z,
        algorithm,
        start_p,
        start_t,
    )?;

    Ok(VhFlashResult {
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
