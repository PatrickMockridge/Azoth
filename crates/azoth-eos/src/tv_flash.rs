//! `eos.tv_flash` - the pressure a mixture reaches at a given temperature and molar volume.
//!
//! Spec: `specs/models/eos/tv_flash.toml`
//!
//! An outer quasi-Newton on pressure, over the molar volume assembled from the phase split
//! at a trial pressure ([`crate::pt_flash`]) and each phase's compressibility. The volume
//! falls monotonically with pressure at a fixed temperature, which is what makes the
//! inversion well posed. The loop lives in [`crate::flash_property`].

use azoth_core::units::{MolarVolume, ThermodynamicTemperature, pascals};
use azoth_core::{Result, apply_checks};

use crate::algorithm_of;
use crate::flash_property::{Property, solve_pressure};
use crate::mixture::Mixture;
use crate::model_gen;
use crate::molar_enthalpy_entropy::IdealGasModel;
use crate::pr_molar_volume::MOLAR_GAS_CONSTANT;
use crate::results::TvFlashResult;

/// The pressure at which a mixture has a given molar volume at a temperature.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `t` or `v` is not positive.
/// * [`azoth_core::AzothError::SolverNotConverged`] if the iteration reaches its cap.
///
/// # Example
/// ```
/// use azoth_core::units::{cubic_meters_per_mole, kelvins};
/// use azoth_eos::{databank, tv_flash};
///
/// let (mixture, ideal_gas) = databank::mixture_of(&["methane", "n-butane"], None)
///     .expect("the pair resolves");
/// let r = tv_flash::tv_flash(
///     &mixture,
///     &ideal_gas,
///     kelvins(300.0),
///     cubic_meters_per_mole(1.9610505645240935e-3),
///     &[0.6, 0.4],
/// )?;
/// assert!((r.pressure.value - 1.0e6).abs() < 10.0);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn tv_flash(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t: ThermodynamicTemperature,
    v: MolarVolume,
    z: &[f64],
) -> Result<TvFlashResult> {
    let spec = &model_gen::TV_FLASH_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "V" => Some(v.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let algorithm = algorithm_of(spec)?;
    // The ideal-gas pressure for this temperature and volume: a first guess within an
    // order of magnitude of the true pressure, for either a gas or a condensed phase.
    let start = MOLAR_GAS_CONSTANT * t.value / v.value;
    let solved = solve_pressure(
        mixture,
        ideal_gas,
        t,
        v.value,
        z,
        Property::Volume,
        algorithm,
        start,
    )?;

    apply_checks(
        spec.derived_checks(),
        |quantity| (quantity == "P").then_some(solved.pressure),
        &mut warnings,
    )?;

    Ok(TvFlashResult {
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
