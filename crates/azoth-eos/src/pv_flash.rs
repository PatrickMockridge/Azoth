//! `eos.pv_flash` - the temperature a mixture reaches at a given pressure and molar volume.
//!
//! Spec: `specs/models/eos/pv_flash.toml`
//!
//! An outer quasi-Newton on temperature, over the molar volume assembled from the phase
//! split at a trial temperature ([`crate::pt_flash`]). The volume rises monotonically with
//! temperature at a fixed pressure. The loop lives in [`crate::flash_property`].

use azoth_core::units::{MolarVolume, Pressure, kelvins};
use azoth_core::{Result, apply_checks};

use crate::algorithm_of;
use crate::flash_property::{Property, solve_temperature};
use crate::mixture::Mixture;
use crate::model_gen;
use crate::molar_enthalpy_entropy::IdealGasModel;
use crate::results::PvFlashResult;

/// The temperature at which a mixture has a given molar volume at a pressure.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `p` or `v` is not positive.
/// * [`azoth_core::AzothError::SolverNotConverged`] if the iteration reaches its cap.
///
/// # Example
/// ```
/// use azoth_eos::Cubic;
/// use azoth_core::units::{cubic_meters_per_mole, pascals};
/// use azoth_eos::{databank, pv_flash};
///
/// let (mixture, ideal_gas) = databank::mixture_of(&["methane", "n-butane"], Cubic::Pr, None)
///     .expect("the pair resolves");
/// let r = pv_flash::pv_flash(
///     &mixture,
///     &ideal_gas,
///     pascals(1.0e6),
///     cubic_meters_per_mole(1.9610505645240935e-3),
///     &[0.6, 0.4],
/// )?;
/// assert!((r.temperature.value - 300.0).abs() < 1.0e-3);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn pv_flash(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    p: Pressure,
    v: MolarVolume,
    z: &[f64],
) -> Result<PvFlashResult> {
    let spec = &model_gen::PV_FLASH_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "P" => Some(p.value),
            "V" => Some(v.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let algorithm = algorithm_of(spec)?;
    let solved = solve_temperature(
        mixture,
        ideal_gas,
        p,
        v.value,
        z,
        Property::Volume,
        algorithm,
    )?;

    apply_checks(
        spec.derived_checks(),
        |quantity| (quantity == "T").then_some(solved.temperature),
        &mut warnings,
    )?;

    Ok(PvFlashResult {
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
