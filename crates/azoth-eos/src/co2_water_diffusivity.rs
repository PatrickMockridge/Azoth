//! `eos.co2_water_diffusivity` - the CO2-in-water binary diffusivity from the
//! Tammi correlation.
//!
//! Spec: `specs/calcs/eos/co2_water_diffusivity.toml`, which records the exponential
//! correlation.

use azoth_core::units::{ThermodynamicTemperature, square_meters_per_second};
use azoth_core::{Result, apply_checks};

use crate::results::Co2WaterDiffusivityResult;
use crate::spec_gen;

/// The CO2-in-water binary diffusion coefficient, from NeqSim's `CO2water`
/// correlation.
///
/// The correlation is temperature-only: it carries no solute or solvent argument,
/// which is what NeqSim's method does - it returns the same value whatever pair it
/// is asked for.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `T` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::kelvins;
/// use azoth_eos::co2_water_diffusivity;
///
/// let r = co2_water_diffusivity(kelvins(298.15))?;
/// assert!((r.d.value - 2.0208212539579613e-9).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T` is the symbol in the published equation
pub fn co2_water_diffusivity(T: ThermodynamicTemperature) -> Result<Co2WaterDiffusivityResult> {
    let spec = &spec_gen::CO2_WATER_DIFFUSIVITY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| (quantity == "T").then_some(T.value),
        &mut warnings,
    )?;

    let d = 0.03389 * (-2213.7 / T.value).exp() * 1.0e-4;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "d").then_some(d),
        &mut warnings,
    )?;

    Ok(Co2WaterDiffusivityResult {
        d: square_meters_per_second(d),
        warnings,
    })
}
