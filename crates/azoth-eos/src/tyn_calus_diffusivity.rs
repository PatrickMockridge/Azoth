//! `eos.tyn_calus_diffusivity` - the liquid binary diffusivity from the Tyn-Calus
//! correlation.
//!
//! Spec: `specs/calcs/eos/tyn_calus_diffusivity.toml`, which records the correlation
//! and the molar-volume and viscosity clamps.

use azoth_core::units::{
    DynamicViscosity, MolarVolume, ThermodynamicTemperature, square_meters_per_second,
};
use azoth_core::{Result, apply_checks};

use crate::results::TynCalusDiffusivityResult;
use crate::spec_gen;

/// The binary diffusion coefficient at infinite dilution, from the Tyn-Calus
/// correlation.
///
/// `VA` and `VB` are the solute and solvent liquid molar volumes at the normal
/// boiling point, and `eta` the solvent viscosity. All three are clamped to
/// NeqSim's `[20, 600]` cm**3/mol and `[0.01, 500]` cP before the correlation is
/// applied.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `T` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{cubic_meters_per_mole, kelvins, pascal_seconds};
/// use azoth_eos::tyn_calus_diffusivity;
///
/// let r = tyn_calus_diffusivity(cubic_meters_per_mole(8.816478555304741e-5),
///     cubic_meters_per_mole(1.0578760045924226e-4), kelvins(298.15),
///     pascal_seconds(9.163501315189954e-4))?;
/// assert!((r.d.value - 1.4502274002446276e-9).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `VA`, `VB`, `T` and `D` are the symbols in the published equation
pub fn tyn_calus_diffusivity(
    VA: MolarVolume,
    VB: MolarVolume,
    T: ThermodynamicTemperature,
    eta: DynamicViscosity,
) -> Result<TynCalusDiffusivityResult> {
    let spec = &spec_gen::TYN_CALUS_DIFFUSIVITY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "VA" => Some(VA.value),
            "VB" => Some(VB.value),
            "T" => Some(T.value),
            "eta" => Some(eta.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // The published constants are tuned to cm**3/mol and cP, and NeqSim clamps both
    // to its reasonable limits before applying the correlation.
    let va_cm3 = (VA.value * 1.0e6).clamp(20.0, 600.0);
    let vb_cm3 = (VB.value * 1.0e6).clamp(20.0, 600.0);
    let eta_cp = (eta.value * 1000.0).clamp(0.01, 500.0);

    let d_cm2s = 8.93e-8 * vb_cm3.powf(0.267) * T.value / (eta_cp * va_cm3.powf(0.433));
    let d = d_cm2s * 1.0e-4;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "d").then_some(d),
        &mut warnings,
    )?;

    Ok(TynCalusDiffusivityResult {
        d: square_meters_per_second(d),
        warnings,
    })
}
