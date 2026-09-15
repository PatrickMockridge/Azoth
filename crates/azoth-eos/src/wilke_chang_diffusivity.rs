//! `eos.wilke_chang_diffusivity` - the liquid binary diffusivity from the
//! Wilke-Chang correlation.
//!
//! Spec: `specs/calcs/eos/wilke_chang_diffusivity.toml`, which records the
//! correlation, the association parameter and the clamps.

use azoth_core::units::{
    DynamicViscosity, MolarMass, MolarVolume, ThermodynamicTemperature, square_meters_per_second,
};
use azoth_core::{Result, apply_checks};

use crate::results::WilkeChangDiffusivityResult;
use crate::spec_gen;

/// The binary diffusion coefficient at infinite dilution, from the Wilke-Chang
/// correlation.
///
/// `phi` is the solvent association parameter (see [`crate::databank::wilke_chang_phi`]),
/// `M` the solvent molar mass and `VA` the solute molar volume at its normal boiling
/// point. `VA` and `eta` are clamped to NeqSim's `[20, 600]` cm**3/mol and
/// `[0.01, 500]` cP before the correlation is applied.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `T` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{cubic_meters_per_mole, kelvins, kilograms_per_mole, pascal_seconds};
/// use azoth_eos::wilke_chang_diffusivity;
///
/// let r = wilke_chang_diffusivity(2.26, kilograms_per_mole(0.018015), kelvins(298.15),
///     pascal_seconds(8.915447896200597e-4), cubic_meters_per_mole(4.0203262233375156e-5))?;
/// assert!((r.d.value - 1.7212261801805907e-9).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `M`, `T` and `VA` are the symbols in the published equation
pub fn wilke_chang_diffusivity(
    phi: f64,
    M: MolarMass,
    T: ThermodynamicTemperature,
    eta: DynamicViscosity,
    VA: MolarVolume,
) -> Result<WilkeChangDiffusivityResult> {
    let spec = &spec_gen::WILKE_CHANG_DIFFUSIVITY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "phi" => Some(phi),
            "M" => Some(M.value),
            "T" => Some(T.value),
            "eta" => Some(eta.value),
            "VA" => Some(VA.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // The published constants are tuned to g/mol, cm**3/mol and cP, and NeqSim
    // clamps the volume and viscosity to its reasonable limits.
    let m_g = M.value * 1000.0;
    let va_cm3 = (VA.value * 1.0e6).clamp(20.0, 600.0);
    let eta_cp = (eta.value * 1000.0).clamp(0.01, 500.0);

    let d_cm2s = 7.4e-8 * (phi * m_g).sqrt() * T.value / (eta_cp * va_cm3.powf(0.6));
    let d = d_cm2s * 1.0e-4;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "d").then_some(d),
        &mut warnings,
    )?;

    Ok(WilkeChangDiffusivityResult {
        d: square_meters_per_second(d),
        warnings,
    })
}
