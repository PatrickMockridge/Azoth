//! `eos.pr_mass_density` - mass density from a molar volume.
//!
//! ```text
//! rho = M/v
//! ```
//!
//! Spec: `specs/calcs/eos/pr_mass_density.toml`, which carries why one division earns
//! a calc and the `g/mol`-versus-`kg/mol` trap that the units layer cannot catch.

use azoth_core::units::{MolarMass, MolarVolume, kilograms_per_cubic_meter};
use azoth_core::{Result, apply_checks};

use crate::results::PrMassDensityResult;
use crate::spec_gen;

/// Mass density from a molar mass and a molar volume.
///
/// `M` is the caller's - this library ships no component data - and `v` comes from
/// [`crate::pr_molar_volume`].
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `M <= 0` or `v <= 0`. Both are
///   strictly positive for any real substance and state, and `v` is a divisor.
///
/// # Example
/// ```
/// use azoth_core::units::{cubic_meters_per_mole, kilograms_per_mole};
/// use azoth_eos::pr_mass_density;
///
/// let r = pr_mass_density(
///     kilograms_per_mole(0.0440956),
///     cubic_meters_per_mole(0.0018317107825229842),
/// )?;
/// assert!((r.rho.value - 24.073451125981286).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `M` is the symbol in the equation
pub fn pr_mass_density(M: MolarMass, v: MolarVolume) -> Result<PrMassDensityResult> {
    let spec = &spec_gen::PR_MASS_DENSITY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "M" => Some(M.value),
            "v" => Some(v.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let rho = M.value / v.value;

    apply_checks(
        spec.derived_checks(),
        |quantity| (quantity == "rho").then_some(rho),
        &mut warnings,
    )?;

    Ok(PrMassDensityResult {
        rho: kilograms_per_cubic_meter(rho),
        warnings,
    })
}
