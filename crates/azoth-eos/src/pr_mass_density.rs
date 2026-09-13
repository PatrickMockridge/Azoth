//! `eos.pr_mass_density` - mass density from a molar volume.
//!
//! ```text
//! rho = M/v
//! ```
//!
//! Spec: `specs/calcs/eos/pr_mass_density.yaml`
//!
//! # Why one division earns a calc
//!
//! It is the last step of the path a caller actually wants - `pr_z_factor` gives a
//! compressibility factor, `pr_molar_volume` turns it into a volume, and this turns
//! that into the density somebody asked for. It is also the step where the units are
//! load-bearing in both directions, and the one that gives `kg/mol` its first
//! consumer: that unit has been in the vocabulary since the vocabulary was made
//! checkable, with a correct conversion path and no calculation using it. A unit
//! whose conversion is tested is not a defect, but it is a unit nothing has ever
//! proved works end to end. This is that proof.
//!
//! # The unit that catches people
//!
//! Molar masses are conventionally quoted in `g/mol` and this takes `kg/mol`, a
//! factor of 1000 apart. The units layer cannot catch it - both are plausible
//! magnitudes of the same dimension - so the spec says so in the input's
//! description, in the derivation and in the assumptions, which is as much as a
//! comment can do.

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
