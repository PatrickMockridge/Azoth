//! `eos.pr_molar_volume` - molar volume from a compressibility factor.
//!
//! ```text
//! v = z*R*T/P
//! ```
//!
//! Spec: `specs/calcs/eos/pr_molar_volume.toml`, which carries why this is the one
//! dimensional calc in the namespace, why the gas constant is computed from its two
//! defining constants, and why the evaluation order matters.

use azoth_core::units::{Pressure, ThermodynamicTemperature, cubic_meters_per_mole};
use azoth_core::{Result, apply_checks};

use crate::results::PrMolarVolumeResult;
use crate::spec_gen;

/// The Avogadro constant, `6.02214076e23 /mol`. Exact by the 2019 SI definition.
pub const AVOGADRO_PER_MOL: f64 = 6.022_140_76e23;

/// The Boltzmann constant, `1.380649e-23 J/K`. Exact by the 2019 SI definition.
pub const BOLTZMANN_J_PER_K: f64 = 1.380_649e-23;

/// The molar gas constant, `N_A * k_B`. Exact, because both factors are.
pub const MOLAR_GAS_CONSTANT: f64 = AVOGADRO_PER_MOL * BOLTZMANN_J_PER_K;

/// Molar volume at a state, from its compressibility factor.
///
/// `z` comes from [`crate::pr_z_factor`]; either admissible root may be passed, the
/// vapour one giving the vapour volume and the liquid one the liquid volume. That
/// is why this takes a `z` rather than a phase - the caller already chose what the
/// answer means.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `z <= 0`, or if `T` or `P` is not
///   positive. A molar volume is positive, so a non-positive input can only be a
///   caller error, and `P` is a divisor besides.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, pascals};
/// use azoth_eos::pr_molar_volume;
///
/// let r = pr_molar_volume(0.7907789662973796, kelvins(295.864), pascals(1_062_000.0))?;
/// assert!((r.v.value - 0.0018317107825229842).abs() < 1e-18);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T` and `P` are the symbols in the published equation
pub fn pr_molar_volume(
    z: f64,
    T: ThermodynamicTemperature,
    P: Pressure,
) -> Result<PrMolarVolumeResult> {
    let spec = &spec_gen::PR_MOLAR_VOLUME_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "z" => Some(z),
            "T" => Some(T.value),
            "P" => Some(P.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // Written in the equation's own order, left to right; the spec's worked example
    // records why the order is stated rather than incidental.
    let v = z * MOLAR_GAS_CONSTANT * T.value / P.value;

    apply_checks(
        spec.derived_checks(),
        |quantity| (quantity == "v").then_some(v),
        &mut warnings,
    )?;

    Ok(PrMolarVolumeResult {
        v: cubic_meters_per_mole(v),
        warnings,
    })
}
