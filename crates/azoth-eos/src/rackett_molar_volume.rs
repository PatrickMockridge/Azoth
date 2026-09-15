//! `eos.rackett_molar_volume` - the saturated liquid molar volume from the Rackett
//! equation.
//!
//! ```text
//! v = (R*Tc/Pc)*Z_RA**(1 + (1 - Tr)**(2/7)),  Z_RA = 0.29056 - 0.08775*omega
//! ```
//!
//! Spec: `specs/calcs/eos/rackett_molar_volume.toml`, which records the Yamada-Gunn
//! `Z_RA` and why the exact gas constant is used.

use azoth_core::units::{Pressure, ThermodynamicTemperature, cubic_meters_per_mole};
use azoth_core::{Result, apply_checks};

use crate::pr_molar_volume::MOLAR_GAS_CONSTANT;
use crate::results::RackettMolarVolumeResult;
use crate::spec_gen;

/// The saturated liquid molar volume of a pure component, from the Spencer-Danner
/// Rackett equation.
///
/// `omega`, `Tc` and `Pc` are the caller's.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tc`, `Pc` or `T` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, pascals};
/// use azoth_eos::rackett_molar_volume;
///
/// let r = rackett_molar_volume(0.152, kelvins(369.83), pascals(4_248_000.0), kelvins(298.15))?;
/// assert!((r.v.value - 8.991467403942754e-05).abs() < 1e-18);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tc`, `Pc` and `T` are the symbols in the published equation
pub fn rackett_molar_volume(
    omega: f64,
    Tc: ThermodynamicTemperature,
    Pc: Pressure,
    T: ThermodynamicTemperature,
) -> Result<RackettMolarVolumeResult> {
    let spec = &spec_gen::RACKETT_MOLAR_VOLUME_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "omega" => Some(omega),
            "Tc" => Some(Tc.value),
            "Pc" => Some(Pc.value),
            "T" => Some(T.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let z_ra = 0.29056 - 0.08775 * omega;
    let tr = T.value / Tc.value;
    let exponent = 1.0 + (1.0 - tr).powf(2.0 / 7.0);
    let v = MOLAR_GAS_CONSTANT * Tc.value / Pc.value * z_ra.powf(exponent);

    apply_checks(
        spec.derived_checks(),
        |name| (name == "v").then_some(v),
        &mut warnings,
    )?;

    Ok(RackettMolarVolumeResult {
        v: cubic_meters_per_mole(v),
        warnings,
    })
}
