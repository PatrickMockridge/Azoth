//! `eos.fuller_schettler_giddings_diffusivity` - the gas binary diffusivity from the
//! Fuller-Schettler-Giddings correlation.
//!
//! Spec: `specs/calcs/eos/fuller_schettler_giddings_diffusivity.toml`, which records the
//! correlation, the diffusion-volume ladder and the two ways the class states its pressure unit.

use azoth_core::units::{
    MolarMass, MolarVolume, Pressure, ThermodynamicTemperature, square_meters_per_second,
};
use azoth_core::{Result, apply_checks};

use crate::results::FullerSchettlerGiddingsDiffusivityResult;
use crate::spec_gen;

/// The binary diffusion coefficient of a gas pair, from the Fuller-Schettler-Giddings
/// correlation.
///
/// `MA`/`MB` are the pair's molar masses and `VA`/`VB` their *diffusion volumes* - the
/// class's own ladder, [`crate::databank::fuller_diffusion_volume`], resolves a name to one.
///
/// **The pressure is in bar and the constant is `1.013e-3`**, which is the pair
/// `calcBinaryDiffusionCoefficient` uses; the class's javadoc writes `1.013e-2` with `P` in
/// atm, ten times the code's constant and inconsistent with its own unit. The port follows
/// the code, and the capture is what says so.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `T` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{cubic_meters_per_mole, kelvins, kilograms_per_mole, pascals};
/// use azoth_eos::fuller_schettler_giddings_diffusivity;
///
/// let r = fuller_schettler_giddings_diffusivity(
///     kilograms_per_mole(0.016043),
///     kilograms_per_mole(0.0280135),
///     cubic_meters_per_mole(2.514e-5),
///     cubic_meters_per_mole(1.85e-5),
///     kelvins(298.15),
///     pascals(101325.0),
/// )?;
/// assert!((r.d.value - 2.155058583907974e-5).abs() < 1e-18);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `MA`, `MB`, `VA`, `VB` and `T` are the symbols in the source
pub fn fuller_schettler_giddings_diffusivity(
    MA: MolarMass,
    MB: MolarMass,
    VA: MolarVolume,
    VB: MolarVolume,
    T: ThermodynamicTemperature,
    P: Pressure,
) -> Result<FullerSchettlerGiddingsDiffusivityResult> {
    let spec = &spec_gen::FULLER_SCHETTLER_GIDDINGS_DIFFUSIVITY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "MA" => Some(MA.value),
            "MB" => Some(MB.value),
            "VA" => Some(VA.value),
            "VB" => Some(VB.value),
            "T" => Some(T.value),
            "P" => Some(P.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // The published constants are tuned to g/mol, cm**3/mol, bar and cm**2/s.
    let ma_g = MA.value * 1000.0;
    let mb_g = MB.value * 1000.0;
    let va_cm3 = VA.value * 1.0e6;
    let vb_cm3 = VB.value * 1.0e6;
    let p_bar = P.value * 1.0e-5;

    let sigma_v = va_cm3.cbrt() + vb_cm3.cbrt();
    let d_cm2s = 1.013e-3 * T.value.powf(1.75) * (1.0 / ma_g + 1.0 / mb_g).sqrt()
        / (p_bar * sigma_v * sigma_v);
    let d = d_cm2s * 1.0e-4;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "d").then_some(d),
        &mut warnings,
    )?;

    Ok(FullerSchettlerGiddingsDiffusivityResult {
        d: square_meters_per_second(d),
        warnings,
    })
}
