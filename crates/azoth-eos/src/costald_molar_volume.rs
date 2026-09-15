//! `eos.costald_molar_volume` - the saturated liquid molar volume from the COSTALD
//! equation.
//!
//! Spec: `specs/calcs/eos/costald_molar_volume.toml`, which records the two dimensionless
//! volume functions, the characteristic-volume back-calculation, and what is not applied
//! (the compressed-liquid and polar corrections).

use azoth_core::units::{
    MassDensity, MolarMass, MolarVolume, ThermodynamicTemperature, cubic_meters_per_mole,
};
use azoth_core::{Result, apply_checks};

use crate::results::CostaldMolarVolumeResult;
use crate::spec_gen;

/// The dimensionless saturated-volume function `V_R^(0)`.
fn vr0(tr: f64) -> f64 {
    let tau = 1.0 - tr;
    if tau <= 0.0 {
        return 1.0;
    }
    let t13 = tau.cbrt();
    1.0 - 1.52816 * t13 + 1.43907 * t13 * t13 - 0.81446 * tau + 0.190454 * t13 * tau
}

/// The dimensionless departure-volume function `V_R^(delta)`.
fn vrdelta(tr: f64) -> f64 {
    (-0.296123 + 0.386914 * tr - 0.0427258 * tr * tr - 0.0480645 * tr * tr * tr) / (tr - 1.00001)
}

/// The saturated liquid molar volume of a pure component, from the Hankinson-Thomson
/// COSTALD equation.
///
/// `Vc` and `rho_normal` feed the characteristic volume: it is back-calculated from
/// the normal liquid density at 288.71 K, falling back to `Vc` when the density is
/// absent or the standard reduced temperature reaches 0.9.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tc`, `Vc` or `T` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, cubic_meters_per_mole, kilograms_per_cubic_meter, kilograms_per_mole};
/// use azoth_eos::costald_molar_volume;
///
/// let r = costald_molar_volume(0.3013, kelvins(507.6), cubic_meters_per_mole(3.7e-4),
///     kilograms_per_mole(0.086177), kilograms_per_cubic_meter(664.0), kelvins(298.15))?;
/// assert!((r.v.value - 0.0001315006475383877).abs() < 1e-18);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tc`, `Vc`, `M` and `T` are the symbols in the published equation
pub fn costald_molar_volume(
    omega: f64,
    Tc: ThermodynamicTemperature,
    Vc: MolarVolume,
    M: MolarMass,
    rho_normal: MassDensity,
    T: ThermodynamicTemperature,
) -> Result<CostaldMolarVolumeResult> {
    let spec = &spec_gen::COSTALD_MOLAR_VOLUME_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "omega" => Some(omega),
            "Tc" => Some(Tc.value),
            "Vc" => Some(Vc.value),
            "M" => Some(M.value),
            "rho_normal" => Some(rho_normal.value),
            "T" => Some(T.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let tc = Tc.value;
    let mut v_star = Vc.value;
    if rho_normal.value > 0.0 {
        let tr_std = 288.71 / tc;
        if tr_std < 0.9 {
            let factor = vr0(tr_std) * (1.0 - omega * vrdelta(tr_std));
            if factor.abs() > 1e-15 {
                let v_molar_std = M.value / rho_normal.value;
                let estimate = v_molar_std / factor;
                if estimate > 0.0 {
                    v_star = estimate;
                }
            }
        }
    }

    let tr = T.value / tc;
    let v = v_star * vr0(tr) * (1.0 - omega * vrdelta(tr));

    apply_checks(
        spec.derived_checks(),
        |name| (name == "v").then_some(v),
        &mut warnings,
    )?;

    Ok(CostaldMolarVolumeResult {
        v: cubic_meters_per_mole(v),
        warnings,
    })
}
