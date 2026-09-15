//! `eos.chung_viscosity` - the gas viscosity from the Chung correlation.
//!
//! Spec: `specs/calcs/eos/chung_viscosity.toml`, which records the ten coefficient
//! rows, the collision integral, and the dense-gas correction.

use azoth_core::units::{MolarMass, MolarVolume, ThermodynamicTemperature, pascal_seconds};
use azoth_core::{Result, apply_checks};

use crate::results::ChungViscosityResult;
use crate::spec_gen;

/// The ten rows of Chung's high-pressure coefficient table (TPoLG Table 9-5).
/// Each row is `[a, b, c, d]` for `E_j = a + b*omega + c*mu_r**4 + d*kappa`.
const CHUNG_HP: [[f64; 4]; 10] = [
    [6.324, 50.412, -51.680, 1189.0],
    [1.210e-3, -1.154e-3, -6.257e-3, 0.03728],
    [5.283, 254.209, -168.48, 3898.0],
    [6.623, 38.096, -8.464, 31.42],
    [19.745, 7.630, -14.354, 31.53],
    [-1.9, -12.537, 4.985, -18.15],
    [24.275, 3.450, -11.291, 69.35],
    [0.7972, 1.117, 0.01235, -4.117],
    [-0.2382, 0.06770, -0.8163, 4.025],
    [0.06863, 0.3479, 0.5926, -0.727],
];

/// The gas dynamic viscosity of a pure component, from Chung's correlation.
///
/// `dipole` is in debye and `kappa` is the dimensionless viscosity correction
/// factor; both come from the component databank. `V` is the gas's molar volume,
/// which the dense-gas correction `y = Vc/(6*V)` consumes.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tc`, `Vc`, `T` or `V` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, cubic_meters_per_mole, kilograms_per_mole};
/// use azoth_eos::chung_viscosity;
///
/// let r = chung_viscosity(0.0115, kelvins(190.56), cubic_meters_per_mole(9.9e-5),
///     kilograms_per_mole(0.016043), 0.0, 0.0, kelvins(300.0),
///     cubic_meters_per_mole(2.4409707154781444e-3))?;
/// assert!((r.mu.value - 1.124009415430133e-5).abs() < 1e-18);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tc`, `Vc`, `M`, `T` and `V` are the symbols in the published equation
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn chung_viscosity(
    omega: f64,
    Tc: ThermodynamicTemperature,
    Vc: MolarVolume,
    M: MolarMass,
    dipole: f64,
    kappa: f64,
    T: ThermodynamicTemperature,
    V: MolarVolume,
) -> Result<ChungViscosityResult> {
    let spec = &spec_gen::CHUNG_VISCOSITY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "omega" => Some(omega),
            "Tc" => Some(Tc.value),
            "Vc" => Some(Vc.value),
            "M" => Some(M.value),
            "dipole" => Some(dipole),
            "kappa" => Some(kappa),
            "T" => Some(T.value),
            "V" => Some(V.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // The published constants are tuned to Vc in cm**3/mol and M in g/mol.
    let vc_cm3 = Vc.value * 1.0e6;
    let m_g = M.value * 1000.0;

    let rel_visc = 131.3 * dipole / (vc_cm3 * Tc.value).sqrt();
    let rel4 = rel_visc * rel_visc * rel_visc * rel_visc;
    let fc = 1.0 - 0.2756 * omega + 0.059035 * rel4 + kappa;
    let e: [f64; 10] = CHUNG_HP.map(|row| row[0] + row[1] * omega + row[2] * rel4 + row[3] * kappa);

    let temp_var = 1.2593 * T.value / Tc.value;
    let omega_visc = 1.16145 * temp_var.powf(-0.14874)
        + 0.52487 * (-0.77320 * temp_var).exp()
        + 2.16178 * (-2.43787 * temp_var).exp();

    let chungy = Vc.value / (6.0 * V.value);
    let g1 = (1.0 - 0.5 * chungy) / (1.0 - chungy).powi(3);
    let g2 = (e[0] * ((1.0 - (-e[3] * chungy).exp()) / chungy)
        + e[1] * g1 * (e[4] * chungy).exp()
        + e[2] * g1)
        / (e[0] * e[3] + e[1] + e[2]);
    let viskstarstar =
        e[6] * chungy * chungy * g2 * (e[7] + e[8] / temp_var + e[9] * temp_var.powi(-2)).exp();
    let viskstar = temp_var.sqrt() / omega_visc * (fc * (1.0 / g2 + e[5] * chungy)) + viskstarstar;

    // microPoise, then SI.
    let mu_micropoise = viskstar * 36.344 * (m_g * Tc.value).sqrt() / vc_cm3.powf(2.0 / 3.0);
    let mu = mu_micropoise * 1.0e-7;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "mu").then_some(mu),
        &mut warnings,
    )?;

    Ok(ChungViscosityResult {
        mu: pascal_seconds(mu),
        warnings,
    })
}
