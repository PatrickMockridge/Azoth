//! `eos.chung_conductivity` - the gas thermal conductivity from the Chung
//! correlation.
//!
//! Spec: `specs/calcs/eos/chung_conductivity.toml`, which records the dilute-gas
//! viscosity, the Eucken correction and the gas constant.

use azoth_core::units::{
    MolarHeatCapacity, MolarMass, MolarVolume, ThermodynamicTemperature, watts_per_meter_kelvin,
};
use azoth_core::{Result, apply_checks};

use crate::results::ChungConductivityResult;
use crate::spec_gen;

/// NeqSim's gas constant, not the exact SI value. The correlation is tuned to this
/// value and the port reproduces NeqSim to the last digit only by using it.
const R: f64 = 8.3144621;

/// The gas thermal conductivity of a pure component, from Chung's correlation.
///
/// `Cv0` is the ideal-gas heat capacity at constant volume at `T` (`Cp0 - R`);
/// `dipole` is in debye and `kappa` the dimensionless viscosity correction factor.
/// The viscosity inside the correlation is the *dilute-gas* Chung viscosity, not
/// `eos.chung_viscosity`'s dense-gas form.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tc`, `Vc` or `T` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, cubic_meters_per_mole, kilograms_per_mole,
///     joules_per_mole_kelvin};
/// use azoth_eos::chung_conductivity;
///
/// let r = chung_conductivity(joules_per_mole_kelvin(27.544151394),
///     kilograms_per_mole(0.016043), 0.0115, kelvins(190.56), cubic_meters_per_mole(9.9e-5),
///     0.0, 0.0, kelvins(300.0))?;
/// assert!((r.k.value - 0.03384023073765177).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Cv0`, `Tc`, `Vc`, `M` and `T` are the symbols in the published equation
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
#[allow(clippy::approx_constant)] // `0.6366` is Chung's coefficient, not 2/pi.
pub fn chung_conductivity(
    Cv0: MolarHeatCapacity,
    M: MolarMass,
    omega: f64,
    Tc: ThermodynamicTemperature,
    Vc: MolarVolume,
    dipole: f64,
    kappa: f64,
    T: ThermodynamicTemperature,
) -> Result<ChungConductivityResult> {
    let spec = &spec_gen::CHUNG_CONDUCTIVITY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "Cv0" => Some(Cv0.value),
            "M" => Some(M.value),
            "omega" => Some(omega),
            "Tc" => Some(Tc.value),
            "Vc" => Some(Vc.value),
            "dipole" => Some(dipole),
            "kappa" => Some(kappa),
            "T" => Some(T.value),
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

    let t_star = 1.2593 * T.value / Tc.value;
    let omega_v = 1.16145 * t_star.powf(-0.14874)
        + 0.52487 * (-0.77320 * t_star).exp()
        + 2.16178 * (-2.43787 * t_star).exp();

    // The dilute-gas Chung viscosity, in microPoise then Pa*s.
    let eta0_micropoise = 40.785 * fc * (m_g * T.value).sqrt() / (vc_cm3.powf(2.0 / 3.0) * omega_v);
    let eta0 = eta0_micropoise * 1.0e-7;

    let alpha = Cv0.value / R - 1.5;
    let beta = 0.7862 - 0.7109 * omega + 1.3168 * omega * omega;
    let z = 2.0 + 10.5 * (T.value / Tc.value).powi(2);
    let psi = 1.0
        + alpha
            * ((0.215 + 0.28288 * alpha - 1.061 * beta + 0.26665 * z)
                / (0.6366 + beta * z + 1.061 * alpha * beta));

    let k = 3.75 * R * psi * eta0 / M.value;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "k").then_some(k),
        &mut warnings,
    )?;

    Ok(ChungConductivityResult {
        k: watts_per_meter_kelvin(k),
        warnings,
    })
}
