//! `eos.friction_theory_conductivity` - NeqSim's `FrictionTheoryConductivityMethod`.
//!
//! The dilute-gas term is Chung's pure-component conductivity mixed by Mason-Saxena, and
//! the dense term converts the phase's viscosity excess into conductivity with the
//! Eucken factor: `lambda = lambda_0 + max(eta - eta_0, 0) (Cv0 + 1.25 R) / M`.

use azoth_core::units::{
    DynamicViscosity, MolarHeatCapacity, MolarMass, MolarVolume, ThermodynamicTemperature,
    watts_per_meter_kelvin,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::chung_conductivity::chung_dilute;
use crate::model_gen;
use crate::results::FrictionTheoryConductivityResult;

/// NeqSim's gas constant, which the correlation is tuned to.
const R: f64 = 8.3144621;

/// The dense-fluid thermal conductivity from a phase's viscosity.
///
/// **The phase's own viscosity is an input.** The class reads `phase.getViscosity()` -
/// whatever model the phase already carries - so this calc does not run one; a caller
/// takes the other id's answer and passes it. `Cv0` is the ideal-gas heat capacity at
/// constant volume (`Cp0 - R`, from `eos.ideal_gas_cp`).
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] for vectors that are not one entry per
///   component, or a composition that is not the length of the rest.
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's range.
///
/// # Example
/// ```
/// use azoth_core::units::{
///     cubic_meters_per_mole, joules_per_mole_kelvin, kelvins, kilograms_per_mole,
///     pascal_seconds,
/// };
/// use azoth_eos::friction_theory_conductivity;
///
/// let r = friction_theory_conductivity(
///     &[joules_per_mole_kelvin(27.544151394)],
///     &[kelvins(190.56)],
///     &[cubic_meters_per_mole(9.9e-5)],
///     &[0.0115],
///     &[0.0],
///     &[kilograms_per_mole(0.016043)],
///     &[1.0],
///     kelvins(300.0),
///     pascal_seconds(1.1403910448332615e-5),
/// )?;
/// assert!((r.conductivity.value - 0.034385915611225566).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T` is the symbol in the published equation
#[allow(clippy::too_many_arguments)]
// The signature is the spec's declared inputs.
// The parameter names are the spec's inputs, which is what lets the binding be generated.
#[allow(clippy::similar_names)]
pub fn friction_theory_conductivity(
    cv0: &[MolarHeatCapacity],
    tc: &[ThermodynamicTemperature],
    critical_volume: &[MolarVolume],
    omega: &[f64],
    dipole: &[f64],
    molar_mass: &[MolarMass],
    x: &[f64],
    T: ThermodynamicTemperature,
    viscosity: DynamicViscosity,
) -> Result<FrictionTheoryConductivityResult> {
    let spec = &model_gen::FRICTION_THEORY_CONDUCTIVITY_SPEC;
    let mut warnings = Vec::new();

    let count = x.len();
    if cv0.len() != count
        || tc.len() != count
        || critical_volume.len() != count
        || omega.len() != count
        || dipole.len() != count
        || molar_mass.len() != count
    {
        return Err(AzothError::invalid_input(
            "x",
            format!(
                "{count} mole fraction(s) against {} Cv0, {} Tc, {} Vc, {} omega, {} dipole \
                 and {} molar mass value(s): every one is one entry per component",
                cv0.len(),
                tc.len(),
                critical_volume.len(),
                omega.len(),
                dipole.len(),
                molar_mass.len()
            ),
        ));
    }
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "x" => x.iter().copied().reduce(f64::min),
            "T" => Some(T.value),
            "viscosity" => Some(viscosity.value),
            "cv0" => cv0.iter().map(|v| v.value).reduce(f64::min),
            "tc" => tc.iter().map(|v| v.value).reduce(f64::min),
            "critical_volume" => critical_volume.iter().map(|v| v.value).reduce(f64::min),
            "molar_mass" => molar_mass.iter().map(|v| v.value).reduce(f64::min),
            _ => None,
        },
        &mut warnings,
    )?;

    // Each component's dilute-gas conductivity and viscosity, from the one Chung block
    // `eos.chung_conductivity` also answers from.
    let mut pure_conductivity = Vec::with_capacity(count);
    let mut pure_viscosity = Vec::with_capacity(count);
    for i in 0..count {
        let (mut lambda, eta0) = chung_dilute(
            cv0[i].value,
            molar_mass[i].value,
            omega[i],
            tc[i].value,
            critical_volume[i].value,
            dipole[i],
            0.0,
            T.value,
        );
        if lambda < 1.0e-50 {
            lambda = 1.0e-50;
        }
        pure_conductivity.push(lambda);
        pure_viscosity.push(eta0);
    }

    // Mason-Saxena mixing. **The class's two mass exponents disagree**: the conductivity
    // rule raises `M_i/M_j` to a quarter and the viscosity rule `M_j/M_i`, which is the
    // class's own asymmetry and is reproduced rather than repaired.
    let mut dilute_conductivity = 0.0;
    for i in 0..count {
        let mut denominator = 0.0;
        for j in 0..count {
            let ratio = molar_mass[i].value / molar_mass[j].value;
            let aij = (1.0
                + (pure_conductivity[i] / pure_conductivity[j]).sqrt() * ratio.powf(0.25))
            .powi(2)
                / (8.0 * (1.0 + ratio)).sqrt();
            denominator += x[j] * aij;
        }
        dilute_conductivity += x[i] * pure_conductivity[i] / denominator;
    }

    let mut dilute_viscosity = 0.0;
    for i in 0..count {
        let mut denominator = 0.0;
        for j in 0..count {
            let ratio = molar_mass[i].value / molar_mass[j].value;
            let phi = (1.0
                + (pure_viscosity[i] / pure_viscosity[j]).sqrt()
                    * (molar_mass[j].value / molar_mass[i].value).powf(0.25))
            .powi(2)
                / (8.0 * (1.0 + ratio)).sqrt();
            denominator += x[j] * phi;
        }
        dilute_viscosity += x[i] * pure_viscosity[i] / denominator;
    }

    let cv0_mix: f64 = (0..count).map(|i| x[i] * cv0[i].value).sum();
    let eucken_factor = cv0_mix + 1.25 * R;
    let total_mass: f64 = (0..count).map(|i| x[i] * molar_mass[i].value).sum();

    let excess = viscosity.value - dilute_viscosity;
    let dense_conductivity = if excess > 0.0 && total_mass > 1.0e-30 {
        excess * eucken_factor / total_mass
    } else {
        0.0
    };

    let mut conductivity = dilute_conductivity + dense_conductivity;
    if conductivity < 1.0e-10 {
        conductivity = 1.0e-10;
    }

    apply_checks(
        spec.derived_checks(),
        |name| (name == "conductivity").then_some(conductivity),
        &mut warnings,
    )?;

    Ok(FrictionTheoryConductivityResult {
        conductivity: watts_per_meter_kelvin(conductivity),
        pure_component_conductivity: pure_conductivity,
        dilute_conductivity: watts_per_meter_kelvin(dilute_conductivity),
        dense_conductivity: watts_per_meter_kelvin(dense_conductivity),
        warnings,
    })
}
