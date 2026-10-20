//! `eos.filippov_conductivity` - a liquid's thermal conductivity from the Filippov
//! mixing rule, with the optional Missenard pressure correction.
//!
//! Spec: `specs/models/eos/filippov_conductivity.toml`, which records
//! NeqSim's `FilippovConductivityMethod.calcConductivity` and its two helpers.

use azoth_core::units::{MolarMass, Pressure, ThermodynamicTemperature, watts_per_meter_kelvin};
use azoth_core::{AzothError, Result, apply_checks};

use crate::model_gen;
use crate::results::FilippovConductivityResult;

/// The floor NeqSim puts under a pure component's polynomial, in W/(m K).
const MIN_PURE_CONDUCTIVITY: f64 = 1.0e-10;

/// A liquid mixture's thermal conductivity, from the components' polynomials and the
/// Filippov mixing rule.
///
/// The pure conductivities are the LIQCOND polynomials at `T`, each floored at
/// [`MIN_PURE_CONDUCTIVITY`]. They are mixed by **mass** fraction:
/// `k = sum_i w_i k_i - c sum_{i<j} w_i w_j |k_i - k_j|`. The Missenard factor
/// `0.98 + 0.0079 Pr sqrt(Tr)` is applied where `use_pressure_correction` is set and the
/// state is inside its window (`Tr < 0.95` and `Pr >= 0.01`); it is never below one.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if the vectors differ in length or the
///   matrix is not one row of three per component.
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's range.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_mole, kelvins, pascals};
/// use azoth_eos::filippov_conductivity::filippov_conductivity;
///
/// // The CO2 absorber's aqueous phase at 313.15 K and 50 bar.
/// let r = filippov_conductivity(
///     &[[0.251502, 5.238919e-4, -3.82111e-6], [-0.384, 0.00525, -6.37e-6]],
///     &[kilograms_per_mole(0.04401), kilograms_per_mole(0.018015)],
///     &[6.691762234084198e-4, 0.9993308237765915],
///     &[kelvins(304.19), kelvins(647.3)],
///     &[pascals(7.3815e6), pascals(22.089e6)],
///     kelvins(313.15),
///     pascals(5.0e6),
///     0.72,
///     true,
/// )?;
/// assert!((r.k.value - 0.6337077403924924).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T` and `P` are the symbols in the published equation
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn filippov_conductivity(
    liquid_conductivity: &[[f64; 3]],
    molar_mass: &[MolarMass],
    z: &[f64],
    tc: &[ThermodynamicTemperature],
    pc: &[Pressure],
    T: ThermodynamicTemperature,
    P: Pressure,
    filippov_coefficient: f64,
    use_pressure_correction: bool,
) -> Result<FilippovConductivityResult> {
    let spec = &model_gen::FILIPPOV_CONDUCTIVITY_SPEC;
    let mut warnings = Vec::new();

    let count = z.len();
    if liquid_conductivity.len() != count
        || molar_mass.len() != count
        || tc.len() != count
        || pc.len() != count
    {
        return Err(AzothError::invalid_input(
            "liquid_conductivity",
            format!(
                "{} coefficient row(s), {} molar mass(es), {} critical temperature(s), {} \
                 critical pressure(s) and {} mole fraction(s): the mean is over the components, \
                 so all five are one entry each",
                liquid_conductivity.len(),
                molar_mass.len(),
                tc.len(),
                pc.len(),
                count
            ),
        ));
    }

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T.value),
            "P" => Some(P.value),
            "z" => z.iter().copied().reduce(f64::min),
            "molar_mass" => molar_mass.iter().map(|m| m.value).reduce(f64::min),
            "tc" => tc.iter().map(|t| t.value).reduce(f64::min),
            "pc" => pc.iter().map(|p| p.value).reduce(f64::min),
            _ => None,
        },
        &mut warnings,
    )?;

    // The mass fractions, which are the mean's weights.
    let masses: Vec<f64> = z
        .iter()
        .zip(molar_mass)
        .map(|(fraction, mass)| fraction * mass.value)
        .collect();
    let total_mass: f64 = masses.iter().sum();

    let pure: Vec<f64> = liquid_conductivity
        .iter()
        .map(|row| {
            (row[0] + row[1] * T.value + row[2] * T.value.powi(2)).max(MIN_PURE_CONDUCTIVITY)
        })
        .collect();

    let mut linear = 0.0;
    for (mass, pure_k) in masses.iter().zip(&pure) {
        let weight = if total_mass > 0.0 {
            mass / total_mass
        } else {
            0.0
        };
        linear += weight * pure_k;
    }
    let mut correction = 0.0;
    for i in 0..count {
        let wi = if total_mass > 0.0 {
            masses[i] / total_mass
        } else {
            0.0
        };
        for j in (i + 1)..count {
            let wj = if total_mass > 0.0 {
                masses[j] / total_mass
            } else {
                0.0
            };
            correction += wi * wj * (pure[i] - pure[j]).abs();
        }
    }
    let mut k = (linear - filippov_coefficient * correction).max(MIN_PURE_CONDUCTIVITY);

    let tc_mix: f64 = z.iter().zip(tc).map(|(x, t)| x * t.value).sum();
    let pc_mix: f64 = z.iter().zip(pc).map(|(x, p)| x * p.value).sum();
    let applied = if use_pressure_correction {
        missenard_pressure_correction(T.value, P.value, tc_mix, pc_mix)
    } else {
        1.0
    };
    k *= applied;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "k").then_some(k),
        &mut warnings,
    )?;

    Ok(FilippovConductivityResult {
        k: watts_per_meter_kelvin(k),
        pressure_correction: applied,
        pure_component_conductivity: pure,
        warnings,
    })
}

/// `getMissenardPressureCorrection`: `0.98 + 0.0079 Pr sqrt(Tr)`, clamped up to one and
/// switched off outside `Tr < 0.95` and `Pr >= 0.01`.
fn missenard_pressure_correction(temperature: f64, pressure: f64, tc_mix: f64, pc_mix: f64) -> f64 {
    if pc_mix <= 0.0 || tc_mix <= 0.0 {
        return 1.0;
    }
    let reduced_temperature = temperature / tc_mix;
    let reduced_pressure = pressure / pc_mix;
    if reduced_temperature >= 0.95 || reduced_pressure < 0.01 {
        return 1.0;
    }
    let factor = 0.98 + 0.0079 * reduced_pressure * reduced_temperature.sqrt();
    factor.max(1.0)
}
