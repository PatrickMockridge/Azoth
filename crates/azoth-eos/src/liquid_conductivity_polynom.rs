//! `eos.liquid_conductivity_polynom` - a liquid's thermal conductivity from the LIQCOND
//! polynomial.
//!
//! Spec: `specs/models/eos/liquid_conductivity_polynom.toml`, which records the mass-fraction
//! mean, the polynomial and the floor a negative one takes.

use azoth_core::units::{MolarMass, ThermodynamicTemperature, watts_per_meter_kelvin};
use azoth_core::{AzothError, Result, apply_checks};

use crate::model_gen;
use crate::results::LiquidConductivityPolynomResult;

/// The floor NeqSim puts under a pure component's polynomial, in W/(m K).
const MIN_PURE_CONDUCTIVITY: f64 = 1.0e-10;

/// A liquid mixture's thermal conductivity, from the components' LIQCOND polynomials.
///
/// **A mass-fraction mean**, which is `calcConductivity`: each component's polynomial is
/// evaluated at `T`, floored at [`MIN_PURE_CONDUCTIVITY`], and weighted by its **mass** fraction -
/// so `molar_mass` and `z` are both read even though the polynomial itself is a function of
/// temperature alone.
///
/// **This is the aqueous phase's conductivity.** A gas or a hydrocarbon liquid takes PFCT
/// (`eos.thermal_conductivity`); the two answer different numbers at the same state, and which
/// one a phase gets is the dispatch this id exists beside.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if the vectors differ in length or the matrix is
///   not one row of three per component.
/// * [`azoth_core::AzothError::OutOfRange`] if `T` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::kelvins;
/// use azoth_eos::liquid_conductivity_polynom::liquid_conductivity_polynom;
///
/// // The CO2 absorber's aqueous phase at 313.15 K.
/// let r = liquid_conductivity_polynom::liquid_conductivity_polynom(
///     &[[0.251502, 0.0005238919, -3.82111e-6], [-0.384, 0.00525, -6.37e-6]],
///     &[0.04401, 0.018015],
///     &[0.0006691762234084198, 0.9993308237765915],
///     kelvins(313.15),
/// )?;
/// assert!((r.k.value - 0.6344057039895419).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T` is the symbol in the polynomial
pub fn liquid_conductivity_polynom(
    liquid_conductivity: &[[f64; 3]],
    molar_mass: &[MolarMass],
    z: &[f64],
    T: ThermodynamicTemperature,
) -> Result<LiquidConductivityPolynomResult> {
    let spec = &model_gen::LIQUID_CONDUCTIVITY_POLYNOM_SPEC;
    let mut warnings = Vec::new();

    let count = z.len();
    if liquid_conductivity.len() != count || molar_mass.len() != count {
        return Err(AzothError::invalid_input(
            "liquid_conductivity",
            format!(
                "{} row(s) of coefficients, {} molar mass(es) and {} mole fraction(s): the mean \
                 is over the components, so the three are one entry each",
                liquid_conductivity.len(),
                molar_mass.len(),
                count
            ),
        ));
    }

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T.value),
            // The class reads a mole fraction without checking it; a negative one would take
            // mass away from the mean.
            "z" => z.iter().copied().reduce(f64::min),
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

    let mut k = 0.0;
    for (row, mass) in liquid_conductivity.iter().zip(&masses) {
        let weight = if total_mass > 0.0 {
            mass / total_mass
        } else {
            0.0
        };
        let pure = row[0] + row[1] * T.value + row[2] * T.value.powi(2);
        k += weight * pure.max(MIN_PURE_CONDUCTIVITY);
    }

    apply_checks(
        spec.derived_checks(),
        |name| (name == "k").then_some(k),
        &mut warnings,
    )?;

    Ok(LiquidConductivityPolynomResult {
        k: watts_per_meter_kelvin(k),
        warnings,
    })
}
