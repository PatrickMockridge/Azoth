//! `standards.iso6976` - the calorific values and density of a natural gas.
//!
//! ```text
//! M    = sum z_i M_i
//! Z    = 1 - (sum z_i sqrt(b_i))**2
//! d    = (sum z_i M_i) / M_air * Z_air / Z
//! Hsup = sum z_i Hsup_i     Hinf = sum z_i Hinf_i
//! rho_ideal = P M / (R T)   rho_real = rho_ideal / Z
//! ```
//!
//! Spec: `specs/models/standards/iso6976.toml`, which carries the standard's citation and
//! the two constants below.
//!
//! # What the class does, and where this follows it
//!
//! `Standard_ISO6976.calculate` fills five z-weighted sums - the molar mass, the compression
//! factor's summation term at each of three reference temperatures, the two relative-density
//! terms and the ten calorific values - and `getValue` then selects among them by the two
//! reference temperatures and shapes the answer by the reference type. **This returns every
//! one of those quantities by name, each in its own unit, instead of selecting one by a
//! flag**: a caller that wants the number NeqSim's `getValue("LCV")` would return composes
//! it from the molar calorific value and the density, and the composition is then visible
//! rather than hidden behind a string argument.
//!
//! # The two constants, and why they are not today's
//!
//! `R` is `8.314510` J/(mol*K) and the reference pressure `1.01325` bar, both as the class
//! declares them. The current CODATA molar gas constant is `8.314462618`, and the
//! difference is in the eighth digit - small, but the port reproduces the class rather than
//! improving it, and a reader comparing this against NeqSim's output should not have to
//! wonder which constant moved.

use azoth_core::units::{
    ThermodynamicTemperature, joules_per_mole, kilograms_per_cubic_meter, kilograms_per_mole,
};
use azoth_core::{AzothError, Result, apply_checks};
use uom::si::f64::Ratio;

use crate::iso6976_constants::{Iso6976Row, row};
use crate::model_gen;
use crate::results::Iso6976Result;

/// The class's own molar gas constant, J/(mol*K). Not today's CODATA value - see the module.
const R: f64 = 8.314510;

/// The standard's reference pressure, Pa: 1.01325 bar.
const REFERENCE_PRESSURE: f64 = 101_325.0;

/// Dry air's molar mass, kg/mol, from `ThermodynamicConstantsInterface.molarMassAir`.
const MOLAR_MASS_AIR: f64 = 0.028_965_46;

/// 273.15 K, so the standard's degrees Celsius are reachable from a kelvin temperature.
const KELVIN_AT_ZERO_CELSIUS: f64 = 273.15;

/// Air's own compression factor at the volumetric reference temperature, from the class's
/// field initialisers: `Zair0`, `Zair15` and `Zair20`.
fn air_compression_factor(celsius_value: f64) -> Result<f64> {
    // The same tolerance the table's own selectors use, and for the same reason: the caller
    // states kelvin, so 288.7 K is 15.550000000000011 C and an exact match would refuse it.
    const TOLERANCE: f64 = 1.0e-9;
    match celsius_value {
        other if (other - 0.0).abs() < TOLERANCE => Ok(0.999_41),
        other if (other - 15.0).abs() < TOLERANCE || (other - 15.55).abs() < TOLERANCE => {
            Ok(0.999_58)
        }
        other if (other - 20.0).abs() < TOLERANCE => Ok(0.999_63),
        other => Err(AzothError::invalid_input(
            "volumetric_reference_temperature",
            format!(
                "the standard's air compression factor is tabulated at 0, 15 and 20 C, and \
                 {other} C is not one of them (15.55 C is the 60 F reference and reads the \
                 15 C value)"
            ),
        )),
    }
}

/// The calorific values and density of a natural gas, from its composition.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the composition and the component list disagree in
///   length, if a component has no row in the standard's table, or if a reference
///   temperature is outside the standard's set - which the class silently corrects and
///   this refuses.
///
/// # Example
/// ```
/// use azoth_core::units::kelvins;
/// use azoth_standards::iso6976;
///
/// let r = iso6976(&["methane".to_string()], &[1.0], kelvins(288.15), kelvins(298.15))?;
/// // The table's methane row, summed over one component.
/// assert!((r.molar_mass.value - 0.016043).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn iso6976(
    components: &[String],
    z: &[f64],
    volumetric_reference_temperature: ThermodynamicTemperature,
    energy_reference_temperature: ThermodynamicTemperature,
) -> Result<Iso6976Result> {
    let spec = &model_gen::ISO6976_SPEC;
    let mut warnings = Vec::new();

    // **The two reference temperatures cross into Celsius here and are never Celsius
    // again**: the standard states its rows in degrees Celsius, this library states a
    // temperature in kelvin, and the offset is exact. The table's own selectors take the
    // Celsius value, which is what its columns are named for.
    let vol_c = volumetric_reference_temperature.value - KELVIN_AT_ZERO_CELSIUS;
    let energy_c = energy_reference_temperature.value - KELVIN_AT_ZERO_CELSIUS;
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            // The declared bounds are in kelvin, which is what a caller states; the Celsius
            // values below are the table's own reference frame.
            "volumetric_reference_temperature" => Some(volumetric_reference_temperature.value),
            "energy_reference_temperature" => Some(energy_reference_temperature.value),
            _ => None,
        },
        &mut warnings,
    )?;

    if components.len() != z.len() {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "{} components and {} mole fractions, which is not a composition",
                components.len(),
                z.len()
            ),
        ));
    }
    if components.is_empty() {
        return Err(AzothError::invalid_input(
            "components",
            "a gas of no components has no molar mass and no calorific value",
        ));
    }
    if let Some(bad) = z.iter().position(|value| *value < 0.0) {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "z[{bad}] is {} but a mole fraction cannot be negative",
                z[bad]
            ),
        ));
    }

    // The rows, in the composition's order, and the sums the class fills.
    let rows: Vec<Iso6976Row> = components
        .iter()
        .map(|name| row(name))
        .collect::<Result<_>>()?;

    let mut molar_mass = 0.0;
    let mut summation = 0.0;
    let mut relative_ideal = 0.0;
    let mut superior = 0.0;
    let mut inferior = 0.0;
    for (row, fraction) in rows.iter().zip(z) {
        molar_mass += fraction * row.molar_mass;
        summation += fraction * row.summation_factor_at(vol_c)?;
        relative_ideal += fraction * row.molar_mass / MOLAR_MASS_AIR;
        superior += fraction * row.calorific_value_at(energy_c, true)?;
        inferior += fraction * row.calorific_value_at(energy_c, false)?;
    }

    // The class's own arithmetic, in its own order: `Zmix = 1 - (sum z sqrt(b))^2`, then the
    // real-state relative density, then the two densities from the ideal-gas law.
    let compression_factor = 1.0 - summation * summation;
    let z_air = air_compression_factor(vol_c)?;
    let relative_density = relative_ideal * z_air / compression_factor;
    let reference_temperature = vol_c + KELVIN_AT_ZERO_CELSIUS;
    let density_ideal = REFERENCE_PRESSURE * molar_mass / (R * reference_temperature);
    let density_real = density_ideal / compression_factor;

    apply_checks(
        spec.derived_checks(),
        |quantity| match quantity {
            "compression_factor" => Some(compression_factor),
            _ => None,
        },
        &mut warnings,
    )?;

    Ok(Iso6976Result {
        molar_mass: kilograms_per_mole(molar_mass),
        compression_factor: Ratio::new::<uom::si::ratio::ratio>(compression_factor),
        relative_density: Ratio::new::<uom::si::ratio::ratio>(relative_density),
        density_ideal: kilograms_per_cubic_meter(density_ideal),
        density_real: kilograms_per_cubic_meter(density_real),
        superior_calorific_value: joules_per_mole(superior),
        inferior_calorific_value: joules_per_mole(inferior),
        warnings,
    })
}
