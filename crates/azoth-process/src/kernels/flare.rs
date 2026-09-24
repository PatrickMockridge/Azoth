//! `unit_ops.flare` - the pass-through whose steady state is its report.

use azoth_core::units::{MassRate, Power, kelvins, kilograms_per_second, watts};
use azoth_core::{AzothError, Result};
use azoth_reactions::databank::element_composition;
use azoth_standards::iso6976;

use crate::stream::Stream;

/// The two numbers `Flare.run` reports beside the record.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlareNumbers {
    /// The heat the flare releases, W - `inStream.LCV() * flowSm3sec`.
    pub heat_duty: Power,
    /// The carbon dioxide the combustion forms, kg/s.
    pub co2_emission: MassRate,
}

/// NeqSim's `ThermodynamicConstantsInterface` values the two conversions use.
///
/// **Two gas constants, and the class really does use both.** `Stream.LCV()` reaches the
/// standard's `molRefm3` through `Standard_ISO6976`'s own `R = 8.314510`, and
/// `getFlowRate("Sm3/sec")` divides by `ThermodynamicConstantsInterface.R = 8.3144621`. They
/// differ in the sixth digit, and a port that used one for both would be a port of neither.
const R_STANDARD: f64 = 8.314510;
const R_THERMO: f64 = 8.3144621;
/// The standard's reference pressure, Pa.
const REFERENCE_PRESSURE: f64 = 101_325.0;
/// The standard's volumetric reference temperature, K - 0 °C, the `Nm3` basis.
const NORMAL_TEMPERATURE: f64 = 273.15;
/// NeqSim's `standardStateTemperature`, K - 15 °C, the `Sm3` basis.
const STANDARD_STATE_TEMPERATURE: f64 = 288.15;
/// The 60 °F reference, K, which `Stream.LCV()` states the energy at.
const SIXTY_FAHRENHEIT: f64 = 288.7;
/// Molar mass of CO2, kg/mol, as `Flare.run` hard-codes it.
const CO2_MOLAR_MASS: f64 = 44.01e-3;

/// **A flare is a pass-through with a report.**
///
/// `Flare.run` clones the inlet into the outlet and computes two numbers:
///
/// ```text
/// heatDuty    = inStream.LCV() * inStream.getFlowRate("Sm3/sec")
/// co2Emission = (sum_i z_i * n * nC_i) * 44.01e-3
/// ```
///
/// `LCV()` is `new Standard_ISO6976(fluid, 0, 15.55, "volume").getValue("InferiorCalorificValue")
/// * 1.0e3`: the inferior calorific value in **joules per normal cubic metre at 0 °C**, from
/// the standard's 60 °F column and the mixture's compression factor at 0 °C. `Sm3/sec` is
/// `n * R * 288.15 / atm` - a volume at **15 °C**, ideal. **So the reported duty multiplies
/// an energy density at one reference by a volumetric flow at another**, and is five and a
/// half per cent above the same gas measured consistently. That is what the class computes
/// and what this reproduces; the case says so.
///
/// The carbon comes from the element table rather than the standard's, which is what
/// `getElements().getNumberOfElements("C")` reads - so a gas the standard has no row for
/// still gets its CO2 number, and only its duty is refused.
///
/// # Errors
/// [`AzothError::InvalidInput`] if a component has no row in ISO 6976's table, or if the
/// element table does not carry one; the outlet is the inlet either way, which is why this
/// returns the pair rather than failing silently with a zero duty.
pub fn flare(inlet: &Stream) -> Result<(Stream, FlareNumbers)> {
    if inlet.n < 0.0 {
        return Err(AzothError::invalid_input(
            "inlet_n",
            format!(
                "a molar flow cannot be negative, and this one is {}",
                inlet.n
            ),
        ));
    }

    // The standard's own two numbers at the reference pair `Stream.LCV()` states: the
    // inferior calorific value at 60 °F, and the compression factor at 0 °C.
    let quality = iso6976(
        &inlet.components,
        &inlet.z,
        kelvins(NORMAL_TEMPERATURE),
        kelvins(SIXTY_FAHRENHEIT),
    )?;

    // J/Nm3: the molar value times the normal-state mole density.
    let moles_per_normal_cubic_metre =
        REFERENCE_PRESSURE / (R_STANDARD * NORMAL_TEMPERATURE * quality.compression_factor.value);
    let volumetric_calorific_value =
        quality.inferior_calorific_value.value * moles_per_normal_cubic_metre;

    // m3/s at 15 °C, ideal - `getFlowRate("Sm3/sec")`.
    let standard_volumetric_flow =
        inlet.n * R_THERMO * STANDARD_STATE_TEMPERATURE / REFERENCE_PRESSURE;

    let heat_duty = volumetric_calorific_value * standard_volumetric_flow;

    let mut carbon = 0.0;
    for (name, fraction) in inlet.components.iter().zip(&inlet.z) {
        let composition = element_composition(name)?.ok_or_else(|| {
            AzothError::invalid_input(
                "components",
                format!(
                    "{name} is not in the element table, so the carbon it carries is unknown \
                     and the CO2 emission cannot be computed"
                ),
            )
        })?;
        let atoms: f64 = composition
            .iter()
            .filter(|(element, _)| element == "C")
            .map(|(_, count)| count)
            .sum();
        carbon += fraction * inlet.n * atoms;
    }

    Ok((
        inlet.clone(),
        FlareNumbers {
            heat_duty: watts(heat_duty),
            co2_emission: kilograms_per_second(carbon * CO2_MOLAR_MASS),
        },
    ))
}
