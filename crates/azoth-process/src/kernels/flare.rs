//! `unit_ops.flare` - the pass-through whose steady state is its report.

use azoth_core::units::{MassRate, Power, kelvins, kilograms_per_second, watts};
use azoth_core::{AzothError, Result};
use azoth_reactions::databank::element_composition;
use azoth_standards::iso6976;

use crate::stream::Stream;

/// The numbers `Flare.run` reports beside the record, and the surface it exposes around them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlareNumbers {
    /// The heat the flare releases, W - `inStream.LCV() * flowSm3sec`.
    pub heat_duty: Power,
    /// The carbon dioxide the combustion forms, kg/s.
    pub co2_emission: MassRate,
    /// `getFlowRate("Sm3/sec")`, m3/s - the surrogate's standard volume at the operating point.
    pub standard_volumetric_flow: f64,
    /// `getThermoSystem().getDensity("kg/m3")`, kg/m3.
    pub density: f64,
    /// `getFlowRate("kg/sec")`, kg/s.
    pub mass_rate: f64,
    /// `getFlowRate("mole/sec")`, mol/s.
    pub molar_rate: f64,
    /// `estimateRadiationHeatFlux`, W/m2; absent without a ground distance.
    pub radiation_flux: Option<f64>,
    /// `radiationDistanceForFlux`, m; absent without a threshold.
    pub radiation_distance: Option<f64>,
    /// `getDispersionSurrogate().getExitVelocityMs()`, m/s.
    pub exit_velocity: f64,
    /// The surrogate's `momentumFlux`, `rho v^2`, Pa.
    pub momentum_flux: f64,
    /// The duty over the design heat release; absent without a positive design.
    pub heat_utilization: Option<f64>,
    /// The operating mass flow over the design mass flow; absent as above.
    pub mass_utilization: Option<f64>,
    /// The operating molar flow over the design molar flow; absent as above.
    pub molar_utilization: Option<f64>,
    /// `isOverloaded`; absent until a utilisation exists.
    pub overloaded: Option<bool>,
}

/// The setters' state the radiation, dispersion and capacity surfaces read.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlareSettings {
    /// `setFlameHeight`: the point source's height, m.
    pub flame_height: f64,
    /// `setRadiantFraction`: the share of the duty radiated.
    pub radiant_fraction: f64,
    /// `setTipDiameter`: the tip's outlet diameter, m.
    pub tip_diameter: f64,
    /// The horizontal distance the flux is read at, m.
    pub ground_distance: Option<f64>,
    /// The flux the distance is solved for, W/m2.
    pub radiation_flux_threshold: Option<f64>,
    /// `setDesignHeatDutyCapacity`, W.
    pub design_heat_duty_capacity: Option<f64>,
    /// `setDesignMassFlowCapacity`, kg/s.
    pub design_mass_flow_capacity: Option<f64>,
    /// `setDesignMolarFlowCapacity`, mol/s.
    pub design_molar_flow_capacity: Option<f64>,
}

impl Default for FlareSettings {
    fn default() -> Self {
        Self {
            flame_height: DEFAULT_FLAME_HEIGHT,
            radiant_fraction: DEFAULT_RADIANT_FRACTION,
            tip_diameter: DEFAULT_TIP_DIAMETER,
            ground_distance: None,
            radiation_flux_threshold: None,
            design_heat_duty_capacity: None,
            design_mass_flow_capacity: None,
            design_molar_flow_capacity: None,
        }
    }
}

/// `Flare.flameHeight`'s initialiser, m.
pub const DEFAULT_FLAME_HEIGHT: f64 = 30.0;
/// `Flare.radiantFraction`'s initialiser.
pub const DEFAULT_RADIANT_FRACTION: f64 = 0.18;
/// `Flare.tipDiameter`'s initialiser, m.
pub const DEFAULT_TIP_DIAMETER: f64 = 0.3;
/// `setTipDiameter`'s floor, m.
pub const MIN_TIP_DIAMETER: f64 = 1.0e-4;
/// The surrogate's density floor, kg/m3.
pub const MIN_DENSITY: f64 = 1.0e-3;
/// The surrogate's tip-area floor, m2.
pub const MIN_TIP_AREA: f64 = 1.0e-6;
/// `CapacityCheckResult.exceeds`: a utilisation above it is an overload.
pub const OVERLOAD_TOLERANCE: f64 = 1.0e-6;

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
/// `LCV()` is `Standard_ISO6976(fluid, 0, 15.55, "volume").getValue("InferiorCalorificValue")`,
/// scaled by `1.0e3`: the inferior calorific value in **joules per normal cubic metre at
/// 0 °C**, from the standard's 60 °F column and the mixture's compression factor at 0 °C.
/// `Sm3/sec` is `n * R * 288.15 / atm`, a volume at **15 °C**, ideal. **So the reported duty
/// multiplies an energy density at one reference by a volumetric flow at another**, and is
/// five and a half per cent above the same gas measured consistently. That is what the class
/// computes and what this reproduces; the case says so.
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
    flare_with_settings(inlet, &FlareSettings::default())
}

/// The same flare, with the setters' state the radiation, dispersion and capacity surfaces read.
///
/// # Errors
/// The same refusals as [`flare`]. A two-phase stream is not one: its dispersion surrogate reads
/// the system's volume-weighted mean density, as `getThermoSystem().getDensity` does.
pub fn flare_with_settings(
    inlet: &Stream,
    settings: &FlareSettings,
) -> Result<(Stream, FlareNumbers)> {
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

    // **The dispersion surrogate's density is the stream's system density**,
    // `getThermoSystem().getDensity`'s reading: the phases' volume-weighted mean, which answers on
    // a two-phase feed where a single phase's cubic volume would not.
    let density = crate::segment::phase::system_mass_density(inlet)?;
    let mass_rate = inlet.mass_flow()?;
    let radiation_flux = settings.ground_distance.map(|distance| {
        radiation_heat_flux(
            heat_duty,
            distance,
            settings.flame_height,
            settings.radiant_fraction,
        )
    });
    let radiation_distance = settings.radiation_flux_threshold.map(|threshold| {
        radiation_distance_for_flux(
            heat_duty,
            threshold,
            settings.flame_height,
            settings.radiant_fraction,
        )
    });
    let (exit_velocity, momentum_flux) = dispersion(mass_rate, density, settings.tip_diameter);
    let (heat_utilization, mass_utilization, molar_utilization, overloaded) = capacity(
        (heat_duty, mass_rate, inlet.n),
        (
            settings.design_heat_duty_capacity,
            settings.design_mass_flow_capacity,
            settings.design_molar_flow_capacity,
        ),
    );

    Ok((
        inlet.clone(),
        FlareNumbers {
            heat_duty: watts(heat_duty),
            co2_emission: kilograms_per_second(carbon * CO2_MOLAR_MASS),
            standard_volumetric_flow,
            density,
            mass_rate,
            molar_rate: inlet.n,
            radiation_flux,
            radiation_distance,
            exit_velocity,
            momentum_flux,
            heat_utilization,
            mass_utilization,
            molar_utilization,
            overloaded,
        },
    ))
}

/// `estimateRadiationHeatFlux`: the point source at the flame height over a sphere.
fn radiation_heat_flux(
    heat_duty: f64,
    ground_distance: f64,
    flame_height: f64,
    radiant_fraction: f64,
) -> f64 {
    let radial = (ground_distance * ground_distance + flame_height * flame_height)
        .max(0.0)
        .sqrt();
    if radial < 1.0e-6 {
        return 0.0;
    }
    radiant_fraction * heat_duty / (4.0 * std::f64::consts::PI * radial * radial)
}

/// `radiationDistanceForFlux`: the inverse, zero where the contour is inside the sphere.
fn radiation_distance_for_flux(
    heat_duty: f64,
    threshold: f64,
    flame_height: f64,
    radiant_fraction: f64,
) -> f64 {
    if threshold <= 0.0 {
        return 0.0;
    }
    let numerator = radiant_fraction * heat_duty;
    if numerator <= 0.0 {
        return 0.0;
    }
    let radial_squared = numerator / (4.0 * std::f64::consts::PI * threshold);
    let horizontal_squared = radial_squared - flame_height * flame_height;
    if horizontal_squared <= 0.0 {
        return 0.0;
    }
    horizontal_squared.sqrt()
}

/// `getDispersionSurrogate`: the tip's exit velocity and momentum flux.
fn dispersion(mass_rate: f64, density: f64, tip_diameter: f64) -> (f64, f64) {
    let density = density.max(MIN_DENSITY);
    let area = (std::f64::consts::PI * tip_diameter * tip_diameter / 4.0).max(MIN_TIP_AREA);
    let velocity = (mass_rate / (density * area)).max(0.0);
    (velocity, density * velocity * velocity)
}

/// `evaluateCapacity`: the three utilisations and the overload flag.
fn capacity(
    actual: (f64, f64, f64),
    design: (Option<f64>, Option<f64>, Option<f64>),
) -> (Option<f64>, Option<f64>, Option<f64>, Option<bool>) {
    fn utilization(actual: f64, design: Option<f64>) -> Option<f64> {
        match design {
            Some(value) if value.is_finite() && value > 0.0 => Some(actual / value),
            _ => None,
        }
    }

    let heat = utilization(actual.0, design.0);
    let mass = utilization(actual.1, design.1);
    let molar = utilization(actual.2, design.2);
    let present: Vec<f64> = [heat, mass, molar].into_iter().flatten().collect();
    let overloaded = if present.is_empty() {
        None
    } else {
        Some(
            present
                .iter()
                .any(|value| *value > 1.0 + OVERLOAD_TOLERANCE),
        )
    };
    (heat, mass, molar, overloaded)
}
