//! `process.gas_dispersion` - the ambient-flash screening as a registered id.
//!
//! Spec: `specs/models/process/gas_dispersion.toml`. The arithmetic is
//! `process/safety/dispersion/GasDispersionAnalyzer.java` at the pin: the ambient flash, the
//! fuel data, the density-ratio model selection, and the flammable and toxic endpoints over
//! [`crate::models::gaussian_plume`]'s plume or the class's Britter-McQuaid fit.

use azoth_core::units::{
    Length, MassDensity, MassRate, Pressure, ThermodynamicTemperature, Velocity, Volume,
    cubic_meters, cubic_meters_per_mole, kelvins, kilograms_per_cubic_meter, kilograms_per_mole,
    kilograms_per_second, meters, meters_per_second, pascals,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use azoth_eos::cubic::Cubic;
use azoth_eos::databank;
use azoth_eos::entropy_state;
use azoth_eos::phase_boundary::check_composition;
use azoth_eos::results::Phase;
use azoth_eos::{pr_mass_density, pr_molar_volume, pt_flash};

use crate::executor::json::{optional_scalar, scalar, warnings as wire_warnings};
use crate::model_gen;
use crate::models::gaussian_plume::{
    Stability, Terrain, concentration, distance_to_threshold, sigma_y, sigma_z,
};

/// `GasDispersionAnalyzer.GAS_CONSTANT`, J/(mol K).
const GAS_CONSTANT: f64 = 8.314_462_618;
/// `ToxicLibrary`'s gas constant, J/(mol K), which the class does not share with the above.
const TOXIC_GAS_CONSTANT: f64 = 8.314;
/// `GasDispersionAnalyzer.AIR_GAS_CONSTANT`, J/(kg K).
const AIR_GAS_CONSTANT: f64 = 287.05;
/// `GasDispersionAnalyzer.DENSE_GAS_RATIO`.
const DENSE_GAS_RATIO: f64 = 1.2;
/// `GasDispersionAnalyzer.DEFAULT_RELEASE_HEIGHT_M`.
const DEFAULT_RELEASE_HEIGHT: f64 = 0.0;
/// The BoundaryConditions defaults the calculation reads.
const DEFAULT_AMBIENT_TEMPERATURE: f64 = 288.15;
const DEFAULT_ATMOSPHERIC_PRESSURE: f64 = 101_325.0;
const DEFAULT_WIND_SPEED: f64 = 5.0;
const DEFAULT_SURFACE_ROUGHNESS: f64 = 0.1;
const DEFAULT_STABILITY: &str = "D";
/// The class's `flammableEndpointFractionOfLfl` default.
const DEFAULT_ENDPOINT_FRACTION: f64 = 0.5;
/// The Gaussian cloud-volume integral's slab count.
const CLOUD_STEPS: usize = 200;

/// The class's lower-flammable-limit table, by normalized name, volume fraction.
const LOWER_FLAMMABLE_LIMITS: [(&str, f64); 13] = [
    ("methane", 0.044),
    ("ethane", 0.030),
    ("propane", 0.021),
    ("nbutane", 0.018),
    ("ibutane", 0.018),
    ("npentane", 0.014),
    ("ipentane", 0.014),
    ("nhexane", 0.012),
    ("hydrogen", 0.040),
    ("h2s", 0.040),
    ("co", 0.125),
    ("methanol", 0.060),
    ("benzene", 0.012),
];

/// The model the screen runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Selection {
    Auto,
    Gaussian,
    Heavy,
}

impl Selection {
    fn parse(value: Option<&str>) -> Result<Self> {
        Ok(match value.unwrap_or("auto") {
            "auto" => Self::Auto,
            "gaussian_plume" => Self::Gaussian,
            "heavy_gas" => Self::Heavy,
            other => {
                return Err(AzothError::invalid_input(
                    "model_selection",
                    format!(
                        "the model is `auto`, `gaussian_plume` or `heavy_gas`, and {other:?} is not"
                    ),
                ));
            }
        })
    }
}

/// Result of `process.gas_dispersion`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct GasDispersionResult {
    /// `gaussian_plume` or `heavy_gas`.
    pub selected_model: String,
    /// The rate the screen used.
    #[serde(serialize_with = "scalar")]
    pub mass_release_rate: MassRate,
    /// The flammable share of the rate.
    #[serde(serialize_with = "scalar")]
    pub flammable_mass_release_rate: MassRate,
    /// The gas phase's density at the ambient flash.
    #[serde(serialize_with = "scalar")]
    pub source_density: MassDensity,
    /// `P / (287.05 T)`.
    #[serde(serialize_with = "scalar")]
    pub air_density: MassDensity,
    /// The flammable components' mole fraction.
    pub fuel_mole_fraction: f64,
    /// The flammable components' mass fraction.
    pub fuel_mass_fraction: f64,
    /// The mixture's Le Chatelier LFL, or the stated override; absent when the phase has
    /// no flammable component.
    pub lower_flammable_limit: Option<f64>,
    /// The endpoint fraction the first distance was asked at.
    pub flammable_endpoint_fraction: f64,
    /// The distance at the endpoint fraction; absent when the arm answers none.
    #[serde(serialize_with = "optional_scalar")]
    pub distance_to_flammable_endpoint: Option<Length>,
    /// The distance at half the LFL; absent when the arm answers none.
    #[serde(serialize_with = "optional_scalar")]
    pub distance_to_half_lfl: Option<Length>,
    /// The distance at the LFL; absent when there is no fuel.
    #[serde(serialize_with = "optional_scalar")]
    pub distance_to_lfl: Option<Length>,
    /// The cloud volume above the LFL.
    #[serde(serialize_with = "scalar")]
    pub flammable_cloud_volume: Volume,
    /// The screened component's name, empty when none was stated.
    pub toxic_component: String,
    /// The toxic threshold by volume fraction; absent when none was asked.
    pub toxic_threshold: Option<f64>,
    /// The distance at the toxic threshold; absent when none was asked or the component is
    /// not in the phase.
    #[serde(serialize_with = "optional_scalar")]
    pub toxic_distance: Option<Length>,
    /// The wind speed the screen ran at.
    #[serde(serialize_with = "scalar")]
    pub wind_speed: Velocity,
    /// The Pasquill class the screen ran at.
    pub stability_class: String,
    /// `gaussian_plume_briggs` or `heavy_gas_britter_mcquaid`.
    pub screening_basis: String,
    /// Whether the LFL distance exists and is above zero.
    pub has_flammable_cloud: bool,
    /// Whether the toxic distance exists and is above zero.
    pub has_toxic_endpoint: bool,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

/// The class's name normalization: lower case, with `-`, `_` and spaces stripped.
fn normalize(name: &str) -> String {
    name.to_lowercase().replace(['-', '_', ' '], "")
}

/// The LFL of a component, by the class's table.
fn lfl_of(name: &str) -> Option<f64> {
    let key = normalize(name);
    LOWER_FLAMMABLE_LIMITS
        .iter()
        .find(|(word, _)| *word == key)
        .map(|(_, value)| *value)
}

/// The gas phase's fuel data.
struct FuelData {
    mole_fraction: f64,
    mass_fraction: f64,
    molar_mass: f64,
    lfl: Option<f64>,
}

impl FuelData {
    fn has_fuel(&self) -> bool {
        self.mole_fraction > 0.0 && self.mass_fraction > 0.0 && self.lfl.is_some()
    }
}

/// The class's `calculateFuelData`, on the phase's composition and molar masses.
fn fuel_data(
    names: &[String],
    phase_z: &[f64],
    molar_masses: &[Option<f64>],
    phase_molar_mass: f64,
    manual: Option<f64>,
) -> Result<FuelData> {
    if let Some(lfl) = manual {
        return Ok(FuelData {
            mole_fraction: 1.0,
            mass_fraction: 1.0,
            molar_mass: phase_molar_mass,
            lfl: Some(lfl),
        });
    }
    let mut total_mass = 0.0;
    let mut fuel_mole = 0.0;
    let mut fuel_mass = 0.0;
    let mut reciprocal = 0.0;
    for (index, (name, fraction)) in names.iter().zip(phase_z).enumerate() {
        let mass = molar_masses[index].ok_or_else(|| {
            AzothError::invalid_input(
                "molar_mass",
                format!("component {index} carries no molar mass"),
            )
        })?;
        let fraction = fraction.max(0.0);
        total_mass += fraction * mass;
        if let Some(lfl) = lfl_of(name) {
            if fraction > 0.0 {
                fuel_mole += fraction;
                fuel_mass += fraction * mass;
                reciprocal += fraction / lfl;
            }
        }
    }
    if fuel_mole <= 0.0 || reciprocal <= 0.0 || total_mass <= 0.0 {
        return Ok(FuelData {
            mole_fraction: 0.0,
            mass_fraction: 0.0,
            molar_mass: f64::NAN,
            lfl: None,
        });
    }
    Ok(FuelData {
        mole_fraction: fuel_mole,
        mass_fraction: fuel_mass / total_mass,
        molar_mass: fuel_mass / fuel_mole,
        lfl: Some(fuel_mole / reciprocal),
    })
}

/// The toxic component's phase fractions.
struct ComponentData {
    mole_fraction: f64,
    mass_fraction: f64,
    molar_mass: f64,
}

/// The class's `componentData`.
fn component_data(
    names: &[String],
    phase_z: &[f64],
    molar_masses: &[Option<f64>],
    target: &str,
) -> Result<ComponentData> {
    let wanted = normalize(target);
    let mut total_mass = 0.0;
    let mut component_mole = 0.0;
    let mut component_mass = 0.0;
    let mut component_molar_mass = f64::NAN;
    for (index, (name, fraction)) in names.iter().zip(phase_z).enumerate() {
        let mass = molar_masses[index].ok_or_else(|| {
            AzothError::invalid_input(
                "molar_mass",
                format!("component {index} carries no molar mass"),
            )
        })?;
        let fraction = fraction.max(0.0);
        total_mass += fraction * mass;
        if normalize(name) == wanted {
            component_mole = fraction;
            component_mass = fraction * mass;
            component_molar_mass = mass;
        }
    }
    if component_mole <= 0.0 || total_mass <= 0.0 {
        return Ok(ComponentData {
            mole_fraction: 0.0,
            mass_fraction: 0.0,
            molar_mass: f64::NAN,
        });
    }
    Ok(ComponentData {
        mole_fraction: component_mole,
        mass_fraction: component_mass / total_mass,
        molar_mass: component_molar_mass,
    })
}

/// The dense-gas distance for a concentration ratio; `None` outside `(0, 1)`.
fn dense_distance(
    volumetric_rate: f64,
    source_density: f64,
    air_density: f64,
    wind_speed: f64,
    ratio: f64,
) -> Option<f64> {
    if ratio.is_nan() || ratio <= 0.0 || ratio >= 1.0 {
        return None;
    }
    let g0 = 9.81 * (source_density - air_density) / air_density;
    let alpha = (g0 * g0 * volumetric_rate / wind_speed.powi(5)).powf(0.2);
    let scale = (volumetric_rate / wind_speed).sqrt();
    let log_c = ratio.log10();
    let b = -0.25 * alpha.max(0.01).log10();
    let log_d_over_l = -0.85 * log_c + b;
    Some(scale * 10.0f64.powf(log_d_over_l))
}

/// The class's Gaussian cloud-volume integral above the LFL.
fn gaussian_cloud_volume(
    plume_rate: f64,
    stack_height: f64,
    wind_speed: f64,
    stability: Stability,
    terrain: Terrain,
    threshold: f64,
    distance: f64,
) -> f64 {
    if !distance.is_finite() || distance <= 1.0 || threshold <= 0.0 {
        return 0.0;
    }
    let steps = CLOUD_STEPS;
    let start = 1.0;
    let dx = (distance - start) / steps as f64;
    let mut volume = 0.0;
    for index in 0..steps {
        let x = start + (index as f64 + 0.5) * dx;
        let centreline = concentration(
            plume_rate,
            stack_height,
            wind_speed,
            stability,
            terrain,
            x,
            0.0,
            0.0,
        );
        if centreline <= threshold {
            continue;
        }
        let factor = (2.0 * (centreline / threshold).ln()).sqrt();
        let semi_width = sigma_y(x, stability, terrain) * factor;
        let semi_height = sigma_z(x, stability, terrain) * factor;
        volume += 0.5 * std::f64::consts::PI * semi_width * semi_height * dx;
    }
    volume
}

/// The class's dense-gas cloud-volume estimate.
fn heavy_cloud_volume(distance: f64) -> f64 {
    if !distance.is_finite() || distance <= 0.0 {
        return 0.0;
    }
    let width = (0.25 * distance).max(1.0);
    let height = (0.05 * distance).clamp(1.0, 5.0);
    0.5 * distance * width * height
}

/// The flammable and toxic endpoints of a continuous gas release at ambient conditions.
///
/// # Errors
/// `InvalidInput` for a non-positive release rate, wind speed or stated manual LFL, an
/// unknown stability, terrain or model word, a composition that is not one entry per
/// component and sums to one, an endpoint fraction outside `(0, 1]`, or whatever the
/// databank or the ambient flash refuses.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn gas_dispersion(
    components: &[String],
    composition: &[f64],
    mass_release_rate: MassRate,
    ambient_temperature: Option<ThermodynamicTemperature>,
    atmospheric_pressure: Option<Pressure>,
    wind_speed: Option<Velocity>,
    pasquill_stability_class: Option<&str>,
    surface_roughness: Option<Length>,
    offshore: Option<bool>,
    release_height: Option<Length>,
    model_selection: Option<&str>,
    lower_flammable_limit: Option<f64>,
    flammable_endpoint_fraction_of_lfl: Option<f64>,
    toxic_component: Option<&str>,
    toxic_threshold: Option<f64>,
) -> Result<GasDispersionResult> {
    let spec = &model_gen::GAS_DISPERSION_SPEC;
    let mut warnings = Vec::new();
    let ambient_t = ambient_temperature.unwrap_or_else(|| kelvins(DEFAULT_AMBIENT_TEMPERATURE));
    let ambient_p = atmospheric_pressure.unwrap_or_else(|| pascals(DEFAULT_ATMOSPHERIC_PRESSURE));
    let wind = wind_speed.unwrap_or_else(|| meters_per_second(DEFAULT_WIND_SPEED));
    let roughness = surface_roughness.unwrap_or_else(|| meters(DEFAULT_SURFACE_ROUGHNESS));
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "mass_release_rate" => Some(mass_release_rate.value),
            "wind_speed" => Some(wind.value),
            "lower_flammable_limit" => lower_flammable_limit,
            _ => None,
        },
        &mut warnings,
    )?;

    let stability = Stability::parse(pasquill_stability_class.unwrap_or(DEFAULT_STABILITY))?;
    let selection = Selection::parse(model_selection)?;
    let endpoint_fraction = flammable_endpoint_fraction_of_lfl.unwrap_or(DEFAULT_ENDPOINT_FRACTION);
    if !(0.0 < endpoint_fraction && endpoint_fraction <= 1.0) {
        return Err(AzothError::invalid_input(
            "flammable_endpoint_fraction_of_lfl",
            format!("the endpoint fraction must be in (0, 1], and {endpoint_fraction} is not"),
        ));
    }
    let release_height = release_height
        .map(|height| height.value.max(0.0))
        .unwrap_or(DEFAULT_RELEASE_HEIGHT);
    let terrain = if offshore.unwrap_or(false) || roughness.value < 0.5 {
        Terrain::Rural
    } else {
        Terrain::Urban
    };

    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, _ideal_gas) = databank::mixture_of(&names, Cubic::Pr, None)?;
    if composition.len() != mixture.len() {
        return Err(AzothError::invalid_input(
            "composition",
            format!(
                "a mixture of {} components needs that many mole fractions, and {composition:?} \
                 has {}",
                mixture.len(),
                composition.len()
            ),
        ));
    }
    check_composition(composition, mixture.len(), "composition")?;

    let flash = pt_flash(&mixture, ambient_t, ambient_p, composition)?;
    let (phase_z, phase_root) = match flash.phase {
        Phase::AllLiquid => (flash.x.clone(), flash.z_liquid),
        _ => (flash.y.clone(), flash.z_vapour),
    };
    let molar_masses: Vec<Option<f64>> = mixture
        .components()
        .iter()
        .map(|component| component.molar_mass)
        .collect();
    let phase_molar_mass = entropy_state::molar_mass_of(&mixture, &phase_z)?;
    let molar_volume = pr_molar_volume(phase_root, ambient_t, ambient_p)?.v.value;
    let shift = mixture.volume_shift(&phase_z);
    let source_density = pr_mass_density(
        kilograms_per_mole(phase_molar_mass),
        cubic_meters_per_mole(molar_volume - shift),
    )?
    .rho
    .value;
    let air_density = ambient_p.value / (AIR_GAS_CONSTANT * ambient_t.value);

    let fuel = fuel_data(
        components,
        &phase_z,
        &molar_masses,
        phase_molar_mass,
        lower_flammable_limit,
    )?;
    let selected = match selection {
        Selection::Auto if source_density / air_density > DENSE_GAS_RATIO => Selection::Heavy,
        Selection::Auto => Selection::Gaussian,
        forced => forced,
    };

    let rate = mass_release_rate.value;
    let mut distance_to_flammable_endpoint = None;
    let mut distance_to_half_lfl = None;
    let mut distance_to_lfl = None;
    let mut cloud_volume = 0.0;
    if fuel.has_fuel() {
        let lfl = fuel.lfl.unwrap_or(f64::NAN);
        match selected {
            Selection::Heavy => {
                let volumetric_rate = rate / source_density.max(1.0e-9);
                let lfl_ratio = lfl / fuel.mole_fraction;
                distance_to_lfl = dense_distance(
                    volumetric_rate,
                    source_density,
                    air_density,
                    wind.value,
                    lfl_ratio,
                )
                .map(meters);
                distance_to_half_lfl = dense_distance(
                    volumetric_rate,
                    source_density,
                    air_density,
                    wind.value,
                    0.5 * lfl_ratio,
                )
                .map(meters);
                distance_to_flammable_endpoint = dense_distance(
                    volumetric_rate,
                    source_density,
                    air_density,
                    wind.value,
                    endpoint_fraction * lfl_ratio,
                )
                .map(meters);
                cloud_volume =
                    distance_to_lfl.map_or(0.0, |distance| heavy_cloud_volume(distance.value));
            }
            _ => {
                let fuel_rate = rate * fuel.mass_fraction;
                let lfl_kg =
                    lfl * fuel.molar_mass * ambient_p.value / (GAS_CONSTANT * ambient_t.value);
                distance_to_lfl = distance_to_threshold(
                    fuel_rate,
                    release_height,
                    wind.value,
                    stability,
                    terrain,
                    lfl_kg,
                )
                .map(meters);
                distance_to_half_lfl = distance_to_threshold(
                    fuel_rate,
                    release_height,
                    wind.value,
                    stability,
                    terrain,
                    0.5 * lfl_kg,
                )
                .map(meters);
                distance_to_flammable_endpoint = distance_to_threshold(
                    fuel_rate,
                    release_height,
                    wind.value,
                    stability,
                    terrain,
                    endpoint_fraction * lfl_kg,
                )
                .map(meters);
                cloud_volume = distance_to_lfl.map_or(0.0, |distance| {
                    gaussian_cloud_volume(
                        fuel_rate,
                        release_height,
                        wind.value,
                        stability,
                        terrain,
                        lfl_kg,
                        distance.value,
                    )
                });
            }
        }
    }

    let asked_toxic = toxic_component
        .map(str::trim)
        .is_some_and(|name| !name.is_empty())
        && toxic_threshold.is_some_and(|threshold| threshold > 0.0);
    let mut toxic_distance = None;
    let mut reported_threshold = None;
    let mut toxic_name = String::new();
    if let (true, Some(name), Some(threshold)) = (asked_toxic, toxic_component, toxic_threshold) {
        toxic_name = name.to_string();
        reported_threshold = Some(threshold);
        let component = component_data(components, &phase_z, &molar_masses, name)?;
        if component.mole_fraction > 0.0
            && component.mass_fraction > 0.0
            && component.molar_mass.is_finite()
        {
            match selected {
                Selection::Heavy => {
                    let volumetric_rate = rate / source_density.max(1.0e-9);
                    let ratio = threshold / component.mole_fraction;
                    toxic_distance = dense_distance(
                        volumetric_rate,
                        source_density,
                        air_density,
                        wind.value,
                        ratio,
                    )
                    .map(meters);
                }
                _ => {
                    let toxic_rate = rate * component.mass_fraction;
                    let threshold_kg = threshold * component.molar_mass * ambient_p.value
                        / (TOXIC_GAS_CONSTANT * ambient_t.value);
                    toxic_distance = distance_to_threshold(
                        toxic_rate,
                        release_height,
                        wind.value,
                        stability,
                        terrain,
                        threshold_kg,
                    )
                    .map(meters);
                }
            }
        }
    }

    let basis = match selected {
        Selection::Heavy => "heavy_gas_britter_mcquaid",
        _ => "gaussian_plume_briggs",
    };
    Ok(GasDispersionResult {
        selected_model: match selected {
            Selection::Heavy => "heavy_gas",
            _ => "gaussian_plume",
        }
        .to_string(),
        mass_release_rate: kilograms_per_second(rate),
        flammable_mass_release_rate: kilograms_per_second(rate * fuel.mass_fraction),
        source_density: kilograms_per_cubic_meter(source_density),
        air_density: kilograms_per_cubic_meter(air_density),
        fuel_mole_fraction: fuel.mole_fraction,
        fuel_mass_fraction: fuel.mass_fraction,
        lower_flammable_limit: fuel.lfl,
        flammable_endpoint_fraction: endpoint_fraction,
        distance_to_flammable_endpoint,
        distance_to_half_lfl,
        distance_to_lfl,
        flammable_cloud_volume: cubic_meters(cloud_volume),
        toxic_component: toxic_name,
        toxic_threshold: reported_threshold,
        toxic_distance,
        wind_speed: wind,
        stability_class: pasquill_stability_class
            .unwrap_or(DEFAULT_STABILITY)
            .to_string(),
        screening_basis: basis.to_string(),
        has_flammable_cloud: distance_to_lfl.is_some_and(|distance| distance.value > 0.0),
        has_toxic_endpoint: toxic_distance.is_some_and(|distance| distance.value > 0.0),
        warnings,
    })
}
