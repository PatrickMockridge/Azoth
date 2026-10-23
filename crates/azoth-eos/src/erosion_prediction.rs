//! `eos.erosion_prediction` - NeqSim's `ErosionPredictionCalculator`.
//!
//! The API RP 14E erosional velocity and the DNV RP O501 sand-erosion rate, with the
//! geometry, material and particle branches the class carries, then the cumulative wall
//! loss and the risk band.

use azoth_core::units::{
    DynamicViscosity, Length, MassDensity, MassFlux, Time, Velocity, meters, meters_per_second,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::{ErosionPredictionResult, ErosionRiskLevel};
use crate::spec_gen;

/// The class's kilograms per cubic metre to pounds per cubic foot factor.
const LB_FT3_PER_KG_M3: f64 = 0.062428;
/// The class's 365.25-day year.
const SECONDS_PER_YEAR: f64 = 31_557_600.0;
/// The reference particle diameter and density the DNV rate is normalised by.
const REFERENCE_PARTICLE_DIAMETER: f64 = 0.00025;
const REFERENCE_PARTICLE_DENSITY: f64 = 2650.0;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Geometry {
    Straight,
    Elbow,
    Tee,
    BlindTee,
    Reducer,
    Choke,
    Weld,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Material {
    CarbonSteel,
    DuplexSteel,
    SuperDuplex,
    ThirteenCr,
    Inconel,
    Titanium,
    TungstenCarbide,
}

impl Geometry {
    fn parse(name: &str) -> Result<Self> {
        match name {
            "straight" => Ok(Self::Straight),
            "elbow" => Ok(Self::Elbow),
            "tee" => Ok(Self::Tee),
            "blind_tee" => Ok(Self::BlindTee),
            "reducer" => Ok(Self::Reducer),
            "choke" => Ok(Self::Choke),
            "weld" => Ok(Self::Weld),
            other => Err(AzothError::invalid_input(
                "geometry",
                format!(
                    "`{other}` is not one of the geometry values: straight, elbow, tee, \
                     blind_tee, reducer, choke or weld"
                ),
            )),
        }
    }

    fn geometry_factor(self, bend_radius_ratio: f64) -> f64 {
        match self {
            Self::Straight => 0.0,
            Self::Elbow => {
                if bend_radius_ratio <= 1.0 {
                    1.0
                } else if bend_radius_ratio <= 1.5 {
                    0.8
                } else if bend_radius_ratio <= 3.0 {
                    0.5
                } else {
                    0.3
                }
            }
            Self::Tee => 0.6,
            Self::BlindTee => 0.3,
            Self::Reducer => 0.4,
            Self::Choke => 2.0,
            Self::Weld => 1.2,
        }
    }

    fn impact_angle_degrees(self, bend_radius_ratio: f64) -> f64 {
        match self {
            Self::Straight => 1.0,
            Self::Elbow => (45.0 / bend_radius_ratio).max(5.0),
            Self::Tee => 45.0,
            Self::BlindTee => 90.0,
            Self::Reducer => 15.0,
            Self::Choke => 30.0,
            Self::Weld => 20.0,
        }
    }

    fn eroded_area(self, diameter: f64) -> f64 {
        let radius = diameter / 2.0;
        match self {
            Self::Straight => std::f64::consts::PI * diameter,
            Self::Elbow => std::f64::consts::PI * diameter * diameter * 0.1,
            Self::Tee | Self::BlindTee => std::f64::consts::PI * radius * radius * 0.5,
            Self::Reducer => std::f64::consts::PI * radius * radius * 0.3,
            Self::Choke => std::f64::consts::PI * radius * radius * 0.2,
            Self::Weld => std::f64::consts::PI * diameter * 0.005,
        }
    }
}

impl Material {
    fn parse(name: &str) -> Result<Self> {
        match name {
            "carbon_steel" => Ok(Self::CarbonSteel),
            "duplex_steel" => Ok(Self::DuplexSteel),
            "super_duplex" => Ok(Self::SuperDuplex),
            "13cr" => Ok(Self::ThirteenCr),
            "inconel" => Ok(Self::Inconel),
            "titanium" => Ok(Self::Titanium),
            "tungsten_carbide" => Ok(Self::TungstenCarbide),
            other => Err(AzothError::invalid_input(
                "pipe_material",
                format!(
                    "`{other}` is not one of the material values: carbon_steel, \
                     duplex_steel, super_duplex, 13cr, inconel, titanium or tungsten_carbide"
                ),
            )),
        }
    }

    fn erosion_constant(self) -> f64 {
        match self {
            Self::Inconel | Self::Titanium => 1.0e-9,
            Self::TungstenCarbide => 5.0e-10,
            _ => 2.0e-9,
        }
    }

    fn velocity_exponent(self) -> f64 {
        match self {
            Self::TungstenCarbide => 2.3,
            _ => 2.6,
        }
    }

    fn density(self) -> f64 {
        match self {
            Self::Titanium => 4500.0,
            Self::Inconel => 8440.0,
            _ => 7850.0,
        }
    }
}

/// The API RP 14E erosional velocity and the DNV RP O501 sand-erosion rate.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's hard range.
/// * [`azoth_core::AzothError::InvalidInput`] if a geometry or material name is unknown.
///
/// # Example
/// ```
/// use azoth_core::units::{
///     kilograms_per_cubic_meter, kilograms_per_square_meter_second, meters,
///     meters_per_second, pascal_seconds, seconds,
/// };
/// use azoth_eos::erosion_prediction::erosion_prediction;
///
/// let r = erosion_prediction(
///     kilograms_per_cubic_meter(100.0),
///     meters_per_second(10.0),
///     meters(0.1524),
///     meters(0.01),
///     kilograms_per_square_meter_second(0.0),
///     meters(0.00025),
///     kilograms_per_cubic_meter(2650.0),
///     pascal_seconds(0.001),
///     "elbow",
///     1.5,
///     "carbon_steel",
///     100.0,
///     meters(0.003),
///     seconds(788_940_000.0),
/// )?;
/// assert_eq!(r.risk_level.as_str(), "medium");
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn erosion_prediction(
    mixture_density: MassDensity,
    mixture_velocity: Velocity,
    pipe_diameter: Length,
    wall_thickness: Length,
    sand_rate: MassFlux,
    sand_particle_diameter: Length,
    sand_particle_density: MassDensity,
    mixture_viscosity: DynamicViscosity,
    geometry: &str,
    bend_radius_ratio: f64,
    pipe_material: &str,
    api_c_factor: f64,
    corrosion_allowance: Length,
    design_life: Time,
) -> Result<ErosionPredictionResult> {
    let spec = &spec_gen::EROSION_PREDICTION_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "mixture_density" => Some(mixture_density.value),
            "mixture_velocity" => Some(mixture_velocity.value),
            "pipe_diameter" => Some(pipe_diameter.value),
            "wall_thickness" => Some(wall_thickness.value),
            "sand_rate" => Some(sand_rate.value),
            "sand_particle_diameter" => Some(sand_particle_diameter.value),
            "sand_particle_density" => Some(sand_particle_density.value),
            "mixture_viscosity" => Some(mixture_viscosity.value),
            "bend_radius_ratio" => Some(bend_radius_ratio),
            "api_c_factor" => Some(api_c_factor),
            "corrosion_allowance" => Some(corrosion_allowance.value),
            "design_life" => Some(design_life.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let geometry = Geometry::parse(geometry)?;
    let material = Material::parse(pipe_material)?;
    let diameter = pipe_diameter.value;
    let velocity = mixture_velocity.value;

    let erosional_velocity = if mixture_density.value > 0.0 {
        let density_lb_ft3 = mixture_density.value * LB_FT3_PER_KG_M3;
        (api_c_factor / density_lb_ft3.sqrt()) * 0.3048
    } else {
        0.0
    };
    let velocity_ratio = if erosional_velocity > 0.0 {
        velocity / erosional_velocity
    } else {
        0.0
    };
    let within_api_limits = velocity_ratio <= 1.0;

    let sand_rate_value = sand_rate.value;
    let particle_diameter = sand_particle_diameter.value;
    let stokes = sand_particle_density.value * particle_diameter * particle_diameter * velocity
        / (18.0 * mixture_viscosity.value * diameter);
    let particle_ratio = if stokes > 10.0 {
        0.9
    } else if stokes > 1.0 {
        0.8
    } else {
        0.7
    };
    let particle_velocity = velocity * particle_ratio;

    let sand_rate_kg_yr = sand_rate_value * SECONDS_PER_YEAR;
    let rate_mm_yr = if sand_rate_value > 0.0 && geometry != Geometry::Straight {
        let impact = geometry
            .impact_angle_degrees(bend_radius_ratio)
            .to_radians()
            .sin();
        let impact_function = impact * (1.0 + 0.6 * (1.0 - impact));
        let erosion_kg_kg = material.erosion_constant()
            * impact_function
            * particle_velocity.powf(material.velocity_exponent())
            * (particle_diameter / REFERENCE_PARTICLE_DIAMETER).powf(0.2)
            * geometry.geometry_factor(bend_radius_ratio)
            * (sand_particle_density.value / REFERENCE_PARTICLE_DENSITY);
        let area = geometry.eroded_area(diameter);
        if area > 0.0 && material.density() > 0.0 {
            (erosion_kg_kg * sand_rate_kg_yr / material.density() / area) * 1000.0
        } else {
            0.0
        }
    } else if sand_rate_value > 0.0 {
        material.erosion_constant()
            * particle_velocity.powf(material.velocity_exponent())
            * sand_rate_kg_yr
            * 0.01
    } else {
        0.0
    };

    let life_years = design_life.value / SECONDS_PER_YEAR;
    let cumulative_mm = rate_mm_yr * life_years;
    let remaining_mm =
        wall_thickness.value * 1000.0 - cumulative_mm - corrosion_allowance.value * 1000.0;

    let risk_level = if rate_mm_yr > 1.0 || velocity_ratio > 1.5 {
        ErosionRiskLevel::Critical
    } else if rate_mm_yr > 0.5 || velocity_ratio > 1.0 {
        ErosionRiskLevel::High
    } else if rate_mm_yr > 0.1 || velocity_ratio > 0.8 {
        ErosionRiskLevel::Medium
    } else {
        ErosionRiskLevel::Low
    };

    Ok(ErosionPredictionResult {
        erosional_velocity: meters_per_second(erosional_velocity),
        velocity_ratio,
        erosion_rate: meters_per_second(rate_mm_yr * 1.0e-3 / SECONDS_PER_YEAR),
        cumulative_erosion: meters(cumulative_mm * 1.0e-3),
        remaining_wall_thickness: meters(remaining_mm * 1.0e-3),
        within_api_limits,
        within_erosion_limits: remaining_mm > 0.0,
        risk_level,
        warnings,
    })
}
