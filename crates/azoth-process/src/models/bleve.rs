//! `process.bleve` - the CCPS/TNO fireball as a registered id.
//!
//! Spec: `specs/models/process/bleve.toml`. The arithmetic is
//! `process/safety/fire/BLEVECalculator.java` at the pin: the fireball's diameter and
//! duration, the surface emissive power, the capped point-source view factor and its
//! bisection inverse. The transmissivity is `process.jet_fire`'s.

use azoth_core::units::{
    HeatFluxDensity, Length, Mass, Time, meters, seconds, watts_per_square_meter,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{optional_scalar, scalar, warnings as wire_warnings};
use crate::model_gen;
use crate::models::jet_fire::wayne_transmissivity;

/// The diameter correlation's coefficient, m per kg^(1/3).
const DIAMETER_COEFFICIENT: f64 = 5.8;
/// The duration correlation's coefficient, s per kg^(1/3).
const DURATION_COEFFICIENT: f64 = 0.45;
/// The bisection's cap.
const BISECTION_STEPS: usize = 60;
/// The bisection's far end, m.
const SEARCH_MAX: f64 = 5000.0;
/// The near bound's clearance beyond the fireball radius, m.
const NEAR_CLEARANCE: f64 = 0.1;
/// The class's default relative humidity.
const DEFAULT_HUMIDITY: f64 = 0.5;

/// Result of `process.bleve`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct BleveResult {
    /// `5.8 M^(1/3)`, the fireball's diameter.
    #[serde(serialize_with = "scalar")]
    pub fireball_diameter: Length,
    /// `0.45 M^(1/3)`, the fireball's lifetime.
    #[serde(serialize_with = "scalar")]
    pub fireball_duration: Time,
    /// `eta M dHc / (pi D^2 t)`, the surface emissive power.
    #[serde(serialize_with = "scalar")]
    pub surface_emissive_power: HeatFluxDensity,
    /// Wayne's transmissivity at the stated distance.
    pub transmissivity: Option<f64>,
    /// The incident heat flux at the stated distance.
    #[serde(serialize_with = "optional_scalar")]
    pub incident_heat_flux: Option<HeatFluxDensity>,
    /// The distance at which the flux falls to the target.
    #[serde(serialize_with = "optional_scalar")]
    pub distance_to_flux: Option<Length>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

/// `5.8 M^(1/3)`, m.
fn fireball_diameter(mass: f64) -> f64 {
    DIAMETER_COEFFICIENT * mass.cbrt()
}

/// `0.45 M^(1/3)`, s.
fn fireball_duration(mass: f64) -> f64 {
    DURATION_COEFFICIENT * mass.cbrt()
}

/// `eta M dHc / (pi D^2 t)`, W/m2.
fn surface_emissive_power(mass: f64, heat_of_combustion: f64, radiative_fraction: f64) -> f64 {
    let diameter = fireball_diameter(mass);
    let duration = fireball_duration(mass);
    radiative_fraction * mass * heat_of_combustion
        / (std::f64::consts::PI * diameter * diameter * duration)
}

/// The incident flux at a distance, W/m2.
fn incident_heat_flux(
    mass: f64,
    heat_of_combustion: f64,
    radiative_fraction: f64,
    relative_humidity: f64,
    distance: f64,
) -> f64 {
    if distance <= 0.0 {
        return f64::INFINITY;
    }
    let diameter = fireball_diameter(mass);
    let view = (diameter * diameter / (4.0 * distance * distance)).min(1.0);
    surface_emissive_power(mass, heat_of_combustion, radiative_fraction)
        * view
        * wayne_transmissivity(distance, relative_humidity)
}

/// The class's bisection, absent when the fireball never reaches the target.
fn distance_to_flux(
    mass: f64,
    heat_of_combustion: f64,
    radiative_fraction: f64,
    relative_humidity: f64,
    target: f64,
) -> Option<f64> {
    let mut lo = fireball_diameter(mass) / 2.0 + NEAR_CLEARANCE;
    let mut hi = SEARCH_MAX;
    if incident_heat_flux(
        mass,
        heat_of_combustion,
        radiative_fraction,
        relative_humidity,
        lo,
    ) < target
    {
        return None;
    }
    for _ in 0..BISECTION_STEPS {
        let mid = 0.5 * (lo + hi);
        if incident_heat_flux(
            mass,
            heat_of_combustion,
            radiative_fraction,
            relative_humidity,
            mid,
        ) > target
        {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

/// Refuse a value the arithmetic cannot carry to a published number.
fn require_finite(name: &'static str, value: f64) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(AzothError::invalid_input(
            name,
            "a non-finite value is not a number this calculation publishes",
        ))
    }
}

/// The fireball radiation at one distance, or the distance to a target flux.
///
/// # Errors
/// `InvalidInput` for a value the class's constructor refuses, a radiative fraction
/// outside `[0, 1]`, a non-positive distance or target, a humidity outside `(0, 1]`, or a
/// non-finite input.
pub fn bleve(
    fuel_mass: Mass,
    heat_of_combustion: f64,
    radiative_fraction: f64,
    relative_humidity: Option<f64>,
    distance: Option<Length>,
    target_flux: Option<HeatFluxDensity>,
) -> Result<BleveResult> {
    let spec = &model_gen::BLEVE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "fuel_mass" => Some(fuel_mass.value),
            "heat_of_combustion" => Some(heat_of_combustion),
            "radiative_fraction" => Some(radiative_fraction),
            "relative_humidity" => relative_humidity,
            "distance" => distance.map(|value| value.value),
            "target_flux" => target_flux.map(|value| value.value),
            _ => None,
        },
        &mut warnings,
    )?;

    require_finite("fuel_mass", fuel_mass.value)?;
    require_finite("heat_of_combustion", heat_of_combustion)?;
    require_finite("radiative_fraction", radiative_fraction)?;
    let humidity = relative_humidity.unwrap_or(DEFAULT_HUMIDITY);
    require_finite("relative_humidity", humidity)?;
    if let Some(value) = distance {
        require_finite("distance", value.value)?;
    }
    if let Some(value) = target_flux {
        require_finite("target_flux", value.value)?;
    }

    let mass = fuel_mass.value;
    let stated_distance = distance.map(|value| value.value);

    Ok(BleveResult {
        fireball_diameter: meters(fireball_diameter(mass)),
        fireball_duration: seconds(fireball_duration(mass)),
        surface_emissive_power: watts_per_square_meter(surface_emissive_power(
            mass,
            heat_of_combustion,
            radiative_fraction,
        )),
        transmissivity: stated_distance.map(|value| wayne_transmissivity(value, humidity)),
        incident_heat_flux: stated_distance.map(|value| {
            watts_per_square_meter(incident_heat_flux(
                mass,
                heat_of_combustion,
                radiative_fraction,
                humidity,
                value,
            ))
        }),
        distance_to_flux: target_flux.and_then(|value| {
            distance_to_flux(
                mass,
                heat_of_combustion,
                radiative_fraction,
                humidity,
                value.value,
            )
            .map(meters)
        }),
        warnings,
    })
}
