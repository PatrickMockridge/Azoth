//! `process.pool_fire` - the API 521 / Mudan pool fire as a registered id.
//!
//! Spec: `specs/models/process/pool_fire.toml`. The arithmetic is
//! `process/safety/fire/PoolFireModel.java` at the pin: Thomas's flame height, Mudan's
//! smoke-corrected SEP, the vertical-cylinder view factor and its bisection inverse.

use azoth_core::units::{
    HeatFluxDensity, Length, MassFlux, MassRate, kilograms_per_second, meters,
    watts_per_square_meter,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{optional_scalar, scalar, warnings as wire_warnings};
use crate::model_gen;

/// Thomas's gravitational constant, m/s2.
const GRAVITY: f64 = 9.81;
/// Mudan's clean-flame emissive power, W/m2.
const SEP_MAX: f64 = 140_000.0;
/// Mudan's smoke-blocked emissive power, W/m2.
const SEP_SMOKE: f64 = 20_000.0;
/// Mudan's smoke-decay constant, 1/m.
const SEP_DECAY: f64 = 0.12;
/// The class's default air density, kg/m3.
const DEFAULT_AIR_DENSITY: f64 = 1.2;
/// The bisection's cap.
const BISECTION_STEPS: usize = 60;
/// The bisection's far end, m.
const SEARCH_MAX: f64 = 1000.0;
/// The near bound's clearance beyond the pool radius, m.
const NEAR_CLEARANCE: f64 = 0.1;
/// The view factor's clearance inside the pool radius, m.
const VIEW_CLEARANCE: f64 = 0.001;

/// Result of `process.pool_fire`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PoolFireResult {
    /// `m'' pi D^2 / 4`, the mass the pool burns per second.
    #[serde(serialize_with = "scalar")]
    pub total_burn_rate: MassRate,
    /// Thomas's flame height above the pool.
    #[serde(serialize_with = "scalar")]
    pub flame_height: Length,
    /// Mudan's smoke-corrected surface emissive power.
    #[serde(serialize_with = "scalar")]
    pub surface_emissive_power: HeatFluxDensity,
    /// Mudan's vertical-cylinder view factor at the stated distance.
    pub view_factor: Option<f64>,
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

/// `m'' pi D^2 / 4`, kg/s.
fn total_burn_rate(diameter: f64, burn_rate: f64) -> f64 {
    burn_rate * 0.25 * std::f64::consts::PI * diameter * diameter
}

/// Thomas's no-wind flame height, m.
fn flame_height(diameter: f64, burn_rate: f64, air_density: f64) -> f64 {
    let term = burn_rate / (air_density * (GRAVITY * diameter).sqrt());
    diameter * 42.0 * term.powf(0.61)
}

/// Mudan's smoke-corrected SEP, W/m2.
fn surface_emissive_power(diameter: f64) -> f64 {
    let clean = (-SEP_DECAY * diameter).exp();
    SEP_MAX * clean + SEP_SMOKE * (1.0 - clean)
}

/// Mudan's vertical-cylinder view factor at a ground-level receiver.
fn view_factor(diameter: f64, burn_rate: f64, air_density: f64, distance: f64) -> f64 {
    let radius = diameter / 2.0;
    let height = flame_height(diameter, burn_rate, air_density) / radius;
    let s = distance.max(radius + VIEW_CLEARANCE) / radius;
    let a = (height * height + s * s + 1.0) / (2.0 * s);
    let f = (1.0 / std::f64::consts::PI)
        * ((height / (s * s - 1.0).sqrt()).atan()
            - (a - (s * s + 1.0) / (2.0 * s)) / (a * a - 1.0).sqrt());
    f.clamp(0.0, 1.0)
}

/// The incident flux at a distance, W/m2.
fn incident_heat_flux(diameter: f64, burn_rate: f64, air_density: f64, distance: f64) -> f64 {
    surface_emissive_power(diameter) * view_factor(diameter, burn_rate, air_density, distance)
}

/// The class's bisection, absent when the fire never reaches the target.
fn distance_to_flux(diameter: f64, burn_rate: f64, air_density: f64, target: f64) -> Option<f64> {
    let mut lo = diameter / 2.0 + NEAR_CLEARANCE;
    let mut hi = SEARCH_MAX;
    if incident_heat_flux(diameter, burn_rate, air_density, lo) < target {
        return None;
    }
    for _ in 0..BISECTION_STEPS {
        let mid = 0.5 * (lo + hi);
        if incident_heat_flux(diameter, burn_rate, air_density, mid) > target {
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

/// The pool-fire radiation at one distance, or the distance to a target flux.
///
/// # Errors
/// `InvalidInput` for a value the class's constructor refuses, a non-positive target, or a
/// non-finite input.
pub fn pool_fire(
    pool_diameter: Length,
    burn_rate: MassFlux,
    air_density: Option<f64>,
    distance: Option<Length>,
    target_flux: Option<HeatFluxDensity>,
) -> Result<PoolFireResult> {
    let spec = &model_gen::POOL_FIRE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "pool_diameter" => Some(pool_diameter.value),
            "burn_rate" => Some(burn_rate.value),
            "air_density" => air_density,
            "distance" => distance.map(|value| value.value),
            "target_flux" => target_flux.map(|value| value.value),
            _ => None,
        },
        &mut warnings,
    )?;

    require_finite("pool_diameter", pool_diameter.value)?;
    require_finite("burn_rate", burn_rate.value)?;
    let air = air_density.unwrap_or(DEFAULT_AIR_DENSITY);
    require_finite("air_density", air)?;
    if let Some(value) = distance {
        require_finite("distance", value.value)?;
    }
    if let Some(value) = target_flux {
        require_finite("target_flux", value.value)?;
    }

    let diameter = pool_diameter.value;
    let burning = burn_rate.value;
    let stated_distance = distance.map(|value| value.value);

    Ok(PoolFireResult {
        total_burn_rate: kilograms_per_second(total_burn_rate(diameter, burning)),
        flame_height: meters(flame_height(diameter, burning, air)),
        surface_emissive_power: watts_per_square_meter(surface_emissive_power(diameter)),
        view_factor: stated_distance.map(|value| view_factor(diameter, burning, air, value)),
        incident_heat_flux: stated_distance
            .map(|value| watts_per_square_meter(incident_heat_flux(diameter, burning, air, value))),
        distance_to_flux: target_flux
            .and_then(|value| distance_to_flux(diameter, burning, air, value.value).map(meters)),
        warnings,
    })
}
