//! `process.vce` - the TNO Multi-Energy overpressure as a registered id.
//!
//! Spec: `specs/models/process/vce.toml`. The arithmetic is
//! `process/safety/fire/VCEModel.java` at the pin: the Sachs-scaled distance, the
//! strength-class plateau and power-law fit, and its bisection inverse.

use azoth_core::units::{Energy, Length, Mass, Pressure, joules, meters, pascals};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{optional_scalar, scalar, warnings as wire_warnings};
use crate::model_gen;

/// The ten strength-class plateaus at `Rbar <= 1`, dimensionless.
const PLATEAU: [f64; 10] = [0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.00, 2.00, 5.00, 10.0];
/// The bisection's cap.
const BISECTION_STEPS: usize = 80;
/// The bisection's near end, m.
const SEARCH_MIN: f64 = 0.1;
/// The bisection's far end, m.
const SEARCH_MAX: f64 = 50_000.0;
/// The class's default ambient pressure, Pa.
const DEFAULT_AMBIENT_PRESSURE: f64 = 101_325.0;

/// Result of `process.vce`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct VceResult {
    /// `m dHc`, the cloud's combustion energy.
    #[serde(serialize_with = "scalar")]
    pub combustion_energy: Energy,
    /// `d / cbrt(E / p0)`, the Sachs-scaled distance.
    pub scaled_distance: Option<f64>,
    /// The strength-class plateau or power-law value.
    pub dimensionless_overpressure: Option<f64>,
    /// The dimensionless overpressure times the ambient pressure.
    #[serde(serialize_with = "optional_scalar")]
    pub overpressure: Option<Pressure>,
    /// The distance at which the overpressure falls to the target.
    #[serde(serialize_with = "optional_scalar")]
    pub distance_to_overpressure: Option<Length>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

/// `m dHc`, J.
fn combustion_energy(mass: f64, heat_of_combustion: f64) -> f64 {
    mass * heat_of_combustion
}

/// `d / cbrt(E / p0)`, the Sachs-scaled distance.
fn scaled_distance(distance: f64, energy: f64, ambient_pressure: f64) -> f64 {
    distance / (energy / ambient_pressure).cbrt()
}

/// The class's plateau or power-law dimensionless overpressure.
fn dimensionless_overpressure(
    distance: f64,
    energy: f64,
    ambient_pressure: f64,
    strength_class: usize,
) -> f64 {
    let rbar = scaled_distance(distance, energy, ambient_pressure);
    if rbar <= 0.0 {
        return f64::INFINITY;
    }
    let plateau = PLATEAU[strength_class - 1];
    if rbar <= 1.0 {
        plateau
    } else {
        let beta = if strength_class <= 5 { 1.0 } else { 1.2 };
        plateau * rbar.powf(-beta)
    }
}

/// The side-on overpressure at a distance, Pa.
fn overpressure(distance: f64, energy: f64, ambient_pressure: f64, strength_class: usize) -> f64 {
    dimensionless_overpressure(distance, energy, ambient_pressure, strength_class)
        * ambient_pressure
}

/// The class's bisection, absent when the cloud never reaches the target.
fn distance_to_overpressure(
    energy: f64,
    ambient_pressure: f64,
    strength_class: usize,
    target: f64,
) -> Option<f64> {
    let mut lo = SEARCH_MIN;
    let mut hi = SEARCH_MAX;
    if overpressure(lo, energy, ambient_pressure, strength_class) < target {
        return None;
    }
    for _ in 0..BISECTION_STEPS {
        let mid = 0.5 * (lo + hi);
        if overpressure(mid, energy, ambient_pressure, strength_class) > target {
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

/// The TNO Multi-Energy overpressure at one distance, or its inverse.
///
/// # Errors
/// `InvalidInput` for a value the class's constructor refuses, a non-positive distance or
/// target, a non-integral strength class, or a non-finite input.
pub fn vce(
    cloud_fuel_mass: Mass,
    heat_of_combustion: f64,
    strength_class: f64,
    ambient_pressure: Option<Pressure>,
    distance: Option<Length>,
    target_overpressure: Option<Pressure>,
) -> Result<VceResult> {
    let spec = &model_gen::VCE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "cloud_fuel_mass" => Some(cloud_fuel_mass.value),
            "heat_of_combustion" => Some(heat_of_combustion),
            "strength_class" => Some(strength_class),
            "ambient_pressure" => ambient_pressure.map(|value| value.value),
            "distance" => distance.map(|value| value.value),
            "target_overpressure" => target_overpressure.map(|value| value.value),
            _ => None,
        },
        &mut warnings,
    )?;

    require_finite("cloud_fuel_mass", cloud_fuel_mass.value)?;
    require_finite("heat_of_combustion", heat_of_combustion)?;
    require_finite("strength_class", strength_class)?;
    let ambient = ambient_pressure.map_or(DEFAULT_AMBIENT_PRESSURE, |value| value.value);
    require_finite("ambient_pressure", ambient)?;
    if let Some(value) = distance {
        require_finite("distance", value.value)?;
    }
    if let Some(value) = target_overpressure {
        require_finite("target_overpressure", value.value)?;
    }
    if strength_class.fract() != 0.0 {
        return Err(AzothError::invalid_input(
            "strength_class",
            "the TNO strength class is a whole number from 1 to 10",
        ));
    }

    let class = strength_class as usize;
    let energy = combustion_energy(cloud_fuel_mass.value, heat_of_combustion);
    let stated_distance = distance.map(|value| value.value);

    Ok(VceResult {
        combustion_energy: joules(energy),
        scaled_distance: stated_distance.map(|value| scaled_distance(value, energy, ambient)),
        dimensionless_overpressure: stated_distance
            .map(|value| dimensionless_overpressure(value, energy, ambient, class)),
        overpressure: stated_distance
            .map(|value| pascals(overpressure(value, energy, ambient, class))),
        distance_to_overpressure: target_overpressure.and_then(|value| {
            distance_to_overpressure(energy, ambient, class, value.value).map(meters)
        }),
        warnings,
    })
}
