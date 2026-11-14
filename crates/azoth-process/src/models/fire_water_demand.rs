//! `process.fire_water_demand` - the fire-water demand as a registered id.
//!
//! Spec: `specs/models/process/fire_water_demand.toml`. The arithmetic is
//! `process/safety/firewater/FireWaterDemandCalculator.java` at the pin: the area and
//! object demand, the simultaneous factor, the hourly rate and the water and foam volumes.

use azoth_core::units::{
    Area, Time, Volume, VolumeRate, cubic_meters, cubic_meters_per_hour, liters_per_minute,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{scalar, warnings as wire_warnings};
use crate::model_gen;

/// The class's NORSOK process-area rate, `L/(min m2)`.
const DEFAULT_APPLICATION_RATE_LPM_M2: f64 = 10.0;
/// The class's default duration, s.
const DEFAULT_DURATION_S: f64 = 1800.0;
/// The class's default foam concentrate percent.
const DEFAULT_FOAM_CONCENTRATE_PERCENT: f64 = 0.0;
/// The class's default simultaneous area factor.
const DEFAULT_SIMULTANEOUS_AREA_FACTOR: f64 = 1.0;

/// One dedicated object: a wetted surface and its own application rate.
#[derive(Debug, Clone, PartialEq)]
pub struct ProtectedObject {
    /// The wetted surface.
    pub surface_area: Area,
    /// The object's application rate in `L/(min m2)`; dimensionless.
    pub application_rate: f64,
}

/// Result of `process.fire_water_demand`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct FireWaterDemandResult {
    /// The general-area demand.
    #[serde(serialize_with = "scalar")]
    pub area_demand: VolumeRate,
    /// The summed object demand.
    #[serde(serialize_with = "scalar")]
    pub object_demand: VolumeRate,
    /// The summed demand after the simultaneous factor.
    #[serde(serialize_with = "scalar")]
    pub total_demand: VolumeRate,
    /// The total demand as an hourly rate.
    #[serde(serialize_with = "scalar")]
    pub total_demand_m3_per_hour: VolumeRate,
    /// The water volume over the duration.
    #[serde(serialize_with = "scalar")]
    pub water_volume: Volume,
    /// The foam concentrate volume.
    #[serde(serialize_with = "scalar")]
    pub foam_concentrate_volume: Volume,
    /// `object_demand / area_demand`; absent where the area demand is zero.
    pub object_to_area_demand_ratio: Option<f64>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
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

/// The area and object demand, the simultaneous factor and the volumes.
///
/// # Errors
/// `InvalidInput` for a negative protected area, rate, row area or row rate, a non-positive
/// duration, a foam percent outside `[0, 10]`, a simultaneous factor below one, or a
/// non-finite input.
pub fn fire_water_demand(
    protected_area: Area,
    application_rate: Option<f64>,
    duration: Option<Time>,
    foam_concentrate_percent: Option<f64>,
    simultaneous_area_factor: Option<f64>,
    protected_objects: Option<&[ProtectedObject]>,
) -> Result<FireWaterDemandResult> {
    let spec = &model_gen::FIRE_WATER_DEMAND_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "protected_area" => Some(protected_area.value),
            "application_rate" => application_rate,
            "duration" => duration.map(|value| value.value),
            "foam_concentrate_percent" => foam_concentrate_percent,
            "simultaneous_area_factor" => simultaneous_area_factor,
            _ => None,
        },
        &mut warnings,
    )?;

    require_finite("protected_area", protected_area.value)?;
    let rate = application_rate.unwrap_or(DEFAULT_APPLICATION_RATE_LPM_M2);
    require_finite("application_rate", rate)?;
    let duration_s = duration
        .map(|value| value.value)
        .unwrap_or(DEFAULT_DURATION_S);
    require_finite("duration", duration_s)?;
    let foam = foam_concentrate_percent.unwrap_or(DEFAULT_FOAM_CONCENTRATE_PERCENT);
    require_finite("foam_concentrate_percent", foam)?;
    let factor = simultaneous_area_factor.unwrap_or(DEFAULT_SIMULTANEOUS_AREA_FACTOR);
    require_finite("simultaneous_area_factor", factor)?;

    let mut object_demand = 0.0;
    for object in protected_objects.unwrap_or(&[]) {
        require_finite("protected_objects.surface_area", object.surface_area.value)?;
        require_finite(
            "protected_objects.application_rate",
            object.application_rate,
        )?;
        if object.surface_area.value < 0.0 || object.application_rate < 0.0 {
            return Err(AzothError::invalid_input(
                "protected_objects",
                "a negative surface area or application rate is not an object the demand can \
                 cover; the class's `addObject` accepts it and publishes a negative demand",
            ));
        }
        object_demand += object.surface_area.value * object.application_rate;
    }

    // The class's operation order, to the letter.
    let area_demand = protected_area.value * rate;
    let total = (area_demand + object_demand) * factor;
    let m3_per_hour = total * 60.0 / 1000.0;
    let water = total * duration_s / 60_000.0;
    let foam_volume = water * foam / 100.0;
    let ratio = if area_demand <= 0.0 {
        None
    } else {
        Some(object_demand / area_demand)
    };

    Ok(FireWaterDemandResult {
        area_demand: liters_per_minute(area_demand),
        object_demand: liters_per_minute(object_demand),
        total_demand: liters_per_minute(total),
        total_demand_m3_per_hour: cubic_meters_per_hour(m3_per_hour),
        water_volume: cubic_meters(water),
        foam_concentrate_volume: cubic_meters(foam_volume),
        object_to_area_demand_ratio: ratio,
        warnings,
    })
}
