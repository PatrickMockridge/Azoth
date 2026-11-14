//! `process.jet_fire` - the API 521 point-source jet fire as a registered id.
//!
//! Spec: `specs/models/process/jet_fire.toml`. The arithmetic is
//! `process/safety/fire/JetFireModel.java` at the pin: the radiated power, Wayne's
//! transmissivity, the incident flux and its fixed-point inverse.

use azoth_core::units::{
    HeatFluxDensity, Length, MassRate, Power, meters, watts, watts_per_square_meter,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{optional_scalar, scalar, warnings as wire_warnings};
use crate::model_gen;

/// The fixed point's iteration cap, the class's own.
const INVERSE_STEPS: usize = 25;
/// The fixed point's stopping width, m.
const INVERSE_WIDTH: f64 = 0.01;
/// The class's default relative humidity.
const DEFAULT_HUMIDITY: f64 = 0.5;

/// Result of `process.jet_fire`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct JetFireResult {
    /// `eta m_dot dHc`, the total heat the jet radiates.
    #[serde(serialize_with = "scalar")]
    pub total_radiative_power: Power,
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

/// Wayne's simplified atmospheric transmissivity; shared with `process.bleve`.
pub(crate) fn wayne_transmissivity(distance: f64, relative_humidity: f64) -> f64 {
    let water = relative_humidity * 1700.0;
    let xpw = water * distance.max(1.0);
    let log = xpw.log10();
    (1.006 - 0.01171 * log - 0.02368 * log * log).clamp(0.0, 1.0)
}

/// The class's fixed-point inverse of the point-source flux.
fn distance_to_flux(power: f64, relative_humidity: f64, target: f64) -> f64 {
    let denominator = 4.0 * std::f64::consts::PI * target;
    let mut radius = (power / denominator).sqrt();
    for _ in 0..INVERSE_STEPS {
        let tau = wayne_transmissivity(radius, relative_humidity);
        let next = (power * tau / denominator).sqrt();
        if (next - radius).abs() < INVERSE_WIDTH {
            radius = next;
            break;
        }
        radius = next;
    }
    radius
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

/// The point-source jet-fire radiation at one distance, or the distance to a target flux.
///
/// # Errors
/// `InvalidInput` for a value the class's constructor refuses, a non-positive distance or
/// target, a humidity outside `(0, 1]`, or a non-finite input.
#[allow(clippy::too_many_arguments)]
pub fn jet_fire(
    mass_flow: MassRate,
    heat_of_combustion: f64,
    radiative_fraction: f64,
    relative_humidity: Option<f64>,
    distance: Option<Length>,
    target_flux: Option<HeatFluxDensity>,
) -> Result<JetFireResult> {
    let spec = &model_gen::JET_FIRE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "mass_flow" => Some(mass_flow.value),
            "heat_of_combustion" => Some(heat_of_combustion),
            "radiative_fraction" => Some(radiative_fraction),
            "relative_humidity" => relative_humidity,
            "distance" => distance.map(|value| value.value),
            "target_flux" => target_flux.map(|value| value.value),
            _ => None,
        },
        &mut warnings,
    )?;

    require_finite("mass_flow", mass_flow.value)?;
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

    let power = radiative_fraction * mass_flow.value * heat_of_combustion;
    let stated_distance = distance.map(|value| value.value);
    let stated_target = target_flux.map(|value| value.value);

    Ok(JetFireResult {
        total_radiative_power: watts(power),
        transmissivity: stated_distance.map(|value| wayne_transmissivity(value, humidity)),
        incident_heat_flux: stated_distance.map(|value| {
            watts_per_square_meter(
                power * wayne_transmissivity(value, humidity)
                    / (4.0 * std::f64::consts::PI * value * value),
            )
        }),
        distance_to_flux: stated_target
            .map(|value| meters(distance_to_flux(power, humidity, value))),
        warnings,
    })
}
