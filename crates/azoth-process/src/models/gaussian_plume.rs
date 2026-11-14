//! `process.gaussian_plume` - the Briggs plume as a registered id.
//!
//! Spec: `specs/models/process/gaussian_plume.toml`. The arithmetic is
//! `process/safety/dispersion/GaussianPlume.java` at the pin: the sigma curves, the two-image
//! concentration, and the log scan with its bisection.

use azoth_core::units::{
    Length, MassDensity, MassRate, Velocity, kilograms_per_cubic_meter, meters,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{optional_scalar, scalar, warnings as wire_warnings};
use crate::model_gen;

/// The scan's logarithmic sample count above the first point.
const SCAN_STEPS: usize = 200;
/// The scan's lower end, m.
const SCAN_MIN_DISTANCE: f64 = 1.0;
/// The scan's upper end, and its answer when the profile never falls below the target, m.
const SCAN_MAX_DISTANCE: f64 = 100_000.0;
/// The bisection cap.
const BISECTION_STEPS: usize = 60;
/// The bisection's stopping width, m.
const BISECTION_WIDTH: f64 = 0.5;

/// One of the six Pasquill classes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Stability {
    A,
    B,
    C,
    D,
    E,
    F,
}

impl Stability {
    pub(crate) fn parse(value: &str) -> Result<Self> {
        Ok(match value {
            "A" => Self::A,
            "B" => Self::B,
            "C" => Self::C,
            "D" => Self::D,
            "E" => Self::E,
            "F" => Self::F,
            other => {
                return Err(AzothError::invalid_input(
                    "stability",
                    format!("a Pasquill class is one of A-F, and {other:?} is not"),
                ));
            }
        })
    }
}

/// The terrain category.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Terrain {
    Rural,
    Urban,
}

impl Terrain {
    pub(crate) fn parse(value: &str) -> Result<Self> {
        Ok(match value {
            "rural" => Self::Rural,
            "urban" => Self::Urban,
            other => {
                return Err(AzothError::invalid_input(
                    "terrain",
                    format!("a terrain category is `rural` or `urban`, and {other:?} is not"),
                ));
            }
        })
    }
}

/// Result of `process.gaussian_plume`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct GaussianPlumeResult {
    /// The horizontal coefficient at the stated distance.
    #[serde(serialize_with = "scalar")]
    pub sigma_y: Length,
    /// The vertical coefficient at the stated distance.
    #[serde(serialize_with = "scalar")]
    pub sigma_z: Length,
    /// The concentration at the stated point.
    #[serde(serialize_with = "scalar")]
    pub concentration: MassDensity,
    /// The concentration at the same distance on the ground centreline.
    #[serde(serialize_with = "scalar")]
    pub centerline_ground_concentration: MassDensity,
    /// The distance to the threshold; absent when the profile never reaches it.
    #[serde(serialize_with = "optional_scalar")]
    pub distance_to_threshold: Option<Length>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

/// `sigma_y = a x / sqrt(1 + b x)`, the class's Briggs table.
pub(crate) fn sigma_y(distance: f64, stability: Stability, terrain: Terrain) -> f64 {
    let (a, b) = match (stability, terrain) {
        (Stability::A, Terrain::Rural) => (0.22, 0.0001),
        (Stability::A, Terrain::Urban) => (0.32, 0.0004),
        (Stability::B, Terrain::Rural) => (0.16, 0.0001),
        (Stability::B, Terrain::Urban) => (0.32, 0.0004),
        (Stability::C, Terrain::Rural) => (0.11, 0.0001),
        (Stability::C, Terrain::Urban) => (0.22, 0.0004),
        (Stability::D, Terrain::Rural) => (0.08, 0.0001),
        (Stability::D, Terrain::Urban) => (0.16, 0.0004),
        (Stability::E, Terrain::Rural) => (0.06, 0.0001),
        (Stability::E, Terrain::Urban) => (0.11, 0.0004),
        (Stability::F, Terrain::Rural) => (0.04, 0.0001),
        (Stability::F, Terrain::Urban) => (0.11, 0.0004),
    };
    a * distance / (1.0 + b * distance).powf(0.5)
}

/// `sigma_z = a x / (1 + b x)^c`, the class's Briggs table.
pub(crate) fn sigma_z(distance: f64, stability: Stability, terrain: Terrain) -> f64 {
    let (a, b, c) = match (stability, terrain) {
        (Stability::A, Terrain::Rural) => (0.20, 0.0, 1.0),
        (Stability::A, Terrain::Urban) => (0.24, 0.001, 0.5),
        (Stability::B, Terrain::Rural) => (0.12, 0.0, 1.0),
        (Stability::B, Terrain::Urban) => (0.24, 0.001, 0.5),
        (Stability::C, Terrain::Rural) => (0.08, 0.0002, 0.5),
        (Stability::C, Terrain::Urban) => (0.20, 0.0, 1.0),
        (Stability::D, Terrain::Rural) => (0.06, 0.0015, 0.5),
        (Stability::D, Terrain::Urban) => (0.14, 0.0003, 0.5),
        (Stability::E, Terrain::Rural) => (0.03, 0.0003, 1.0),
        (Stability::E, Terrain::Urban) => (0.08, 0.00015, 0.5),
        (Stability::F, Terrain::Rural) => (0.016, 0.0003, 1.0),
        (Stability::F, Terrain::Urban) => (0.08, 0.00015, 0.5),
    };
    a * distance / (1.0 + b * distance).powf(c)
}

/// The two-image Gaussian concentration, kg/m3.
#[allow(clippy::too_many_arguments)] // the class's own five state terms
pub(crate) fn concentration(
    emission: f64,
    stack_height: f64,
    wind_speed: f64,
    stability: Stability,
    terrain: Terrain,
    x: f64,
    y: f64,
    z: f64,
) -> f64 {
    if x <= 0.0 || emission == 0.0 {
        return 0.0;
    }
    let sy = sigma_y(x, stability, terrain);
    let sz = sigma_z(x, stability, terrain);
    let h = stack_height;
    let term1 = emission / (2.0 * std::f64::consts::PI * wind_speed * sy * sz);
    let y_term = (-0.5 * (y * y) / (sy * sy)).exp();
    let z_term1 = (-0.5 * (z - h).powi(2) / (sz * sz)).exp();
    let z_term2 = (-0.5 * (z + h).powi(2) / (sz * sz)).exp();
    term1 * y_term * (z_term1 + z_term2)
}

/// The class's scan and bisection for the ground centreline distance to a threshold.
#[allow(clippy::too_many_arguments)]
pub(crate) fn distance_to_threshold(
    emission: f64,
    stack_height: f64,
    wind_speed: f64,
    stability: Stability,
    terrain: Terrain,
    threshold: f64,
) -> Option<f64> {
    let centreline = |x: f64| {
        concentration(
            emission,
            stack_height,
            wind_speed,
            stability,
            terrain,
            x,
            0.0,
            0.0,
        )
    };
    let mut previous_distance = SCAN_MIN_DISTANCE;
    let mut previous_concentration = centreline(previous_distance);
    let mut above = (previous_concentration >= threshold).then_some(previous_distance);
    let mut below: Option<f64> = None;
    for step in 1..=SCAN_STEPS {
        let fraction = step as f64 / SCAN_STEPS as f64;
        let distance = SCAN_MIN_DISTANCE * (SCAN_MAX_DISTANCE / SCAN_MIN_DISTANCE).powf(fraction);
        let value = centreline(distance);
        if value >= threshold {
            above = Some(distance);
        } else if previous_concentration >= threshold {
            above = Some(previous_distance);
            below = Some(distance);
        }
        previous_distance = distance;
        previous_concentration = value;
    }
    let mut lo = above?;
    let below = match below {
        Some(value) => value,
        None => return Some(SCAN_MAX_DISTANCE),
    };
    let mut hi = below;
    for _ in 0..BISECTION_STEPS {
        let mid = 0.5 * (lo + hi);
        if centreline(mid) > threshold {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < BISECTION_WIDTH {
            break;
        }
    }
    Some(0.5 * (lo + hi))
}

/// The concentration of a continuous Gaussian plume at one point, and the distance to a
/// threshold on its ground centreline.
///
/// # Errors
/// `InvalidInput` for a value the class's constructor refuses, a non-positive threshold, a
/// distance the spec's range closes, or an unknown stability or terrain word.
#[allow(clippy::too_many_arguments)]
pub fn gaussian_plume(
    emission_rate: MassRate,
    effective_stack_height: Length,
    wind_speed: Velocity,
    stability: &str,
    terrain: &str,
    distance: Length,
    crosswind: Length,
    height: Length,
    threshold: MassDensity,
) -> Result<GaussianPlumeResult> {
    let spec = &model_gen::GAUSSIAN_PLUME_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "emission_rate" => Some(emission_rate.value),
            "wind_speed" => Some(wind_speed.value),
            "distance" => Some(distance.value),
            "threshold" => Some(threshold.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let stability = Stability::parse(stability)?;
    let terrain = Terrain::parse(terrain)?;
    let emission = emission_rate.value;
    let stack = effective_stack_height.value;
    let wind = wind_speed.value;
    let x = distance.value;

    Ok(GaussianPlumeResult {
        sigma_y: meters(sigma_y(x, stability, terrain)),
        sigma_z: meters(sigma_z(x, stability, terrain)),
        concentration: kilograms_per_cubic_meter(concentration(
            emission,
            stack,
            wind,
            stability,
            terrain,
            x,
            crosswind.value,
            height.value,
        )),
        centerline_ground_concentration: kilograms_per_cubic_meter(concentration(
            emission, stack, wind, stability, terrain, x, 0.0, 0.0,
        )),
        distance_to_threshold: distance_to_threshold(
            emission,
            stack,
            wind,
            stability,
            terrain,
            threshold.value,
        )
        .map(meters),
        warnings,
    })
}
