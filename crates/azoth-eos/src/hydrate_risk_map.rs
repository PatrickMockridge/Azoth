//! `eos.hydrate_risk_map` - NeqSim's `HydrateRiskMapper`.
//!
//! The hydrate formation temperature at each profile point, its subcooling against the
//! actual temperature, and the risk band the two decide.

use azoth_core::units::{TemperatureInterval, kelvin_intervals, kelvins, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::Cubic;
use crate::hydrate::{HydrateModel, hydrate_mixture_of};
use crate::hydrate_formation_temperature::hydrate_formation_temperature;
use crate::model_gen;
use crate::results::{HydrateRiskLevel, HydrateRiskMapResult, HydrateRiskPoint};

/// The class's defaults for the two subcooling thresholds.
const DEFAULT_HIGH_RISK_SUBCOOLING: f64 = 3.0;
const DEFAULT_MEDIUM_RISK_SUBCOOLING: f64 = 6.0;

/// The hydrate risk along a profile.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if the three point vectors differ in length,
///   or a component name is unknown.
///
/// # Example
/// ```
/// use azoth_eos::hydrate_risk_map::hydrate_risk_map;
///
/// let r = hydrate_risk_map(
///     &[
///         "methane".to_owned(),
///         "ethane".to_owned(),
///         "propane".to_owned(),
///         "water".to_owned(),
///     ],
///     &[0.82, 0.06, 0.03, 0.09],
///     &[0.0],
///     &[1.0e7],
///     &[300.0],
///     Some("srk"),
///     Some("pvtsim"),
///     None,
///     None,
/// )?;
/// assert_eq!(r.total_points, 1);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn hydrate_risk_map(
    components: &[String],
    z: &[f64],
    distances: &[f64],
    pressures: &[f64],
    temperatures: &[f64],
    eos: Option<&str>,
    hydrate_model: Option<&str>,
    high_risk_subcooling: Option<TemperatureInterval>,
    medium_risk_subcooling: Option<TemperatureInterval>,
) -> Result<HydrateRiskMapResult> {
    let spec = &model_gen::HYDRATE_RISK_MAP_SPEC;
    let mut warnings = Vec::new();

    apply_checks(spec.input_checks(), |_| None, &mut warnings)?;

    let count = distances.len();
    if pressures.len() != count || temperatures.len() != count {
        return Err(AzothError::invalid_input(
            "distances",
            format!(
                "the profile needs one distance, pressure and temperature per point; got \
                 {}, {} and {}",
                distances.len(),
                pressures.len(),
                temperatures.len()
            ),
        ));
    }

    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    for name in &names {
        crate::databank::entry(name, None)?;
    }
    let cubic: Cubic = eos
        .unwrap_or("srk")
        .parse()
        .map_err(|_| AzothError::invalid_input("eos", "the cubic is not one of pr, srk or rk"))?;
    let hydration_model: HydrateModel = hydrate_model
        .unwrap_or("pvtsim")
        .parse()
        .unwrap_or_default();

    let high =
        high_risk_subcooling.unwrap_or_else(|| kelvin_intervals(DEFAULT_HIGH_RISK_SUBCOOLING));
    let medium =
        medium_risk_subcooling.unwrap_or_else(|| kelvin_intervals(DEFAULT_MEDIUM_RISK_SUBCOOLING));

    let mixture = hydrate_mixture_of(&names, cubic, None, hydration_model)
        .map(|(mixture, _ideal_gas)| mixture);

    let mut points = Vec::with_capacity(count);
    let mut critical_count = 0_u32;
    let mut unknown_count = 0_u32;
    let mut minimum_subcooling: Option<f64> = None;
    let mut worst = HydrateRiskLevel::Low;

    for index in 0..count {
        let distance = distances[index];
        let pressure = pressures[index];
        let actual = temperatures[index];
        let (hydrate_temperature, subcooling, risk) = match &mixture {
            Ok(mixture) => match hydrate_formation_temperature(mixture, pascals(pressure), z) {
                Ok(result) => {
                    let subcooling = actual - result.temperature.value;
                    let risk = if subcooling < 0.0 {
                        HydrateRiskLevel::Critical
                    } else if subcooling < high.value {
                        HydrateRiskLevel::High
                    } else if subcooling < medium.value {
                        HydrateRiskLevel::Medium
                    } else {
                        HydrateRiskLevel::Low
                    };
                    (Some(result.temperature), Some(subcooling), risk)
                }
                Err(_) => (None, None, HydrateRiskLevel::Unknown),
            },
            Err(_) => (None, None, HydrateRiskLevel::Unknown),
        };

        if risk == HydrateRiskLevel::Critical {
            critical_count += 1;
        }
        if risk == HydrateRiskLevel::Unknown {
            unknown_count += 1;
        } else {
            if minimum_subcooling.is_none_or(|current| subcooling.unwrap_or(f64::NAN) < current) {
                minimum_subcooling = subcooling;
            }
            if rank(risk) < rank(worst) {
                worst = risk;
            }
        }

        points.push(HydrateRiskPoint {
            distance: azoth_core::units::meters(distance),
            pressure: pascals(pressure),
            actual_temperature: kelvins(actual),
            hydrate_temperature,
            subcooling: subcooling.map(kelvin_intervals),
            risk_level: risk.as_str().to_owned(),
        });
    }

    // A failed point can hide a hydrate region, so only a confirmed CRITICAL outranks UNKNOWN.
    let overall_risk = if unknown_count > 0 && worst != HydrateRiskLevel::Critical {
        HydrateRiskLevel::Unknown
    } else {
        worst
    };

    Ok(HydrateRiskMapResult {
        points,
        overall_risk,
        minimum_subcooling: minimum_subcooling.map(kelvin_intervals),
        critical_point_count: critical_count,
        unknown_point_count: unknown_count,
        total_points: u32::try_from(count).unwrap_or(u32::MAX),
        all_points_evaluated: unknown_count == 0,
        warnings,
    })
}

fn rank(risk: HydrateRiskLevel) -> u8 {
    match risk {
        HydrateRiskLevel::Critical => 0,
        HydrateRiskLevel::High => 1,
        HydrateRiskLevel::Medium => 2,
        HydrateRiskLevel::Low => 3,
        HydrateRiskLevel::Unknown => 4,
    }
}
