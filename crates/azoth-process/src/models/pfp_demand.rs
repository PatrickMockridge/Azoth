//! `process.pfp_demand` - the passive fire protection demand as a registered id.
//!
//! Spec: `specs/models/process/pfp_demand.toml`. The arithmetic is
//! `process/safety/fire/PfpDemandCalculator.java` at the pin: the lumped bare-steel time to
//! critical temperature, the strict PFP test, the conduction-limited thickness and the
//! H/J rating.

use azoth_core::units::{
    HeatFluxDensity, Length, ThermalConductivity, ThermodynamicTemperature, Time, meters,
    millimeters, seconds,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{scalar, warnings as wire_warnings};
use crate::model_gen;

/// Carbon steel density, kg/m3.
const STEEL_DENSITY: f64 = 7850.0;
/// Carbon steel specific heat, J/(kg K).
const STEEL_CP: f64 = 490.0;
/// The class's default absorptivity.
const DEFAULT_ABSORPTIVITY: f64 = 0.9;
/// The class's default initial temperature, K.
const DEFAULT_INITIAL_TEMPERATURE_K: f64 = 288.15;
/// The class's default critical temperature, K.
const DEFAULT_CRITICAL_TEMPERATURE_K: f64 = 673.15;
/// The class's default fire temperature, K.
const DEFAULT_FIRE_TEMPERATURE_K: f64 = 1173.15;
/// The class's default PFP conductivity, W/(m K).
const DEFAULT_PFP_CONDUCTIVITY: f64 = 0.2;
/// The rating's 60-minute boundary, s.
const LONG_DURATION_S: f64 = 3600.0;
/// The boundary's epsilon, s: exactly 3600 s is short.
const LONG_DURATION_EPSILON_S: f64 = 1.0e-6;
/// The floor on the fire-to-steel driving temperature, K.
const MIN_DRIVING_TEMPERATURE_K: f64 = 1.0;

/// The class's two rating families.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FireType {
    /// H-class, UL 1709.
    Pool,
    /// J-class, ISO 22899.
    Jet,
}

impl FireType {
    /// The stated type, the class's `pool` default when absent.
    fn parse(value: Option<&str>) -> Result<Self> {
        match value {
            None | Some("pool") => Ok(Self::Pool),
            Some("jet") => Ok(Self::Jet),
            Some(other) => Err(AzothError::invalid_input(
                "fire_type",
                format!("the fire type is `pool` or `jet`, and {other:?} is not"),
            )),
        }
    }
}

/// Result of `process.pfp_demand`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PfpDemandResult {
    /// The steel mass per exposed area in kg/m2; dimensionless.
    pub steel_mass_per_area: f64,
    /// The lumped time to the critical temperature.
    #[serde(serialize_with = "scalar")]
    pub bare_steel_time_to_critical: Time,
    /// The stated survival time, echoed.
    #[serde(serialize_with = "scalar")]
    pub required_survival_time: Time,
    /// Whether the bare time is shorter than the survival time.
    pub pfp_required: bool,
    /// The conduction-limited protection thickness; zero when not required.
    #[serde(serialize_with = "scalar")]
    pub required_pfp_thickness: Length,
    /// The same thickness in millimetres, the class's second getter.
    #[serde(serialize_with = "scalar")]
    pub required_pfp_thickness_mm: Length,
    /// `none`, `h60`, `h120`, `j60` or `j120`.
    pub rating: String,
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

/// The bare-steel time, the strict PFP verdict, the conduction-limited thickness and rating.
///
/// # Errors
/// `InvalidInput` for a non-positive flux, thickness or survival time, an absorptivity
/// outside `(0, 1]`, a non-positive conductivity or temperature, a critical temperature at
/// or below the initial one, an unknown fire type, or a non-finite input.
#[allow(clippy::too_many_arguments)]
pub fn pfp_demand(
    incident_heat_flux: HeatFluxDensity,
    steel_thickness: Length,
    required_survival_time: Time,
    absorptivity: Option<f64>,
    initial_temperature: Option<ThermodynamicTemperature>,
    critical_temperature: Option<ThermodynamicTemperature>,
    fire_temperature: Option<ThermodynamicTemperature>,
    pfp_conductivity: Option<ThermalConductivity>,
    fire_type: Option<&str>,
) -> Result<PfpDemandResult> {
    let spec = &model_gen::PFP_DEMAND_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "incident_heat_flux" => Some(incident_heat_flux.value),
            "steel_thickness" => Some(steel_thickness.value),
            "required_survival_time" => Some(required_survival_time.value),
            "absorptivity" => absorptivity,
            "initial_temperature" => initial_temperature.map(|value| value.value),
            "critical_temperature" => critical_temperature.map(|value| value.value),
            "fire_temperature" => fire_temperature.map(|value| value.value),
            "pfp_conductivity" => pfp_conductivity.map(|value| value.value),
            _ => None,
        },
        &mut warnings,
    )?;

    require_finite("incident_heat_flux", incident_heat_flux.value)?;
    require_finite("steel_thickness", steel_thickness.value)?;
    require_finite("required_survival_time", required_survival_time.value)?;
    let alpha = absorptivity.unwrap_or(DEFAULT_ABSORPTIVITY);
    require_finite("absorptivity", alpha)?;
    let initial = initial_temperature
        .map(|value| value.value)
        .unwrap_or(DEFAULT_INITIAL_TEMPERATURE_K);
    require_finite("initial_temperature", initial)?;
    let critical = critical_temperature
        .map(|value| value.value)
        .unwrap_or(DEFAULT_CRITICAL_TEMPERATURE_K);
    require_finite("critical_temperature", critical)?;
    let fire = fire_temperature
        .map(|value| value.value)
        .unwrap_or(DEFAULT_FIRE_TEMPERATURE_K);
    require_finite("fire_temperature", fire)?;
    let conductivity = pfp_conductivity
        .map(|value| value.value)
        .unwrap_or(DEFAULT_PFP_CONDUCTIVITY);
    require_finite("pfp_conductivity", conductivity)?;
    let rating_family = FireType::parse(fire_type)?;

    if critical <= initial {
        return Err(AzothError::invalid_input(
            "critical_temperature",
            "the critical temperature must exceed the initial one; at or below it the class \
             publishes a zero or negative bare-steel time with an infinite or negative thickness",
        ));
    }

    // The class's operation order, to the letter.
    let mass_per_area = STEEL_DENSITY * steel_thickness.value;
    let delta_t = critical - initial;
    let q_abs = alpha * incident_heat_flux.value;
    let bare = mass_per_area * STEEL_CP * delta_t / q_abs;
    let survival = required_survival_time.value;
    let pfp_required = bare < survival;
    let thickness = if pfp_required {
        let q_allow = mass_per_area * STEEL_CP * delta_t / survival;
        let steel_average = 0.5 * (initial + critical);
        let driving = (fire - steel_average).max(MIN_DRIVING_TEMPERATURE_K);
        conductivity * driving / q_allow
    } else {
        0.0
    };
    let long = survival > LONG_DURATION_S + LONG_DURATION_EPSILON_S;
    let rating = match (pfp_required, long, rating_family) {
        (false, _, _) => "none",
        (true, true, FireType::Pool) => "h120",
        (true, false, FireType::Pool) => "h60",
        (true, true, FireType::Jet) => "j120",
        (true, false, FireType::Jet) => "j60",
    };

    Ok(PfpDemandResult {
        steel_mass_per_area: mass_per_area,
        bare_steel_time_to_critical: seconds(bare),
        required_survival_time: seconds(survival),
        pfp_required,
        required_pfp_thickness: meters(thickness),
        required_pfp_thickness_mm: millimeters(thickness * 1000.0),
        rating: rating.to_string(),
        warnings,
    })
}
