//! `process.cooler` - `Heater.run` reached through `Cooler`, as a registered id.
//!
//! Spec: `specs/models/process/cooler.toml`. The arithmetic is
//! [`crate::kernels::cooler`], which is [`crate::kernels::heater`]: `Cooler` overrides
//! `runTransient` and some getters and not `run`, and the probe's two captures are
//! byte-identical.

use azoth_core::units::{
    MassRate, MolarEnergy, Power, Pressure, ThermodynamicTemperature, Time, joules_per_mole, watts,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{optional_scalar, scalar, warnings as wire_warnings};
use crate::kernels::cooler as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.cooler`.
///
/// **The same fields as [`crate::models::HeaterResult`], and its own type.** A result class
/// is the registry's handle on an id - `result_fields` answers by id - so one class shared
/// by two ids would be one handle for two things, which is the shape every other pair in
/// this registry refuses too.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct CoolerResult {
    /// Molar flow out, mol/s: the inlet's.
    pub outlet_n: f64,
    /// Outlet composition, the inlet's.
    pub outlet_z: Vec<f64>,
    /// Outlet pressure: `inlet_p - pressure_drop`.
    #[serde(serialize_with = "scalar")]
    pub outlet_p: Pressure,
    /// Outlet temperature, which is the stated one or the one the shifted enthalpy reaches.
    #[serde(serialize_with = "scalar")]
    pub outlet_t: ThermodynamicTemperature,
    /// Outlet molar enthalpy.
    #[serde(serialize_with = "scalar")]
    pub outlet_h: MolarEnergy,
    /// The duty moved, W - negative when heat was removed, which is this entry's usual case.
    #[serde(serialize_with = "scalar")]
    pub outlet_duty: Power,
    /// The post-step valve opening, in percent; present only on the dynamic branch.
    pub cooling_valve_opening: Option<f64>,
    /// The equilibrium the thermal lag moved toward, K; present only on the dynamic branch.
    #[serde(serialize_with = "optional_scalar")]
    pub equilibrium_outlet_temperature: Option<ThermodynamicTemperature>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

impl CoolerResult {
    /// The result of one kernel call, for the reason [`crate::models::HeaterResult::of`] gives
    /// about the warnings.
    #[must_use]
    pub fn of(outcome: &crate::kernels::heater::HeaterOutcome, warnings: Vec<Warning>) -> Self {
        Self {
            outlet_n: outcome.outlet.n,
            outlet_z: outcome.outlet.z.clone(),
            outlet_p: outcome.outlet.p,
            outlet_t: outcome.outlet.t,
            outlet_h: joules_per_mole(outcome.outlet.h.value),
            outlet_duty: watts(outcome.duty.value),
            cooling_valve_opening: None,
            equilibrium_outlet_temperature: None,
            warnings,
        }
    }
}

/// Cool a stream to a stated temperature, or by a stated duty.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the shapes disagree, if a temperature and a
/// duty are both given, if the pressure drop leaves a non-positive pressure, or if a duty is
/// asked of a stream carrying no flow, and whatever the databank or the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input
pub fn cooler(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
    outlet_temperature: Option<ThermodynamicTemperature>,
    duty: Option<Power>,
    pressure_drop: Option<Pressure>,
    time_step: Option<Time>,
    dynamic_temperature_control_enabled: Option<bool>,
    design_mass_flow: Option<MassRate>,
    design_inlet_temperature: Option<ThermodynamicTemperature>,
    design_outlet_temperature: Option<ThermodynamicTemperature>,
    cooling_medium_temperature: Option<ThermodynamicTemperature>,
    design_valve_opening: Option<f64>,
    actuator_time_constant: Option<Time>,
    thermal_time_constant: Option<Time>,
    cooling_valve_opening: Option<f64>,
    requested_valve_opening: Option<f64>,
    dynamic_outlet_temperature: Option<ThermodynamicTemperature>,
) -> Result<CoolerResult> {
    let spec = &model_gen::COOLER_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "inlet_t" => Some(inlet_t.value),
            "outlet_temperature" => outlet_temperature.map(|t| t.value),
            "time_step" => time_step.map(|t| t.value),
            "design_mass_flow" => design_mass_flow.map(|m| m.value),
            "design_inlet_temperature" => design_inlet_temperature.map(|t| t.value),
            "design_outlet_temperature" => design_outlet_temperature.map(|t| t.value),
            "cooling_medium_temperature" => cooling_medium_temperature.map(|t| t.value),
            "design_valve_opening" => design_valve_opening,
            "actuator_time_constant" => actuator_time_constant.map(|t| t.value),
            "thermal_time_constant" => thermal_time_constant.map(|t| t.value),
            "cooling_valve_opening" => cooling_valve_opening,
            "requested_valve_opening" => requested_valve_opening,
            "dynamic_outlet_temperature" => dynamic_outlet_temperature.map(|t| t.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let transient_requested = time_step.is_some();
    if !transient_requested {
        for (name, stated) in [
            (
                "dynamic_temperature_control_enabled",
                dynamic_temperature_control_enabled.is_some(),
            ),
            ("design_mass_flow", design_mass_flow.is_some()),
            (
                "design_inlet_temperature",
                design_inlet_temperature.is_some(),
            ),
            (
                "design_outlet_temperature",
                design_outlet_temperature.is_some(),
            ),
            (
                "cooling_medium_temperature",
                cooling_medium_temperature.is_some(),
            ),
            ("design_valve_opening", design_valve_opening.is_some()),
            ("actuator_time_constant", actuator_time_constant.is_some()),
            ("thermal_time_constant", thermal_time_constant.is_some()),
            ("cooling_valve_opening", cooling_valve_opening.is_some()),
            ("requested_valve_opening", requested_valve_opening.is_some()),
            (
                "dynamic_outlet_temperature",
                dynamic_outlet_temperature.is_some(),
            ),
        ] {
            if stated {
                return Err(AzothError::invalid_input(
                    name,
                    "a transient input is stated without `time_step`",
                ));
            }
        }
    }
    for (name, value) in [
        ("time_step", time_step.map(|v| v.value)),
        ("design_mass_flow", design_mass_flow.map(|v| v.value)),
        (
            "design_inlet_temperature",
            design_inlet_temperature.map(|v| v.value),
        ),
        (
            "design_outlet_temperature",
            design_outlet_temperature.map(|v| v.value),
        ),
        (
            "cooling_medium_temperature",
            cooling_medium_temperature.map(|v| v.value),
        ),
        ("design_valve_opening", design_valve_opening),
        (
            "actuator_time_constant",
            actuator_time_constant.map(|v| v.value),
        ),
        (
            "thermal_time_constant",
            thermal_time_constant.map(|v| v.value),
        ),
        ("cooling_valve_opening", cooling_valve_opening),
        ("requested_valve_opening", requested_valve_opening),
    ] {
        if let Some(value) = value {
            if !value.is_finite() {
                return Err(AzothError::invalid_input(
                    name,
                    format!("a finite value is required, and {value} is not"),
                ));
            }
        }
    }
    if let Some(value) = dynamic_outlet_temperature.map(|v| v.value) {
        // **NaN is the class's first-step sentinel**; an infinity is not a state.
        if value.is_infinite() {
            return Err(AzothError::invalid_input(
                "dynamic_outlet_temperature",
                format!("a carried temperature is finite or NaN, and {value} is neither"),
            ));
        }
    }

    let active_transient =
        transient_requested && dynamic_temperature_control_enabled.unwrap_or(false);
    if active_transient {
        for (name, stated) in [
            ("design_mass_flow", design_mass_flow.is_some()),
            (
                "design_inlet_temperature",
                design_inlet_temperature.is_some(),
            ),
            (
                "design_outlet_temperature",
                design_outlet_temperature.is_some(),
            ),
            (
                "cooling_medium_temperature",
                cooling_medium_temperature.is_some(),
            ),
        ] {
            if !stated {
                return Err(AzothError::invalid_input(
                    name,
                    "the configured branch needs the design point \
                     `configureDynamicTemperatureControl` takes",
                ));
            }
        }
        for (name, stated) in [
            ("outlet_temperature", outlet_temperature.is_some()),
            ("duty", duty.is_some()),
        ] {
            if stated {
                return Err(AzothError::invalid_input(
                    name,
                    "the dynamic branch sets the outlet itself, so a steady specification \
                     is ambiguous",
                ));
            }
        }
    }

    let feed = Stream::from_pt(
        components.to_vec(),
        inlet_z.to_vec(),
        inlet_n,
        inlet_p,
        inlet_t,
    )?;

    if active_transient {
        let design_opening = design_valve_opening.unwrap_or(50.0);
        let carried_opening = cooling_valve_opening.unwrap_or(design_opening);
        let settings = kernel::CoolerTransientSettings {
            time_step_s: time_step.map_or(0.0, |v| v.value),
            design_mass_flow_kg_per_hour: design_mass_flow.map_or(0.0, |v| v.value * 3600.0),
            design_inlet_temperature_k: design_inlet_temperature.map_or(0.0, |v| v.value),
            design_outlet_temperature_k: design_outlet_temperature.map_or(0.0, |v| v.value),
            cooling_medium_temperature_k: cooling_medium_temperature.map_or(0.0, |v| v.value),
            design_valve_opening_percent: design_opening,
            actuator_time_constant_s: actuator_time_constant.map_or(10.0, |v| v.value),
            thermal_time_constant_s: thermal_time_constant.map_or(20.0, |v| v.value),
            cooling_valve_opening_percent: carried_opening,
            requested_valve_opening_percent: requested_valve_opening,
            dynamic_outlet_temperature_k: dynamic_outlet_temperature.map(|v| v.value),
            pressure_drop_pa: pressure_drop.map(|v| v.value),
        };
        let outcome = kernel::cooler_transient(&feed, &settings)?;
        return Ok(CoolerResult {
            outlet_n: outcome.outlet.n,
            outlet_z: outcome.outlet.z.clone(),
            outlet_p: outcome.outlet.p,
            outlet_t: outcome.outlet.t,
            outlet_h: joules_per_mole(outcome.outlet.h.value),
            outlet_duty: watts(outcome.duty.value),
            cooling_valve_opening: Some(outcome.cooling_valve_opening_percent),
            equilibrium_outlet_temperature: Some(outcome.equilibrium_outlet_temperature),
            warnings,
        });
    }

    let outcome = kernel(&feed, outlet_temperature, duty, pressure_drop)?;
    Ok(CoolerResult::of(&outcome, warnings))
}
