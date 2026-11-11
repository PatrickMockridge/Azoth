//! `eos.surf_cooldown` - NeqSim's `SurfCooldownAnalyzer`.
//!
//! The hydrate equilibrium temperature at the operating pressure, the fluid's density and
//! specific heat at the operating point, and the lumped cooldown that follows: the no-touch
//! time is when the fluid reaches the hydrate temperature plus the margin.

use azoth_core::units::{
    HeatTransfer, Length, Pressure, TemperatureInterval, ThermalConductivity,
    ThermodynamicTemperature, Time, kelvin_intervals, kelvins, seconds,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::Cubic;
use crate::databank;
use crate::hydrate::{HydrateModel, hydrate_mixture_of};
use crate::hydrate_formation_temperature::hydrate_formation_temperature;
use crate::model_gen;
use crate::molar_enthalpy_entropy::molar_enthalpy_entropy;
use crate::pipe_beggs_and_brills::{corrected_molar_volume, molar_mass};
use crate::pipeline_cooldown::pipeline_cooldown;
use crate::results::{CooldownVerdict, SurfCooldownResult};
use crate::tp_multiflash::tp_multiflash;

/// The class's defaults: a 72-hour horizon at a five-minute step, a 3 K hydrate margin and
/// the steel properties the cooldown carries.
const DEFAULT_TOTAL_TIME: f64 = 259_200.0;
const DEFAULT_TIME_STEP: f64 = 300.0;
const DEFAULT_HYDRATE_MARGIN: f64 = 3.0;
const DEFAULT_STEEL_DENSITY: f64 = 7850.0;
const DEFAULT_STEEL_SPECIFIC_HEAT: f64 = 500.0;
/// The class's water-presence floor.
const WATER_PRESENCE_TOLERANCE: f64 = 1.0e-12;

/// The surface-cooldown no-touch-time analysis.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's hard range.
/// * [`azoth_core::AzothError::InvalidInput`] if the fluid carries no water and
///   `require_water` is set, or the flash fails.
///
/// # Example
/// ```
/// use azoth_core::units::{
///     kelvins, meters, pascals, seconds, watts_per_meter_kelvin,
///     watts_per_square_meter_kelvin,
/// };
/// use azoth_eos::surf_cooldown::surf_cooldown;
///
/// let r = surf_cooldown(
///     &[
///         "methane".to_owned(),
///         "ethane".to_owned(),
///         "propane".to_owned(),
///         "water".to_owned(),
///     ],
///     &[0.82, 0.06, 0.03, 0.09],
///     pascals(1.2e7),
///     kelvins(338.15),
///     meters(0.254),
///     meters(0.0159),
///     meters(0.060),
///     meters(0.006),
///     watts_per_meter_kelvin(0.22),
///     kelvins(277.15),
///     None,
///     Some(watts_per_square_meter_kelvin(2.5)),
///     None,
///     None,
///     None,
///     None,
///     None,
///     None,
///     Some("srk"),
///     Some("pvtsim"),
/// )?;
/// assert!(r.water_present);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn surf_cooldown(
    components: &[String],
    z: &[f64],
    operating_pressure: Pressure,
    flash_temperature: ThermodynamicTemperature,
    internal_diameter: Length,
    wall_thickness: Length,
    insulation_thickness: Length,
    coating_thickness: Length,
    insulation_conductivity: ThermalConductivity,
    seabed_temperature: ThermodynamicTemperature,
    operating_temperature: Option<ThermodynamicTemperature>,
    overall_u_value: Option<HeatTransfer>,
    external_htc: Option<HeatTransfer>,
    hydrate_margin: Option<TemperatureInterval>,
    required_no_touch_time: Option<Time>,
    total_time: Option<Time>,
    time_step: Option<Time>,
    require_water: Option<bool>,
    eos: Option<&str>,
    hydrate_model: Option<&str>,
) -> Result<SurfCooldownResult> {
    let spec = &model_gen::SURF_COOLDOWN_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "operating_pressure" => Some(operating_pressure.value),
            "flash_temperature" => Some(flash_temperature.value),
            "seabed_temperature" => Some(seabed_temperature.value),
            "operating_temperature" => operating_temperature.map(|value| value.value),
            "overall_u_value" => overall_u_value.map(|value| value.value),
            "external_htc" => external_htc.map(|value| value.value),
            "required_no_touch_time" => required_no_touch_time.map(|value| value.value),
            "total_time" => total_time.map(|value| value.value),
            "time_step" => time_step.map(|value| value.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    if z.len() != names.len() {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "the composition has {} entries for {} components",
                z.len(),
                names.len()
            ),
        ));
    }
    for name in &names {
        databank::entry(name, None)?;
    }
    let cubic: Cubic = eos
        .unwrap_or("srk")
        .parse()
        .map_err(|_| AzothError::invalid_input("eos", "the cubic is not one of pr, srk or rk"))?;
    let hydration_model: HydrateModel = hydrate_model
        .unwrap_or("pvtsim")
        .parse()
        .unwrap_or_default();

    let water_index = names
        .iter()
        .position(|name| name.eq_ignore_ascii_case("water"));
    let water_mole_fraction = water_index.map_or(0.0, |index| z[index]);
    let water_present = water_mole_fraction > WATER_PRESENCE_TOLERANCE;
    if !water_present && require_water.unwrap_or(false) {
        return Err(AzothError::invalid_input(
            "components",
            "the fluid carries no water, so no hydrate equilibrium temperature can be \
             computed and the no-touch time is not meaningful; add water to the mixture or \
             leave `require_water` unset",
        ));
    }

    let (mixture, ideal_gas) = databank::mixture_of(&names, cubic, None)?;
    let srk_shift = matches!(cubic, Cubic::Srk);
    let flash = tp_multiflash(&mixture, flash_temperature, operating_pressure, z)?;

    let components_resolved = mixture.components();
    let mut total_mass = 0.0;
    let mut total_volume = 0.0;
    let mut total_heat_capacity = 0.0;
    for (index, fraction) in flash.phase_fractions.iter().enumerate() {
        if *fraction <= 0.0 {
            continue;
        }
        let composition = &flash.x[index];
        let z_factor = flash.z_factor[index];
        let volume = corrected_molar_volume(
            z_factor,
            flash_temperature.value,
            operating_pressure.value,
            composition,
            components_resolved,
            srk_shift,
        );
        let mass = molar_mass(composition, components_resolved, components)?;
        let cp = molar_enthalpy_entropy(
            &mixture,
            &ideal_gas,
            flash_temperature,
            operating_pressure,
            composition,
            z_factor,
        )?
        .cp
        .value;
        total_mass += fraction * mass;
        total_volume += fraction * volume;
        total_heat_capacity += fraction * cp;
    }
    if total_mass <= 0.0 || total_volume <= 0.0 {
        return Err(AzothError::invalid_input(
            "components",
            "the flash left no phase with a positive mass and volume",
        ));
    }
    let fluid_density = total_mass / total_volume;
    let fluid_specific_heat = total_heat_capacity / total_mass;

    let hydrate_temperature = hydrate_mixture_of(&names, cubic, None, hydration_model)
        .ok()
        .and_then(|(hydrated, _ideal_gas)| {
            hydrate_formation_temperature(&hydrated, operating_pressure, z)
                .ok()
                .map(|result| result.temperature)
        });
    let margin = hydrate_margin.unwrap_or_else(|| kelvin_intervals(DEFAULT_HYDRATE_MARGIN));

    let initial_temperature = operating_temperature.unwrap_or(flash_temperature);
    let target = hydrate_temperature.map(|temperature| kelvins(temperature.value + margin.value));
    let cooldown = pipeline_cooldown(
        internal_diameter,
        wall_thickness,
        insulation_thickness,
        coating_thickness,
        azoth_core::units::kilograms_per_cubic_meter(DEFAULT_STEEL_DENSITY),
        azoth_core::units::joules_per_kilogram_kelvin(DEFAULT_STEEL_SPECIFIC_HEAT),
        insulation_conductivity,
        azoth_core::units::kilograms_per_cubic_meter(fluid_density),
        azoth_core::units::joules_per_kilogram_kelvin(fluid_specific_heat),
        initial_temperature,
        seabed_temperature,
        time_step.unwrap_or_else(|| seconds(DEFAULT_TIME_STEP)),
        total_time.unwrap_or_else(|| seconds(DEFAULT_TOTAL_TIME)),
        overall_u_value,
        Some(overall_u_value.is_none()),
        external_htc,
        target,
    )?;

    let classify = |duration: f64| match required_no_touch_time {
        Some(required) if required.value > 0.0 => {
            if duration >= required.value {
                CooldownVerdict::Ok
            } else if duration >= 0.75 * required.value {
                CooldownVerdict::Marginal
            } else {
                CooldownVerdict::Critical
            }
        }
        _ => {
            if duration >= 12.0 * 3600.0 {
                CooldownVerdict::Ok
            } else if duration >= 6.0 * 3600.0 {
                CooldownVerdict::Marginal
            } else {
                CooldownVerdict::Critical
            }
        }
    };

    // The class reports an absent or seabed-warm hydrate temperature as +infinity with its
    // no-risk verdict, and a target the march never reaches as +infinity with a clean one.
    let (no_touch_time, verdict) = match hydrate_temperature {
        None => (None, CooldownVerdict::NoHydrateRisk),
        Some(temperature) if temperature.value <= seabed_temperature.value => {
            (None, CooldownVerdict::NoHydrateRisk)
        }
        Some(temperature) => {
            let target_temperature = temperature.value + margin.value;
            if initial_temperature.value <= target_temperature {
                (Some(0.0), classify(0.0))
            } else {
                match cooldown.time_to_reach_target {
                    Some(time) => (Some(time.value), classify(time.value)),
                    None => (None, CooldownVerdict::Ok),
                }
            }
        }
    };

    Ok(SurfCooldownResult {
        verdict,
        water_present,
        water_mole_fraction,
        no_touch_time: no_touch_time.map(seconds),
        initial_fluid_temperature: initial_temperature,
        hydrate_equilibrium_temperature: hydrate_temperature,
        hydrate_margin: margin,
        seabed_temperature,
        overall_u_value: cooldown.effective_u_value,
        time_constant: cooldown.time_constant,
        fluid_density: azoth_core::units::kilograms_per_cubic_meter(fluid_density),
        fluid_specific_heat: azoth_core::units::joules_per_kilogram_kelvin(fluid_specific_heat),
        internal_diameter,
        insulation_thickness,
        warnings,
    })
}
