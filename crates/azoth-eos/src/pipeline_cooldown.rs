//! `eos.pipeline_cooldown` - NeqSim's `PipelineCooldownCalculator`.
//!
//! The lumped thermal model the class marches with explicit Euler: the fluid, steel and
//! insulation are one thermal mass losing heat through the layers (or a supplied U) to the
//! ambient temperature.

use azoth_core::units::{
    HeatTransfer, Length, MassDensity, SpecificHeatCapacity, ThermalConductivity,
    ThermodynamicTemperature, Time, seconds, watts_per_square_meter_kelvin,
};
use azoth_core::{Result, apply_checks};

use crate::model_gen;
use crate::results::PipelineCooldownResult;

/// The class's steel conductivity, insulation density and insulation specific heat, none
/// of which it exposes a setter for.
const STEEL_CONDUCTIVITY: f64 = 50.0;
const INSULATION_DENSITY: f64 = 600.0;
const INSULATION_SPECIFIC_HEAT: f64 = 1700.0;
/// The class's defaults for the supplied U value and the external coefficient.
const DEFAULT_U_VALUE: f64 = 3.0;
const DEFAULT_EXTERNAL_HTC: f64 = 300.0;
/// The class ignores an insulation thinner than a millimetre in the layer resistance.
const MINIMUM_INSULATION: f64 = 0.001;

/// The explicit-Euler pipeline cooldown march.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's hard range.
///
/// # Example
/// ```
/// use azoth_core::units::{
///     joules_per_kilogram_kelvin, kelvins, kilograms_per_cubic_meter, meters, seconds,
///     watts_per_meter_kelvin, watts_per_square_meter_kelvin,
/// };
/// use azoth_eos::pipeline_cooldown::pipeline_cooldown;
///
/// let r = pipeline_cooldown(
///     meters(0.254),
///     meters(0.0127),
///     meters(0.050),
///     meters(0.006),
///     kilograms_per_cubic_meter(7850.0),
///     joules_per_kilogram_kelvin(500.0),
///     watts_per_meter_kelvin(0.17),
///     kilograms_per_cubic_meter(750.0),
///     joules_per_kilogram_kelvin(2200.0),
///     kelvins(353.15),
///     kelvins(277.15),
///     seconds(3600.0),
///     seconds(36000.0),
///     Some(watts_per_square_meter_kelvin(3.0)),
///     None,
///     None,
///     None,
/// )?;
/// assert_eq!(r.time.len(), 11);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn pipeline_cooldown(
    internal_diameter: Length,
    wall_thickness: Length,
    insulation_thickness: Length,
    coating_thickness: Length,
    steel_density: MassDensity,
    steel_specific_heat: SpecificHeatCapacity,
    insulation_conductivity: ThermalConductivity,
    fluid_density: MassDensity,
    fluid_specific_heat: SpecificHeatCapacity,
    initial_fluid_temperature: ThermodynamicTemperature,
    ambient_temperature: ThermodynamicTemperature,
    time_step: Time,
    total_time: Time,
    overall_u_value: Option<HeatTransfer>,
    use_layer_calculation: Option<bool>,
    external_htc: Option<HeatTransfer>,
    target_temperature: Option<ThermodynamicTemperature>,
) -> Result<PipelineCooldownResult> {
    let spec = &model_gen::PIPELINE_COOLDOWN_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "internal_diameter" => Some(internal_diameter.value),
            "wall_thickness" => Some(wall_thickness.value),
            "insulation_thickness" => Some(insulation_thickness.value),
            "coating_thickness" => Some(coating_thickness.value),
            "steel_density" => Some(steel_density.value),
            "steel_specific_heat" => Some(steel_specific_heat.value),
            "insulation_conductivity" => Some(insulation_conductivity.value),
            "fluid_density" => Some(fluid_density.value),
            "fluid_specific_heat" => Some(fluid_specific_heat.value),
            "initial_fluid_temperature" => Some(initial_fluid_temperature.value),
            "ambient_temperature" => Some(ambient_temperature.value),
            "time_step" => Some(time_step.value),
            "total_time" => Some(total_time.value),
            "overall_u_value" => overall_u_value.map(|value| value.value),
            "external_htc" => external_htc.map(|value| value.value),
            "target_temperature" => target_temperature.map(|value| value.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let dt = time_step.value;
    let step_count = (total_time.value / dt).ceil() as usize;

    let radius_inner = internal_diameter.value / 2.0;
    let radius_outer = radius_inner + wall_thickness.value;
    let radius_insulation = radius_outer + insulation_thickness.value;
    let radius_external = radius_insulation + coating_thickness.value;

    let steel_resistance =
        (radius_outer / radius_inner).ln() / (2.0 * std::f64::consts::PI * STEEL_CONDUCTIVITY);
    let insulation_resistance = if insulation_thickness.value > MINIMUM_INSULATION {
        (radius_insulation / radius_outer).ln()
            / (2.0 * std::f64::consts::PI * insulation_conductivity.value)
    } else {
        0.0
    };
    let external_resistance = 1.0
        / (2.0
            * std::f64::consts::PI
            * radius_external
            * external_htc.map_or(DEFAULT_EXTERNAL_HTC, |value| value.value));
    let calculated_u = 1.0
        / (2.0
            * std::f64::consts::PI
            * radius_external
            * (steel_resistance + insulation_resistance + external_resistance));

    let supplied_u = overall_u_value.map_or(DEFAULT_U_VALUE, |value| value.value);
    let use_layers = use_layer_calculation.unwrap_or(false);
    let effective_u = if use_layers { calculated_u } else { supplied_u };
    let conductance_per_length = effective_u * 2.0 * std::f64::consts::PI * radius_external;

    let fluid_area = std::f64::consts::PI * radius_inner * radius_inner;
    let steel_area =
        std::f64::consts::PI * (radius_outer * radius_outer - radius_inner * radius_inner);
    let insulation_area = std::f64::consts::PI
        * (radius_insulation * radius_insulation - radius_outer * radius_outer);
    let heat_capacity_per_length = fluid_density.value * fluid_area * fluid_specific_heat.value
        + steel_density.value * steel_area * steel_specific_heat.value
        + INSULATION_DENSITY * insulation_area * INSULATION_SPECIFIC_HEAT;

    let mut time = Vec::with_capacity(step_count + 1);
    let mut fluid_temperature = Vec::with_capacity(step_count + 1);
    time.push(0.0);
    fluid_temperature.push(initial_fluid_temperature.value);
    for index in 0..step_count {
        let previous = fluid_temperature[index];
        let heat_loss = conductance_per_length * (previous - ambient_temperature.value);
        fluid_temperature.push(previous - heat_loss / heat_capacity_per_length * dt);
        time.push((index + 1) as f64 * dt);
    }

    let time_to_reach_target = target_temperature.and_then(|target| {
        fluid_temperature
            .iter()
            .position(|temperature| *temperature <= target.value)
            .map(|index| {
                if index == 0 {
                    0.0
                } else {
                    let t1 = time[index - 1];
                    let t2 = time[index];
                    let temperature_1 = fluid_temperature[index - 1];
                    let temperature_2 = fluid_temperature[index];
                    t1 + (t2 - t1) * (temperature_1 - target.value)
                        / (temperature_1 - temperature_2)
                }
            })
    });

    Ok(PipelineCooldownResult {
        time,
        fluid_temperature,
        time_constant: seconds(heat_capacity_per_length / conductance_per_length),
        calculated_u_value: watts_per_square_meter_kelvin(calculated_u),
        effective_u_value: watts_per_square_meter_kelvin(effective_u),
        time_to_reach_target: time_to_reach_target.map(seconds),
        warnings,
    })
}
