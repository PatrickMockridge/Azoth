//! `process.heater` - the heater's kernel as a registered id.
//!
//! Spec: `specs/models/process/heater.toml`. The arithmetic is
//! [`crate::kernels::heater`]; what is here is the boundary a case and a cross-impl test
//! address.

use azoth_core::units::{
    Area, HeatTransfer, Length, Mass, MolarEnergy, Power, Pressure, ThermalConductance,
    ThermodynamicTemperature, joules_per_mole, kelvins, kilograms, meters, square_meters, watts,
    watts_per_kelvin, watts_per_square_meter_kelvin,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{optional_scalar, scalar, warnings as wire_warnings};
use crate::kernels::heater as kernel;
use crate::kernels::heater::HeaterCapacity;
use crate::kernels::heater_electrical;
use crate::kernels::heater_mechanical;
use crate::model_gen;
use crate::stream::Stream;

/// One row of `Heater`'s instrument design.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct HeaterInstrument {
    /// The tag the list assigns: `SYMBOL-(1000 + position)`.
    pub tag: String,
    /// The ISA symbol.
    pub isa_symbol: String,
    /// The service description.
    pub service: String,
    /// `AI`, `AO`, `DI` or `DO`.
    pub io_type: String,
    /// The range's unit; absent for a digital trip.
    pub range_unit: Option<String>,
    /// The range's lower bound, in `range_unit`.
    pub range_min: f64,
    /// The range's upper bound, in `range_unit`.
    pub range_max: f64,
    /// The SIL rating; zero for a non-safety instrument.
    pub sil_rating: f64,
    /// Whether the row is safety-related.
    pub safety_related: bool,
}

/// One row of the mechanical design's candidate table.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct HeaterMechanicalType {
    pub mechanical_type: String,
    #[serde(serialize_with = "scalar")]
    pub required_area: Area,
    #[serde(serialize_with = "scalar")]
    pub required_ua: ThermalConductance,
    #[serde(serialize_with = "scalar")]
    pub overall_heat_transfer_coefficient: HeatTransfer,
    #[serde(serialize_with = "scalar")]
    pub approach_temperature: ThermodynamicTemperature,
    #[serde(serialize_with = "scalar")]
    pub estimated_length: Length,
    #[serde(serialize_with = "scalar")]
    pub inner_diameter: Length,
    #[serde(serialize_with = "scalar")]
    pub outer_diameter: Length,
    #[serde(serialize_with = "scalar")]
    pub wall_thickness: Length,
    pub tube_count: f64,
    pub tube_passes: f64,
    #[serde(serialize_with = "scalar")]
    pub fin_surface_area: Area,
    pub estimated_pressure_drop: f64,
    #[serde(serialize_with = "scalar")]
    pub estimated_weight: Mass,
    #[serde(serialize_with = "scalar")]
    pub module_length: Length,
    #[serde(serialize_with = "scalar")]
    pub module_width: Length,
    #[serde(serialize_with = "scalar")]
    pub module_height: Length,
}

/// Result of `process.heater`.
///
/// **Serialisable, so that a flowsheet's run publishes it.** The wire writes each quantity
/// through [`crate::executor::json::scalar`], which reads the unit's name from the field's own
/// dimension - and the keys are this struct's field names, which are the ones `FIELDS` declares
/// and the Python dataclass carries.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct HeaterResult {
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
    /// The duty moved, W: `inlet_n * (outlet_h - inlet_h)`, in every branch.
    #[serde(serialize_with = "scalar")]
    pub outlet_duty: Power,
    /// `maxOutletTemperature`'s utilization, when the maximum limit is set.
    pub max_outlet_temperature_utilization: Option<f64>,
    /// `minOutletTemperature`'s utilization, when the minimum limit is positive.
    pub min_outlet_temperature_utilization: Option<f64>,
    /// `getMaxUtilization`, 0.0 on an empty map.
    pub max_utilization: f64,
    /// `isCapacityExceeded`.
    pub capacity_exceeded: bool,
    /// `isHardLimitExceeded`.
    pub hard_limit_exceeded: bool,
    /// `getBottleneckConstraint`'s class-side name.
    pub bottleneck_constraint: Option<String>,
    /// The instruments in the design's list.
    pub instrument_count: Option<f64>,
    /// The `AI` rows.
    pub analog_input_count: Option<f64>,
    /// The `AO` rows.
    pub analog_output_count: Option<f64>,
    /// The `DI` rows.
    pub digital_input_count: Option<f64>,
    /// The `DO` rows.
    pub digital_output_count: Option<f64>,
    /// The four counts' sum.
    pub total_io_count: Option<f64>,
    /// The safety-related rows.
    pub safety_instrument_count: Option<f64>,
    /// The design's instrument list, in the class's insertion order; empty when it did not run.
    pub instruments: Vec<HeaterInstrument>,
    /// The shaft power the electrical type implies.
    #[serde(serialize_with = "optional_scalar")]
    pub shaft_power: Option<Power>,
    /// The electrical input.
    #[serde(serialize_with = "optional_scalar")]
    pub electrical_input_power: Option<Power>,
    /// The apparent power in VA.
    pub apparent_power: Option<f64>,
    /// The reactive power in VAR.
    pub reactive_power: Option<f64>,
    /// The power factor.
    pub power_factor: Option<f64>,
    /// The selected voltage in volts.
    pub rated_voltage: Option<f64>,
    /// The stated supply frequency in Hz.
    pub frequency_hz: Option<f64>,
    /// The full-load current in amperes.
    pub full_load_current: Option<f64>,
    /// The starting current in amperes.
    pub starting_current: Option<f64>,
    /// The electrical losses.
    #[serde(serialize_with = "optional_scalar")]
    pub total_electrical_losses: Option<Power>,
    /// The selected motor's rated power.
    #[serde(serialize_with = "optional_scalar")]
    pub motor_rated_power: Option<Power>,
    /// The motor's efficiency as a fraction.
    pub motor_efficiency: Option<f64>,
    /// The motor's power factor.
    pub motor_power_factor: Option<f64>,
    /// The motor's rated speed in rpm.
    pub motor_rated_speed: Option<f64>,
    /// The locked-rotor current multiplier.
    pub motor_locked_rotor_multiplier: Option<f64>,
    /// The auxiliary load.
    #[serde(serialize_with = "optional_scalar")]
    pub total_auxiliary_power: Option<Power>,
    /// The input plus the auxiliary load.
    #[serde(serialize_with = "optional_scalar")]
    pub total_connected_power: Option<Power>,
    /// The power cable's route length.
    #[serde(serialize_with = "optional_scalar")]
    pub cable_length: Option<Length>,
    /// The selected cable cross-section, in mm2.
    pub cable_cross_section: Option<f64>,
    /// The derated cable ampacity, in amperes.
    pub cable_ampacity: Option<f64>,
    /// The cable's voltage drop, in percent.
    pub cable_voltage_drop: Option<f64>,
    /// The cable's 1-second short-circuit withstand, in kA.
    pub cable_short_circuit_withstand: Option<f64>,
    /// The switchgear's starter type.
    pub switchgear_starter: Option<String>,
    /// The switchgear's type.
    pub switchgear_type: Option<String>,
    /// The circuit-breaker rating, in amperes.
    pub switchgear_circuit_breaker: Option<f64>,
    /// The fuse rating, in amperes.
    pub switchgear_fuse: Option<f64>,
    /// The switchgear's rated current, in amperes.
    pub switchgear_rated_current: Option<f64>,
    /// The switchgear's short-circuit rating, in kA.
    pub switchgear_short_circuit: Option<f64>,
    /// The switchgear's estimated weight.
    #[serde(serialize_with = "optional_scalar")]
    pub switchgear_weight: Option<Mass>,
    /// Whether a utility supply temperature is stated.
    pub has_supply_temperature: Option<bool>,
    /// Whether a utility return temperature is stated.
    pub has_return_temperature: Option<bool>,
    /// Whether an approach temperature is stated.
    pub has_approach_temperature: Option<bool>,
    /// Whether a heat-capacity rate is stated and positive.
    pub has_heat_capacity_rate: Option<bool>,
    /// Whether an overall heat-transfer coefficient is stated and positive.
    pub has_overall_heat_transfer_coefficient: Option<bool>,
    /// The selected mechanical configuration.
    pub mechanical_selected_type: Option<String>,
    /// The selected configuration's required area.
    #[serde(serialize_with = "optional_scalar")]
    pub mechanical_selected_area: Option<Area>,
    /// The log-mean temperature difference.
    #[serde(serialize_with = "optional_scalar")]
    pub log_mean_temperature_difference: Option<ThermodynamicTemperature>,
    /// The approach temperature.
    #[serde(serialize_with = "optional_scalar")]
    pub mechanical_approach_temperature: Option<ThermodynamicTemperature>,
    /// The overall coefficient the estimate used.
    #[serde(serialize_with = "optional_scalar")]
    pub used_overall_heat_transfer_coefficient: Option<HeatTransfer>,
    /// The calculated UA.
    #[serde(serialize_with = "optional_scalar")]
    pub calculated_ua: Option<ThermalConductance>,
    /// The five candidates, in the class's enum order.
    pub mechanical_types: Vec<HeaterMechanicalType>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

impl HeaterResult {
    /// The result of one kernel call.
    ///
    /// **The warnings are the caller's, because the two callers have different ones.** A case
    /// runs the spec's input checks through [`apply_checks`]; a flowsheet's equivalent is the
    /// *checker's*, which reports through the envelope's diagnostics rather than through a result.
    /// So an empty list here means "this kernel raised none" and not "nothing checked this call".
    #[must_use]
    pub fn of(
        outcome: &kernel::HeaterOutcome,
        capacity: HeaterCapacity,
        instruments: Option<Vec<HeaterInstrument>>,
        electrical: Option<&heater_electrical::ElectricalDesign>,
        utility: Option<&heater_mechanical::UtilitySpecification>,
        mechanical: Option<&heater_mechanical::MechanicalDesign>,
        warnings: Vec<Warning>,
    ) -> Self {
        let counts = instruments.as_ref().map(|rows| {
            let count =
                |io_type: &str| rows.iter().filter(|row| row.io_type == io_type).count() as f64;
            let analog_input = count("AI");
            let analog_output = count("AO");
            let digital_input = count("DI");
            let digital_output = count("DO");
            let safety = rows.iter().filter(|row| row.safety_related).count() as f64;
            (
                rows.len() as f64,
                analog_input,
                analog_output,
                digital_input,
                digital_output,
                analog_input + analog_output + digital_input + digital_output,
                safety,
            )
        });
        Self {
            outlet_n: outcome.outlet.n,
            outlet_z: outcome.outlet.z.clone(),
            outlet_p: outcome.outlet.p,
            outlet_t: outcome.outlet.t,
            outlet_h: joules_per_mole(outcome.outlet.h.value),
            outlet_duty: watts(outcome.duty.value),
            max_outlet_temperature_utilization: capacity.max_outlet_temperature_utilization,
            min_outlet_temperature_utilization: capacity.min_outlet_temperature_utilization,
            max_utilization: capacity.max_utilization,
            capacity_exceeded: capacity.capacity_exceeded,
            hard_limit_exceeded: capacity.hard_limit_exceeded,
            bottleneck_constraint: capacity.bottleneck_constraint.map(str::to_owned),
            instrument_count: counts.map(|values| values.0),
            analog_input_count: counts.map(|values| values.1),
            analog_output_count: counts.map(|values| values.2),
            digital_input_count: counts.map(|values| values.3),
            digital_output_count: counts.map(|values| values.4),
            total_io_count: counts.map(|values| values.5),
            safety_instrument_count: counts.map(|values| values.6),
            instruments: instruments.unwrap_or_default(),
            shaft_power: electrical.map(|value| watts(value.shaft_power_w)),
            electrical_input_power: electrical.map(|value| watts(value.electrical_input_w)),
            apparent_power: electrical.map(|value| value.apparent_power_va),
            reactive_power: electrical.map(|value| value.reactive_power_var),
            power_factor: electrical.map(|value| value.power_factor),
            rated_voltage: electrical.map(|value| value.rated_voltage_v),
            frequency_hz: electrical.map(|value| value.frequency_hz),
            full_load_current: electrical.map(|value| value.full_load_current_a),
            starting_current: electrical.map(|value| value.starting_current_a),
            total_electrical_losses: electrical.map(|value| watts(value.total_electrical_losses_w)),
            motor_rated_power: electrical.map(|value| watts(value.motor_rated_power_w)),
            motor_efficiency: electrical.map(|value| value.motor_efficiency),
            motor_power_factor: electrical.map(|value| value.motor_power_factor),
            motor_rated_speed: electrical.map(|value| value.motor_rated_speed_rpm),
            motor_locked_rotor_multiplier: electrical
                .map(|value| value.motor_locked_rotor_multiplier),
            total_auxiliary_power: electrical.map(|value| watts(value.total_auxiliary_w)),
            total_connected_power: electrical.map(|value| watts(value.total_connected_w)),
            cable_length: electrical.map(|value| meters(value.cable_length_m)),
            cable_cross_section: electrical.map(|value| value.cable_cross_section_mm2),
            cable_ampacity: electrical.map(|value| value.cable_ampacity_a),
            cable_voltage_drop: electrical.map(|value| value.cable_voltage_drop_percent),
            cable_short_circuit_withstand: electrical
                .map(|value| value.cable_short_circuit_withstand_ka),
            switchgear_starter: electrical.map(|value| value.switchgear_starter.to_string()),
            switchgear_type: electrical.map(|value| value.switchgear_type.to_string()),
            switchgear_circuit_breaker: electrical.map(|value| value.switchgear_circuit_breaker_a),
            switchgear_fuse: electrical.map(|value| value.switchgear_fuse_a),
            switchgear_rated_current: electrical.map(|value| value.switchgear_rated_current_a),
            switchgear_short_circuit: electrical.map(|value| value.switchgear_short_circuit_ka),
            switchgear_weight: electrical.map(|value| kilograms(value.switchgear_weight_kg)),
            has_supply_temperature: utility
                .map(heater_mechanical::UtilitySpecification::has_supply_temperature),
            has_return_temperature: utility
                .map(heater_mechanical::UtilitySpecification::has_return_temperature),
            has_approach_temperature: utility
                .map(heater_mechanical::UtilitySpecification::has_approach_temperature),
            has_heat_capacity_rate: utility
                .map(heater_mechanical::UtilitySpecification::has_heat_capacity_rate),
            has_overall_heat_transfer_coefficient: utility.map(
                heater_mechanical::UtilitySpecification::has_overall_heat_transfer_coefficient,
            ),
            mechanical_selected_type: mechanical
                .map(|value| value.selected.kind.spelling().to_string()),
            mechanical_selected_area: mechanical
                .map(|value| square_meters(value.selected.required_area_m2)),
            log_mean_temperature_difference: mechanical
                .map(|value| kelvins(value.log_mean_temperature_difference)),
            mechanical_approach_temperature: mechanical
                .map(|value| kelvins(value.approach_temperature)),
            used_overall_heat_transfer_coefficient: mechanical.map(|value| {
                watts_per_square_meter_kelvin(value.used_overall_heat_transfer_coefficient)
            }),
            calculated_ua: mechanical.map(|value| watts_per_kelvin(value.calculated_ua)),
            mechanical_types: mechanical
                .map(|value| {
                    value
                        .types
                        .iter()
                        .map(|row| HeaterMechanicalType {
                            mechanical_type: row.kind.spelling().to_string(),
                            required_area: square_meters(row.required_area_m2),
                            required_ua: watts_per_kelvin(row.required_ua_w_per_k),
                            overall_heat_transfer_coefficient: watts_per_square_meter_kelvin(
                                row.overall_heat_transfer_coefficient,
                            ),
                            approach_temperature: kelvins(row.approach_temperature_k),
                            estimated_length: meters(row.estimated_length_m),
                            inner_diameter: meters(row.inner_diameter_m),
                            outer_diameter: meters(row.outer_diameter_m),
                            wall_thickness: meters(row.wall_thickness_m),
                            tube_count: row.tube_count,
                            tube_passes: row.tube_passes,
                            fin_surface_area: square_meters(row.fin_surface_area_m2),
                            estimated_pressure_drop: row.estimated_pressure_drop_bar,
                            estimated_weight: kilograms(row.estimated_weight_kg),
                            module_length: meters(row.module_length_m),
                            module_width: meters(row.module_width_m),
                            module_height: meters(row.module_height_m),
                        })
                        .collect()
                })
                .unwrap_or_default(),
            warnings,
        }
    }
}

/// Heat or cool a stream to a stated temperature, or by a stated duty.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the shapes disagree, if the pressure drop
/// leaves a non-positive pressure, or if a duty is asked of a stream carrying no flow, and
/// whatever the databank or the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are fourteen
pub fn heater(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
    outlet_temperature: Option<ThermodynamicTemperature>,
    duty: Option<Power>,
    pressure_drop: Option<Pressure>,
    max_outlet_temperature: Option<ThermodynamicTemperature>,
    min_outlet_temperature: Option<ThermodynamicTemperature>,
    instrument_design: Option<bool>,
    include_safety_instruments: Option<bool>,
    default_sil_level: Option<f64>,
    heat_exchanger_type: Option<&str>,
    electrical_design: Option<bool>,
    motor_sizing_margin: Option<f64>,
    motor_poles: Option<f64>,
    efficiency_class: Option<&str>,
    frequency_hz: Option<f64>,
    phases: Option<f64>,
    number_of_fans: Option<f64>,
    instrumentation_power: Option<Power>,
    cooling_water_pump_power: Option<Power>,
    cable_length: Option<Length>,
    utility_specification: Option<bool>,
    utility_supply_temperature: Option<ThermodynamicTemperature>,
    utility_return_temperature: Option<ThermodynamicTemperature>,
    utility_approach_temperature: Option<ThermodynamicTemperature>,
    utility_heat_capacity_rate: Option<f64>,
    utility_overall_heat_transfer_coefficient: Option<f64>,
    mechanical_design: Option<bool>,
    mechanical_type: Option<&str>,
    mechanical_selection: Option<&str>,
) -> Result<HeaterResult> {
    let spec = &model_gen::HEATER_SPEC;
    let mut warnings = Vec::new();
    // One call over every declared input, for the reason `reactions.chemical_equilibrium`
    // gives. The three optional inputs report `None` where they are omitted, which is a
    // skipped check and not a pass.
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "inlet_t" => Some(inlet_t.value),
            "outlet_temperature" => outlet_temperature.map(|t| t.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let instrument_requested = instrument_design.unwrap_or(false);
    let electrical_requested = electrical_design.unwrap_or(false);
    for (name, stated) in [
        (
            "include_safety_instruments",
            include_safety_instruments.is_some(),
        ),
        ("default_sil_level", default_sil_level.is_some()),
    ] {
        if stated && !instrument_requested {
            return Err(AzothError::invalid_input(
                name,
                "an instrument-design input is stated without `instrument_design`",
            ));
        }
    }
    if heat_exchanger_type.is_some() && !instrument_requested && !electrical_requested {
        return Err(AzothError::invalid_input(
            "heat_exchanger_type",
            "the type is shared by the two designs, and neither flag is set",
        ));
    }
    for (name, stated) in [
        ("motor_sizing_margin", motor_sizing_margin.is_some()),
        ("motor_poles", motor_poles.is_some()),
        ("efficiency_class", efficiency_class.is_some()),
        ("frequency_hz", frequency_hz.is_some()),
        ("phases", phases.is_some()),
        ("number_of_fans", number_of_fans.is_some()),
        ("instrumentation_power", instrumentation_power.is_some()),
        (
            "cooling_water_pump_power",
            cooling_water_pump_power.is_some(),
        ),
        ("cable_length", cable_length.is_some()),
    ] {
        if stated && !electrical_requested {
            return Err(AzothError::invalid_input(
                name,
                "an electrical-design input is stated without `electrical_design`",
            ));
        }
    }
    for (name, value) in [
        ("motor_sizing_margin", motor_sizing_margin),
        ("motor_poles", motor_poles),
        ("frequency_hz", frequency_hz),
        ("phases", phases),
        ("number_of_fans", number_of_fans),
        (
            "instrumentation_power",
            instrumentation_power.map(|value| value.value),
        ),
        (
            "cooling_water_pump_power",
            cooling_water_pump_power.map(|value| value.value),
        ),
        ("cable_length", cable_length.map(|value| value.value)),
    ] {
        if let Some(value) = value {
            if !value.is_finite() {
                return Err(AzothError::invalid_input(
                    name,
                    format!("{value} is not a finite number"),
                ));
            }
        }
    }
    let efficiency_class = match efficiency_class {
        Some(spelling) => heater_electrical::EfficiencyClass::parse(spelling).ok_or_else(|| {
            AzothError::invalid_input(
                "efficiency_class",
                format!("{spelling} is not one of ie1, ie2, ie3 or ie4"),
            )
        })?,
        None => heater_electrical::EfficiencyClass::Ie3,
    };
    let motor_poles = motor_poles.unwrap_or(4.0);
    let phases = phases.unwrap_or(3.0);
    let number_of_fans = number_of_fans.unwrap_or(2.0);
    for (name, value) in [
        ("motor_poles", motor_poles),
        ("phases", phases),
        ("number_of_fans", number_of_fans),
    ] {
        if value < 1.0 || value.fract() != 0.0 {
            return Err(AzothError::invalid_input(
                name,
                format!("{value} is not an integer of at least one"),
            ));
        }
    }
    let sil_level = default_sil_level.unwrap_or(2.0);
    if !sil_level.is_finite()
        || sil_level < 0.0
        || sil_level.fract() != 0.0
        || sil_level > f64::from(u32::MAX)
    {
        return Err(AzothError::invalid_input(
            "default_sil_level",
            format!("the SIL level is {sil_level}, and the class's setter takes an integer"),
        ));
    }
    let instrument_type = match heat_exchanger_type {
        Some(spelling) => {
            heater_electrical::HeatExchangerType::parse(spelling).ok_or_else(|| {
                AzothError::invalid_input(
                    "heat_exchanger_type",
                    format!(
                        "{spelling} is not one of electric_heater, air_cooler or shell_and_tube"
                    ),
                )
            })?
        }
        None => heater_electrical::HeatExchangerType::ElectricHeater,
    };

    let mechanical_requested = mechanical_design.unwrap_or(false);
    for (name, stated) in [
        ("mechanical_type", mechanical_type.is_some()),
        ("mechanical_selection", mechanical_selection.is_some()),
    ] {
        if stated && !mechanical_requested {
            return Err(AzothError::invalid_input(
                name,
                "a mechanical-design input is stated without `mechanical_design`",
            ));
        }
    }
    let manual_selection = match mechanical_type {
        Some(spelling) => Some(
            heater_mechanical::MechanicalType::parse(spelling).ok_or_else(|| {
                AzothError::invalid_input(
                    "mechanical_type",
                    format!("{spelling} is not one of the five mechanical configurations"),
                )
            })?,
        ),
        None => None,
    };
    let selection_criterion = match mechanical_selection {
        Some(spelling) => Some(
            heater_mechanical::SelectionCriterion::parse(spelling).ok_or_else(|| {
                AzothError::invalid_input(
                    "mechanical_selection",
                    format!("{spelling} is not one of min_area, min_weight or min_pressure_drop"),
                )
            })?,
        ),
        None => None,
    };
    let mechanical_request = heater_mechanical::MechanicalRequest {
        manual_selection,
        selection_criterion,
    };

    let utility_requested = utility_specification.unwrap_or(false);
    for (name, stated) in [
        (
            "utility_supply_temperature",
            utility_supply_temperature.is_some(),
        ),
        (
            "utility_return_temperature",
            utility_return_temperature.is_some(),
        ),
        (
            "utility_approach_temperature",
            utility_approach_temperature.is_some(),
        ),
        (
            "utility_heat_capacity_rate",
            utility_heat_capacity_rate.is_some(),
        ),
        (
            "utility_overall_heat_transfer_coefficient",
            utility_overall_heat_transfer_coefficient.is_some(),
        ),
    ] {
        if stated && !utility_requested {
            return Err(AzothError::invalid_input(
                name,
                "a utility-specification input is stated without `utility_specification`",
            ));
        }
    }
    for (name, value) in [
        (
            "utility_supply_temperature",
            utility_supply_temperature.map(|value| value.value),
        ),
        (
            "utility_return_temperature",
            utility_return_temperature.map(|value| value.value),
        ),
        (
            "utility_approach_temperature",
            utility_approach_temperature.map(|value| value.value),
        ),
        ("utility_heat_capacity_rate", utility_heat_capacity_rate),
        (
            "utility_overall_heat_transfer_coefficient",
            utility_overall_heat_transfer_coefficient,
        ),
    ] {
        if let Some(value) = value {
            if !value.is_finite() {
                return Err(AzothError::invalid_input(
                    name,
                    format!("{value} is not a finite number"),
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
    let outcome = kernel(&feed, outlet_temperature, duty, pressure_drop)?;
    let capacity = kernel::capacity(
        outcome.outlet.t,
        max_outlet_temperature,
        min_outlet_temperature,
    );
    let instruments = if instrument_requested {
        let rows = kernel::instruments(
            instrument_type,
            include_safety_instruments.unwrap_or(true),
            sil_level as u32,
        );
        Some(
            rows.into_iter()
                .enumerate()
                .map(|(index, row)| HeaterInstrument {
                    tag: format!("{}-{}", row.isa_symbol, 1000 + index + 1),
                    isa_symbol: row.isa_symbol,
                    service: row.service,
                    io_type: row.io_type,
                    range_unit: row.range_unit,
                    range_min: row.range_min,
                    range_max: row.range_max,
                    sil_rating: f64::from(row.sil_rating),
                    safety_related: row.safety_related,
                })
                .collect(),
        )
    } else {
        None
    };
    let electrical = if electrical_requested {
        let request = heater_electrical::ElectricalRequest {
            motor_sizing_margin: motor_sizing_margin.unwrap_or(1.10),
            motor_poles: motor_poles as u32,
            efficiency_class,
            frequency_hz: frequency_hz.unwrap_or(50.0),
            phases: phases as u32,
            number_of_fans: number_of_fans as u32,
            instrumentation_power_w: instrumentation_power.map_or(1500.0, |value| value.value),
            cooling_water_pump_power_w: cooling_water_pump_power.map_or(0.0, |value| value.value),
            cable_length_m: cable_length.map_or(50.0, |value| value.value),
        };
        Some(heater_electrical::design(
            instrument_type,
            outcome.duty.value,
            &request,
        ))
    } else {
        None
    };

    let utility_specification = heater_mechanical::UtilitySpecification {
        supply_temperature: utility_supply_temperature.map(|value| value.value),
        return_temperature: utility_return_temperature.map(|value| value.value),
        approach_temperature: utility_approach_temperature.map(|value| value.value),
        heat_capacity_rate: utility_heat_capacity_rate,
        overall_heat_transfer_coefficient: utility_overall_heat_transfer_coefficient,
    };
    let utility = if utility_requested {
        Some(utility_specification)
    } else {
        None
    };
    let mechanical = if mechanical_requested {
        Some(heater_mechanical::design(
            outcome.duty.value,
            feed.t.value,
            outcome.outlet.t.value,
            &utility_specification,
            &mechanical_request,
        ))
    } else {
        None
    };

    Ok(HeaterResult::of(
        &outcome,
        capacity,
        instruments,
        electrical.as_ref(),
        utility.as_ref(),
        mechanical.as_ref(),
        warnings,
    ))
}
