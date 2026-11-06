//! `process.heater` - the heater's kernel as a registered id.
//!
//! Spec: `specs/models/process/heater.toml`. The arithmetic is
//! [`crate::kernels::heater`]; what is here is the boundary a case and a cross-impl test
//! address.

use azoth_core::units::{
    MolarEnergy, Power, Pressure, ThermodynamicTemperature, joules_per_mole, watts,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{optional_scalar, scalar, warnings as wire_warnings};
use crate::kernels::heater as kernel;
use crate::kernels::heater::HeaterCapacity;
use crate::kernels::heater_electrical;
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
            total_electrical_losses: electrical
                .map(|value| watts(value.total_electrical_losses_w)),
            motor_rated_power: electrical.map(|value| watts(value.motor_rated_power_w)),
            motor_efficiency: electrical.map(|value| value.motor_efficiency),
            motor_power_factor: electrical.map(|value| value.motor_power_factor),
            motor_rated_speed: electrical.map(|value| value.motor_rated_speed_rpm),
            motor_locked_rotor_multiplier: electrical
                .map(|value| value.motor_locked_rotor_multiplier),
            total_auxiliary_power: electrical.map(|value| watts(value.total_auxiliary_w)),
            total_connected_power: electrical.map(|value| watts(value.total_connected_w)),
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
        Some(spelling) => {
            heater_electrical::EfficiencyClass::parse(spelling).ok_or_else(|| {
                AzothError::invalid_input(
                    "efficiency_class",
                    format!("{spelling} is not one of ie1, ie2, ie3 or ie4"),
                )
            })?
        }
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
        Some(spelling) => kernel::InstrumentType::parse(spelling).ok_or_else(|| {
            AzothError::invalid_input(
                "heat_exchanger_type",
                format!("{spelling} is not one of electric_heater, air_cooler or shell_and_tube"),
            )
        })?,
        None => kernel::InstrumentType::ElectricHeater,
    };

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
            instrumentation_power_w: instrumentation_power
                .map_or(1500.0, |value| value.value),
            cooling_water_pump_power_w: cooling_water_pump_power
                .map_or(0.0, |value| value.value),
        };
        Some(heater_electrical::design(
            instrument_type,
            outcome.duty.value,
            &request,
        ))
    } else {
        None
    };

    Ok(HeaterResult::of(
        &outcome,
        capacity,
        instruments,
        electrical.as_ref(),
        warnings,
    ))
}
