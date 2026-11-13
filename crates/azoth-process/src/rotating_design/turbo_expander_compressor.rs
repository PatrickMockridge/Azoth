//! `TurboExpanderCompressorMechanicalDesign.calcDesign`, the coupled machine's sizing.
//!
//! The design reads both ends: the expander's states size its wheel and casing, the
//! compressor's suction and discharge size its impeller and casing, and one shaft takes the
//! two powers. The class's auxiliary surfaces - shear pins, the seal-gas envelope, the
//! thrust balance and the anti-surge report - are configured fields and are not published.
//!
//! **The equipment's design field is null until `initTECMechanicalDesign()`**, and nothing
//! in `run` constructs it; the probe calls it, and this function is the object it builds.

use azoth_core::Result;
use azoth_core::units::{
    Length, Mass, Power, Pressure, ThermodynamicTemperature, Velocity, kelvins, kilograms, meters,
    meters_per_second, pascals, watts,
};

use crate::segment::phase::{molar_mass_of, system_volumetric_flow_from_feed};
use crate::stream::Stream;

/// The expander's and compressor's tip-speed caps, m/s.
const MAX_EXPANDER_TIP_SPEED: f64 = 450.0;
const MAX_COMPRESSOR_TIP_SPEED: f64 = 350.0;
/// The casing design margins.
const DESIGN_PRESSURE_MARGIN: f64 = 1.10;
const DESIGN_TEMPERATURE_MARGIN: f64 = 30.0;
/// Steel's density, kg/m3.
const STEEL_DENSITY: f64 = 7850.0;
/// The radial inflow turbine's velocity ratio.
const VELOCITY_RATIO: f64 = 0.7;
/// The compressor's stage head cap, kJ/kg.
const MAX_HEAD_PER_STAGE: f64 = 30.0;
/// The shaft's allowable shear stress, Pa.
const ALLOWABLE_SHEAR_STRESS: f64 = 50.0e6;
/// The casing's base allowable stress, MPa.
const BASE_ALLOWABLE_STRESS_MPA: f64 = 137.9;

/// The expander's type, by inlet volume flow.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExpanderType {
    RadialInflow,
    Axial,
}

impl ExpanderType {
    /// The class's spelling.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::RadialInflow => "RADIAL_INFLOW",
            Self::Axial => "AXIAL",
        }
    }
}

/// The compressor casing's type, by design pressure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompressorCasingType {
    VerticallySplit,
    HorizontallySplit,
    Barrel,
}

impl CompressorCasingType {
    /// The class's spelling.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::VerticallySplit => "VERTICALLY_SPLIT",
            Self::HorizontallySplit => "HORIZONTALLY_SPLIT",
            Self::Barrel => "BARREL",
        }
    }
}

/// `TurboExpanderCompressorMechanicalDesign.calcDesign`'s stored result, in SI.
#[derive(Debug, Clone, PartialEq)]
pub struct TurboExpanderCompressorDesign {
    /// The expander wheel's diameter.
    pub expander_wheel_diameter: Length,
    /// The expander wheel's tip speed.
    pub expander_tip_speed: Velocity,
    /// The expander type.
    pub expander_type: ExpanderType,
    /// The expander's stage count.
    pub expander_stages: u32,
    /// The expander's design inlet pressure.
    pub expander_design_inlet_pressure: Pressure,
    /// The expander's design outlet pressure.
    pub expander_design_outlet_pressure: Pressure,
    /// The expander casing's design pressure.
    pub expander_casing_pressure: Pressure,
    /// The expander casing's design temperature.
    pub expander_casing_temperature: ThermodynamicTemperature,
    /// The expander's design power.
    pub expander_power: Power,
    /// The expander's enthalpy drop.
    pub expander_enthalpy_drop: f64,
    /// The compressor impeller's diameter.
    pub compressor_impeller_diameter: Length,
    /// The compressor impeller's tip speed.
    pub compressor_tip_speed: Velocity,
    /// The compressor's stage count.
    pub compressor_stages: u32,
    /// The compressor's polytropic head.
    pub compressor_polytropic_head: f64,
    /// The compressor's head per stage.
    pub compressor_head_per_stage: f64,
    /// The compressor casing's type.
    pub compressor_casing_type: CompressorCasingType,
    /// The compressor's design suction pressure.
    pub compressor_design_suction_pressure: Pressure,
    /// The compressor's design discharge pressure.
    pub compressor_design_discharge_pressure: Pressure,
    /// The compressor casing's design pressure.
    pub compressor_casing_pressure: Pressure,
    /// The compressor casing's design temperature.
    pub compressor_casing_temperature: ThermodynamicTemperature,
    /// The compressor's design power.
    pub compressor_power: Power,
    /// The shared shaft's diameter.
    pub shaft_diameter: Length,
    /// The shared shaft's weight.
    pub shaft_weight: Mass,
    /// The operating speed.
    pub operating_speed: f64,
    /// The maximum continuous speed.
    pub max_continuous_speed: f64,
    /// The trip speed.
    pub trip_speed: f64,
    /// The first critical speed.
    pub first_critical_speed: f64,
    /// The second critical speed.
    pub second_critical_speed: f64,
    /// The bearing span.
    pub bearing_span: Length,
    /// The bearing type.
    pub bearing_type: &'static str,
    /// The expander seal type.
    pub expander_seal_type: &'static str,
    /// The compressor seal type.
    pub compressor_seal_type: &'static str,
    /// The gear ratio the class's field holds.
    pub gear_ratio: f64,
    /// The design torque.
    pub design_torque: f64,
    /// The shaft stress at the design torque.
    pub design_shaft_stress: f64,
    /// The compressor's design surge flow.
    pub design_surge_flow: f64,
    /// The compressor's design stonewall flow.
    pub design_stonewall_flow: f64,
    /// The expander's nozzle throat area.
    pub expander_nozzle_throat_area: f64,
    /// The expander wheel's weight.
    pub expander_wheel_weight: Mass,
    /// The expander casing's weight.
    pub expander_casing_weight: Mass,
    /// The compressor impellers' weight.
    pub compressor_impeller_weight: Mass,
    /// The compressor casing's weight.
    pub compressor_casing_weight: Mass,
    /// The base class's vessel weight.
    pub weight_vessel: Mass,
    /// The base class's internals weight.
    pub weight_internals: Mass,
    /// The piping weight.
    pub weight_piping: Mass,
    /// The electrical and instrument weight.
    pub weight_electro_instrument: Mass,
    /// The structural steel weight.
    pub weight_structural_steel: Mass,
    /// The total weight.
    pub weight_total: Mass,
    /// The module length.
    pub module_length: Length,
    /// The module width.
    pub module_width: Length,
    /// The module height.
    pub module_height: Length,
}

/// `calcDesign` on the coupled machine's run state.
///
/// `expander_power` and `compressor_power` are magnitudes; `compressor_head` is the class's
/// `getCompressorPolytropicHead`, kJ/kg; `gear_ratio` is the equipment's field.
///
/// # Errors
/// Whatever the flash, molar mass or the machine's own degenerate state refuses.
#[allow(clippy::too_many_arguments)] // both ends and the shared state, as the class reads them
pub fn design(
    expander_feed: &Stream,
    expander_outlet: &Stream,
    compressor_feed: &Stream,
    compressor_outlet: &Stream,
    expander_power: Power,
    compressor_power: Power,
    compressor_head: f64,
    operating_speed: f64,
    gear_ratio: f64,
) -> Result<TurboExpanderCompressorDesign> {
    let (expander_mixture, _) = expander_feed.mixture()?;
    let expander_molar_mass = molar_mass_of(&expander_mixture, &expander_feed.z)?;
    let expander_mass_flow = expander_feed.n * expander_molar_mass * 3600.0;
    let expander_volume = system_volumetric_flow_from_feed(expander_feed)? * 3600.0;
    let compressor_volume = system_volumetric_flow_from_feed(compressor_feed)? * 3600.0;
    let expander_power_kw = expander_power.value.abs() / 1000.0;
    let compressor_power_kw = compressor_power.value.abs() / 1000.0;

    let expander_inlet_pressure = expander_feed.p.value / 1.0e5;
    let expander_outlet_pressure = expander_outlet.p.value / 1.0e5;
    let expander_inlet_temperature_c = expander_feed.t.value - 273.15;
    let expander_outlet_temperature_c = expander_outlet.t.value - 273.15;
    let compressor_inlet_pressure = compressor_feed.p.value / 1.0e5;
    let compressor_outlet_pressure = compressor_outlet.p.value / 1.0e5;
    let compressor_inlet_temperature_c = compressor_feed.t.value - 273.15;
    let compressor_outlet_temperature_c = compressor_outlet.t.value - 273.15;

    let expander_casing_pressure = expander_inlet_pressure * DESIGN_PRESSURE_MARGIN;
    let expander_casing_temperature =
        expander_inlet_temperature_c.max(expander_outlet_temperature_c) + DESIGN_TEMPERATURE_MARGIN;
    let compressor_casing_pressure = compressor_outlet_pressure * DESIGN_PRESSURE_MARGIN;
    let compressor_casing_temperature = compressor_inlet_temperature_c
        .max(compressor_outlet_temperature_c)
        + DESIGN_TEMPERATURE_MARGIN;

    let expander_enthalpy_drop = if expander_mass_flow > 0.0 {
        expander_power_kw * 3600.0 / expander_mass_flow
    } else {
        0.0
    };

    let expander = size_expander(expander_volume, expander_enthalpy_drop, operating_speed);
    let compressor = size_compressor(
        compressor_volume,
        compressor_head,
        operating_speed,
        compressor_casing_pressure,
    );
    let shaft = size_shaft(
        expander_power_kw + compressor_power_kw,
        operating_speed,
        expander.wheel_diameter,
        compressor.impeller_diameter,
        compressor.stages,
    );
    let (bearing, expander_seal, compressor_seal) = select_bearings_and_seals(
        operating_speed,
        expander_casing_pressure,
        compressor_casing_pressure,
    );
    let expander_casing_wall = casing_wall(
        expander.wheel_diameter,
        expander_casing_pressure,
        expander_casing_temperature,
    );
    let compressor_casing_wall = casing_wall(
        compressor.impeller_diameter,
        compressor_casing_pressure,
        compressor_casing_temperature,
    );
    let nozzle_area = nozzle_throat_area(expander_volume, expander_enthalpy_drop);
    let design_surge_flow = compressor_volume * 0.65;
    let design_stonewall_flow = compressor_volume * 1.25;
    let weights = weights(WeightsInput {
        expander_wheel_diameter: expander.wheel_diameter,
        compressor_stages: compressor.stages,
        compressor_impeller_diameter: compressor.impeller_diameter,
        shaft_diameter: shaft.diameter,
        bearing_span: shaft.bearing_span,
        expander_casing_wall,
        compressor_casing_wall,
        expander_power_kw,
        compressor_power_kw,
        compressor_seal,
    });
    let module = module_dimensions(
        expander.wheel_diameter,
        expander_casing_wall,
        compressor.stages,
        compressor.impeller_diameter,
        compressor_casing_wall,
    );

    Ok(TurboExpanderCompressorDesign {
        expander_wheel_diameter: meters(expander.wheel_diameter / 1000.0),
        expander_tip_speed: meters_per_second(expander.tip_speed),
        expander_type: expander.kind,
        expander_stages: expander.stages,
        expander_design_inlet_pressure: pascals(expander_inlet_pressure * 1.0e5),
        expander_design_outlet_pressure: pascals(expander_outlet_pressure * 1.0e5),
        expander_casing_pressure: pascals(expander_casing_pressure * 1.0e5),
        expander_casing_temperature: kelvins(expander_casing_temperature + 273.15),
        expander_power: watts(expander_power_kw * 1000.0),
        expander_enthalpy_drop,
        compressor_impeller_diameter: meters(compressor.impeller_diameter / 1000.0),
        compressor_tip_speed: meters_per_second(compressor.tip_speed),
        compressor_stages: compressor.stages,
        compressor_polytropic_head: compressor_head,
        compressor_head_per_stage: compressor.head_per_stage,
        compressor_casing_type: compressor.casing,
        compressor_design_suction_pressure: pascals(compressor_inlet_pressure * 1.0e5),
        compressor_design_discharge_pressure: pascals(compressor_outlet_pressure * 1.0e5),
        compressor_casing_pressure: pascals(compressor_casing_pressure * 1.0e5),
        compressor_casing_temperature: kelvins(compressor_casing_temperature + 273.15),
        compressor_power: watts(compressor_power_kw * 1000.0),
        shaft_diameter: meters(shaft.diameter / 1000.0),
        shaft_weight: kilograms(weights.shaft),
        operating_speed,
        max_continuous_speed: shaft.max_continuous_speed,
        trip_speed: shaft.trip_speed,
        first_critical_speed: shaft.first_critical_speed,
        second_critical_speed: shaft.second_critical_speed,
        bearing_span: meters(shaft.bearing_span / 1000.0),
        bearing_type: bearing,
        expander_seal_type: expander_seal,
        compressor_seal_type: compressor_seal,
        gear_ratio,
        design_torque: shaft.torque,
        design_shaft_stress: shaft.stress,
        design_surge_flow,
        design_stonewall_flow,
        expander_nozzle_throat_area: nozzle_area,
        expander_wheel_weight: kilograms(weights.expander_wheel),
        expander_casing_weight: kilograms(weights.expander_casing),
        compressor_impeller_weight: kilograms(weights.compressor_impeller),
        compressor_casing_weight: kilograms(weights.compressor_casing),
        weight_vessel: kilograms(weights.vessel),
        weight_internals: kilograms(weights.internals),
        weight_piping: kilograms(weights.piping),
        weight_electro_instrument: kilograms(weights.electrical),
        weight_structural_steel: kilograms(weights.structural),
        weight_total: kilograms(weights.total),
        module_length: meters(module.0),
        module_width: meters(module.1),
        module_height: meters(module.2),
    })
}

struct Expander {
    wheel_diameter: f64,
    tip_speed: f64,
    kind: ExpanderType,
    stages: u32,
}

/// `sizeExpanderWheel`, in millimetres.
fn size_expander(
    volume_flow_m3_hr: f64,
    enthalpy_drop_kj_kg: f64,
    operating_speed: f64,
) -> Expander {
    if enthalpy_drop_kj_kg <= 0.0 || volume_flow_m3_hr <= 0.0 {
        return Expander {
            wheel_diameter: 300.0,
            tip_speed: 300.0,
            kind: ExpanderType::RadialInflow,
            stages: 1,
        };
    }
    let (kind, stages) = if volume_flow_m3_hr > 100_000.0 {
        (
            ExpanderType::Axial,
            (enthalpy_drop_kj_kg / 50.0).ceil().max(1.0) as u32,
        )
    } else {
        (ExpanderType::RadialInflow, 1)
    };
    let spouting = (2.0 * enthalpy_drop_kj_kg * 1000.0).sqrt();
    let mut tip_speed = (VELOCITY_RATIO * spouting).min(MAX_EXPANDER_TIP_SPEED);
    let volume_m3_s = volume_flow_m3_hr / 3600.0;
    let enthalpy_per_stage = enthalpy_drop_kj_kg / f64::from(stages);
    let mut calc_speed = 70.0 * enthalpy_per_stage.powf(0.75) / (volume_m3_s * 60.0).sqrt();
    calc_speed = calc_speed.clamp(3000.0, 100_000.0);
    let speed = if operating_speed > 0.0 {
        operating_speed
    } else {
        calc_speed
    };
    let mut wheel_diameter = tip_speed * 60.0 / (std::f64::consts::PI * speed) * 1000.0;
    wheel_diameter = wheel_diameter.clamp(100.0, 1500.0);
    tip_speed = std::f64::consts::PI * (wheel_diameter / 1000.0) * speed / 60.0;
    Expander {
        wheel_diameter,
        tip_speed,
        kind,
        stages,
    }
}

struct Compressor {
    impeller_diameter: f64,
    tip_speed: f64,
    stages: u32,
    head_per_stage: f64,
    casing: CompressorCasingType,
}

/// `sizeCompressorImpeller`, in millimetres.
fn size_compressor(
    volume_flow_m3_hr: f64,
    polytropic_head: f64,
    operating_speed: f64,
    casing_design_pressure: f64,
) -> Compressor {
    let casing = if casing_design_pressure > 100.0 {
        CompressorCasingType::Barrel
    } else if casing_design_pressure > 40.0 {
        CompressorCasingType::HorizontallySplit
    } else {
        CompressorCasingType::VerticallySplit
    };
    if polytropic_head <= 0.0 || volume_flow_m3_hr <= 0.0 {
        return Compressor {
            impeller_diameter: 300.0,
            tip_speed: 0.0,
            stages: 1,
            head_per_stage: 0.0,
            casing,
        };
    }
    let stages = ((polytropic_head / MAX_HEAD_PER_STAGE).ceil() as u32).max(1);
    let head_per_stage = polytropic_head / f64::from(stages);
    let mut tip_speed = (head_per_stage * 1000.0 / 0.50).sqrt();
    tip_speed = tip_speed.min(MAX_COMPRESSOR_TIP_SPEED);
    let speed = if operating_speed > 0.0 {
        operating_speed
    } else {
        10_000.0
    };
    let mut impeller_diameter = tip_speed * 60.0 / (std::f64::consts::PI * speed) * 1000.0;
    impeller_diameter = impeller_diameter.clamp(100.0, 1500.0);
    tip_speed = std::f64::consts::PI * (impeller_diameter / 1000.0) * speed / 60.0;
    Compressor {
        impeller_diameter,
        tip_speed,
        stages,
        head_per_stage,
        casing,
    }
}

struct Shaft {
    diameter: f64,
    bearing_span: f64,
    torque: f64,
    stress: f64,
    max_continuous_speed: f64,
    trip_speed: f64,
    first_critical_speed: f64,
    second_critical_speed: f64,
}

/// `sizeShaft` and `calculateRotorDynamics`, in millimetres, rpm and Nm.
fn size_shaft(
    total_power_kw: f64,
    speed_rpm: f64,
    expander_wheel_diameter: f64,
    compressor_impeller_diameter: f64,
    compressor_stages: u32,
) -> Shaft {
    let torque = if speed_rpm > 0.0 {
        total_power_kw * 1000.0 * 60.0 / (2.0 * std::f64::consts::PI * speed_rpm)
    } else {
        0.0
    };
    let diameter = if speed_rpm <= 0.0 || total_power_kw <= 0.0 {
        50.0
    } else {
        let d = (16.0 * torque / (std::f64::consts::PI * ALLOWABLE_SHEAR_STRESS)).powf(1.0 / 3.0)
            * 1000.0
            * 1.5;
        ((d / 5.0).ceil() * 5.0).max(30.0)
    };
    let d_m = diameter / 1000.0;
    let polar = std::f64::consts::PI * d_m.powi(4) / 32.0;
    let stress = if polar > 0.0 {
        (torque * d_m / 2.0) / polar / 1.0e6
    } else {
        0.0
    };
    let expander_end = expander_wheel_diameter * 1.3;
    let compressor_end = f64::from(compressor_stages) * compressor_impeller_diameter * 0.8
        + compressor_impeller_diameter;
    let bearing_span = expander_end + compressor_end + 200.0;
    let max_continuous_speed = speed_rpm * 1.05;
    let trip_speed = max_continuous_speed * 1.05;
    let l_span = bearing_span / 1000.0;
    let mass_per_length = STEEL_DENSITY * std::f64::consts::PI * d_m * d_m / 4.0;
    let first_critical_speed = if l_span > 0.0 && mass_per_length > 0.0 {
        let stiffness = 48.0 * 2.1e11 * std::f64::consts::PI * d_m.powi(4) / 64.0;
        (946.0 / (2.0 * std::f64::consts::PI))
            * (stiffness / (mass_per_length * l_span.powi(3))).sqrt()
            * 60.0
    } else {
        0.0
    };
    Shaft {
        diameter,
        bearing_span,
        torque,
        stress,
        max_continuous_speed,
        trip_speed,
        first_critical_speed,
        second_critical_speed: first_critical_speed * 4.0,
    }
}

/// `selectBearingsAndSeals`, by speed and the two casing pressures.
fn select_bearings_and_seals(
    speed_rpm: f64,
    expander_casing_pressure: f64,
    compressor_casing_pressure: f64,
) -> (&'static str, &'static str, &'static str) {
    let bearing = if speed_rpm > 50_000.0 {
        "Active Magnetic"
    } else if speed_rpm > 20_000.0 {
        "Tilting Pad"
    } else {
        "Sleeve"
    };
    let expander_seal = if expander_casing_pressure > 50.0 {
        "Dry Gas"
    } else {
        "Labyrinth"
    };
    let compressor_seal = if compressor_casing_pressure > 20.0 {
        "Dry Gas"
    } else {
        "Labyrinth"
    };
    (bearing, expander_seal, compressor_seal)
}

/// `getTemperatureDeratedStress`, MPa.
fn derated_stress(temperature_c: f64) -> f64 {
    if temperature_c > 200.0 {
        BASE_ALLOWABLE_STRESS_MPA * 0.85
    } else if temperature_c < -50.0 {
        BASE_ALLOWABLE_STRESS_MPA * 0.90
    } else {
        BASE_ALLOWABLE_STRESS_MPA
    }
}

/// `calculateExpanderCasing`/`calculateCompressorCasing`, in millimetres.
fn casing_wall(diameter_mm: f64, pressure_bara: f64, temperature_c: f64) -> f64 {
    let inner = diameter_mm * 1.3;
    let pressure_mpa = pressure_bara * 0.1;
    let allowable = derated_stress(temperature_c);
    let joint = 0.85;
    let thickness = (pressure_mpa * inner) / (2.0 * allowable * joint - 0.6 * pressure_mpa) + 3.0;
    thickness.max(10.0)
}

/// `estimateExpanderNozzleThroatArea`, m2.
fn nozzle_throat_area(volume_flow_m3_hr: f64, enthalpy_drop_kj_kg: f64) -> f64 {
    if volume_flow_m3_hr <= 0.0 || enthalpy_drop_kj_kg <= 0.0 {
        return 0.01;
    }
    let velocity = 0.95 * (2.0 * enthalpy_drop_kj_kg * 1000.0).sqrt();
    let mut area = (volume_flow_m3_hr / 3600.0) / velocity.max(1.0);
    area *= 1.2;
    area
}

struct WeightsInput {
    expander_wheel_diameter: f64,
    compressor_stages: u32,
    compressor_impeller_diameter: f64,
    shaft_diameter: f64,
    bearing_span: f64,
    expander_casing_wall: f64,
    compressor_casing_wall: f64,
    expander_power_kw: f64,
    compressor_power_kw: f64,
    compressor_seal: &'static str,
}

struct Weights {
    expander_wheel: f64,
    expander_casing: f64,
    compressor_impeller: f64,
    compressor_casing: f64,
    shaft: f64,
    vessel: f64,
    internals: f64,
    piping: f64,
    electrical: f64,
    structural: f64,
    total: f64,
}

/// `calculateWeights`, in kilograms.
fn weights(input: WeightsInput) -> Weights {
    let wheel_volume = std::f64::consts::PI
        * (input.expander_wheel_diameter / 2000.0).powi(2)
        * (input.expander_wheel_diameter / 3000.0);
    let expander_wheel = wheel_volume * STEEL_DENSITY * 0.7;

    let impeller_volume = f64::from(input.compressor_stages)
        * std::f64::consts::PI
        * (input.compressor_impeller_diameter / 2000.0).powi(2)
        * (input.compressor_impeller_diameter / 3000.0);
    let compressor_impeller = impeller_volume * STEEL_DENSITY * 0.7;

    let shaft_length = input.bearing_span * 1.5;
    let shaft_volume =
        std::f64::consts::PI * (input.shaft_diameter / 2000.0).powi(2) * (shaft_length / 1000.0);
    let shaft = shaft_volume * STEEL_DENSITY;

    let exp_out = input.expander_wheel_diameter * 1.3 + 2.0 * input.expander_casing_wall;
    let exp_in = input.expander_wheel_diameter * 1.3;
    let exp_length = input.expander_wheel_diameter * 2.0;
    let exp_volume = std::f64::consts::PI / 4.0
        * ((exp_out / 1000.0).powi(2) - (exp_in / 1000.0).powi(2))
        * (exp_length / 1000.0);
    let expander_casing = exp_volume * STEEL_DENSITY * 1.5;

    let comp_out = input.compressor_impeller_diameter * 1.3 + 2.0 * input.compressor_casing_wall;
    let comp_in = input.compressor_impeller_diameter * 1.3;
    let comp_length = f64::from(input.compressor_stages) * input.compressor_impeller_diameter * 1.2;
    let comp_volume = std::f64::consts::PI / 4.0
        * ((comp_out / 1000.0).powi(2) - (comp_in / 1000.0).powi(2))
        * (comp_length / 1000.0);
    let compressor_casing = comp_volume * STEEL_DENSITY * 1.5;

    let mut internals = expander_wheel + compressor_impeller + shaft;
    internals *= 1.3;

    let lube_system = 200.0 + (input.expander_power_kw + input.compressor_power_kw) * 0.15;
    let seal_system = if input.compressor_seal == "Dry Gas" {
        500.0
    } else {
        100.0
    };
    let baseplate = (expander_casing + compressor_casing + internals) * 0.2;
    let piping_weight = (expander_casing + compressor_casing) * 0.15;
    let electrical = 50.0;
    let structural = (expander_casing + compressor_casing + internals) * 0.1;
    let vessel = expander_casing + compressor_casing;
    let piping = piping_weight + lube_system + seal_system;
    let structural_steel = structural + baseplate;
    let total = vessel + internals + piping + electrical + structural_steel;

    Weights {
        expander_wheel,
        expander_casing,
        compressor_impeller,
        compressor_casing,
        shaft,
        vessel,
        internals,
        piping,
        electrical,
        structural: structural_steel,
        total,
    }
}

/// `calculateModuleDimensions`, in metres.
fn module_dimensions(
    expander_wheel_diameter: f64,
    expander_casing_wall: f64,
    compressor_stages: u32,
    compressor_impeller_diameter: f64,
    compressor_casing_wall: f64,
) -> (f64, f64, f64) {
    let expander_len = expander_wheel_diameter / 1000.0 * 2.5;
    let compressor_len =
        f64::from(compressor_stages) * compressor_impeller_diameter / 1000.0 * 1.2 + 1.0;
    let module_length = expander_len + compressor_len + 3.0;
    let max_casing_od = (expander_wheel_diameter * 1.3 + 2.0 * expander_casing_wall)
        .max(compressor_impeller_diameter * 1.3 + 2.0 * compressor_casing_wall);
    let module_width = 2.5_f64.max(max_casing_od / 1000.0 * 2.0 + 1.5);
    let module_height = 2.0_f64.max(max_casing_od / 1000.0 + 1.5);
    (module_length, module_width, module_height)
}
