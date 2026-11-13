//! `ExpanderMechanicalDesign.calcDesign`, the API 617 surface.
//!
//! The class exposes no setter: every quantity is the expander's run state and the
//! class's constants. The arithmetic is transcribed in the order `calcDesign` runs it,
//! and the struct converts the class's millimetres, bar and degrees Celsius at the
//! boundary.

use azoth_core::Result;
use azoth_core::units::{
    Length, Mass, Power, Pressure, ThermodynamicTemperature, Velocity, kelvins, kilograms, meters,
    meters_per_second, pascals,
};

use crate::segment::phase::{molar_mass_of, system_volumetric_flow_from_feed};
use crate::stream::Stream;

/// The material limit for the wheel tip speed, m/s.
const MAX_TIP_SPEED: f64 = 450.0;
/// The casing design pressure's margin.
const DESIGN_PRESSURE_MARGIN: f64 = 1.10;
/// The casing design temperature's margin, degC.
const DESIGN_TEMPERATURE_MARGIN: f64 = 30.0;
/// Steel's density, kg/m3.
const STEEL_DENSITY: f64 = 7850.0;
/// The radial inflow turbine's velocity ratio.
const VELOCITY_RATIO: f64 = 0.7;
/// The casing's allowable stress, MPa.
const CASING_ALLOWABLE_STRESS_MPA: f64 = 137.9;
/// The shaft's allowable shear stress, Pa.
const SHAFT_ALLOWABLE_SHEAR_PA: f64 = 50.0e6;

/// `ExpanderMechanicalDesign.calcDesign`'s stored result, in SI.
///
/// The base class's geometry is millimetres in this class, unlike the pump's metres; the
/// boundary conversions are the field constructors below.
#[derive(Debug, Clone, PartialEq)]
pub struct ExpanderDesign {
    /// `ExpanderType`: `radial_inflow` or `axial`. `mixed_flow` is never selected.
    pub expander_type: &'static str,
    /// `LoadType`: `generator` or `brake`. `compressor` is never selected.
    pub load_type: &'static str,
    /// The number of stages `selectExpanderType` assigns.
    pub number_of_stages: u32,
    /// The wheel's outer diameter.
    pub wheel_diameter: Length,
    /// The wheel's tip speed.
    pub tip_speed: Velocity,
    /// The shaft diameter, after the 1.5 factor and the 5 mm rounding.
    pub shaft_diameter: Length,
    /// The rated speed.
    pub rated_speed: f64,
    /// The maximum continuous speed, 1.05 times the rated.
    pub max_continuous_speed: f64,
    /// The trip speed, 1.10 times the rated.
    pub trip_speed: f64,
    /// The first critical speed, floored at half the rated.
    pub first_critical_speed: f64,
    /// The casing's design pressure.
    pub casing_pressure: Pressure,
    /// The casing's design temperature.
    pub casing_temperature: ThermodynamicTemperature,
    /// The casing wall thickness.
    pub casing_wall_thickness: Length,
    /// The bearing span, 2.5 wheel diameters.
    pub bearing_span: Length,
    /// The wheel weight, a hollowed disc.
    pub wheel_weight: Mass,
    /// The shaft weight.
    pub shaft_weight: Mass,
    /// The generator weight; zero for a brake load.
    pub generator_weight: Mass,
    /// The seal type, by casing pressure.
    pub seal_type: &'static str,
    /// The bearing type, by speed.
    pub bearing_type: &'static str,
    /// The casing's inner diameter.
    pub inner_diameter: Length,
    /// The casing's outer diameter.
    pub outer_diameter: Length,
    /// The tangent-to-tangent length.
    pub tantan_length: Length,
    /// The weight the base class calls the vessel: the casing.
    pub weight_vessel: Mass,
    /// The weight the base class calls the internals: the rotor bundle.
    pub weight_internals: Mass,
    /// The piping, lube and seal-system weight.
    pub weight_piping: Mass,
    /// The electrical and instrument weight.
    pub weight_electro_instrument: Mass,
    /// The structural steel and baseplate weight.
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

/// `ExpanderMechanicalDesign.calcDesign`, on the expander's run state.
///
/// `outlet_t` is the run's outlet temperature and `power` its shaft power in watts
/// (negative for an expansion). The class stores the run's isentropic efficiency but no
/// sizing arithmetic reads it, so it is not a parameter here.
///
/// # Errors
/// Whatever the density, molar mass or flash refuses.
pub fn design(
    feed: &Stream,
    outlet_t: ThermodynamicTemperature,
    power: Power,
) -> Result<ExpanderDesign> {
    let (mixture, _) = feed.mixture()?;
    let molar_mass = molar_mass_of(&mixture, &feed.z)?;

    // The class reads the inlet's own `getFlowRate` values, in kg/hr and m3/hr.
    // **The volume flow is the feed's root, not the flash's.** A negative-flash answer
    // carries the trial composition's root, the other side of the cubic; NeqSim's
    // `getDensity` reads the system's, which is what the separator family's capacity
    // surfaces use `from_feed` for.
    let mass_flow_kg_hr = feed.n * molar_mass * 3600.0;
    let volume_flow_m3_hr = system_volumetric_flow_from_feed(feed)? * 3600.0;
    let inlet_pressure_bara = feed.p.value / 1.0e5;
    let inlet_temperature_c = feed.t.value - 273.15;
    let outlet_temperature_c = outlet_t.value - 273.15;

    // `getPower("kW")` is negative for an expansion; the design takes its magnitude.
    let recovered_power_kw = (power.value / 1000.0).abs();
    let enthalpy_drop_kj_kg = if mass_flow_kg_hr > 0.0 {
        recovered_power_kw * 3600.0 / mass_flow_kg_hr
    } else {
        0.0
    };

    let casing_pressure_bara = inlet_pressure_bara * DESIGN_PRESSURE_MARGIN;
    let casing_temperature_c =
        inlet_temperature_c.max(outlet_temperature_c) + DESIGN_TEMPERATURE_MARGIN;

    let (expander_type, load_type, number_of_stages) =
        select_expander_type(volume_flow_m3_hr, enthalpy_drop_kj_kg, recovered_power_kw);
    let wheel = wheel_sizing(volume_flow_m3_hr, enthalpy_drop_kj_kg, number_of_stages);
    let shaft_diameter = shaft_sizing(wheel.rated_speed, recovered_power_kw);
    let rotor = rotor_dynamics(
        wheel.wheel_diameter,
        shaft_diameter,
        wheel.rated_speed,
        casing_pressure_bara,
    );
    let casing = casing_design(
        wheel.wheel_diameter,
        casing_pressure_bara,
        casing_temperature_c,
    );
    let generator_weight = generator_sizing(recovered_power_kw, load_type);
    let weights = weights(
        wheel.wheel_diameter,
        shaft_diameter,
        rotor.bearing_span,
        &casing,
        generator_weight,
        load_type,
        rotor.seal_type,
        recovered_power_kw,
    );
    let (module_length, module_width, module_height) = module_dimensions(
        rotor.bearing_span,
        casing.outer_diameter,
        load_type,
        recovered_power_kw,
    );

    Ok(ExpanderDesign {
        expander_type,
        load_type,
        number_of_stages,
        wheel_diameter: meters(wheel.wheel_diameter / 1000.0),
        tip_speed: meters_per_second(wheel.tip_speed),
        shaft_diameter: meters(shaft_diameter / 1000.0),
        rated_speed: wheel.rated_speed,
        max_continuous_speed: wheel.max_continuous_speed,
        trip_speed: wheel.trip_speed,
        first_critical_speed: rotor.first_critical_speed,
        casing_pressure: pascals(casing_pressure_bara * 1.0e5),
        casing_temperature: kelvins(casing_temperature_c + 273.15),
        casing_wall_thickness: meters(casing.wall_thickness / 1000.0),
        bearing_span: meters(rotor.bearing_span / 1000.0),
        wheel_weight: kilograms(weights.wheel),
        shaft_weight: kilograms(weights.shaft),
        generator_weight: kilograms(generator_weight),
        seal_type: rotor.seal_type,
        bearing_type: rotor.bearing_type,
        inner_diameter: meters(casing.inner_diameter / 1000.0),
        outer_diameter: meters(casing.outer_diameter / 1000.0),
        tantan_length: meters(rotor.bearing_span * 2.0 / 1000.0),
        weight_vessel: kilograms(weights.vessel),
        weight_internals: kilograms(weights.internals),
        weight_piping: kilograms(weights.piping),
        weight_electro_instrument: kilograms(weights.electro_instrument),
        weight_structural_steel: kilograms(weights.structural_steel),
        weight_total: kilograms(weights.total),
        module_length: meters(module_length),
        module_width: meters(module_width),
        module_height: meters(module_height),
    })
}

/// `selectExpanderType`: axial above the volume-flow threshold, radial otherwise; the
/// generator above 100 kW, the oil brake below.
fn select_expander_type(
    volume_flow_m3_hr: f64,
    enthalpy_drop_kj_kg: f64,
    recovered_power_kw: f64,
) -> (&'static str, &'static str, u32) {
    let (expander_type, stages) = if volume_flow_m3_hr > 100_000.0 {
        ("axial", (enthalpy_drop_kj_kg / 50.0).ceil().max(0.0) as u32)
    } else {
        ("radial_inflow", 1)
    };
    let load_type = if recovered_power_kw > 100.0 {
        "generator"
    } else {
        "brake"
    };
    (expander_type, load_type, stages)
}

struct Wheel {
    wheel_diameter: f64,
    tip_speed: f64,
    rated_speed: f64,
    max_continuous_speed: f64,
    trip_speed: f64,
}

/// `calculateWheelSizing`, in the class's millimetres and rpm.
fn wheel_sizing(volume_flow_m3_hr: f64, enthalpy_drop_kj_kg: f64, stages: u32) -> Wheel {
    if enthalpy_drop_kj_kg <= 0.0 || volume_flow_m3_hr <= 0.0 {
        return Wheel {
            wheel_diameter: 300.0,
            tip_speed: 300.0,
            rated_speed: 20_000.0,
            max_continuous_speed: 21_000.0,
            trip_speed: 23_100.0,
        };
    }

    let spouting_velocity = (2.0 * enthalpy_drop_kj_kg * 1000.0).sqrt();
    let tip_speed = (VELOCITY_RATIO * spouting_velocity).min(MAX_TIP_SPEED);
    let volume_flow_m3_s = volume_flow_m3_hr / 3600.0;

    let enthalpy_per_stage = enthalpy_drop_kj_kg / f64::from(stages.max(1));
    let mut rated_speed = 70.0 * enthalpy_per_stage.powf(0.75) / (volume_flow_m3_s * 60.0).sqrt();
    rated_speed = rated_speed.clamp(3000.0, 100_000.0);

    let mut wheel_diameter = tip_speed * 60.0 / (std::f64::consts::PI * rated_speed) * 1000.0;
    wheel_diameter = wheel_diameter.clamp(100.0, 1000.0);
    rated_speed = tip_speed * 60.0 / (std::f64::consts::PI * wheel_diameter / 1000.0);

    Wheel {
        wheel_diameter,
        tip_speed,
        rated_speed,
        max_continuous_speed: rated_speed * 1.05,
        trip_speed: rated_speed * 1.10,
    }
}

/// `calculateShaftSizing`, in millimetres.
fn shaft_sizing(rated_speed: f64, power_kw: f64) -> f64 {
    if rated_speed <= 0.0 || power_kw <= 0.0 {
        return 50.0;
    }
    let torque_knm = power_kw * 60.0 / (2.0 * std::f64::consts::PI * rated_speed);
    let torque_nm = torque_knm * 1000.0;
    let diameter_m =
        (16.0 * torque_nm / (std::f64::consts::PI * SHAFT_ALLOWABLE_SHEAR_PA)).powf(1.0 / 3.0);
    let diameter_mm = diameter_m * 1000.0 * 1.5;
    ((diameter_mm / 5.0).ceil() * 5.0).max(30.0)
}

struct Rotor {
    bearing_span: f64,
    first_critical_speed: f64,
    seal_type: &'static str,
    bearing_type: &'static str,
}

/// `calculateRotorDynamics`: the bearing span, the Rayleigh-Ritz critical speed and the
/// two selection bands.
fn rotor_dynamics(
    wheel_diameter: f64,
    shaft_diameter: f64,
    rated_speed: f64,
    casing_pressure_bara: f64,
) -> Rotor {
    let bearing_span = wheel_diameter * 2.5;
    let d = shaft_diameter / 1000.0;
    let length = bearing_span / 1000.0;
    let mass_per_length = STEEL_DENSITY * std::f64::consts::PI * d * d / 4.0;
    let mut first_critical_speed = 946.0 * (d.powi(4) / (mass_per_length * length.powi(3))).sqrt();
    first_critical_speed = first_critical_speed.max(rated_speed * 0.5);

    let bearing_type = if rated_speed > 50_000.0 {
        "active_magnetic"
    } else if rated_speed > 20_000.0 {
        "tilting_pad"
    } else {
        "sleeve"
    };
    let seal_type = if casing_pressure_bara > 50.0 {
        "dry_gas"
    } else {
        "labyrinth"
    };

    Rotor {
        bearing_span,
        first_critical_speed,
        seal_type,
        bearing_type,
    }
}

struct Casing {
    inner_diameter: f64,
    outer_diameter: f64,
    wall_thickness: f64,
}

/// `calculateCasingDesign`, in millimetres.
fn casing_design(
    wheel_diameter: f64,
    casing_pressure_bara: f64,
    casing_temperature_c: f64,
) -> Casing {
    let inner_diameter = wheel_diameter * 1.3;
    let pressure_mpa = casing_pressure_bara * 0.1;
    let mut allowable_stress = CASING_ALLOWABLE_STRESS_MPA;
    if casing_temperature_c > 200.0 {
        allowable_stress *= 0.85;
    } else if casing_temperature_c < -50.0 {
        allowable_stress *= 0.9;
    }
    let joint_efficiency = 0.85;
    let wall_thickness = (pressure_mpa * inner_diameter)
        / (2.0 * allowable_stress * joint_efficiency - 0.6 * pressure_mpa)
        + 3.0;
    let wall_thickness = wall_thickness.max(10.0);
    Casing {
        inner_diameter,
        outer_diameter: inner_diameter + 2.0 * wall_thickness,
        wall_thickness,
    }
}

/// `calculateGeneratorSizing`, in kilograms.
fn generator_sizing(recovered_power_kw: f64, load_type: &str) -> f64 {
    if load_type != "generator" {
        return 0.0;
    }
    if recovered_power_kw <= 100.0 {
        100.0 + recovered_power_kw * 3.0
    } else if recovered_power_kw <= 1000.0 {
        400.0 + recovered_power_kw * 2.0
    } else {
        2400.0 + recovered_power_kw
    }
}

struct Weights {
    wheel: f64,
    shaft: f64,
    vessel: f64,
    internals: f64,
    piping: f64,
    electro_instrument: f64,
    structural_steel: f64,
    total: f64,
}

/// `calculateWeights`, in kilograms.
#[allow(clippy::too_many_arguments)]
fn weights(
    wheel_diameter: f64,
    shaft_diameter: f64,
    bearing_span: f64,
    casing: &Casing,
    generator_weight: f64,
    load_type: &str,
    seal_type: &str,
    recovered_power_kw: f64,
) -> Weights {
    let wheel_volume =
        std::f64::consts::PI * (wheel_diameter / 2000.0).powi(2) * (wheel_diameter / 3000.0);
    let wheel = wheel_volume * STEEL_DENSITY * 0.7;

    let shaft_length = bearing_span * 1.5;
    let shaft_volume =
        std::f64::consts::PI * (shaft_diameter / 2000.0).powi(2) * (shaft_length / 1000.0);
    let shaft = shaft_volume * STEEL_DENSITY;

    let casing_length = bearing_span * 2.0;
    let casing_volume = std::f64::consts::PI / 4.0
        * ((casing.outer_diameter / 1000.0).powi(2) - (casing.inner_diameter / 1000.0).powi(2))
        * (casing_length / 1000.0);
    let vessel = casing_volume * STEEL_DENSITY * 1.5;

    let internals = wheel + shaft + wheel * 0.3;
    let equipment = vessel + internals + generator_weight;
    let lube_system = 200.0 + recovered_power_kw * 0.2;
    let seal_system = if seal_type == "dry_gas" { 500.0 } else { 100.0 };
    let baseplate = equipment * 0.2;
    let piping = equipment * 0.15;
    let electrical = if load_type == "generator" {
        recovered_power_kw * 0.3
    } else {
        50.0
    };
    let structural = equipment * 0.1;

    Weights {
        wheel,
        shaft,
        vessel,
        internals,
        piping: piping + lube_system + seal_system,
        electro_instrument: electrical,
        structural_steel: structural + baseplate,
        total: equipment + lube_system + seal_system + baseplate + piping + electrical + structural,
    }
}

/// `calculateModuleDimensions`, in metres.
fn module_dimensions(
    bearing_span: f64,
    outer_diameter: f64,
    load_type: &str,
    recovered_power_kw: f64,
) -> (f64, f64, f64) {
    let expander_length = bearing_span / 1000.0 * 2.5;
    let generator_length = if load_type == "generator" {
        1.0 + recovered_power_kw.sqrt() * 0.03
    } else {
        0.5
    };
    (
        expander_length + generator_length + 2.0,
        2.0_f64.max(outer_diameter / 1000.0 * 2.0 + 1.0),
        1.5_f64.max(outer_diameter / 1000.0 + 1.0),
    )
}
