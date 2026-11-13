//! `CompressorMechanicalDesign.calcDesign`, the API 617 sizing surface.
//!
//! The stored sizing, weights and module are transcribed in the order `calcDesign` runs
//! them. The class's `CompressorCasingDesignCalculator` audit - API 617/ASME/NACE checks
//! and their prose - is a separate report and is not published here; its inputs
//! (`casing_material_grade`, `casing_corrosion_allowance`, `nace_compliance`,
//! `h2s_partial_pressure`) change no field below.
//!
//! **The impeller screen is not clipped.** `setImpellerCandidate` computes the diameter
//! and the flow coefficient without clamping them, and the violations are reported as
//! issues; a stage count that no candidate satisfies leaves the head-based candidate's
//! numbers and the issue list.

use azoth_core::Result;
use azoth_core::units::{
    Length, Mass, Power, Pressure, ThermodynamicTemperature, Velocity, kelvins, kilograms, meters,
    meters_per_second, pascals, watts,
};

use crate::segment::phase::{molar_mass_of, system_volumetric_flow_from_feed};
use crate::stream::Stream;

/// The stage head cap, kJ/kg.
const MAX_HEAD_PER_STAGE: f64 = 30.0;
/// The work coefficient.
const WORK_COEFFICIENT: f64 = 0.50;
/// The flow-coefficient band.
const FLOW_COEFFICIENT_MIN: f64 = 0.01;
const FLOW_COEFFICIENT_MAX: f64 = 0.15;
/// The impeller-diameter band, mm.
const MIN_IMPELLER_DIAMETER: f64 = 100.0;
const MAX_IMPELLER_DIAMETER: f64 = 1500.0;
/// The driver margins: below 150 kW, below 750 kW, above.
const DRIVER_MARGIN_SMALL: f64 = 1.25;
const DRIVER_MARGIN_MEDIUM: f64 = 1.15;
const DRIVER_MARGIN_LARGE: f64 = 1.10;
/// The casing design margins.
const DESIGN_PRESSURE_MARGIN: f64 = 1.10;
const DESIGN_TEMPERATURE_MARGIN: f64 = 30.0;
/// Steel's density, kg/m3.
const STEEL_DENSITY: f64 = 7850.0;
/// The shaft's allowable shear stress, Pa.
const SHAFT_ALLOWABLE_SHEAR_PA: f64 = 50.0e6;
/// The gravity the class converts a chart head to metres with - `9.81`, not standard gravity.
const HEAD_TO_METRES_GRAVITY: f64 = 9.81;

/// `CompressorMechanicalDesign.CasingType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CasingType {
    VerticallySplit,
    HorizontallySplit,
    Barrel,
}

impl CasingType {
    /// The enum's spelling.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::VerticallySplit => "vertically_split",
            Self::HorizontallySplit => "horizontally_split",
            Self::Barrel => "barrel",
        }
    }
}

/// The caller's design inputs, in SI.
pub struct CompressorDesignInputs {
    /// `maxStagesPerCasing`, the impeller screen's cap.
    pub max_stages_per_casing: u32,
}

/// `CompressorMechanicalDesign.calcDesign`'s stored result, in SI.
#[derive(Debug, Clone, PartialEq)]
pub struct CompressorDesign {
    /// The stage count the impeller screen settled on.
    pub number_of_stages: u32,
    /// The impeller's outer diameter.
    pub impeller_diameter: Length,
    /// The shaft diameter.
    pub shaft_diameter: Length,
    /// The tip speed.
    pub tip_speed: Velocity,
    /// The inlet flow coefficient.
    pub flow_coefficient: f64,
    /// The speed the impeller screen used.
    pub impeller_sizing_speed: f64,
    /// The runner's inlet volume flow.
    pub impeller_sizing_flow: f64,
    /// The head the screen divided.
    pub impeller_sizing_head: f64,
    /// Whether the screen's limits are met.
    pub impeller_sizing_ok: bool,
    /// The driver power.
    pub driver_power: Power,
    /// The driver margin.
    pub driver_margin: f64,
    /// The casing design pressure.
    pub design_pressure: Pressure,
    /// The casing design temperature.
    pub design_temperature: ThermodynamicTemperature,
    /// The maximum continuous speed.
    pub max_continuous_speed: f64,
    /// The trip speed.
    pub trip_speed: f64,
    /// The first critical speed.
    pub first_critical_speed: f64,
    /// The casing type, by design pressure.
    pub casing_type: CasingType,
    /// The head per stage, kJ/kg.
    pub head_per_stage: f64,
    /// The bearing span, mm - this class stores it in millimetres.
    pub bearing_span: Length,
    /// The casing weight, kg.
    pub casing_weight: Mass,
    /// The rotor bundle weight, kg.
    pub bundle_weight: Mass,
    /// The rotor weight, kg.
    pub rotor_weight: Mass,
    /// The casing's outer diameter, m.
    pub outer_diameter: Length,
    /// The casing's inner diameter, m.
    pub inner_diameter: Length,
    /// The casing's wall thickness, m.
    pub wall_thickness: Length,
    /// The tangent-to-tangent length, m.
    pub tantan_length: Length,
    /// The base class's vessel-shell weight, kg.
    pub weight_vessel_shell: Mass,
    /// The base class's internals weight, kg.
    pub weight_internals: Mass,
    /// The base class's nozzle weight, kg - zero by the class's own comment.
    pub weight_nozzle: Mass,
    /// The piping weight, kg.
    pub weight_piping: Mass,
    /// The electrical and instrument weight, kg.
    pub weight_electro_instrument: Mass,
    /// The structural steel weight, kg.
    pub weight_structural_steel: Mass,
    /// The total weight, kg.
    pub weight_total: Mass,
    /// The module length, m.
    pub module_length: Length,
    /// The module width, m.
    pub module_width: Length,
    /// The module height, m.
    pub module_height: Length,
}

/// `calcDesign` on the compressor's run state.
///
/// `polytropic_head` is the class's `getPolytropicFluidHead`, kJ/kg. The design reads the
/// inlet and outlet states, the shaft power and the speed; the mass flow is derived from
/// the inlet.
///
/// # Errors
/// Whatever the density, molar mass or flash refuses.
pub fn design(
    feed: &Stream,
    outlet: &Stream,
    power: Power,
    polytropic_head: f64,
    speed: f64,
    inputs: &CompressorDesignInputs,
) -> Result<CompressorDesign> {
    let (mixture, _) = feed.mixture()?;
    let molar_mass = molar_mass_of(&mixture, &feed.z)?;
    let mass_flow_kg_hr = feed.n * molar_mass * 3600.0;
    // **The feed's root, not a negative flash's trial root.** NeqSim's `getFlowRate("m3/hr")`
    // reads the system density at the feed composition; the flash's trial root is the other
    // side of the cubic on a near-dew-point suction.
    let volume_flow_m3_hr = system_volumetric_flow_from_feed(feed)? * 3600.0;
    let discharge_pressure_bara = outlet.p.value / 1.0e5;
    let discharge_temperature_c = outlet.t.value - 273.15;

    let design_pressure_bara = discharge_pressure_bara * DESIGN_PRESSURE_MARGIN;
    let design_temperature_c = discharge_temperature_c + DESIGN_TEMPERATURE_MARGIN;
    let casing_type = select_casing_type(design_pressure_bara);

    let (head_based_stages, _) = number_of_stages(polytropic_head);
    let impeller = impeller_sizing(
        volume_flow_m3_hr,
        polytropic_head,
        speed,
        head_based_stages,
        inputs.max_stages_per_casing,
    );
    if !impeller.diameter.is_finite() {
        // The class returns before the rest of the design when the screen cannot size.
        return Err(azoth_core::AzothError::invalid_input(
            "mechanical_design",
            "the impeller screen could not size the machine",
        ));
    }

    let shaft_diameter = shaft_diameter(power.value / 1000.0, speed);
    let driver_margin = driver_margin(power.value / 1000.0);
    let driver_power_kw = power.value / 1000.0 * driver_margin;
    let rotor = rotor_dynamics(
        impeller.number_of_stages,
        impeller.diameter,
        shaft_diameter,
        speed,
    );
    let weights = weights(
        impeller.number_of_stages,
        impeller.diameter,
        shaft_diameter,
        rotor.bearing_span,
        casing_type,
        design_pressure_bara,
        driver_power_kw,
        mass_flow_kg_hr,
    );
    let outer_diameter = impeller.diameter * 1.5 / 1000.0;
    let inner_diameter = impeller.diameter / 1000.0;
    let module = module_dimensions(rotor.bearing_span, outer_diameter, driver_power_kw);

    Ok(CompressorDesign {
        number_of_stages: impeller.number_of_stages,
        impeller_diameter: meters(impeller.diameter / 1000.0),
        shaft_diameter: meters(shaft_diameter / 1000.0),
        tip_speed: meters_per_second(impeller.tip_speed),
        flow_coefficient: impeller.flow_coefficient,
        impeller_sizing_speed: speed,
        impeller_sizing_flow: volume_flow_m3_hr,
        impeller_sizing_head: polytropic_head,
        impeller_sizing_ok: impeller.issues.is_empty(),
        driver_power: watts(driver_power_kw * 1000.0),
        driver_margin,
        design_pressure: pascals(design_pressure_bara * 1.0e5),
        design_temperature: kelvins(design_temperature_c + 273.15),
        max_continuous_speed: rotor.max_continuous_speed,
        trip_speed: rotor.trip_speed,
        first_critical_speed: rotor.first_critical_speed,
        casing_type,
        head_per_stage: impeller.head_per_stage,
        bearing_span: meters(rotor.bearing_span / 1000.0),
        casing_weight: kilograms(weights.casing),
        bundle_weight: kilograms(weights.bundle),
        rotor_weight: kilograms(weights.rotor),
        outer_diameter: meters(outer_diameter),
        inner_diameter: meters(inner_diameter),
        wall_thickness: meters((outer_diameter - inner_diameter) / 2.0),
        tantan_length: meters(rotor.bearing_span * 1.3 / 1000.0),
        weight_vessel_shell: kilograms(weights.casing),
        weight_internals: kilograms(weights.bundle),
        weight_nozzle: kilograms(0.0),
        weight_piping: kilograms(weights.piping),
        weight_electro_instrument: kilograms(weights.electrical),
        weight_structural_steel: kilograms(weights.structural),
        weight_total: kilograms(weights.total),
        module_length: meters(module.0),
        module_width: meters(module.1),
        module_height: meters(module.2),
    })
}

/// `selectCasingType`: barrel above 100 bara, horizontal above 40, vertical otherwise.
fn select_casing_type(pressure_bara: f64) -> CasingType {
    if pressure_bara > 100.0 {
        CasingType::Barrel
    } else if pressure_bara > 40.0 {
        CasingType::HorizontallySplit
    } else {
        CasingType::VerticallySplit
    }
}

/// `calculateNumberOfStages`, in kJ/kg.
fn number_of_stages(total_head: f64) -> (u32, f64) {
    if !total_head.is_finite() || total_head <= 0.0 {
        return (1, 0.0);
    }
    let stages = ((total_head / MAX_HEAD_PER_STAGE).ceil() as u32).max(1);
    (stages, total_head / f64::from(stages))
}

struct Impeller {
    number_of_stages: u32,
    head_per_stage: f64,
    diameter: f64,
    tip_speed: f64,
    flow_coefficient: f64,
    issues: Vec<&'static str>,
}

/// `calculateImpellerSizing` and its `setImpellerCandidate`, in millimetres.
fn impeller_sizing(
    volume_flow_m3_hr: f64,
    total_head: f64,
    speed_rpm: f64,
    head_based_stages: u32,
    max_stages: u32,
) -> Impeller {
    if !speed_rpm.is_finite()
        || speed_rpm <= 0.0
        || !volume_flow_m3_hr.is_finite()
        || volume_flow_m3_hr <= 0.0
        || !total_head.is_finite()
        || total_head <= 0.0
    {
        return Impeller {
            number_of_stages: 0,
            head_per_stage: f64::NAN,
            diameter: f64::NAN,
            tip_speed: f64::NAN,
            flow_coefficient: f64::NAN,
            issues: vec!["the inlet state is not positive"],
        };
    }
    let max_stages = max_stages.max(1);
    let mut stages = head_based_stages.max(1);
    while stages > 0 && stages <= max_stages {
        let candidate = impeller_candidate(volume_flow_m3_hr, total_head, speed_rpm, stages);
        if candidate.issues.is_empty() {
            return candidate;
        }
        if candidate.diameter < MIN_IMPELLER_DIAMETER
            || candidate.flow_coefficient > FLOW_COEFFICIENT_MAX
        {
            break;
        }
        stages += 1;
    }
    impeller_candidate(
        volume_flow_m3_hr,
        total_head,
        speed_rpm,
        head_based_stages.max(1),
    )
}

/// `setImpellerCandidate`: the coupled geometry and the issue list, unclipped.
fn impeller_candidate(
    volume_flow_m3_hr: f64,
    total_head: f64,
    speed_rpm: f64,
    stages: u32,
) -> Impeller {
    let head_per_stage = total_head / f64::from(stages);
    let required_tip_speed = (head_per_stage * 1000.0 / WORK_COEFFICIENT).sqrt();
    let diameter = required_tip_speed * 60.0 / (std::f64::consts::PI * speed_rpm) * 1000.0;
    let tip_speed = std::f64::consts::PI * (diameter / 1000.0) * speed_rpm / 60.0;
    let flow_coefficient = (volume_flow_m3_hr / 3600.0) / ((diameter / 1000.0).powi(2) * tip_speed);
    let mut issues = Vec::new();
    if !(MIN_IMPELLER_DIAMETER..=MAX_IMPELLER_DIAMETER).contains(&diameter) {
        issues.push("the impeller diameter is outside 100-1500 mm");
    }
    if !(FLOW_COEFFICIENT_MIN..=FLOW_COEFFICIENT_MAX).contains(&flow_coefficient) {
        issues.push("the inlet flow coefficient is outside 0.01-0.15");
    }
    Impeller {
        number_of_stages: stages,
        head_per_stage,
        diameter,
        tip_speed,
        flow_coefficient,
        issues,
    }
}

/// `calculateShaftDiameter`, in millimetres.
fn shaft_diameter(power_kw: f64, speed_rpm: f64) -> f64 {
    if speed_rpm <= 0.0 || power_kw <= 0.0 {
        return 50.0;
    }
    let torque_nm = power_kw * 1000.0 * 60.0 / (2.0 * std::f64::consts::PI * speed_rpm);
    let diameter = (16.0 * torque_nm / (std::f64::consts::PI * SHAFT_ALLOWABLE_SHEAR_PA))
        .powf(1.0 / 3.0)
        * 1000.0
        * 1.5;
    diameter.clamp(30.0, 300.0)
}

/// `calculateDriverSizing`'s margin.
fn driver_margin(shaft_power_kw: f64) -> f64 {
    if shaft_power_kw < 150.0 {
        DRIVER_MARGIN_SMALL
    } else if shaft_power_kw < 750.0 {
        DRIVER_MARGIN_MEDIUM
    } else {
        DRIVER_MARGIN_LARGE
    }
}

struct Rotor {
    max_continuous_speed: f64,
    trip_speed: f64,
    bearing_span: f64,
    first_critical_speed: f64,
}

/// `calculateRotorDynamics`, in rpm and millimetres.
fn rotor_dynamics(
    stages: u32,
    impeller_diameter: f64,
    shaft_diameter: f64,
    operating_speed: f64,
) -> Rotor {
    let max_continuous_speed = operating_speed * 1.05;
    let trip_speed = max_continuous_speed * 1.05;
    let bearing_span = f64::from(stages) * (impeller_diameter * 0.8) + impeller_diameter;
    let stiffness = 48.0 * 2.1e11 * std::f64::consts::PI * (shaft_diameter / 2000.0).powi(4) / 64.0;
    let mass_per_length = STEEL_DENSITY * std::f64::consts::PI * (shaft_diameter / 2000.0).powi(2);
    let first_critical_speed = if bearing_span > 0.0 {
        (946.0 / (2.0 * std::f64::consts::PI))
            * (stiffness / (mass_per_length * (bearing_span / 1000.0).powi(3))).sqrt()
            * 60.0
    } else {
        f64::NAN
    };
    Rotor {
        max_continuous_speed,
        trip_speed,
        bearing_span,
        first_critical_speed,
    }
}

struct Weights {
    rotor: f64,
    casing: f64,
    bundle: f64,
    piping: f64,
    electrical: f64,
    structural: f64,
    total: f64,
}

/// `calculateWeights`, in kilograms.
#[allow(clippy::too_many_arguments)]
fn weights(
    stages: u32,
    impeller_diameter: f64,
    shaft_diameter: f64,
    bearing_span: f64,
    casing_type: CasingType,
    design_pressure_bara: f64,
    driver_power_kw: f64,
    _mass_flow_kg_hr: f64,
) -> Weights {
    let impeller_weight = f64::from(stages) * 0.5 * (impeller_diameter / 100.0).powf(2.5);
    let shaft_weight = bearing_span / 1000.0
        * STEEL_DENSITY
        * std::f64::consts::PI
        * (shaft_diameter / 2000.0).powi(2);
    let rotor = impeller_weight + shaft_weight;

    let casing_thickness = (design_pressure_bara * impeller_diameter / (2.0 * 150.0)).max(10.0);
    let casing_length = bearing_span * 1.3;
    let casing_outer_diameter = impeller_diameter * 1.5;
    let mut casing = std::f64::consts::PI
        * (casing_outer_diameter / 1000.0)
        * (casing_length / 1000.0)
        * (casing_thickness / 1000.0)
        * STEEL_DENSITY;
    casing *= 1.2;
    if casing_type == CasingType::Barrel {
        casing *= 1.3;
    }

    let internals = f64::from(stages) * 20.0 * (impeller_diameter / 300.0).powi(2);
    let bundle = rotor + internals;
    let seal_system = 100.0 * (shaft_diameter / 100.0);
    let lube_system = 200.0 + driver_power_kw * 0.1;
    let base = casing * 0.3;
    let empty_vessel = casing + base;
    let piping = empty_vessel * 0.2;
    let electrical = driver_power_kw * 0.5;
    let structural = empty_vessel * 0.15;
    let total =
        empty_vessel + bundle + seal_system + lube_system + piping + electrical + structural;

    Weights {
        rotor,
        casing,
        bundle,
        piping,
        electrical,
        structural,
        total,
    }
}

/// `calculateModuleDimensions`, in metres.
fn module_dimensions(
    bearing_span: f64,
    outer_diameter: f64,
    driver_power_kw: f64,
) -> (f64, f64, f64) {
    let compressor_length = bearing_span * 1.5 / 1000.0;
    let driver_length = 0.5 + driver_power_kw * 0.001;
    let coupling_space = 0.5;
    let auxiliary_space = 2.0;
    let module_length =
        (compressor_length + driver_length + coupling_space + auxiliary_space).max(4.0);
    let module_width = (outer_diameter + 3.0).max(3.0);
    let module_height = (outer_diameter + 2.0).max(3.0);
    (module_length, module_width, module_height)
}

/// The class's fluid head, kJ/kg, from its metre head. The chart branch stores the head in
/// metres and converts it back with this gravity, not with standard gravity.
pub fn head_from_metres(head_m: f64) -> f64 {
    head_m * HEAD_TO_METRES_GRAVITY / 1000.0
}

/// The class's fluid head, kJ/kg, from the shaft power in the no-chart polytropic branch:
/// `power / kg_per_sec / 1000 * efficiency`.
pub fn head_from_power(power_w: f64, mass_flow_kg_s: f64, efficiency: f64) -> f64 {
    power_w / mass_flow_kg_s / 1000.0 * efficiency
}
