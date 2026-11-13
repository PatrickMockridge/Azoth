//! `PumpMechanicalDesign.calcDesign`, the API 610 sizing surface.
//!
//! The arithmetic is transcribed in the order `calcDesign` runs it. The class's own
//! gravity is `9.81`, which is not the chart kernel's `9.80665` and not the API 610
//! calculator's; that is why the design's head is a distinct quantity from the model's.
//!
//! The API 610 assessment's driver selection is here because its selected rating
//! overwrites the design's driver power, and every weight and module dimension is built
//! on that value. The assessment's checks are a separate surface.

use azoth_core::Result;
use azoth_core::units::{
    Length, Mass, Power, Pressure, ThermodynamicTemperature, VolumeRate, cubic_meter_per_hour,
    cubic_meters_per_hour, kelvins, kilograms, meters, pascals, watts,
};

use crate::pump_chart::PumpChart;
use crate::segment::phase::{system_mass_density, system_volumetric_flow};
use crate::stream::Stream;

/// The class's own gravity, `PumpMechanicalDesign.GRAVITY`.
const GRAVITY: f64 = 9.81;
/// The API 610 calculator's gravity, which the assessment's pressure basis uses. It is not
/// the design's 9.81, and that is why the design pressure and the required casing pressure
/// differ by a few parts in a thousand.
const API610_GRAVITY: f64 = 9.80665;
/// The driver margins, from absorbed power: below 22 kW, at or below 55 kW, above.
const DRIVER_MARGIN_SMALL: f64 = 1.25;
const DRIVER_MARGIN_MEDIUM: f64 = 1.15;
const DRIVER_MARGIN_LARGE: f64 = 1.10;
/// The design temperature's margin, degC.
const DESIGN_TEMPERATURE_MARGIN: f64 = 30.0;
/// The class's default pump efficiency.
const DEFAULT_PUMP_EFFICIENCY: f64 = 0.75;
/// Steel's density, kg/m3.
const STEEL_DENSITY: f64 = 7850.0;
/// The casing's allowable stress, MPa.
const ALLOWABLE_STRESS_MPA: f64 = 137.9;
/// The shaft's allowable shear stress, Pa.
const SHAFT_ALLOWABLE_SHEAR_PA: f64 = 40.0e6;

/// The API 610 construction type, `PumpApi610DesignCalculator.Api610PumpType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Api610PumpType {
    Unspecified,
    Oh1,
    Oh2,
    Oh3,
    Oh4,
    Oh5,
    Oh6,
    Bb1,
    Bb2,
    Bb3,
    Bb4,
    Bb5,
    Vs1,
    Vs2,
    Vs3,
    Vs4,
    Vs5,
    Vs6,
    Vs7,
}

impl Api610PumpType {
    /// The enum's spelling, as the class prints it.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Unspecified => "unspecified",
            Self::Oh1 => "oh1",
            Self::Oh2 => "oh2",
            Self::Oh3 => "oh3",
            Self::Oh4 => "oh4",
            Self::Oh5 => "oh5",
            Self::Oh6 => "oh6",
            Self::Bb1 => "bb1",
            Self::Bb2 => "bb2",
            Self::Bb3 => "bb3",
            Self::Bb4 => "bb4",
            Self::Bb5 => "bb5",
            Self::Vs1 => "vs1",
            Self::Vs2 => "vs2",
            Self::Vs3 => "vs3",
            Self::Vs4 => "vs4",
            Self::Vs5 => "vs5",
            Self::Vs6 => "vs6",
            Self::Vs7 => "vs7",
        }
    }

    /// The spec's string as the enum, or `None`.
    pub fn parse(value: &str) -> Option<Self> {
        Some(match value {
            "unspecified" => Self::Unspecified,
            "oh1" => Self::Oh1,
            "oh2" => Self::Oh2,
            "oh3" => Self::Oh3,
            "oh4" => Self::Oh4,
            "oh5" => Self::Oh5,
            "oh6" => Self::Oh6,
            "bb1" => Self::Bb1,
            "bb2" => Self::Bb2,
            "bb3" => Self::Bb3,
            "bb4" => Self::Bb4,
            "bb5" => Self::Bb5,
            "vs1" => Self::Vs1,
            "vs2" => Self::Vs2,
            "vs3" => Self::Vs3,
            "vs4" => Self::Vs4,
            "vs5" => Self::Vs5,
            "vs6" => Self::Vs6,
            "vs7" => Self::Vs7,
            _ => return None,
        })
    }

    /// The legacy family the class maps the API 610 code to.
    fn legacy_family(self) -> &'static str {
        match self {
            Self::Bb1 | Self::Bb2 | Self::Bb3 | Self::Bb4 | Self::Bb5 => "between_bearings",
            Self::Vs1 | Self::Vs2 | Self::Vs3 | Self::Vs4 | Self::Vs5 | Self::Vs6 | Self::Vs7 => {
                "vertically_suspended"
            }
            _ => "overhung",
        }
    }
}

/// The seal type, `PumpMechanicalDesign.SealType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SealType {
    Packed,
    SingleMechanical,
    DoubleMechanical,
    MagneticDrive,
    CannedMotor,
}

impl SealType {
    /// The enum's spelling.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Packed => "packed",
            Self::SingleMechanical => "single_mechanical",
            Self::DoubleMechanical => "double_mechanical",
            Self::MagneticDrive => "magnetic_drive",
            Self::CannedMotor => "canned_motor",
        }
    }

    /// The spec's string as the enum, or `None`.
    pub fn parse(value: &str) -> Option<Self> {
        Some(match value {
            "packed" => Self::Packed,
            "single_mechanical" => Self::SingleMechanical,
            "double_mechanical" => Self::DoubleMechanical,
            "magnetic_drive" => Self::MagneticDrive,
            "canned_motor" => Self::CannedMotor,
            _ => return None,
        })
    }

    /// Whether the class treats it as sealless.
    fn is_sealless(self) -> bool {
        matches!(self, Self::MagneticDrive | Self::CannedMotor)
    }
}

/// The API 610 check vocabulary, in the class's order of severity.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Status {
    Pass,
    Warning,
    Fail,
    NotEvaluated,
}

/// `updateAssessmentStatus`: a failure anywhere is a failure; otherwise a warning or a
/// gap beside at least one pass is `PASS_WITH_WARNINGS`.
fn assessment_status(statuses: &[Status]) -> &'static str {
    let mut evaluated = false;
    let mut warning_or_gap = false;
    for status in statuses {
        match status {
            Status::Fail => return "FAIL",
            Status::Pass => evaluated = true,
            Status::Warning | Status::NotEvaluated => warning_or_gap = true,
        }
    }
    if !evaluated {
        "NOT_EVALUATED"
    } else if warning_or_gap {
        "PASS_WITH_WARNINGS"
    } else {
        "PASS"
    }
}

/// The caller's design inputs, in SI.
pub struct PumpDesignInputs<'a> {
    /// The purchaser's API 610 type, or `None` for the head-driven screening selection.
    pub api610_pump_type: Option<Api610PumpType>,
    /// The purchaser's seal type.
    pub seal_type: SealType,
    /// The maximum suction pressure; absent is the inlet's.
    pub maximum_suction_pressure: Option<Pressure>,
    /// The furnished casing MAWP; absent is no evidence.
    pub furnished_casing_mawp: Option<Pressure>,
    /// The purchaser's shutoff head; absent is the chart's, or none.
    pub shutoff_head: Option<Length>,
    /// The NPSH margin factor, effectively `1.20` from the class's database.
    pub npsh_margin_factor: f64,
    /// The minimum NPSH margin.
    pub minimum_npsh_margin: Length,
    /// The preferred operating region's low and high fractions.
    pub por_low: f64,
    pub por_high: f64,
    /// The allowable operating region's low and high fractions.
    pub aor_low: f64,
    pub aor_high: f64,
    /// The rated-point band's low and high fractions.
    pub rated_low: f64,
    pub rated_high: f64,
    /// The pump's speed, rpm.
    pub speed: f64,
    /// The vendor curve, when one is stated.
    pub chart: Option<&'a PumpChart>,
}

/// `PumpMechanicalDesign.calcDesign`'s stored result, in SI.
#[derive(Debug, Clone, PartialEq)]
pub struct PumpDesign {
    /// The legacy family the API 610 type maps to.
    pub pump_type: &'static str,
    /// The API 610 construction type.
    pub api610_pump_type: &'static str,
    /// The seal type.
    pub seal_type: &'static str,
    /// The rated volumetric flow, m3/h.
    pub rated_flow: VolumeRate,
    /// The design's head, m: the differential pressure over `rho * 9.81`.
    pub rated_head: Length,
    /// The hydraulic power, W.
    pub hydraulic_power: Power,
    /// The input's normalized efficiency.
    pub pump_efficiency: f64,
    /// The rated speed, rpm.
    pub rated_speed: f64,
    /// The SI specific speed.
    pub specific_speed: f64,
    /// The US suction specific speed.
    pub suction_specific_speed: f64,
    /// The stage count.
    pub number_of_stages: u32,
    /// The vendor best-efficiency flow, m3/h; absent when no curve states one.
    pub bep_flow: Option<VolumeRate>,
    /// The vendor best-efficiency head, m.
    pub bep_head: Option<Length>,
    /// The BEP's provenance.
    pub bep_source: &'static str,
    /// The NPSHr provenance.
    pub npshr_source: &'static str,
    /// The shutoff head's provenance.
    pub shutoff_source: &'static str,
    /// The effective shutoff head, m; absent when none is known.
    pub shutoff_head: Option<Length>,
    /// The casing design pressure, Pa.
    pub design_pressure: Pressure,
    /// The casing design temperature, K.
    pub design_temperature: ThermodynamicTemperature,
    /// The impeller diameter, m.
    pub impeller_diameter: Length,
    /// The impeller width, m.
    pub impeller_width: Length,
    /// The shaft diameter, m.
    pub shaft_diameter: Length,
    /// The casing wall thickness, m.
    pub casing_wall_thickness: Length,
    /// The maximum allowable working pressure, Pa.
    pub max_allowable_working_pressure: Pressure,
    /// The final driver power, W: the API 610 selection applied after the IEC rounding.
    pub driver_power: Power,
    /// The driver margin the absorbed power selects.
    pub driver_margin: f64,
    /// The suction nozzle size, m.
    pub suction_nozzle_size: Length,
    /// The discharge nozzle size, m.
    pub discharge_nozzle_size: Length,
    /// The minimum continuous flow, m3/h.
    pub minimum_flow: VolumeRate,
    /// The maximum continuous flow, m3/h.
    pub maximum_flow: VolumeRate,
    /// The casing weight, kg.
    pub casing_weight: Mass,
    /// The motor weight, kg.
    pub motor_weight: Mass,
    /// The baseplate weight, kg.
    pub baseplate_weight: Mass,
    /// The casing's inner diameter, m.
    pub inner_diameter: Length,
    /// The casing's outer diameter, m.
    pub outer_diameter: Length,
    /// The tangent-to-tangent length, m.
    pub tantan_length: Length,
    /// The base class's wall thickness, mm - the class stores this field in millimetres.
    pub wall_thickness: Length,
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
    /// The NPSH margin, m.
    pub npsh_margin: Length,
    /// The NPSH margin ratio.
    pub npsh_margin_ratio: f64,
    /// The `assessmentStatus` aggregation.
    pub assessment_status: &'static str,
    /// `getAssessmentBasis`.
    pub assessment_basis: &'static str,
    /// The operating flow over the BEP flow; absent without a BEP.
    pub operating_flow_ratio: f64,
    /// `POR`, `AOR`, `OUTSIDE_AOR` or `NOT_EVALUATED`.
    pub operating_region: &'static str,
    /// The shutoff head the pressure basis uses.
    pub governing_shutoff_head: Length,
    /// The pressure the casing must hold, Pa: suction plus governing shutoff.
    pub required_casing_pressure: Pressure,
    /// The preliminary hydrostatic test pressure, Pa.
    pub preliminary_hydrostatic_test_pressure: Pressure,
    /// The head rise from rated flow to shutoff, as a fraction.
    pub head_rise_to_shutoff: f64,
    /// The driver power the absorbed power and margin require, W.
    pub required_driver_power: Power,
    /// `purchaser_input` or `screening_estimate`.
    pub pump_type_source: &'static str,
}

/// `calcDesign` on the pump's run state.
///
/// `outlet` is the run's outlet, `power` its shaft power, `efficiency` the efficiency the
/// run used, and the NPSH pair the model's own values.
///
/// # Errors
/// Whatever the flash, density or chart refuses.
#[allow(clippy::too_many_arguments)] // the class's own fields, one argument per group
pub fn design(
    feed: &Stream,
    outlet: &Stream,
    power: Power,
    efficiency: f64,
    npsh_required: Length,
    npsh_available: Length,
    inputs: &PumpDesignInputs<'_>,
) -> Result<PumpDesign> {
    let density = system_mass_density(feed)?;
    let volume_flow_m3_hr = system_volumetric_flow(feed)? * 3600.0;
    let suction_pressure_bara = feed.p.value / 1.0e5;
    let discharge_pressure_bara = outlet.p.value / 1.0e5;
    let suction_temperature_c = feed.t.value - 273.15;
    let discharge_temperature_c = outlet.t.value - 273.15;

    let differential_pressure = discharge_pressure_bara - suction_pressure_bara;
    let pump_head = if density > 0.0 {
        differential_pressure * 1.0e5 / (density * GRAVITY)
    } else {
        f64::NAN
    };
    let hydraulic_power_kw = if density > 0.0 && volume_flow_m3_hr > 0.0 && pump_head > 0.0 {
        density * GRAVITY * volume_flow_m3_hr / 3600.0 * pump_head / 1000.0
    } else {
        f64::NAN
    };

    let pump_efficiency = {
        let normalized = normalize_efficiency(efficiency);
        if normalized > 0.0 && normalized <= 1.0 {
            normalized
        } else {
            DEFAULT_PUMP_EFFICIENCY
        }
    };
    let mut absorbed_power_kw = power.value / 1000.0;
    if !(absorbed_power_kw > 0.0) && hydraulic_power_kw > 0.0 && pump_efficiency > 0.0 {
        absorbed_power_kw = hydraulic_power_kw / pump_efficiency;
    }

    let (bep_flow, bep_head, bep_source) = vendor_bep(inputs.chart, inputs.speed, density);
    let (shutoff_head, shutoff_source) =
        vendor_shutoff(inputs.chart, inputs.speed, density, inputs.shutoff_head);
    let npshr_source = if inputs.chart.is_some_and(|chart| chart.has_npsh_curve()) {
        "vendor_curve"
    } else if npsh_required.value > 0.0 {
        "screening_estimate"
    } else {
        "not_available"
    };

    // `bep_flow` is a `VolumeRate` whose `.value` is m3/s; the class's arithmetic is in
    // m3/h, so the conversion is explicit at every read.
    let bep_flow_m3_hr = bep_flow.map(|flow| flow.get::<cubic_meter_per_hour>());
    let specific_speed_flow = if bep_flow_m3_hr.is_some_and(|flow| flow > 0.0) {
        bep_flow_m3_hr.unwrap()
    } else {
        volume_flow_m3_hr
    };
    let specific_speed_head = if bep_head.is_some_and(|head| head.value > 0.0) {
        bep_head.unwrap().value
    } else {
        pump_head
    };
    let (specific_speed, suction_specific_speed) = specific_speeds(
        specific_speed_flow,
        specific_speed_head,
        inputs.speed,
        npsh_required.value,
        inputs.api610_pump_type,
    );

    let (pump_type, api610, number_of_stages) =
        select_pump_type(inputs.api610_pump_type, pump_head);

    let maximum_suction_bara = inputs
        .maximum_suction_pressure
        .map(|pressure| pressure.value / 1.0e5)
        .unwrap_or(suction_pressure_bara);
    let governing_shutoff = shutoff_head
        .filter(|head| head.value > 0.0)
        .map(|head| head.value)
        .unwrap_or(pump_head * 1.10);
    let design_pressure_bara = (maximum_suction_bara
        + density * GRAVITY * governing_shutoff / 1.0e5)
        .max(discharge_pressure_bara);
    let design_temperature_c =
        suction_temperature_c.max(discharge_temperature_c) + DESIGN_TEMPERATURE_MARGIN;

    let (impeller_diameter, impeller_width) =
        impeller_sizing(volume_flow_m3_hr, pump_head, inputs.speed, number_of_stages);
    let shaft_diameter = shaft_sizing(absorbed_power_kw, inputs.speed);
    let casing = casing_design(
        impeller_diameter,
        design_pressure_bara,
        design_temperature_c,
    );
    let driver_margin = driver_margin(absorbed_power_kw);
    let (suction_nozzle, discharge_nozzle) = nozzle_sizes(volume_flow_m3_hr);
    let intermediate_driver_power = calculate_driver_sizing(absorbed_power_kw);
    let weights = weights(
        impeller_diameter,
        impeller_width,
        shaft_diameter,
        casing.wall_thickness,
        pump_type,
        number_of_stages,
        intermediate_driver_power,
    );
    let (module_length, module_width, module_height) =
        module_dimensions(impeller_diameter, intermediate_driver_power);

    // **The API 610 assessment's rating is the design's driver power.** The class's own
    // `calculateDriverSizing` is the intermediate; the assessment selects from its longer
    // candidate list and overwrites the field.
    let driver_power_kw = api610_driver_rating(absorbed_power_kw);

    let npsh_margin = npsh_available.value - npsh_required.value;
    let npsh_margin_ratio = if npsh_required.value > 0.0 {
        npsh_available.value / npsh_required.value
    } else {
        f64::NAN
    };

    // **The assessment's per-check list is a presentation surface and is not published**;
    // `assessment_status` is the class's aggregation of it, and the scalars below are the
    // numeric consequences the check list produces.
    let bep_flow_value = bep_flow_m3_hr.unwrap_or(f64::NAN);
    let operating_flow_ratio = if volume_flow_m3_hr > 0.0 && bep_flow_value > 0.0 {
        volume_flow_m3_hr / bep_flow_value
    } else {
        f64::NAN
    };
    let (operating_region, operating_status) = if operating_flow_ratio.is_nan() {
        ("NOT_EVALUATED", Status::NotEvaluated)
    } else if operating_flow_ratio >= inputs.rated_low && operating_flow_ratio <= inputs.rated_high
    {
        ("POR", Status::Pass)
    } else if operating_flow_ratio >= inputs.por_low && operating_flow_ratio <= inputs.por_high {
        ("POR", Status::Warning)
    } else if operating_flow_ratio >= inputs.aor_low && operating_flow_ratio <= inputs.aor_high {
        ("AOR", Status::Warning)
    } else {
        ("OUTSIDE_AOR", Status::Fail)
    };
    let npsh_status = if npsh_available.value.is_finite() && npsh_required.value > 0.0 {
        let ratio = npsh_available.value / npsh_required.value;
        let margin = npsh_available.value - npsh_required.value;
        if ratio >= inputs.npsh_margin_factor && margin >= inputs.minimum_npsh_margin.value {
            Status::Pass
        } else {
            Status::Fail
        }
    } else {
        Status::NotEvaluated
    };
    let required_casing_pressure_bara =
        maximum_suction_bara + density * API610_GRAVITY * governing_shutoff / 1.0e5;
    let pressure_status = if maximum_suction_bara > 0.0 && pump_head > 0.0 && density > 0.0 {
        match inputs.furnished_casing_mawp {
            None => Status::NotEvaluated,
            Some(mawp) => {
                if mawp.value / 1.0e5 >= required_casing_pressure_bara {
                    Status::Pass
                } else {
                    Status::Fail
                }
            }
        }
    } else {
        Status::NotEvaluated
    };
    let hydrostatic_basis_bara = inputs
        .furnished_casing_mawp
        .filter(|mawp| mawp.value > 0.0)
        .map(|mawp| mawp.value / 1.0e5)
        .unwrap_or(required_casing_pressure_bara);
    let preliminary_hydrostatic_test_pressure_bara = hydrostatic_basis_bara * 1.50;
    let head_rise_to_shutoff = if matches!(shutoff_source, "vendor_curve" | "purchaser_input")
        && shutoff_head.is_some_and(|head| head.value > 0.0)
        && pump_head > 0.0
    {
        shutoff_head.unwrap().value / pump_head - 1.0
    } else {
        f64::NAN
    };
    let head_rise_status = if head_rise_to_shutoff.is_nan() {
        Status::NotEvaluated
    } else if head_rise_to_shutoff >= 0.10 {
        Status::Pass
    } else {
        Status::Fail
    };
    let required_driver_power_kw = absorbed_power_kw * driver_margin;
    let driver_status = if absorbed_power_kw > 0.0 {
        if required_driver_power_kw <= 10_000.0 {
            Status::Pass
        } else {
            Status::Fail
        }
    } else {
        Status::NotEvaluated
    };
    let pump_type_status = match inputs.api610_pump_type {
        None => Status::Warning,
        Some(Api610PumpType::Unspecified) => Status::NotEvaluated,
        Some(_) => Status::Pass,
    };
    let mut statuses = vec![
        operating_status,
        npsh_status,
        pressure_status,
        Status::NotEvaluated,
        head_rise_status,
        driver_status,
        Status::NotEvaluated,
        Status::NotEvaluated,
        Status::NotEvaluated,
        Status::NotEvaluated,
        Status::NotEvaluated,
        pump_type_status,
    ];
    if inputs.seal_type.is_sealless() {
        statuses.push(Status::Warning);
    }
    let status = assessment_status(&statuses);
    let pump_type_source = if inputs.api610_pump_type.is_some() {
        "purchaser_input"
    } else {
        "screening_estimate"
    };

    Ok(PumpDesign {
        pump_type,
        api610_pump_type: api610.as_str(),
        seal_type: inputs.seal_type.as_str(),
        rated_flow: cubic_meters_per_hour(volume_flow_m3_hr),
        rated_head: meters(pump_head),
        hydraulic_power: watts(hydraulic_power_kw * 1000.0),
        pump_efficiency,
        rated_speed: inputs.speed,
        specific_speed,
        suction_specific_speed,
        number_of_stages,
        bep_flow,
        bep_head,
        bep_source,
        npshr_source,
        shutoff_source,
        shutoff_head,
        design_pressure: pascals(design_pressure_bara * 1.0e5),
        design_temperature: kelvins(design_temperature_c + 273.15),
        impeller_diameter: meters(impeller_diameter / 1000.0),
        impeller_width: meters(impeller_width / 1000.0),
        shaft_diameter: meters(shaft_diameter / 1000.0),
        casing_wall_thickness: meters(casing.wall_thickness / 1000.0),
        max_allowable_working_pressure: pascals(casing.mawp_bara * 1.0e5),
        driver_power: watts(driver_power_kw * 1000.0),
        driver_margin,
        suction_nozzle_size: meters(suction_nozzle * 0.0254),
        discharge_nozzle_size: meters(discharge_nozzle * 0.0254),
        minimum_flow: cubic_meters_per_hour(
            (if bep_flow_m3_hr.is_some_and(|flow| flow > 0.0) {
                bep_flow_m3_hr.unwrap()
            } else {
                volume_flow_m3_hr
            }) * 0.65,
        ),
        maximum_flow: cubic_meters_per_hour(
            (if bep_flow_m3_hr.is_some_and(|flow| flow > 0.0) {
                bep_flow_m3_hr.unwrap()
            } else {
                volume_flow_m3_hr
            }) * 1.15,
        ),
        casing_weight: kilograms(weights.casing),
        motor_weight: kilograms(weights.motor),
        baseplate_weight: kilograms(weights.baseplate),
        inner_diameter: meters(casing.inner_diameter / 1000.0),
        outer_diameter: meters(casing.outer_diameter / 1000.0),
        tantan_length: meters(impeller_width * f64::from(number_of_stages) / 1000.0),
        wall_thickness: meters(casing.wall_thickness / 1000.0),
        weight_piping: kilograms(weights.piping),
        weight_electro_instrument: kilograms(weights.electro_instrument),
        weight_structural_steel: kilograms(weights.structural_steel),
        weight_total: kilograms(weights.total),
        module_length: meters(module_length),
        module_width: meters(module_width),
        module_height: meters(module_height),
        npsh_margin: meters(npsh_margin),
        npsh_margin_ratio,
        assessment_status: status,
        assessment_basis: "API 610 screening; vendor and purchaser verification required",
        operating_flow_ratio,
        operating_region,
        governing_shutoff_head: meters(governing_shutoff),
        required_casing_pressure: pascals(required_casing_pressure_bara * 1.0e5),
        preliminary_hydrostatic_test_pressure: pascals(
            preliminary_hydrostatic_test_pressure_bara * 1.0e5,
        ),
        head_rise_to_shutoff,
        required_driver_power: watts(required_driver_power_kw * 1000.0),
        pump_type_source,
    })
}

/// The class's `normalizeEfficiency`: a percentage becomes a fraction.
fn normalize_efficiency(efficiency: f64) -> f64 {
    if efficiency > 1.0 && efficiency <= 100.0 {
        efficiency / 100.0
    } else {
        efficiency
    }
}

/// `populateVendorCurveData`'s BEP pair and provenance.
fn vendor_bep(
    chart: Option<&PumpChart>,
    speed: f64,
    _density: f64,
) -> (Option<VolumeRate>, Option<Length>, &'static str) {
    let Some(chart) = chart else {
        return (None, None, "not_available");
    };
    if !(speed > 0.0) {
        return (None, None, "not_available");
    }
    let flow = chart.best_efficiency_flow_rate(speed);
    if !(flow > 0.0) {
        return (None, None, "not_available");
    }
    // The class's converter returns a metre head unchanged; the model's curve is metres.
    let head = chart.head(flow, speed);
    if !(head > 0.0) {
        return (None, None, "not_available");
    }
    (
        Some(cubic_meters_per_hour(flow)),
        Some(meters(head)),
        "vendor_curve",
    )
}

/// `populateVendorCurveData`'s shutoff head and provenance.
fn vendor_shutoff(
    chart: Option<&PumpChart>,
    speed: f64,
    _density: f64,
    purchaser: Option<Length>,
) -> (Option<Length>, &'static str) {
    if let Some(head) = purchaser {
        return (Some(head), "purchaser_input");
    }
    let Some(chart) = chart else {
        return (None, "not_available");
    };
    if !(speed > 0.0) {
        return (None, "not_available");
    }
    let head = chart.head(0.0, speed);
    if head > 0.0 {
        (Some(meters(head)), "vendor_curve")
    } else {
        (None, "not_available")
    }
}

/// `calculateSpecificSpeeds`, in the class's two conventions.
fn specific_speeds(
    flow_m3_hr: f64,
    head_m: f64,
    speed_rpm: f64,
    npsh_required_m: f64,
    api610: Option<Api610PumpType>,
) -> (f64, f64) {
    let flow_m3_s = flow_m3_hr / 3600.0;
    let specific = if flow_m3_s > 0.0 && head_m > 0.0 && speed_rpm > 0.0 {
        speed_rpm * flow_m3_s.sqrt() / head_m.powf(0.75)
    } else {
        f64::NAN
    };
    let suction = if npsh_required_m > 0.0 {
        let mut flow_gpm = flow_m3_hr * 4.402867;
        if api610 == Some(Api610PumpType::Bb1) {
            flow_gpm /= 2.0;
        }
        let npsh_feet = npsh_required_m * 3.28084;
        speed_rpm * flow_gpm.sqrt() / npsh_feet.powf(0.75)
    } else {
        f64::NAN
    };
    (specific, suction)
}

/// `selectPumpType`: the purchaser's API 610 code, or the head-driven screening choice.
fn select_pump_type(
    api610: Option<Api610PumpType>,
    head_m: f64,
) -> (&'static str, Api610PumpType, u32) {
    if let Some(api610) = api610 {
        let family = api610.legacy_family();
        let mut stages = 1;
        if family == "between_bearings"
            && matches!(
                api610,
                Api610PumpType::Bb3 | Api610PumpType::Bb4 | Api610PumpType::Bb5
            )
        {
            stages = stages.max(2);
        }
        return (family, api610, stages);
    }
    if head_m > 250.0 {
        let mut stages = (head_m / 200.0).ceil() as u32;
        stages = stages.max(2);
        let api610 = if stages > 2 {
            Api610PumpType::Bb3
        } else {
            Api610PumpType::Bb2
        };
        ("between_bearings", api610, stages)
    } else {
        ("overhung", Api610PumpType::Oh2, 1)
    }
}

/// `calculateImpellerSizing`, in millimetres.
fn impeller_sizing(flow_m3_hr: f64, head_m: f64, speed_rpm: f64, stages: u32) -> (f64, f64) {
    if speed_rpm <= 0.0 || flow_m3_hr <= 0.0 || head_m <= 0.0 {
        // The class sets the diameter and returns; the width keeps its initialiser.
        return (200.0, 30.0);
    }
    let head_per_stage = head_m / f64::from(stages.max(1));
    let tip_speed = (2.0 * GRAVITY * head_per_stage).sqrt();
    let diameter =
        (60.0 * tip_speed / (std::f64::consts::PI * speed_rpm) * 1000.0).clamp(50.0, 1000.0);
    let flow_m3_s = flow_m3_hr / 3600.0;
    let d2 = diameter / 1000.0;
    let cm2 = 0.1 * tip_speed;
    let width = if cm2 > 0.0 && d2 > 0.0 {
        (flow_m3_s / (std::f64::consts::PI * d2 * cm2) * 1000.0).clamp(10.0, 200.0)
    } else {
        10.0
    };
    (diameter, width)
}

/// `calculateShaftSizing`, in millimetres.
fn shaft_sizing(power_kw: f64, speed_rpm: f64) -> f64 {
    if speed_rpm <= 0.0 {
        return 40.0;
    }
    let torque_nm = power_kw * 60.0 / (2.0 * std::f64::consts::PI * speed_rpm) * 1000.0;
    let diameter_m =
        (16.0 * torque_nm / (std::f64::consts::PI * SHAFT_ALLOWABLE_SHEAR_PA)).powf(1.0 / 3.0);
    let diameter_mm = diameter_m * 1000.0 * 1.2;
    ((diameter_mm / 5.0).ceil() * 5.0).max(20.0)
}

struct Casing {
    inner_diameter: f64,
    outer_diameter: f64,
    wall_thickness: f64,
    mawp_bara: f64,
}

/// `calculateCasingDesign`, in millimetres and bar.
fn casing_design(
    impeller_diameter: f64,
    design_pressure_bara: f64,
    design_temperature_c: f64,
) -> Casing {
    let inner_diameter = impeller_diameter * 1.15;
    let pressure_mpa = design_pressure_bara * 0.1;
    let mut allowable_stress = ALLOWABLE_STRESS_MPA;
    if design_temperature_c > 200.0 {
        allowable_stress *= 0.85;
    }
    let joint_efficiency = 1.0;
    let wall_thickness = (pressure_mpa * inner_diameter
        / (2.0 * allowable_stress * joint_efficiency - 0.6 * pressure_mpa)
        + 3.0)
        .max(6.0);
    let mawp_bara = 2.0 * allowable_stress * joint_efficiency * (wall_thickness - 3.0)
        / (inner_diameter + 0.6 * (wall_thickness - 3.0))
        * 10.0;
    Casing {
        inner_diameter,
        outer_diameter: inner_diameter + 2.0 * wall_thickness,
        wall_thickness,
        mawp_bara,
    }
}

/// The driver margin's three bands.
fn driver_margin(absorbed_power_kw: f64) -> f64 {
    if absorbed_power_kw < 22.0 {
        DRIVER_MARGIN_SMALL
    } else if absorbed_power_kw <= 55.0 {
        DRIVER_MARGIN_MEDIUM
    } else {
        DRIVER_MARGIN_LARGE
    }
}

/// `selectDriverRating`: the smallest candidate covering the required power.
fn select_driver_rating(required_kw: f64, candidates: &[f64]) -> Option<f64> {
    if !(required_kw.is_finite() && required_kw > 0.0) {
        return None;
    }
    candidates
        .iter()
        .copied()
        .filter(|candidate| *candidate >= required_kw)
        .reduce(f64::min)
}

/// `calculateDriverSizing`'s intermediate: its own 25-rating IEC list, which the weights
/// and the module read before the API 610 assessment overwrites the field.
fn calculate_driver_sizing(absorbed_power_kw: f64) -> f64 {
    const CANDIDATES: [f64; 25] = [
        0.75, 1.1, 1.5, 2.2, 3.0, 4.0, 5.5, 7.5, 11.0, 15.0, 18.5, 22.0, 30.0, 37.0, 45.0, 55.0,
        75.0, 90.0, 110.0, 132.0, 160.0, 200.0, 250.0, 315.0, 400.0,
    ];
    let required = absorbed_power_kw * driver_margin(absorbed_power_kw);
    select_driver_rating(required, &CANDIDATES).unwrap_or(required)
}

/// The API 610 assessment's selected rating, from its longer list.
fn api610_driver_rating(absorbed_power_kw: f64) -> f64 {
    const CANDIDATES: [f64; 39] = [
        0.75, 1.1, 1.5, 2.2, 3.0, 4.0, 5.5, 7.5, 11.0, 15.0, 18.5, 22.0, 30.0, 37.0, 45.0, 55.0,
        75.0, 90.0, 110.0, 132.0, 160.0, 200.0, 250.0, 315.0, 400.0, 500.0, 630.0, 800.0, 1000.0,
        1250.0, 1600.0, 2000.0, 2500.0, 3150.0, 4000.0, 5000.0, 6300.0, 8000.0, 10_000.0,
    ];
    let required = absorbed_power_kw * driver_margin(absorbed_power_kw);
    select_driver_rating(required, &CANDIDATES).unwrap_or(required)
}

/// `calculateNozzleSizes`, in the class's inches.
fn nozzle_sizes(flow_m3_hr: f64) -> (f64, f64) {
    const STANDARD: [f64; 11] = [1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0];
    let flow_m3_s = flow_m3_hr / 3600.0;
    let suction_diameter_mm = 2.0 * (flow_m3_s / 2.5 / std::f64::consts::PI).sqrt() * 1000.0;
    let discharge_diameter_mm = 2.0 * (flow_m3_s / 5.0 / std::f64::consts::PI).sqrt() * 1000.0;
    let mut suction = 4.0;
    for size in STANDARD {
        if size * 25.4 >= suction_diameter_mm {
            suction = size;
            break;
        }
    }
    let mut discharge = 3.0;
    for size in STANDARD {
        if size * 25.4 >= discharge_diameter_mm {
            discharge = size;
            break;
        }
    }
    if discharge >= suction {
        if let Some(index) = STANDARD.iter().position(|size| *size == suction) {
            if index > 0 {
                discharge = STANDARD[index - 1];
            }
        }
    }
    (suction, discharge)
}

struct Weights {
    casing: f64,
    motor: f64,
    baseplate: f64,
    piping: f64,
    electro_instrument: f64,
    structural_steel: f64,
    total: f64,
}

/// `calculateWeights`, in kilograms.
#[allow(clippy::too_many_arguments)]
fn weights(
    impeller_diameter: f64,
    impeller_width: f64,
    shaft_diameter: f64,
    casing_wall_thickness: f64,
    pump_type: &str,
    number_of_stages: u32,
    driver_power_kw: f64,
) -> Weights {
    let impeller_volume = std::f64::consts::PI
        * (impeller_diameter / 2000.0).powi(2)
        * (impeller_width / 1000.0)
        * 0.3;
    let impeller_weight = impeller_volume * STEEL_DENSITY;

    let shaft_length = if pump_type == "overhung" {
        500.0
    } else {
        1000.0
    };
    let shaft_volume =
        std::f64::consts::PI * (shaft_diameter / 2000.0).powi(2) * (shaft_length / 1000.0);
    let shaft_weight = shaft_volume * STEEL_DENSITY;

    let casing_od = (impeller_diameter * 1.15 + 2.0 * casing_wall_thickness) / 1000.0;
    let casing_id = impeller_diameter * 1.15 / 1000.0;
    let casing_length = (impeller_width + 100.0) / 1000.0 * f64::from(number_of_stages);
    let casing_volume =
        std::f64::consts::PI / 4.0 * (casing_od.powi(2) - casing_id.powi(2)) * casing_length;
    let casing_weight = casing_volume * STEEL_DENSITY * 2.0;

    let driver_power = driver_power_kw;
    let motor_weight = if driver_power <= 0.0 {
        50.0
    } else if driver_power <= 10.0 {
        20.0 + driver_power * 3.0
    } else if driver_power <= 100.0 {
        50.0 + driver_power * 2.5
    } else {
        300.0 + driver_power * 1.5
    };
    let baseplate_weight = (casing_weight + motor_weight) * 0.6;
    let seal_weight = 5.0 + shaft_diameter * 0.5;
    let coupling_weight = 10.0 + driver_power * 0.1;

    let empty_weight = casing_weight
        + impeller_weight
        + shaft_weight
        + motor_weight
        + baseplate_weight
        + seal_weight
        + coupling_weight;
    let piping = empty_weight * 0.15;
    let electro_instrument = empty_weight * 0.10;
    let structural_steel = empty_weight * 0.05;

    Weights {
        casing: casing_weight,
        motor: motor_weight,
        baseplate: baseplate_weight,
        piping,
        electro_instrument,
        structural_steel,
        total: empty_weight + piping + electro_instrument + structural_steel,
    }
}

/// `calculateModuleDimensions`, in metres.
fn module_dimensions(impeller_diameter: f64, driver_power_kw: f64) -> (f64, f64, f64) {
    let module_width = 1.0_f64.max(impeller_diameter / 500.0 + 0.5);
    let pump_length = impeller_diameter / 1000.0 * 2.0;
    let motor_length = 0.3 + driver_power_kw.sqrt() * 0.05;
    let module_length = pump_length + motor_length + 1.0;
    let module_height = 0.5_f64.max(impeller_diameter / 1000.0 + 0.3);
    (module_length, module_width, module_height)
}
