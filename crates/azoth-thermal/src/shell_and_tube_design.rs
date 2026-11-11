//! `thermal.shell_and_tube_design` - the detailed TEMA/ASME VIII/NACE design loop.
//!
//! Transcribed from `ShellAndTubeDesignCalculator.java` at NeqSim `f0c7436c69`. The pipeline
//! is single-pass: TEMA configuration, material lookup, tube bundle, shell dimensions,
//! baffles, UHX-13, UG-37, UG-27, UG-99, the optional NACE assessment, the optional thermal
//! and vibration block, weights and costs. Nothing iterates.
//!
//! **The geometry is millimetres**; the ASME formulas convert to MPa internally and the
//! thermal core is fed metres. **Two upstream defects are reproduced**: the NACE test is
//! `"Yes".equalsIgnoreCase` against a `true`/`false` column, and the price test looks for
//! `STAINLESS` in a specification string that only ever holds `ASME-II-A`.
//!
//! Spec: `specs/calcs/thermal/shell_and_tube_design.toml`.

use azoth_core::units::{
    Area, DynamicViscosity, Frequency, HeatTransfer, Length, MassDensity, MassRate, Pressure,
    SpecificHeatCapacity, ThermalConductivity, ThermodynamicTemperature, Velocity, millimeters,
    pascals, square_meters, watts_per_meter_kelvin,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::bell_delaware;
use crate::results::{AppliedStandardRow, NaceIssueRow, ShellAndTubeDesignResult};
use crate::shell_and_tube_rating::{RatingRequest, ShellSideMethod, evaluate};
use crate::shell_and_tube_vibration::{VibrationRequest, screen};
use crate::spec_gen;
use crate::tema_standard::{
    FrontHead, PitchPattern, RearHead, Shell, TemaClass, TubeSize, class_limits,
    configuration_cost_factor, max_unsupported_span, min_baffle_spacing, tube_outside_diameter,
};
use crate::tube_materials;

/// The class's fixed tube wall, mm.
const TUBE_WALL_MM: f64 = 2.108;
/// The class's fixed shell corrosion allowance, mm.
const CORROSION_SHELL_MM: f64 = 3.175;
/// The class's fixed tube corrosion allowance, mm.
const CORROSION_TUBE_MM: f64 = 0.0;
/// The class's fixed pitch ratio.
const PITCH_RATIO: f64 = 1.25;
/// The class's fixed baffle cut.
const BAFFLE_CUT: f64 = 0.25;
/// The class's fixed shell joint efficiency inside the shell-wall formula.
const SHELL_WALL_JOINT_EFFICIENCY: f64 = 0.85;
/// Steel density, kg/m3.
const STEEL_DENSITY: f64 = 7850.0;
/// Water density, kg/m3.
const WATER_DENSITY: f64 = 1000.0;
/// The NACE sour-service threshold, bar.
const SOUR_THRESHOLD_BAR: f64 = 0.003;
/// The class's fixed hardness, HV.
const HARDNESS_HV: f64 = 248.0;
/// Carbon-steel price, USD/kg.
const CARBON_PRICE: f64 = 1.5;
/// Stainless price, USD/kg - unreachable, as the class reads the specification string.
const STAINLESS_PRICE: f64 = 5.0;
/// The standard shell sizes, mm.
const SHELL_SIZES_MM: [f64; 38] = [
    168.3, 219.1, 273.0, 323.9, 355.6, 406.4, 457.2, 508.0, 558.8, 609.6, 660.4, 711.2, 762.0,
    812.8, 863.6, 914.4, 965.2, 1016.0, 1066.8, 1117.6, 1168.4, 1219.2, 1270.0, 1320.8, 1371.6,
    1422.4, 1473.2, 1524.0, 1574.8, 1625.6, 1676.4, 1727.2, 1778.0, 1828.8, 1879.6, 1930.4, 1981.2,
    2032.0,
];
/// The standard plate thicknesses, mm.
const PLATE_THICKNESSES_MM: [f64; 14] = [
    3.175, 4.763, 6.35, 7.938, 9.525, 12.7, 15.875, 19.05, 22.225, 25.4, 31.75, 38.1, 44.45, 50.8,
];

/// The design inputs, in SI magnitudes for the geometry (`mm`) and bar/MPa where the class
/// speaks them.
#[derive(Debug, Clone)]
pub struct ShellAndTubeDesignRequest {
    pub tema_designation: String,
    pub tema_class: TemaClass,
    pub required_area: f64,
    pub shell_side_pressure: f64,
    pub tube_side_pressure: f64,
    pub design_temperature: f64,
    pub shell_material_grade: String,
    pub tube_material_grade: String,
    pub tube_size: TubeSize,
    pub tube_length: f64,
    pub tube_passes: f64,
    pub shell_joint_efficiency: f64,
    pub nozzle_diameter: f64,
    pub nozzle_wall_thickness: f64,
    pub sour_service_assessment: bool,
    pub h2s_partial_pressure: Option<f64>,
    pub fouling_tube: Option<f64>,
    pub fouling_shell: Option<f64>,
    pub shell_side_method: Option<ShellSideMethod>,
    pub tube_side_density: Option<f64>,
    pub tube_side_viscosity: Option<f64>,
    pub tube_side_cp: Option<f64>,
    pub tube_side_conductivity: Option<f64>,
    pub tube_side_mass_flow: Option<f64>,
    pub shell_side_density: Option<f64>,
    pub shell_side_viscosity: Option<f64>,
    pub shell_side_cp: Option<f64>,
    pub shell_side_conductivity: Option<f64>,
    pub shell_side_mass_flow: Option<f64>,
}

impl Default for ShellAndTubeDesignRequest {
    fn default() -> Self {
        Self {
            tema_designation: "AES".to_string(),
            tema_class: TemaClass::R,
            required_area: 100.0,
            shell_side_pressure: 10.0,
            tube_side_pressure: 10.0,
            design_temperature: 423.15,
            shell_material_grade: "SA-516-70".to_string(),
            tube_material_grade: "SA-179".to_string(),
            tube_size: TubeSize::Tube3_4Inch,
            tube_length: 6096.0,
            tube_passes: 2.0,
            shell_joint_efficiency: 0.85,
            nozzle_diameter: 150.0,
            nozzle_wall_thickness: 7.11,
            sour_service_assessment: false,
            h2s_partial_pressure: None,
            fouling_tube: None,
            fouling_shell: None,
            shell_side_method: None,
            tube_side_density: None,
            tube_side_viscosity: None,
            tube_side_cp: None,
            tube_side_conductivity: None,
            tube_side_mass_flow: None,
            shell_side_density: None,
            shell_side_viscosity: None,
            shell_side_cp: None,
            shell_side_conductivity: None,
            shell_side_mass_flow: None,
        }
    }
}

/// One material's properties, as the class loads them.
#[derive(Debug, Clone)]
struct Material {
    name: String,
    allowable_mpa: f64,
    smys_mpa: f64,
    conductivity: f64,
}

fn load_material(grade: &str, is_shell: bool) -> Material {
    match tube_materials::lookup(grade) {
        Some(row) => Material {
            name: format!("{} {grade}", row.specification),
            allowable_mpa: row.allowable_stress_mpa,
            smys_mpa: row.smys_mpa,
            conductivity: row.thermal_conductivity_wpm_k,
        },
        None => Material {
            name: if is_shell {
                "Carbon Steel SA-516-70".to_string()
            } else {
                "Carbon Steel SA-179".to_string()
            },
            allowable_mpa: if is_shell { 138.0 } else { 117.0 },
            smys_mpa: if is_shell { 260.0 } else { 228.0 },
            conductivity: 52.0,
        },
    }
}

fn parse_designation(designation: &str) -> Result<(FrontHead, Shell, RearHead)> {
    let letters: Vec<char> = designation.chars().collect();
    if letters.len() != 3 {
        return Err(AzothError::invalid_input(
            "tema_designation",
            format!("{designation} is not three letters"),
        ));
    }
    let front: FrontHead = letters[0].to_string().parse()?;
    let shell: Shell = letters[1].to_string().parse()?;
    let rear: RearHead = letters[2].to_string().parse()?;
    Ok((front, shell, rear))
}

fn rear_clearance_mm(rear: RearHead) -> f64 {
    match rear {
        RearHead::T => 90.0,
        RearHead::S => 35.0,
        RearHead::U => 15.0,
        _ => 10.0,
    }
}

fn round_shell_size(diameter: f64) -> f64 {
    SHELL_SIZES_MM
        .iter()
        .copied()
        .find(|size| *size >= diameter)
        .unwrap_or(SHELL_SIZES_MM[SHELL_SIZES_MM.len() - 1])
}

fn round_plate(thickness: f64) -> f64 {
    PLATE_THICKNESSES_MM
        .iter()
        .copied()
        .find(|size| *size >= thickness)
        .unwrap_or(PLATE_THICKNESSES_MM[PLATE_THICKNESSES_MM.len() - 1])
}

fn allowable_stress(shell: &Material, tube: &Material, material: &str, temperature_c: f64) -> f64 {
    let mut base = if material == shell.name {
        shell.allowable_mpa
    } else if material == tube.name {
        tube.allowable_mpa
    } else {
        137.9
    };
    if temperature_c > 400.0 {
        base *= 0.55;
    } else if temperature_c > 300.0 {
        base *= 0.70;
    } else if temperature_c > 200.0 {
        base *= 0.85;
    }
    base
}

fn whole(value: f64, field: &str) -> Result<i64> {
    if !value.is_finite() || value <= 0.0 || (value - value.round()).abs() > 0.0 {
        return Err(AzothError::invalid_input(
            field,
            format!("{value} is not a positive whole count"),
        ));
    }
    Ok(value.round() as i64)
}

fn group_present(values: &[Option<f64>], side: &str) -> Result<bool> {
    let present = values.iter().filter(|value| value.is_some()).count();
    if present == 0 {
        return Ok(false);
    }
    if present != values.len() {
        return Err(AzothError::invalid_input(
            side,
            format!("the {side} thermal properties are all-or-none"),
        ));
    }
    Ok(true)
}

/// The detailed design.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a designation that is not three declared letters, a
/// non-integral pass count, a partial thermal side, or a state the ASME formulas cannot take.
#[allow(clippy::too_many_lines)] // one transcription of the class's pipeline
pub fn rate(request: &ShellAndTubeDesignRequest) -> Result<ShellAndTubeDesignResult> {
    let spec = &spec_gen::SHELL_AND_TUBE_DESIGN_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "required_area" => Some(request.required_area),
            "shell_side_pressure" => Some(request.shell_side_pressure),
            "tube_side_pressure" => Some(request.tube_side_pressure),
            "design_temperature" => Some(request.design_temperature),
            "tube_length" => Some(request.tube_length),
            "tube_passes" => Some(request.tube_passes),
            "shell_joint_efficiency" => Some(request.shell_joint_efficiency),
            "nozzle_diameter" => Some(request.nozzle_diameter),
            "nozzle_wall_thickness" => Some(request.nozzle_wall_thickness),
            "h2s_partial_pressure" => request.h2s_partial_pressure,
            "fouling_tube" => request.fouling_tube,
            "fouling_shell" => request.fouling_shell,
            _ => None,
        },
        &mut warnings,
    )?;

    let (front, shell_letter, rear) = parse_designation(&request.tema_designation)?;
    let passes = whole(request.tube_passes, "tube_passes")?;
    if passes < 1 {
        return Err(AzothError::invalid_input(
            "tube_passes",
            "at least one pass",
        ));
    }
    if !request.required_area.is_finite() || request.required_area <= 0.0 {
        return Err(AzothError::invalid_input(
            "required_area",
            "the required area must be positive",
        ));
    }
    let design_temperature_c = request.design_temperature - 273.15;

    let shell_material = load_material(&request.shell_material_grade, true);
    let tube_material = load_material(&request.tube_material_grade, false);

    // Tube bundle.
    let tube_od = tube_outside_diameter(request.tube_size) * 1000.0;
    let tube_id = tube_od - 2.0 * TUBE_WALL_MM;
    let tube_pitch = tube_od * PITCH_RATIO;
    let area_per_tube = std::f64::consts::PI * tube_od / 1000.0 * request.tube_length / 1000.0;
    let mut tube_count = (request.required_area / area_per_tube).ceil() as i64;
    if tube_count % passes != 0 {
        tube_count = (tube_count / passes + 1) * passes;
    }
    let actual_area = tube_count as f64 * area_per_tube;

    // Shell dimensions.
    let layout_angle = PitchPattern::Triangular30;
    let area_per_tube_layout = if matches!(
        layout_angle,
        PitchPattern::Triangular30 | PitchPattern::Triangular60
    ) {
        tube_pitch * tube_pitch * 3.0_f64.sqrt() / 2.0
    } else {
        tube_pitch * tube_pitch
    };
    let pass_lane_factor = 1.0 - 0.05 * (passes as f64 - 1.0);
    let bundle_area = tube_count as f64 * area_per_tube_layout / pass_lane_factor;
    let bundle_diameter = 2.0 * (bundle_area / std::f64::consts::PI).sqrt();
    let shell_inside_diameter =
        round_shell_size(bundle_diameter + 2.0 * tube_od + rear_clearance_mm(rear));
    let design_pressure_mpa = request.shell_side_pressure.max(1.0) * 1.1 * 0.1;
    let shell_allowable = allowable_stress(
        &shell_material,
        &tube_material,
        &shell_material.name,
        design_temperature_c,
    );
    let shell_radius = shell_inside_diameter / 2.0;
    let mut shell_wall = (design_pressure_mpa * shell_radius)
        / (shell_allowable * SHELL_WALL_JOINT_EFFICIENCY - 0.6 * design_pressure_mpa)
        + CORROSION_SHELL_MM;
    let (_, min_wall_m, _) = class_limits(request.tema_class);
    shell_wall = shell_wall.max(min_wall_m * 1000.0);
    shell_wall = round_plate(shell_wall);
    let shell_outside_diameter = shell_inside_diameter + 2.0 * shell_wall;

    // Baffles.
    let minimum_spacing =
        min_baffle_spacing(shell_inside_diameter / 1000.0, request.tema_class) * 1000.0;
    let maximum_spacing = shell_inside_diameter;
    let maximum_span =
        max_unsupported_span(tube_od / 1000.0, Some(tube_material.name.as_str())) * 1000.0;
    let target_spacing = shell_inside_diameter * 0.4;
    let mut baffle_spacing = minimum_spacing.max(target_spacing.min(maximum_spacing));
    baffle_spacing = baffle_spacing.min(maximum_span);
    let effective_length = request.tube_length - 2.0 * baffle_spacing;
    let baffle_count = ((effective_length / baffle_spacing).floor() as i64).max(1);

    // UHX-13 tubesheet.
    let pressure_mpa = request.shell_side_pressure.max(request.tube_side_pressure) * 0.1;
    let ligament = (1.0 - tube_od / (tube_od * PITCH_RATIO)).max(0.25);
    let f_factor = 0.785 / ligament;
    let mut tubesheet = shell_inside_diameter * (f_factor * pressure_mpa / shell_allowable).sqrt()
        + CORROSION_SHELL_MM;
    let min_tubesheet = match request.tema_class {
        TemaClass::C => 15.875,
        TemaClass::B => 12.7,
        TemaClass::R => 19.05,
    };
    tubesheet = tubesheet.max(min_tubesheet);
    let mut standards = vec![AppliedStandardRow {
        name: "ASME VIII Div.1 UHX-13 (Tubesheet)".to_string(),
    }];

    // UG-37 nozzle reinforcement.
    let joint = request.shell_joint_efficiency;
    let t_r = (pressure_mpa * shell_inside_diameter / 2.0)
        / (shell_allowable * joint - 0.6 * pressure_mpa);
    let required_area = request.nozzle_diameter * t_r;
    let t_shell_net = shell_wall - CORROSION_SHELL_MM;
    let a1 = (t_shell_net - t_r).max(0.0) * request.nozzle_diameter * 2.0;
    let t_rn = (pressure_mpa * request.nozzle_diameter / 2.0)
        / (shell_allowable * joint - 0.6 * pressure_mpa);
    let a2 = (request.nozzle_wall_thickness - t_rn).max(0.0) * 2.0 * 2.5 * t_shell_net;
    let available_area = a1 + a2;
    let nozzle_adequate = available_area >= required_area;
    standards.push(AppliedStandardRow {
        name: "ASME VIII Div.1 UG-37 (Nozzle Reinforcement)".to_string(),
    });

    // UG-27 MAWP.
    let mawp_shell_mpa =
        (shell_allowable * joint * t_shell_net) / (shell_radius + 0.6 * t_shell_net);
    let tube_allowable = allowable_stress(
        &shell_material,
        &tube_material,
        &tube_material.name,
        design_temperature_c,
    );
    let t_tube = TUBE_WALL_MM - CORROSION_TUBE_MM;
    let tube_radius = (tube_od - 2.0 * TUBE_WALL_MM) / 2.0;
    let mawp_tube_mpa = (tube_allowable * t_tube) / (tube_radius + 0.6 * t_tube);
    let mawp_shell = mawp_shell_mpa * 10.0;
    let mawp_tube = mawp_tube_mpa * 10.0;
    standards.push(AppliedStandardRow {
        name: "ASME VIII Div.1 UG-27 (MAWP)".to_string(),
    });

    // UG-99 hydrotest.
    let stress_ratio = if design_temperature_c > 100.0 {
        allowable_stress(&shell_material, &tube_material, &shell_material.name, 20.0)
            / allowable_stress(
                &shell_material,
                &tube_material,
                &shell_material.name,
                design_temperature_c,
            )
    } else {
        1.0
    };
    let hydro_shell = 1.3 * mawp_shell * stress_ratio;
    let hydro_tube = 1.3 * mawp_tube * stress_ratio;
    standards.push(AppliedStandardRow {
        name: "ASME VIII Div.1 UG-99 (Hydrostatic Test)".to_string(),
    });

    // NACE assessment.
    let mut nace_issues = Vec::new();
    let mut sour_service_required = false;
    let mut h2s_output = None;
    let mut shell_compliant = false;
    let mut tube_compliant = false;
    if request.sour_service_assessment {
        let h2s = request.h2s_partial_pressure.unwrap_or(0.0);
        h2s_output = Some(h2s);
        sour_service_required = h2s > SOUR_THRESHOLD_BAR;
        shell_compliant = tube_materials::lookup(&request.shell_material_grade)
            .is_some_and(tube_materials::TubeMaterial::is_nace_compliant);
        tube_compliant = tube_materials::lookup(&request.tube_material_grade)
            .is_some_and(tube_materials::TubeMaterial::is_nace_compliant);
        if !sour_service_required {
            standards.push(AppliedStandardRow {
                name: "NACE MR0175/ISO 15156 (Not Required - H2S pp below threshold)".to_string(),
            });
        } else {
            if !shell_compliant {
                nace_issues.push(NaceIssueRow {
                    message: format!(
                        "Shell material {} is not NACE MR0175 compliant for sour service",
                        request.shell_material_grade
                    ),
                });
            }
            if !tube_compliant {
                nace_issues.push(NaceIssueRow {
                    message: format!(
                        "Tube material {} is not NACE MR0175 compliant for sour service",
                        request.tube_material_grade
                    ),
                });
            }
            if HARDNESS_HV > 248.0 {
                nace_issues.push(NaceIssueRow {
                    message: format!("Shell hardness {HARDNESS_HV} HV exceeds 248 HV limit"),
                });
            }
            if HARDNESS_HV > 248.0 {
                nace_issues.push(NaceIssueRow {
                    message: format!("Tube hardness {HARDNESS_HV} HV exceeds 248 HV limit"),
                });
            }
            standards.push(AppliedStandardRow {
                name: "NACE MR0175/ISO 15156 (Sour Service Assessment)".to_string(),
            });
            if !nace_issues.is_empty() {
                standards.push(AppliedStandardRow {
                    name: "NORSOK M-001 Rev.6 (Material Hardness Limits)".to_string(),
                });
            }
        }
    }

    // The optional thermal and vibration block.
    let tube_present = group_present(
        &[
            request.tube_side_density,
            request.tube_side_viscosity,
            request.tube_side_cp,
            request.tube_side_conductivity,
            request.tube_side_mass_flow,
        ],
        "tube_side",
    )?;
    let shell_present = group_present(
        &[
            request.shell_side_density,
            request.shell_side_viscosity,
            request.shell_side_cp,
            request.shell_side_conductivity,
            request.shell_side_mass_flow,
        ],
        "shell_side",
    )?;
    if tube_present != shell_present {
        return Err(AzothError::invalid_input(
            "thermal",
            "the thermal block needs both sides stated",
        ));
    }
    let mut thermal_fields: Option<ThermalBlock> = None;
    let mut vibration_fields: Option<VibrationBlock> = None;
    if tube_present && shell_present {
        let fouling_tube = request.fouling_tube.unwrap_or(1.76e-4);
        let fouling_shell = request.fouling_shell.unwrap_or(1.76e-4);
        let method = request.shell_side_method.unwrap_or(ShellSideMethod::Kern);
        let rating_request = RatingRequest {
            tube_outside_diameter: tube_od / 1000.0,
            tube_inside_diameter: tube_id / 1000.0,
            tube_length: request.tube_length / 1000.0,
            tube_count: tube_count as f64,
            tube_passes: passes as f64,
            tube_pitch: tube_pitch / 1000.0,
            triangular_pitch: true,
            shell_inside_diameter: shell_inside_diameter / 1000.0,
            baffle_spacing: baffle_spacing / 1000.0,
            baffle_count: baffle_count as f64,
            baffle_cut: BAFFLE_CUT,
            tube_wall_conductivity: tube_material.conductivity,
            tube_density: request.tube_side_density.unwrap_or(0.0),
            tube_viscosity: request.tube_side_viscosity.unwrap_or(0.0),
            tube_cp: request.tube_side_cp.unwrap_or(0.0),
            tube_conductivity: request.tube_side_conductivity.unwrap_or(0.0),
            tube_mass_flow: request.tube_side_mass_flow.unwrap_or(0.0),
            shell_density: request.shell_side_density.unwrap_or(0.0),
            shell_viscosity: request.shell_side_viscosity.unwrap_or(0.0),
            shell_cp: request.shell_side_cp.unwrap_or(0.0),
            shell_conductivity: request.shell_side_conductivity.unwrap_or(0.0),
            shell_mass_flow: request.shell_side_mass_flow.unwrap_or(0.0),
            fouling_tube,
            fouling_shell,
            shell_side_method: method,
            ..RatingRequest::default()
        };
        let rating = evaluate(&rating_request)?;
        standards.push(AppliedStandardRow {
            name: "Gnielinski correlation (Tube-side HTC)".to_string(),
        });
        standards.push(AppliedStandardRow {
            name: if method == ShellSideMethod::BellDelaware {
                "Bell-Delaware method (Shell-side HTC)"
            } else {
                "Kern method (Shell-side HTC)"
            }
            .to_string(),
        });
        thermal_fields = Some(ThermalBlock::from(&rating));

        let area = bell_delaware::crossflow_area(
            shell_inside_diameter / 1000.0,
            baffle_spacing / 1000.0,
            tube_od / 1000.0,
            tube_pitch / 1000.0,
        );
        let shell_density = request.shell_side_density.unwrap_or(0.0);
        let velocity = if area > 0.0 {
            request.shell_side_mass_flow.unwrap_or(0.0) / (shell_density * area)
        } else {
            0.0
        };
        let sonic: f64 = (1.4_f64 * 8314.0 * 300.0 / 29.0).sqrt();
        let damping = if shell_density > 100.0 { 0.03 } else { 0.01 };
        let vibration = screen(&VibrationRequest {
            tube_outside_diameter: tube_od / 1000.0,
            tube_inside_diameter: tube_id / 1000.0,
            unsupported_span: baffle_spacing / 1000.0,
            tube_pitch: tube_pitch / 1000.0,
            tube_material_e: 200e9,
            tube_density: STEEL_DENSITY,
            crossflow_velocity: velocity,
            shell_fluid_density: shell_density,
            tube_fluid_density: request.tube_side_density.unwrap_or(0.0),
            shell_inside_diameter: shell_inside_diameter / 1000.0,
            sonic_velocity: sonic,
            damping_ratio: damping,
            triangular_pitch: true,
        });
        standards.push(AppliedStandardRow {
            name: if vibration.passed {
                "TEMA RCB-4.6 (Vibration Screening - PASS)"
            } else {
                "TEMA RCB-4.6 (Vibration Screening - FAIL)"
            }
            .to_string(),
        });
        vibration_fields = Some(VibrationBlock::from(&vibration));
    }

    // Weights.
    let shell_volume = std::f64::consts::PI / 4.0
        * ((shell_outside_diameter / 1000.0).powi(2) - (shell_inside_diameter / 1000.0).powi(2))
        * request.tube_length
        / 1000.0;
    let shell_weight = shell_volume * STEEL_DENSITY;
    let tube_volume = std::f64::consts::PI / 4.0
        * ((tube_od / 1000.0).powi(2) - (tube_id / 1000.0).powi(2))
        * request.tube_length
        / 1000.0;
    let tube_weight = tube_count as f64 * tube_volume * STEEL_DENSITY;
    let estimated_tubesheet = estimate_tubesheet_thickness(
        shell_inside_diameter,
        design_temperature_c,
        request,
        &shell_material,
        &tube_material,
    );
    let tubesheet_volume = std::f64::consts::PI / 4.0
        * (shell_outside_diameter / 1000.0).powi(2)
        * estimated_tubesheet
        / 1000.0;
    let tubesheet_count = if rear == RearHead::U { 1.0 } else { 2.0 };
    let tubesheet_weight = tubesheet_count * tubesheet_volume * STEEL_DENSITY;
    let head_thickness = shell_wall * 0.5;
    let head_volume = 2.0 * std::f64::consts::PI / 3.0
        * (((shell_outside_diameter / 2.0 + head_thickness) / 1000.0).powi(3)
            - ((shell_outside_diameter / 2.0) / 1000.0).powi(3));
    let head_weight = 2.0 * head_volume * STEEL_DENSITY;
    let baffle_area =
        std::f64::consts::PI / 4.0 * (shell_inside_diameter / 1000.0).powi(2) * (1.0 - BAFFLE_CUT);
    let baffle_volume = baffle_area * 6.35 / 1000.0;
    let baffle_weight = baffle_count as f64 * baffle_volume * STEEL_DENSITY;
    let total_dry = shell_weight + tube_weight + tubesheet_weight + head_weight + baffle_weight;
    let shell_side_volume =
        std::f64::consts::PI / 4.0 * (shell_inside_diameter / 1000.0).powi(2) * request.tube_length
            / 1000.0
            - tube_count as f64 * std::f64::consts::PI / 4.0
                * (tube_od / 1000.0).powi(2)
                * request.tube_length
                / 1000.0;
    let tube_side_volume = tube_count as f64 * std::f64::consts::PI / 4.0
        * (tube_id / 1000.0).powi(2)
        * request.tube_length
        / 1000.0;
    let operating_weight = total_dry + (shell_side_volume + tube_side_volume) * WATER_DENSITY;

    // Costs.
    let shell_price = if shell_material.name.to_uppercase().contains("STAINLESS") {
        STAINLESS_PRICE
    } else {
        CARBON_PRICE
    };
    let tube_price = if tube_material.name.to_uppercase().contains("STAINLESS") {
        STAINLESS_PRICE
    } else {
        CARBON_PRICE
    };
    let material_cost = shell_weight * shell_price
        + tube_weight * tube_price
        + (tubesheet_weight + head_weight + baffle_weight) * shell_price;
    let (class_cost_factor, _, _) = class_limits(request.tema_class);
    let fabrication_factor =
        2.5 * configuration_cost_factor(front, shell_letter, rear) * class_cost_factor;
    let fabrication_cost = material_cost * fabrication_factor;
    let total_cost = material_cost + fabrication_cost;

    let thermal = thermal_fields;
    let vibration = vibration_fields;

    Ok(ShellAndTubeDesignResult {
        tema_designation: request.tema_designation.clone(),
        tema_class: class_name(request.tema_class).to_string(),
        shell_inside_diameter: millimeters(shell_inside_diameter),
        shell_outside_diameter: millimeters(shell_outside_diameter),
        shell_wall_thickness: millimeters(shell_wall),
        shell_material: shell_material.name.clone(),
        shell_material_grade: request.shell_material_grade.clone(),
        shell_allowable_stress: pascals(shell_allowable * 1e6),
        shell_smys: pascals(shell_material.smys_mpa * 1e6),
        shell_joint_efficiency: joint,
        tube_count: tube_count as f64,
        tube_outer_diameter: millimeters(tube_od),
        tube_wall_thickness: millimeters(TUBE_WALL_MM),
        tube_length: millimeters(request.tube_length),
        tube_passes: passes as f64,
        tube_pitch: millimeters(tube_pitch),
        tube_pitch_pattern: "Triangular 30°".to_string(),
        tube_material: tube_material.name.clone(),
        tube_material_grade: request.tube_material_grade.clone(),
        tube_allowable_stress: pascals(tube_allowable * 1e6),
        tube_thermal_conductivity: watts_per_meter_kelvin(tube_material.conductivity),
        baffle_type: "Single Segmental".to_string(),
        baffle_count: baffle_count as f64,
        baffle_spacing: millimeters(baffle_spacing),
        baffle_cut: BAFFLE_CUT,
        required_area: square_meters(request.required_area),
        actual_area: square_meters(actual_area),
        area_margin: (actual_area - request.required_area) / request.required_area,
        tubesheet_thickness_uhx: millimeters(tubesheet),
        mawp_shell_side: pascals(mawp_shell * 1e5),
        mawp_tube_side: pascals(mawp_tube * 1e5),
        hydro_test_pressure_shell: pascals(hydro_shell * 1e5),
        hydro_test_pressure_tube: pascals(hydro_tube * 1e5),
        nozzle_reinforcement_area: square_meters(available_area / 1e6),
        nozzle_reinforcement_adequate: nozzle_adequate,
        sour_service_required,
        h2s_partial_pressure: h2s_output.map(|value| pascals(value * 1e5)),
        shell_nace_compliant: shell_compliant,
        tube_nace_compliant: tube_compliant,
        shell_max_hardness: HARDNESS_HV,
        tube_max_hardness: HARDNESS_HV,
        nace_issues,
        applied_standards: standards,
        tube_side_heat_transfer_coefficient: thermal.as_ref().map(|t| t.tube_htc),
        shell_side_heat_transfer_coefficient: thermal.as_ref().map(|t| t.shell_htc),
        overall_heat_transfer_coefficient: thermal.as_ref().map(|t| t.overall_u),
        tube_side_pressure_drop: thermal.as_ref().map(|t| t.tube_dp),
        shell_side_pressure_drop: thermal.as_ref().map(|t| t.shell_dp),
        tube_side_velocity: thermal.as_ref().map(|t| t.tube_velocity),
        shell_side_velocity: thermal.as_ref().map(|t| t.shell_velocity),
        tube_side_reynolds: thermal.as_ref().map(|t| t.tube_re),
        shell_side_reynolds: thermal.as_ref().map(|t| t.shell_re),
        effective_tube_inside_diameter: thermal.as_ref().map(|t| t.effective_id),
        outside_heat_transfer_area: thermal.as_ref().map(|t| t.outside_area),
        tube_fouling_layer_resistance_outside: thermal.as_ref().map(|t| t.deposit_resistance),
        jc: thermal.as_ref().map(|t| t.jc),
        jl: thermal.as_ref().map(|t| t.jl),
        jb: thermal.as_ref().map(|t| t.jb),
        js: thermal.as_ref().map(|t| t.js),
        jr: thermal.as_ref().map(|t| t.jr),
        passed: vibration.as_ref().map(|v| v.passed),
        natural_frequency: vibration.as_ref().map(|v| v.natural_frequency),
        vortex_shedding_frequency: vibration.as_ref().map(|v| v.vortex_frequency),
        vortex_shedding_ratio: vibration.as_ref().map(|v| v.vortex_ratio),
        vortex_shedding_critical: vibration.as_ref().map(|v| v.vortex_critical),
        critical_velocity: vibration.as_ref().map(|v| v.critical_velocity),
        velocity_ratio: vibration.as_ref().map(|v| v.velocity_ratio),
        fluid_elastic_critical: vibration.as_ref().map(|v| v.fluid_elastic_critical),
        acoustic_frequency: vibration.as_ref().map(|v| v.acoustic_frequency),
        acoustic_ratio: vibration.as_ref().map(|v| v.acoustic_ratio),
        acoustic_critical: vibration.as_ref().map(|v| v.acoustic_critical),
        shell_weight,
        tube_weight,
        tubesheet_weight,
        head_weight,
        baffle_weight,
        total_dry_weight: total_dry,
        operating_weight,
        material_cost,
        fabrication_cost,
        total_cost,
        warnings,
    })
}

fn estimate_tubesheet_thickness(
    shell_inside_diameter: f64,
    design_temperature_c: f64,
    request: &ShellAndTubeDesignRequest,
    shell_material: &Material,
    tube_material: &Material,
) -> f64 {
    let pressure_mpa = request.shell_side_pressure.max(request.tube_side_pressure) * 0.1;
    let allowable = allowable_stress(
        shell_material,
        tube_material,
        &shell_material.name,
        design_temperature_c,
    );
    let thickness = shell_inside_diameter * (0.785 * pressure_mpa / allowable).sqrt() / 0.7
        + CORROSION_SHELL_MM;
    let minimum = match request.tema_class {
        TemaClass::C => 15.875,
        TemaClass::B => 12.7,
        TemaClass::R => 19.05,
    };
    thickness.max(minimum)
}

fn class_name(class: TemaClass) -> &'static str {
    match class {
        TemaClass::R => "R",
        TemaClass::C => "C",
        TemaClass::B => "B",
    }
}

struct ThermalBlock {
    tube_htc: HeatTransfer,
    shell_htc: HeatTransfer,
    overall_u: HeatTransfer,
    tube_dp: Pressure,
    shell_dp: Pressure,
    tube_velocity: Velocity,
    shell_velocity: Velocity,
    tube_re: f64,
    shell_re: f64,
    effective_id: Length,
    outside_area: Area,
    deposit_resistance: f64,
    jc: f64,
    jl: f64,
    jb: f64,
    js: f64,
    jr: f64,
}

impl From<&crate::results::ShellAndTubeRatingResult> for ThermalBlock {
    fn from(value: &crate::results::ShellAndTubeRatingResult) -> Self {
        Self {
            tube_htc: value.tube_side_heat_transfer_coefficient,
            shell_htc: value.shell_side_heat_transfer_coefficient,
            overall_u: value.overall_heat_transfer_coefficient,
            tube_dp: value.tube_side_pressure_drop,
            shell_dp: value.shell_side_pressure_drop,
            tube_velocity: value.tube_side_velocity,
            shell_velocity: value.shell_side_velocity,
            tube_re: value.tube_side_reynolds,
            shell_re: value.shell_side_reynolds,
            effective_id: value.effective_tube_inside_diameter,
            outside_area: value.outside_heat_transfer_area,
            deposit_resistance: value.tube_fouling_layer_resistance_outside,
            jc: value.jc,
            jl: value.jl,
            jb: value.jb,
            js: value.js,
            jr: value.jr,
        }
    }
}

struct VibrationBlock {
    passed: bool,
    natural_frequency: Frequency,
    vortex_frequency: Frequency,
    vortex_ratio: f64,
    vortex_critical: bool,
    critical_velocity: Velocity,
    velocity_ratio: f64,
    fluid_elastic_critical: bool,
    acoustic_frequency: Frequency,
    acoustic_ratio: f64,
    acoustic_critical: bool,
}

impl From<&crate::results::ShellAndTubeVibrationResult> for VibrationBlock {
    fn from(value: &crate::results::ShellAndTubeVibrationResult) -> Self {
        Self {
            passed: value.passed,
            natural_frequency: value.natural_frequency,
            vortex_frequency: value.vortex_shedding_frequency,
            vortex_ratio: value.vortex_shedding_ratio,
            vortex_critical: value.vortex_shedding_critical,
            critical_velocity: value.critical_velocity,
            velocity_ratio: value.velocity_ratio,
            fluid_elastic_critical: value.fluid_elastic_critical,
            acoustic_frequency: value.acoustic_frequency,
            acoustic_ratio: value.acoustic_ratio,
            acoustic_critical: value.acoustic_critical,
        }
    }
}

/// The flat entry point the binding and the spec's case runner call.
///
/// # Errors
/// As [`rate`].
#[allow(clippy::too_many_arguments)]
pub fn shell_and_tube_design(
    tema_designation: &str,
    tema_class: &str,
    required_area: f64,
    shell_side_pressure: f64,
    tube_side_pressure: f64,
    design_temperature: ThermodynamicTemperature,
    shell_material_grade: &str,
    tube_material_grade: &str,
    tube_size: &str,
    tube_length: Length,
    tube_passes: f64,
    shell_joint_efficiency: f64,
    nozzle_diameter: Length,
    nozzle_wall_thickness: Length,
    sour_service_assessment: bool,
    h2s_partial_pressure: Option<Pressure>,
    fouling_tube: Option<f64>,
    fouling_shell: Option<f64>,
    shell_side_method: Option<&str>,
    tube_side_density: Option<MassDensity>,
    tube_side_viscosity: Option<DynamicViscosity>,
    tube_side_cp: Option<SpecificHeatCapacity>,
    tube_side_conductivity: Option<ThermalConductivity>,
    tube_side_mass_flow: Option<MassRate>,
    shell_side_density: Option<MassDensity>,
    shell_side_viscosity: Option<DynamicViscosity>,
    shell_side_cp: Option<SpecificHeatCapacity>,
    shell_side_conductivity: Option<ThermalConductivity>,
    shell_side_mass_flow: Option<MassRate>,
) -> Result<ShellAndTubeDesignResult> {
    rate(&ShellAndTubeDesignRequest {
        tema_designation: tema_designation.to_string(),
        tema_class: tema_class.parse()?,
        required_area,
        shell_side_pressure: shell_side_pressure / 1e5,
        tube_side_pressure: tube_side_pressure / 1e5,
        design_temperature: design_temperature.value,
        shell_material_grade: shell_material_grade.to_string(),
        tube_material_grade: tube_material_grade.to_string(),
        tube_size: tube_size.parse()?,
        tube_length: tube_length.value * 1000.0,
        tube_passes,
        shell_joint_efficiency,
        nozzle_diameter: nozzle_diameter.value * 1000.0,
        nozzle_wall_thickness: nozzle_wall_thickness.value * 1000.0,
        sour_service_assessment,
        h2s_partial_pressure: h2s_partial_pressure.map(|value| value.value / 1e5),
        fouling_tube,
        fouling_shell,
        shell_side_method: match shell_side_method {
            Some(value) => Some(value.parse()?),
            None => None,
        },
        tube_side_density: tube_side_density.map(|value| value.value),
        tube_side_viscosity: tube_side_viscosity.map(|value| value.value),
        tube_side_cp: tube_side_cp.map(|value| value.value),
        tube_side_conductivity: tube_side_conductivity.map(|value| value.value),
        tube_side_mass_flow: tube_side_mass_flow.map(|value| value.value),
        shell_side_density: shell_side_density.map(|value| value.value),
        shell_side_viscosity: shell_side_viscosity.map(|value| value.value),
        shell_side_cp: shell_side_cp.map(|value| value.value),
        shell_side_conductivity: shell_side_conductivity.map(|value| value.value),
        shell_side_mass_flow: shell_side_mass_flow.map(|value| value.value),
    })
}
