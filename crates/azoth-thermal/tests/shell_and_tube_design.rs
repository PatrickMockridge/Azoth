//! `thermal.shell_and_tube_design`, against
//! `validation/neqsim/captures/shell_and_tube_design_probe.tsv`.
//!
//! The capture drives `ShellAndTubeDesignCalculator` on its defaults, every rear-head branch,
//! the size table ends, all tube sizes, the derating boundaries, five material pairs, the
//! NACE threshold, joint efficiencies, nozzles and the thermal block. This file holds the
//! defects the port reproduces and the refusals the class's permissive path does not make.

use azoth_test_support as common;
use azoth_thermal::shell_and_tube_design::{ShellAndTubeDesignRequest, rate};
use azoth_thermal::spec_gen;

const CALC_ID: &str = "thermal.shell_and_tube_design";

fn request_for(case: &azoth_core::spec::TestCase) -> ShellAndTubeDesignRequest {
    ShellAndTubeDesignRequest {
        tema_designation: common::input_str(case, "tema_designation").to_string(),
        tema_class: common::input_str(case, "tema_class")
            .parse()
            .expect("a declared class"),
        required_area: common::input(case, "required_area"),
        shell_side_pressure: common::input(case, "shell_side_pressure"),
        tube_side_pressure: common::input(case, "tube_side_pressure"),
        design_temperature: common::input(case, "design_temperature"),
        shell_material_grade: common::input_str(case, "shell_material_grade").to_string(),
        tube_material_grade: common::input_str(case, "tube_material_grade").to_string(),
        tube_size: common::input_str(case, "tube_size")
            .parse()
            .expect("a declared tube size"),
        tube_length: common::input(case, "tube_length"),
        tube_passes: common::input(case, "tube_passes"),
        shell_joint_efficiency: common::input(case, "shell_joint_efficiency"),
        nozzle_diameter: common::input(case, "nozzle_diameter"),
        nozzle_wall_thickness: common::input(case, "nozzle_wall_thickness"),
        sour_service_assessment: case.flag("sour_service_assessment").expect("a stated flag"),
        h2s_partial_pressure: case.input("h2s_partial_pressure"),
        fouling_tube: case.input("fouling_tube"),
        fouling_shell: case.input("fouling_shell"),
        shell_side_method: case
            .string("shell_side_method")
            .map(|value| value.parse().expect("a declared shell-side method")),
        tube_side_density: case.input("tube_side_density"),
        tube_side_viscosity: case.input("tube_side_viscosity"),
        tube_side_cp: case.input("tube_side_cp"),
        tube_side_conductivity: case.input("tube_side_conductivity"),
        tube_side_mass_flow: case.input("tube_side_mass_flow"),
        shell_side_density: case.input("shell_side_density"),
        shell_side_viscosity: case.input("shell_side_viscosity"),
        shell_side_cp: case.input("shell_side_cp"),
        shell_side_conductivity: case.input("shell_side_conductivity"),
        shell_side_mass_flow: case.input("shell_side_mass_flow"),
    }
}

#[test]
fn every_case_in_the_spec() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    common::assert_skips_are_explained(spec);

    let mut executed = 0;
    for case in spec.all_tests() {
        if !case.is_active() {
            continue;
        }
        if case.kind != "worked_example" && case.kind != "reference" {
            continue;
        }
        let out = rate(&request_for(case))
            .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id));
        let context = format!("{}::{}", spec.id, case.id);
        let scalars = [
            (
                "shell_inside_diameter",
                out.shell_inside_diameter.value * 1000.0,
            ),
            (
                "shell_outside_diameter",
                out.shell_outside_diameter.value * 1000.0,
            ),
            (
                "shell_wall_thickness",
                out.shell_wall_thickness.value * 1000.0,
            ),
            (
                "shell_allowable_stress",
                out.shell_allowable_stress.value / 1e6,
            ),
            ("shell_smys", out.shell_smys.value / 1e6),
            ("shell_joint_efficiency", out.shell_joint_efficiency),
            ("tube_count", out.tube_count),
            (
                "tube_outer_diameter",
                out.tube_outer_diameter.value * 1000.0,
            ),
            (
                "tube_wall_thickness",
                out.tube_wall_thickness.value * 1000.0,
            ),
            ("tube_length", out.tube_length.value * 1000.0),
            ("tube_passes", out.tube_passes),
            ("tube_pitch", out.tube_pitch.value * 1000.0),
            (
                "tube_allowable_stress",
                out.tube_allowable_stress.value / 1e6,
            ),
            (
                "tube_thermal_conductivity",
                out.tube_thermal_conductivity.value,
            ),
            ("baffle_count", out.baffle_count),
            ("baffle_spacing", out.baffle_spacing.value * 1000.0),
            ("baffle_cut", out.baffle_cut),
            ("required_area", out.required_area.value),
            ("actual_area", out.actual_area.value),
            ("area_margin", out.area_margin),
            (
                "tubesheet_thickness_uhx",
                out.tubesheet_thickness_uhx.value * 1000.0,
            ),
            ("mawp_shell_side", out.mawp_shell_side.value / 1e5),
            ("mawp_tube_side", out.mawp_tube_side.value / 1e5),
            (
                "hydro_test_pressure_shell",
                out.hydro_test_pressure_shell.value / 1e5,
            ),
            (
                "hydro_test_pressure_tube",
                out.hydro_test_pressure_tube.value / 1e5,
            ),
            (
                "nozzle_reinforcement_area",
                out.nozzle_reinforcement_area.value,
            ),
            ("shell_max_hardness", out.shell_max_hardness),
            ("tube_max_hardness", out.tube_max_hardness),
            ("shell_weight", out.shell_weight),
            ("tube_weight", out.tube_weight),
            ("tubesheet_weight", out.tubesheet_weight),
            ("head_weight", out.head_weight),
            ("baffle_weight", out.baffle_weight),
            ("total_dry_weight", out.total_dry_weight),
            ("operating_weight", out.operating_weight),
            ("material_cost", out.material_cost),
            ("fabrication_cost", out.fabrication_cost),
            ("total_cost", out.total_cost),
            (
                "tube_side_heat_transfer_coefficient",
                out.tube_side_heat_transfer_coefficient
                    .map_or(f64::NAN, |value| value.value),
            ),
            (
                "shell_side_heat_transfer_coefficient",
                out.shell_side_heat_transfer_coefficient
                    .map_or(f64::NAN, |value| value.value),
            ),
            (
                "overall_heat_transfer_coefficient",
                out.overall_heat_transfer_coefficient
                    .map_or(f64::NAN, |value| value.value),
            ),
            (
                "tube_side_pressure_drop",
                out.tube_side_pressure_drop
                    .map_or(f64::NAN, |value| value.value),
            ),
            (
                "shell_side_pressure_drop",
                out.shell_side_pressure_drop
                    .map_or(f64::NAN, |value| value.value),
            ),
            (
                "tube_side_velocity",
                out.tube_side_velocity.map_or(f64::NAN, |value| value.value),
            ),
            (
                "shell_side_velocity",
                out.shell_side_velocity
                    .map_or(f64::NAN, |value| value.value),
            ),
            (
                "tube_side_reynolds",
                out.tube_side_reynolds.unwrap_or(f64::NAN),
            ),
            (
                "shell_side_reynolds",
                out.shell_side_reynolds.unwrap_or(f64::NAN),
            ),
            (
                "effective_tube_inside_diameter",
                out.effective_tube_inside_diameter
                    .map_or(f64::NAN, |value| value.value),
            ),
            (
                "outside_heat_transfer_area",
                out.outside_heat_transfer_area
                    .map_or(f64::NAN, |value| value.value),
            ),
            (
                "tube_fouling_layer_resistance_outside",
                out.tube_fouling_layer_resistance_outside
                    .unwrap_or(f64::NAN),
            ),
            ("jc", out.jc.unwrap_or(f64::NAN)),
            ("jl", out.jl.unwrap_or(f64::NAN)),
            ("jb", out.jb.unwrap_or(f64::NAN)),
            ("js", out.js.unwrap_or(f64::NAN)),
            ("jr", out.jr.unwrap_or(f64::NAN)),
            (
                "natural_frequency",
                out.natural_frequency.map_or(f64::NAN, |value| value.value),
            ),
            (
                "vortex_shedding_frequency",
                out.vortex_shedding_frequency
                    .map_or(f64::NAN, |value| value.value),
            ),
            (
                "vortex_shedding_ratio",
                out.vortex_shedding_ratio.unwrap_or(f64::NAN),
            ),
            (
                "critical_velocity",
                out.critical_velocity.map_or(f64::NAN, |value| value.value),
            ),
            ("velocity_ratio", out.velocity_ratio.unwrap_or(f64::NAN)),
            (
                "acoustic_frequency",
                out.acoustic_frequency.map_or(f64::NAN, |value| value.value),
            ),
            ("acoustic_ratio", out.acoustic_ratio.unwrap_or(f64::NAN)),
        ];
        for (field, value) in scalars {
            if let Some(expected) = case.expected_value(field) {
                common::assert_close(
                    value,
                    expected,
                    case.tolerance,
                    &format!("{context}::{field}"),
                );
            }
        }
        for (field, value) in [
            (
                "nozzle_reinforcement_adequate",
                out.nozzle_reinforcement_adequate,
            ),
            ("sour_service_required", out.sour_service_required),
            ("shell_nace_compliant", out.shell_nace_compliant),
            ("tube_nace_compliant", out.tube_nace_compliant),
            ("passed", out.passed.unwrap_or(false)),
            (
                "vortex_shedding_critical",
                out.vortex_shedding_critical.unwrap_or(false),
            ),
            (
                "fluid_elastic_critical",
                out.fluid_elastic_critical.unwrap_or(false),
            ),
            ("acoustic_critical", out.acoustic_critical.unwrap_or(false)),
        ] {
            if let Some(expected) = case.flag(field) {
                assert_eq!(value, expected, "{context}::{field}");
            }
        }
        for (field, value) in [
            ("tema_designation", out.tema_designation.as_str()),
            ("tema_class", out.tema_class.as_str()),
            ("shell_material", out.shell_material.as_str()),
            ("shell_material_grade", out.shell_material_grade.as_str()),
            ("tube_material", out.tube_material.as_str()),
            ("tube_material_grade", out.tube_material_grade.as_str()),
            ("tube_pitch_pattern", out.tube_pitch_pattern.as_str()),
            ("baffle_type", out.baffle_type.as_str()),
        ] {
            if let Some(expected) = case.expected_string(field) {
                assert_eq!(value, expected, "{context}::{field}");
            }
        }
        if let Some(expected) = case.expected_record("applied_standards") {
            assert_eq!(out.applied_standards.len(), expected.len(), "{context}");
            for (row, claims) in out.applied_standards.iter().zip(expected) {
                assert_eq!(
                    row.name,
                    common::record_name(claims, "name").expect("a stated standard"),
                    "{context}::applied_standards"
                );
            }
        }
        if let Some(expected) = case.expected_record("nace_issues") {
            assert_eq!(out.nace_issues.len(), expected.len(), "{context}");
            for (row, claims) in out.nace_issues.iter().zip(expected) {
                assert_eq!(
                    row.message,
                    common::record_name(claims, "message").expect("a stated issue"),
                    "{context}::nace_issues"
                );
            }
        }
        executed += 1;
    }
    assert!(executed > 0, "no case ran");
}

#[test]
fn a_bad_designation_is_refused_by_name() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    let case = spec
        .all_tests()
        .find(|case| case.is_active() && case.kind == "worked_example")
        .expect("the worked example");
    let mut request = request_for(case);
    request.tema_designation = "AE".to_string();
    let error = rate(&request).expect_err("two letters are not a designation");
    assert!(error.to_string().contains("tema_designation"), "{error}");

    request.tema_designation = "AZE".to_string();
    assert!(rate(&request).is_err(), "an undeclared shell letter");
}

#[test]
fn a_partial_thermal_side_is_refused() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    let case = spec
        .all_tests()
        .find(|case| case.is_active() && case.kind == "worked_example")
        .expect("the worked example");
    let mut request = request_for(case);
    request.tube_side_density = Some(998.0);
    let error = rate(&request).expect_err("one of five is not a side");
    assert!(error.to_string().contains("tube_side"), "{error}");
}

#[test]
fn the_material_table_reads_the_two_upstream_defects() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    let case = spec
        .all_tests()
        .find(|case| case.id == "the_nace_assessment_reads_every_grade_as_not_compliant")
        .expect("the NACE case");
    let out = rate(&request_for(case)).expect("the NACE case computes");
    assert!(out.sour_service_required);
    assert!(!out.shell_nace_compliant);
    assert!(!out.tube_nace_compliant);
    assert_eq!(out.nace_issues.len(), 2);
    // The price test never sees STAINLESS, so a stainless tube costs as carbon.
    let base = spec
        .all_tests()
        .find(|case| case.id == "the_defaults_worked_example")
        .expect("the defaults");
    let stainless = spec
        .all_tests()
        .find(|case| case.id == "a_stainless_tube_moves_the_pressure_design_and_not_the_cost")
        .expect("the stainless case");
    let base_out = rate(&request_for(base)).expect("defaults");
    let stainless_out = rate(&request_for(stainless)).expect("stainless");
    assert_eq!(base_out.material_cost, stainless_out.material_cost);
    assert_eq!(base_out.total_dry_weight, stainless_out.total_dry_weight);
}
