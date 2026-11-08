//! `thermal.tema_standard`, against `validation/neqsim/captures/tema_standard_probe.tsv`.
//!
//! The capture drives `TEMAStandard` on five common and two custom configurations, every enum
//! row, all eight recommendation branches, the tube-count estimates and the span bands.

use azoth_test_support as common;
use azoth_thermal::spec_gen;
use azoth_thermal::tema_standard::{TemaStandardRequest, rate};

const CALC_ID: &str = "thermal.tema_standard";

fn request(case: &azoth_core::spec::TestCase) -> TemaStandardRequest {
    TemaStandardRequest {
        front_head: common::input_str(case, "front_head")
            .parse()
            .expect("a declared head"),
        shell: common::input_str(case, "shell")
            .parse()
            .expect("a declared shell"),
        rear_head: common::input_str(case, "rear_head")
            .parse()
            .expect("a declared rear head"),
        tema_class: common::input_str(case, "tema_class")
            .parse()
            .expect("a declared class"),
        tube_size: common::input_str(case, "tube_size")
            .parse()
            .expect("a declared size"),
        pitch_pattern: common::input_str(case, "pitch_pattern")
            .parse()
            .expect("a declared pattern"),
        baffle_type: common::input_str(case, "baffle_type")
            .parse()
            .expect("a declared baffle"),
        shell_inside_diameter: common::input(case, "shell_inside_diameter"),
        tube_pitch: common::input(case, "tube_pitch"),
        tube_passes: common::input(case, "tube_passes"),
        tube_material: case.string("tube_material").map(str::to_owned),
        needs_cleanable: case.flag("needs_cleanable").expect("a stated flag"),
        large_temperature_difference: case
            .flag("large_temperature_difference")
            .expect("a stated flag"),
        high_pressure: case.flag("high_pressure").expect("a stated flag"),
        hazardous: case.flag("hazardous").expect("a stated flag"),
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
        let out = rate(&request(case))
            .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id));
        let context = &format!("{}::{}", spec.id, case.id);
        for (field, value) in [
            ("configuration_cost_factor", out.configuration_cost_factor),
            ("shell_pressure_drop_factor", out.shell_pressure_drop_factor),
            ("class_cost_factor", out.class_cost_factor),
            ("min_tube_wall_thickness", out.min_tube_wall_thickness.value),
            ("min_corrosion_allowance", out.min_corrosion_allowance.value),
            ("tube_outer_diameter", out.tube_outer_diameter.value),
            (
                "tube_outer_diameter_inch",
                out.tube_outer_diameter_inch.value / 0.0254,
            ),
            ("layout_angle", out.layout_angle),
            ("min_pitch_ratio", out.min_pitch_ratio),
            ("pitch_heat_transfer_factor", out.pitch_heat_transfer_factor),
            (
                "baffle_heat_transfer_factor",
                out.baffle_heat_transfer_factor,
            ),
            (
                "baffle_pressure_drop_factor",
                out.baffle_pressure_drop_factor,
            ),
            ("min_tube_pitch", out.min_tube_pitch.value),
            ("estimated_tube_count", out.estimated_tube_count),
            ("max_unsupported_span", out.max_unsupported_span.value),
            ("min_baffle_spacing", out.min_baffle_spacing.value),
            ("max_baffle_spacing", out.max_baffle_spacing.value),
        ] {
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
            ("bundle_removable", out.bundle_removable),
            ("good_thermal_expansion", out.good_thermal_expansion),
        ] {
            if let Some(expected) = case.flag(field) {
                assert_eq!(value, expected, "{context}::{field}");
            }
        }
        if let Some(expected) = case.expected_string("recommended_designation") {
            assert_eq!(
                out.recommended_designation, expected,
                "{context}::recommended_designation"
            );
        }
        executed += 1;
    }
    assert!(executed > 0, "no case ran");
}

#[test]
fn the_tube_count_truncates_and_the_material_matches_by_substring() {
    use azoth_thermal::tema_standard::{PitchPattern, estimate_tube_count, material_factor};
    common::assert_close(
        estimate_tube_count(0.5, 0.01905, 0.02381, PitchPattern::Triangular30, 2),
        324.0,
        1e-12,
        "two passes",
    );
    assert_eq!(material_factor(Some("STAINLESS_STEEL")), 0.95);
    assert_eq!(material_factor(Some("CUPRONICKEL")), 0.75);
    assert_eq!(material_factor(Some("TITANIUM")), 0.9);
    assert_eq!(material_factor(Some("CARBON_STEEL")), 1.0);
    assert_eq!(material_factor(None), 1.0);
}

#[test]
fn a_non_integral_pass_count_is_refused() {
    let base = TemaStandardRequest {
        front_head: "A".parse().expect("A"),
        shell: "E".parse().expect("E"),
        rear_head: "S".parse().expect("S"),
        tema_class: "R".parse().expect("R"),
        tube_size: "tube_3_4_inch".parse().expect("3/4"),
        pitch_pattern: "triangular_30".parse().expect("30"),
        baffle_type: "single_segmental".parse().expect("single"),
        shell_inside_diameter: 0.5,
        tube_pitch: 0.02381,
        tube_passes: 2.0,
        tube_material: Some("CARBON_STEEL".to_string()),
        needs_cleanable: false,
        large_temperature_difference: false,
        high_pressure: false,
        hazardous: false,
    };
    assert!(rate(&base).is_ok());
    for bad in [
        TemaStandardRequest {
            tube_passes: 2.5,
            ..base.clone()
        },
        TemaStandardRequest {
            tube_passes: 0.0,
            ..base.clone()
        },
        TemaStandardRequest {
            shell_inside_diameter: 0.0,
            ..base.clone()
        },
        TemaStandardRequest {
            tube_pitch: 0.0,
            ..base
        },
    ] {
        assert!(rate(&bad).is_err());
    }
}
