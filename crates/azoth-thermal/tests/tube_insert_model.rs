//! `thermal.tube_insert_model`, against `validation/neqsim/captures/tube_insert_probe.tsv`.
//!
//! The capture drives `TubeInsertModel` on its four insert types across every regime boundary
//! and the twist-ratio clamp, with the enhancement applied to a fixed plain state.

use azoth_test_support as common;
use azoth_thermal::spec_gen;
use azoth_thermal::tube_insert_model::{
    InsertType, TubeInsertRequest, rate, tape_blockage, twisted_tape_heat,
};

const CALC_ID: &str = "thermal.tube_insert_model";

fn request(case: &azoth_core::spec::TestCase) -> TubeInsertRequest {
    TubeInsertRequest {
        insert: common::input_str(case, "insert")
            .parse()
            .expect("a declared insert"),
        reynolds_number: common::input(case, "reynolds_number"),
        prandtl_number: common::input(case, "prandtl_number"),
        twist_ratio: case.input("twist_ratio"),
        tape_thickness: case.input("tape_thickness"),
        tape_width_fraction: case.input("tape_width_fraction"),
        wire_diameter: case.input("wire_diameter"),
        wire_coil_pitch: case.input("wire_coil_pitch"),
        matrix_density: case.input("matrix_density"),
        roughness_ratio: case.input("roughness_ratio"),
        helix_angle: case.input("helix_angle"),
        plain_heat_transfer_coefficient: common::input(case, "plain_heat_transfer_coefficient"),
        plain_pressure_drop: common::input(case, "plain_pressure_drop"),
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
            (
                "heat_transfer_enhancement_ratio",
                out.heat_transfer_enhancement_ratio,
            ),
            (
                "pressure_drop_penalty_ratio",
                out.pressure_drop_penalty_ratio,
            ),
            (
                "performance_evaluation_criteria",
                out.performance_evaluation_criteria,
            ),
            (
                "enhanced_heat_transfer_coefficient",
                out.enhanced_heat_transfer_coefficient.value,
            ),
            ("enhanced_pressure_drop", out.enhanced_pressure_drop.value),
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
        executed += 1;
    }
    assert!(executed > 0, "no case ran");
}

#[test]
fn the_blockage_factor_cancels_the_tape_thickness() {
    common::assert_close(
        tape_blockage(0.001, 0.95),
        tape_blockage(0.003, 0.95),
        1e-15,
        "cancels",
    );
    common::assert_close(
        twisted_tape_heat(10000.0, 5.0, 0.001, 0.95),
        1.4578573849941363,
        1e-12,
        "the worked example",
    );
}

#[test]
fn a_boundary_reynolds_number_takes_the_upper_branch() {
    let at_2300 = rate(&TubeInsertRequest {
        insert: InsertType::WireMatrix,
        reynolds_number: 2300.0,
        prandtl_number: 5.0,
        twist_ratio: None,
        tape_thickness: None,
        tape_width_fraction: None,
        wire_diameter: None,
        wire_coil_pitch: None,
        matrix_density: Some(0.5),
        roughness_ratio: None,
        helix_angle: None,
        plain_heat_transfer_coefficient: 1000.0,
        plain_pressure_drop: 5000.0,
    })
    .expect("2300 computes");
    common::assert_close(
        at_2300.heat_transfer_enhancement_ratio,
        2.677121188565862,
        1e-12,
        "turbulent at 2300",
    );
}

#[test]
fn the_bounds_and_an_unknown_insert_are_refused() {
    let base = TubeInsertRequest {
        insert: InsertType::TwistedTape,
        reynolds_number: 10000.0,
        prandtl_number: 5.0,
        twist_ratio: None,
        tape_thickness: None,
        tape_width_fraction: None,
        wire_diameter: None,
        wire_coil_pitch: None,
        matrix_density: None,
        roughness_ratio: None,
        helix_angle: None,
        plain_heat_transfer_coefficient: 1000.0,
        plain_pressure_drop: 5000.0,
    };
    assert!(rate(&base).is_ok());

    for bad in [
        TubeInsertRequest {
            reynolds_number: -1.0,
            ..base.clone()
        },
        TubeInsertRequest {
            prandtl_number: 0.0,
            ..base.clone()
        },
        TubeInsertRequest {
            tape_thickness: Some(0.0),
            ..base.clone()
        },
        TubeInsertRequest {
            helix_angle: Some(120.0),
            ..base.clone()
        },
        TubeInsertRequest {
            matrix_density: Some(1.5),
            ..base
        },
    ] {
        assert!(rate(&bad).is_err());
    }
    assert!("twisted".parse::<InsertType>().is_err());
}
