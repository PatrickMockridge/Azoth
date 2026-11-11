//! `thermal.lmtd_correction_factor`, against `validation/neqsim/captures/lmtd_correction_probe.tsv`.
//!
//! The capture drives `LMTDcorrectionFactor` on the class's own R/P state and on the edge
//! states: R = 1, no heat transfer, equal inlets, a cross, a four-pass need and an
//! unachievable state.

use azoth_test_support as common;
use azoth_thermal::lmtd_correction_factor::{LmtdRequest, rate};
use azoth_thermal::spec_gen;

const CALC_ID: &str = "thermal.lmtd_correction_factor";

fn request(case: &azoth_core::spec::TestCase) -> LmtdRequest {
    LmtdRequest {
        hot_inlet_temperature: common::input(case, "hot_inlet_temperature"),
        hot_outlet_temperature: common::input(case, "hot_outlet_temperature"),
        cold_inlet_temperature: common::input(case, "cold_inlet_temperature"),
        cold_outlet_temperature: common::input(case, "cold_outlet_temperature"),
        shell_passes: case.input("shell_passes"),
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
            ("r", out.r),
            ("p", out.p),
            ("ft", out.ft),
            ("ft_1_shell_pass", out.ft_1_shell_pass),
            ("ft_2_shell_passes", out.ft_2_shell_passes),
            ("required_shell_passes", out.required_shell_passes),
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
fn the_split_shell_pass_count_is_refused() {
    let base = LmtdRequest {
        hot_inlet_temperature: 150.0,
        hot_outlet_temperature: 90.0,
        cold_inlet_temperature: 30.0,
        cold_outlet_temperature: 80.0,
        shell_passes: Some(1.5),
    };
    assert!(rate(&base).is_err());

    let zero_temperature = LmtdRequest {
        hot_inlet_temperature: 0.0,
        ..base
    };
    assert!(rate(&zero_temperature).is_err());
}
