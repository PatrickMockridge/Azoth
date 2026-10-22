//! Spec-driven tests for `eos.barite_celestite_solid_solution`.

use azoth_eos::barite_celestite_solid_solution::barite_celestite_solid_solution;
use azoth_eos::model_gen;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.barite_celestite_solid_solution";

fn call(
    case: &azoth_core::spec::TestCase,
) -> azoth_eos::results::BariteCelestiteSolidSolutionResult {
    let input = |name: &str| common::input(case, name);
    barite_celestite_solid_solution(
        input("margules"),
        input("a_ba"),
        input("a_sr"),
        input("a_so4"),
        input("ksp_baso4"),
        input("ksp_srso4"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        for (field, value) in [("x_baso4", result.x_baso4), ("x_srso4", result.x_srso4)] {
            common::assert_close(
                value,
                case.expected_value(field)
                    .unwrap_or_else(|| panic!("{context}: the case states {field}")),
                case.tolerance,
                &format!("{context} ({field})"),
            );
        }
        match (
            case.expected_value("total_saturation_index"),
            result.total_saturation_index,
        ) {
            (Some(expected), Some(actual)) => common::assert_close(
                actual,
                expected,
                case.tolerance,
                &format!("{context} (total_saturation_index)"),
            ),
            (None, None) => {}
            (expected, actual) => panic!(
                "{context}: the index is expected {expected:?} and the result carries {actual:?}"
            ),
        }
        assert_eq!(
            result.supersaturated,
            result
                .total_saturation_index
                .is_some_and(|value| value > 0.0),
            "{context}: the flag follows the index"
        );
        common::assert_consistent(&result, context);
    }
}

/// **The class's temperature setter is inert.** The Ksp are the caller's, and there is no
/// temperature input to this model; this test says so by construction rather than by
/// comparing two calls.
#[test]
fn the_end_member_ksp_are_the_callers() {
    let with_class_defaults =
        barite_celestite_solid_solution(2.3, 0.01, 0.0, 0.05, 1.08e-10, 3.44e-7)
            .expect("the state runs");
    let with_custom = barite_celestite_solid_solution(2.3, 0.01, 0.0, 0.05, 1.0e-2, 1.0e-1)
        .expect("the state runs");
    assert!(with_class_defaults.total_saturation_index.expect("index") > 0.0);
    assert!(with_custom.total_saturation_index.expect("index") < 0.0);
}
