//! Spec-driven tests for `eos.brine_mixing_scale_evaluator`.

use azoth_core::spec::{RecordField, TestCase};
use azoth_core::units::{kelvins, kilograms_per_cubic_meter, pascals};
use azoth_eos::brine_mixing_scale_evaluator::{BrineAnalysis, brine_mixing_scale_evaluator};
use azoth_eos::model_gen;
use azoth_eos::results::{BrineMixingPointRow, BrineMixingScaleEvaluatorResult};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.brine_mixing_scale_evaluator";

fn analysis(row: &[RecordField]) -> BrineAnalysis {
    let number = |name: &str| common::record_number(row, name).expect("every field is stated");
    BrineAnalysis {
        calcium: kilograms_per_cubic_meter(number("calcium")),
        barium: kilograms_per_cubic_meter(number("barium")),
        strontium: kilograms_per_cubic_meter(number("strontium")),
        magnesium: kilograms_per_cubic_meter(number("magnesium")),
        sodium: kilograms_per_cubic_meter(number("sodium")),
        potassium: kilograms_per_cubic_meter(number("potassium")),
        chloride: kilograms_per_cubic_meter(number("chloride")),
        sulphate: kilograms_per_cubic_meter(number("sulphate")),
        bicarbonate: kilograms_per_cubic_meter(number("bicarbonate")),
    }
}

fn call(case: &TestCase) -> BrineMixingScaleEvaluatorResult {
    let brine_a: Vec<BrineAnalysis> = common::record_input(case, "brine_a")
        .iter()
        .map(|row| analysis(row))
        .collect();
    let brine_b: Vec<BrineAnalysis> = common::record_input(case, "brine_b")
        .iter()
        .map(|row| analysis(row))
        .collect();
    brine_mixing_scale_evaluator(
        &brine_a,
        &brine_b,
        kelvins(common::input(case, "temperature")),
        pascals(common::input(case, "pressure")),
        common::input(case, "ph"),
        common::input(case, "steps") as usize,
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

fn assert_row(expected: &[RecordField], row: &BrineMixingPointRow, at: &str, tol: f64) {
    for (label, got) in [
        ("fraction_a", row.fraction_a),
        ("si_caco3", row.si_caco3),
        ("si_baso4", row.si_baso4),
        ("si_caso4", row.si_caso4),
        ("si_srso4", row.si_srso4),
    ] {
        common::assert_close(
            got,
            common::record_number(expected, label).expect("every field is stated"),
            tol,
            &format!("{at}::{label}"),
        );
    }
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        for (field, actual) in [
            ("worst_fraction_a", result.worst_fraction_a),
            ("worst_saturation_index", result.worst_saturation_index),
        ] {
            if let Some(expected) = case.expected_value(field) {
                common::assert_close(
                    actual,
                    expected,
                    case.tolerance,
                    &format!("{context} ({field})"),
                );
            }
        }
        if let Some(expected) = case.expected_string("worst_mineral") {
            assert_eq!(
                result.worst_mineral, expected,
                "{context}: the worst-case mineral"
            );
        }
        let expected = common::expected_record(case, "points");
        assert_eq!(
            result.points.len(),
            expected.len(),
            "{context}: the point count"
        );
        for (index, (row, declared)) in result.points.iter().zip(expected).enumerate() {
            assert_row(
                declared,
                row,
                &format!("{context}::points[{index}]"),
                case.tolerance,
            );
        }
        common::assert_consistent(&result, context);
    }
}

/// **A brine is exactly one row.** A list of two is refused rather than read for its first
/// entry, because the choice of which one to read would be this port's and not the class's.
#[test]
fn a_two_row_brine_is_refused() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    let case = &spec.cases[0];
    let brine_a: Vec<BrineAnalysis> = common::record_input(case, "brine_a")
        .iter()
        .map(|row| analysis(row))
        .collect();
    let brine_b: Vec<BrineAnalysis> = common::record_input(case, "brine_b")
        .iter()
        .map(|row| analysis(row))
        .collect();
    let doubled = [brine_a[0], brine_a[0]];
    let error =
        brine_mixing_scale_evaluator(&doubled, &brine_b, kelvins(333.15), pascals(1.0e7), 6.0, 11)
            .expect_err("two brine A rows are not one");
    assert!(error.to_string().contains("brine_a"));
}

/// **`steps` is floored at two**, as the class's `setSteps` floors it; one point is not a
/// mixing range.
#[test]
fn steps_are_floored_at_two() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    let case = &spec.cases[0];
    let brine_a: Vec<BrineAnalysis> = common::record_input(case, "brine_a")
        .iter()
        .map(|row| analysis(row))
        .collect();
    let brine_b: Vec<BrineAnalysis> = common::record_input(case, "brine_b")
        .iter()
        .map(|row| analysis(row))
        .collect();
    let result =
        brine_mixing_scale_evaluator(&brine_a, &brine_b, kelvins(333.15), pascals(1.0e7), 6.0, 1)
            .expect("the state runs");
    assert_eq!(result.points.len(), 2);
    assert_eq!(result.points[0].fraction_a, 0.0);
    assert_eq!(result.points[1].fraction_a, 1.0);
}

/// **A negative ion is refused**, on either brine and under its own name.
#[test]
fn a_negative_ion_is_refused() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    let case = &spec.cases[0];
    let brine_a: Vec<BrineAnalysis> = common::record_input(case, "brine_a")
        .iter()
        .map(|row| analysis(row))
        .collect();
    let brine_b: Vec<BrineAnalysis> = common::record_input(case, "brine_b")
        .iter()
        .map(|row| analysis(row))
        .collect();
    let mut bad = brine_b[0];
    bad.barium = kilograms_per_cubic_meter(-1.0);
    let error =
        brine_mixing_scale_evaluator(&brine_a, &[bad], kelvins(333.15), pascals(1.0e7), 6.0, 11)
            .expect_err("a negative concentration is not a brine");
    assert!(error.to_string().contains("brine_b.barium"));
}
