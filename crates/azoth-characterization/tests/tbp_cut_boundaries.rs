//! Spec-driven tests for the `characterization.tbp_cut_boundaries` model.
//!
//! Every expected value is a row of `validation/neqsim/captures/tbp_cut_boundaries_probe.tsv`.

use azoth_characterization::{TbpCutBoundariesResult, model_gen, tbp_cut_boundaries};
use azoth_core::{AzothError, spec::TestCase};
use azoth_test_support as common;

const MODEL_ID: &str = "characterization.tbp_cut_boundaries";

fn call(case: &TestCase) -> TbpCutBoundariesResult {
    tbp_cut_boundaries(
        case.string("name_prefix").expect("the case names a prefix"),
        case.vector("cumulative_volume_percent")
            .expect("the case states a curve"),
        case.vector("boiling_point").expect("the case states a curve"),
        case.vector("specific_gravity")
            .expect("the case states a curve"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = format!("{}::{}", spec.id, case.id);
        for (name, actual, expected) in [
            (
                "volume_fraction",
                result.volume_fraction.clone(),
                case.expected_vector("volume_fraction"),
            ),
            (
                "specific_gravity",
                result.specific_gravity.clone(),
                case.expected_vector("specific_gravity"),
            ),
            (
                "lower_boiling_point",
                result.lower_boiling_point.iter().map(|q| q.value).collect(),
                case.expected_vector("lower_boiling_point"),
            ),
            (
                "upper_boiling_point",
                result.upper_boiling_point.iter().map(|q| q.value).collect(),
                case.expected_vector("upper_boiling_point"),
            ),
            (
                "average_boiling_point",
                result.average_boiling_point.iter().map(|q| q.value).collect(),
                case.expected_vector("average_boiling_point"),
            ),
        ] {
            let Some(expected) = expected else { continue };
            assert_eq!(actual.len(), expected.len(), "{context}::{name}: the cut count");
            for (index, (got, want)) in actual.iter().zip(expected).enumerate() {
                common::assert_close(
                    *got,
                    *want,
                    case.tolerance,
                    &format!("{context}::{name}[{index}]"),
                );
            }
        }
        // **The names are asserted here and not in the case**, because the case generator carries
        // no string-vector expectation - the same reason `reactive_tp_flash` asserts `phase_type`
        // in its own test.
        assert_eq!(result.cut_name, expected_names(case), "{context}: the cut names");
        common::assert_consistent(&result, &context);
    }
}

/// The names the capture's cut_count and prefix imply: the prefix and the position from one.
fn expected_names(case: &TestCase) -> Vec<String> {
    let prefix = case.string("name_prefix").expect("a prefix");
    let count = case
        .expected_vector("volume_fraction")
        .expect("the case states cuts")
        .len();
    (1..=count).map(|index| format!("{prefix}{index}")).collect()
}

#[test]
fn the_average_is_the_midpoint_of_the_boundaries_that_bracket_it() {
    // The capture's first curve: 320 and 380 bound the first cut, and its average is 350 - a
    // number no boundary carries. A port that took the lower or the upper boundary would give
    // 320 or 380 and agree with `volume_fraction` either way.
    let cuts = tbp_cut_boundaries(
        "TBP",
        &[0.0, 20.0, 45.0, 75.0, 100.0],
        &[320.0, 380.0, 440.0, 520.0, 620.0],
        &[0.68, 0.76, 0.82, 0.88],
    )
    .expect("computes");
    assert_eq!(cuts.cut_name, vec!["TBP1", "TBP2", "TBP3", "TBP4"]);
    assert!((cuts.lower_boiling_point[0].value - 320.0).abs() < 1e-12);
    assert!((cuts.upper_boiling_point[0].value - 380.0).abs() < 1e-12);
    assert!((cuts.average_boiling_point[0].value - 350.0).abs() < 1e-12);
    // And the yields are the interval widths, which sum to the curve's own span.
    assert!((cuts.volume_fraction.iter().sum::<f64>() - 1.0).abs() < 1e-15);
}

#[test]
fn the_class_s_curve_validation_is_a_refusal_for_each_rule() {
    // Each of these is one of the capture's five refused curves. The message names the field, so
    // the check is that the refusal is about the right thing and not a later one.
    let cases: [(&str, Vec<f64>, Vec<f64>, Vec<f64>); 5] = [
        (
            "N-1 gravities",
            vec![0.0, 50.0, 100.0],
            vec![320.0, 420.0, 520.0],
            vec![0.70, 0.80, 0.90],
        ),
        (
            "a curve starting above zero",
            vec![5.0, 50.0, 100.0],
            vec![320.0, 420.0, 520.0],
            vec![0.70, 0.80],
        ),
        (
            "yields that do not increase",
            vec![0.0, 60.0, 50.0, 100.0],
            vec![320.0, 380.0, 440.0, 520.0],
            vec![0.70, 0.80, 0.90],
        ),
        (
            "a temperature that does not increase",
            vec![0.0, 50.0, 100.0],
            vec![320.0, 300.0, 520.0],
            vec![0.70, 0.80],
        ),
        (
            "a gravity that is not positive",
            vec![0.0, 50.0, 100.0],
            vec![320.0, 420.0, 520.0],
            vec![0.70, 0.0],
        ),
    ];
    for (label, yields, boiling, gravity) in cases {
        let refused = tbp_cut_boundaries("TBP", &yields, &boiling, &gravity);
        assert!(
            matches!(refused, Err(AzothError::InvalidInput { .. })),
            "{label} should be refused"
        );
    }
    // An empty prefix names nothing, which the class refuses before it looks at the curve.
    assert!(matches!(
        tbp_cut_boundaries("  ", &[0.0, 100.0], &[320.0, 520.0], &[0.70]),
        Err(AzothError::InvalidInput { .. })
    ));
}
