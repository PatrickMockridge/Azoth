//! Spec-driven tests for the `characterization.tbp_grouping` model.
//!
//! Every expected value is a row of `validation/neqsim/captures/tbp_grouping_probe.tsv`.

use azoth_characterization::{TbpGroupingResult, model_gen, tbp_grouping};
use azoth_core::units::kelvins;
use azoth_core::{AzothError, spec::TestCase};
use azoth_test_support as common;

const MODEL_ID: &str = "characterization.tbp_grouping";

fn call(case: &TestCase) -> TbpGroupingResult {
    let temperatures: Vec<_> = case
        .vector("boiling_point")
        .expect("the case states components")
        .iter()
        .map(|value| kelvins(*value))
        .collect();
    tbp_grouping(
        &temperatures,
        case.vector("mole_fraction")
            .expect("the case states components"),
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
        let expected = case
            .expected_vector("group_fraction")
            .expect("the case expects twenty bins");
        assert_eq!(result.group_fraction.len(), expected.len(), "{context}: the bin count");
        for (index, (got, want)) in result.group_fraction.iter().zip(expected).enumerate() {
            common::assert_close(*got, *want, case.tolerance, &format!("{context}::bin[{index}]"));
        }
        common::assert_consistent(&result, &context);
    }
}

#[test]
fn the_bins_are_always_twenty_and_the_first_six_are_never_filled() {
    // The class returns a fixed twenty-entry array and its lowest threshold opens index six, so
    // bins 0 to 5 are zeros in every case rather than zeros by coincidence on these three.
    let spec = model_gen::model(MODEL_ID).expect("in its own table");
    for case in spec.cases {
        let result = call(case);
        assert_eq!(result.group_fraction.len(), 20, "{}", case.id);
        assert_eq!(
            &result.group_fraction[..6],
            &[0.0; 6],
            "{}: bins 0 to 5 are below the lowest threshold",
            case.id
        );
    }
}

#[test]
fn a_component_below_the_lowest_threshold_is_in_no_bin() {
    // 66.85 C is under the 69.2 the lowest bin opens at, so it is not in bin six and not in bin
    // zero - it is absent. That is what makes the twenty entries sum to less than one, and the
    // capture's first row sums to 0.493 for exactly this reason.
    let result = tbp_grouping(&[kelvins(340.0), kelvins(400.0)], &[0.3, 0.7]).expect("computes");
    assert_eq!(result.group_fraction[6], 0.0);
    assert!((result.group_fraction[8] - 0.7).abs() < 1e-15);
    assert!((result.group_fraction.iter().sum::<f64>() - 0.7).abs() < 1e-15);
    // A fluid with nothing above the threshold returns zeros rather than refusing.
    let light = tbp_grouping(&[kelvins(111.6), kelvins(231.1)], &[0.5, 0.5]).expect("computes");
    assert!(light.group_fraction.iter().all(|value| *value == 0.0));
}

#[test]
fn two_components_between_one_pair_of_thresholds_share_a_bin() {
    // 100.0 C and 120.0 C are both above 98.9 and below 126.1, so the class has nothing to tell
    // them apart by - it carries a boiling point and no carbon number.
    let result = tbp_grouping(&[kelvins(373.15), kelvins(393.15)], &[0.25, 0.75]).expect("computes");
    assert!((result.group_fraction[7] - 1.0).abs() < 1e-15);
    for (index, value) in result.group_fraction.iter().enumerate() {
        if index != 7 {
            assert_eq!(*value, 0.0, "bin[{index}] should be empty");
        }
    }
}

#[test]
fn the_two_vectors_have_to_be_one_component_list() {
    let refused = tbp_grouping(&[kelvins(400.0), kelvins(450.0)], &[0.5]);
    assert!(matches!(refused, Err(AzothError::InvalidInput { .. })));
}
