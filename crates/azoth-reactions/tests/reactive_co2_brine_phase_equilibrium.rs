//! Spec-driven tests for `reactions.reactive_co2_brine_phase_equilibrium`.
//!
//! The oracle is `validation/neqsim/ReactiveBrineProbe.java`, whose capture is
//! `captures/reactive_brine_probe.tsv`. The case carries the scalars the generator can
//! express - the phase fractions and the minimum trial distance - and **the phase
//! compositions are asserted here**, because a matrix is not an expectation the generator
//! carries.

use azoth_core::spec::TestCase;
use azoth_core::units::{kelvins, pascals};
use azoth_reactions::model_gen;
use azoth_reactions::reactive_co2_brine_phase_equilibrium::{
    ReactiveCo2BrinePhaseEquilibriumResult, reactive_co2_brine_phase_equilibrium,
};
use azoth_test_support as common;

const MODEL_ID: &str = "reactions.reactive_co2_brine_phase_equilibrium";

fn call(case: &TestCase) -> ReactiveCo2BrinePhaseEquilibriumResult {
    let components: Vec<String> = case
        .list("components")
        .expect("the case states its components")
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    let moles = case.vector("moles").expect("the case states the amounts");
    reactive_co2_brine_phase_equilibrium(
        &components,
        moles,
        kelvins(common::input(case, "T")),
        pascals(common::input(case, "P")),
    )
    .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model has no cases");
    for case in spec.cases {
        let result = call(case);
        let tolerance = case.tolerance;
        assert_eq!(
            f64::from(result.phase_count),
            common::expected(case, "phase_count"),
            "{}: phase_count",
            case.id
        );
        let expected = case
            .expected_vector("phase_fractions")
            .expect("the case states the phase fractions");
        assert_eq!(result.phase_fractions.len(), expected.len(), "{}", case.id);
        for (actual, wanted) in result.phase_fractions.iter().zip(expected) {
            common::assert_close(
                *actual,
                *wanted,
                tolerance,
                &format!("{}: phase_fractions", case.id),
            );
        }
        common::assert_close(
            result.minimum_trial_distance,
            common::expected(case, "minimum_trial_distance"),
            tolerance,
            &format!("{}: minimum_trial_distance", case.id),
        );
        // The two gates the class accepts on, at the answer.
        assert!(
            result.charge_moles.abs() <= 1.0e-8,
            "{}: net charge {}",
            case.id,
            result.charge_moles
        );
        assert!(
            result.element_balance_residual <= 1.0e-8,
            "{}: element residual {}",
            case.id,
            result.element_balance_residual
        );
    }
}

/// The two captured aqueous compositions, which the case cannot carry.
#[test]
fn the_captured_compositions() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    let captures: &[(&str, &[[f64; 9]; 2])] = &[
        (
            "the_50_bara_brine_at_rest",
            &[
                [
                    0.997906891370804,
                    0.00209310862919582,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                ],
                [
                    0.0333389622860752,
                    0.938518891763245,
                    0.00914161812928211,
                    0.00477757361492435,
                    0.0139191917442065,
                    1.92817079957177e-14,
                    0.000151881228597570,
                    0.000151881233646759,
                    3.13212710035524e-15,
                ],
            ],
        ),
        (
            "the_40_bara_brine_at_ninety_percent_salt",
            &[
                [
                    0.999475012492229,
                    0.000524987507771661,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                ],
                [
                    0.0296605476092939,
                    0.944836013817328,
                    0.00828080144465402,
                    0.00432769537437606,
                    0.0126084968190301,
                    1.95651138791582e-14,
                    0.000143222465106599,
                    0.000143222470188509,
                    3.08143289652555e-15,
                ],
            ],
        ),
    ];
    for (id, rows) in captures {
        let case = spec
            .cases
            .iter()
            .find(|case| case.id == *id)
            .expect("the capture names a case");
        let result = call(case);
        assert_eq!(result.phase_composition.len(), rows.len(), "{id}");
        for (row, wanted) in result.phase_composition.iter().zip(rows.iter()) {
            assert_eq!(row.len(), wanted.len(), "{id}");
            for (index, (actual, expected)) in row.iter().zip(wanted).enumerate() {
                common::assert_close(
                    *actual,
                    *expected,
                    case.tolerance,
                    &format!("{id}: phase_composition[{index}]"),
                );
            }
        }
    }
}
