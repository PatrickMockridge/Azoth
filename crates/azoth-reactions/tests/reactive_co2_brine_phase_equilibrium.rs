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
                    0.997908280273748,
                    0.00209171972625182,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                ],
                [
                    0.0245744601647729,
                    0.947049404880286,
                    0.00922468131144892,
                    0.00482098392389587,
                    0.0140456652353448,
                    1.78562206380913e-14,
                    0.000142402239566919,
                    0.000142402244662883,
                    3.45450282775326e-15,
                ],
            ],
        ),
        (
            "the_40_bara_brine_at_ninety_percent_salt",
            &[
                [
                    0.999475231950557,
                    0.000524768049443163,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                    1.0e-50,
                ],
                [
                    0.0225905405549093,
                    0.951736621393577,
                    0.00834116875419316,
                    0.00435924441319823,
                    0.0127004131673914,
                    1.83943920993111e-14,
                    0.000136005855794769,
                    0.000136005860914645,
                    3.33924548800970e-15,
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
