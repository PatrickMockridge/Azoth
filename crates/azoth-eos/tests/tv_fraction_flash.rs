//! Spec-driven tests for the `eos.tv_fraction_flash` model.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::databank;
use azoth_eos::tv_fraction_flash::tv_fraction_flash;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.tv_fraction_flash";

#[test]
fn every_case_in_the_spec() {
    let spec = azoth_eos::model_gen::model(MODEL_ID).expect("the model");
    assert!(!spec.cases.is_empty(), "the model should have cases");
    for case in spec.cases {
        let (mixture, _) = databank::mixture_of(case.list("components").expect("components"), None)
            .expect("the case's fluid resolves");
        let context = &format!("{}::{}", spec.id, case.id);
        let result = tv_fraction_flash(
            &mixture,
            kelvins(common::input(case, "T")),
            common::input(case, "fraction"),
            pascals(common::input(case, "P")),
            case.vector("z").expect("z"),
        )
        .unwrap_or_else(|e| panic!("{context} should compute but failed: {e}"));
        common::assert_close(
            result.pressure.value,
            common::expected(case, "P"),
            case.tolerance,
            &format!("{context} (P)"),
        );
        common::assert_consistent(&result, context);
    }
}

/// The answer reproduces the fraction it was asked for, across the whole range.
///
/// The two cases sit at 0.5 and 0.9, so a model that answered one pressure for every
/// request near them would pass. This walks the interior, and asserts the *reported*
/// fraction rather than the pressure - the answer is defined by the fraction, and a
/// pressure is only the way there.
#[test]
fn the_answer_reproduces_the_fraction_it_was_asked_for() {
    let (mixture, _) = databank::mixture_of(&["methane", "n-butane"], None).expect("the pair");
    let z = [0.6, 0.4];
    for fraction in [0.05, 0.2, 0.5, 0.8, 0.95] {
        let result = tv_fraction_flash(&mixture, kelvins(330.0), fraction, pascals(2.5e6), &z)
            .unwrap_or_else(|e| panic!("a fraction of {fraction} should have a pressure: {e}"));
        assert!(
            (result.volume_fraction - fraction).abs() < 1e-6,
            "asked for {fraction}, the answer's own fraction is {}",
            result.volume_fraction
        );
        assert_eq!(result.phase, azoth_eos::Phase::TwoPhase);
    }
}

/// The pressure falls as the gas takes more of the volume.
///
/// The monotonicity the search relies on, asserted rather than assumed: a bisection or a
/// Newton on a non-monotone function finds whichever root it happens to approach.
#[test]
fn the_pressure_falls_as_the_gas_fraction_rises() {
    let (mixture, _) = databank::mixture_of(&["methane", "n-butane"], None).expect("the pair");
    let z = [0.6, 0.4];
    let mut previous = f64::INFINITY;
    for fraction in [0.1, 0.3, 0.5, 0.7, 0.9] {
        let p = tv_fraction_flash(&mixture, kelvins(330.0), fraction, pascals(2.5e6), &z)
            .expect("a pressure")
            .pressure
            .value;
        assert!(
            p < previous,
            "a volume fraction of {fraction} needs {p} Pa, and the smaller fraction before it needed {previous}"
        );
        previous = p;
    }
}
