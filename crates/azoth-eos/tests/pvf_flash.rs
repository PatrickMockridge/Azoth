//! Spec-driven tests for the `eos.pvf_flash` model.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::Cubic;
use azoth_eos::databank;
use azoth_eos::pvf_flash::pvf_flash;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.pvf_flash";

#[test]
fn every_case_in_the_spec() {
    let spec = azoth_eos::model_gen::model(MODEL_ID).expect("the model");
    assert!(!spec.cases.is_empty(), "the model should have cases");
    for case in spec.cases {
        let (mixture, _) = databank::mixture_of(
            case.list("components").expect("components"),
            Cubic::Pr,
            None,
        )
        .expect("the case's fluid resolves");
        let context = &format!("{}::{}", spec.id, case.id);
        let result = pvf_flash(
            &mixture,
            pascals(common::input(case, "P")),
            common::input(case, "beta"),
            kelvins(common::input(case, "temperature")),
            case.vector("z").expect("z"),
        )
        .unwrap_or_else(|e| panic!("{context} should compute but failed: {e}"));
        common::assert_close(
            result.t.value,
            common::expected(case, "T"),
            case.tolerance,
            &format!("{context} (T)"),
        );
        common::assert_close(
            result.beta,
            common::input(case, "beta"),
            case.tolerance,
            &format!("{context} (beta)"),
        );
        common::assert_consistent(&result, context);
    }
}

/// The endpoints belong to the bubble- and dew-point models, and are refused by name.
///
/// A vapour fraction of exactly zero or one is a saturation temperature, which
/// `eos.bubble_temperature` and `eos.dew_temperature` already are. Answering them here
/// would be a second implementation of a calculation that has one.
#[test]
fn the_endpoints_are_refused_by_name() {
    let (mixture, _) =
        databank::mixture_of(&["methane", "n-butane"], Cubic::Pr, None).expect("the pair");
    for (beta, which) in [(0.0, "bubble"), (1.0, "dew")] {
        let error = pvf_flash(&mixture, pascals(2.5e6), beta, kelvins(330.0), &[0.6, 0.4])
            .expect_err("an endpoint is another model's");
        let message = format!("{error}");
        assert!(
            message.contains(which),
            "the refusal should name the {which} point: {message}"
        );
    }
}

/// The fraction comes back as what was asked for, at every fraction across the region.
///
/// The two spec cases sit at 0.84 and 0.99 - both near the dew end - so a model that
/// answered a temperature for *any* fraction would pass them. This walks the interior.
#[test]
fn the_answer_reproduces_the_fraction_it_was_asked_for() {
    let (mixture, _) =
        databank::mixture_of(&["methane", "n-butane"], Cubic::Pr, None).expect("the pair");
    for fraction in [0.05, 0.2, 0.5, 0.8, 0.95] {
        let result = pvf_flash(
            &mixture,
            pascals(2.5e6),
            fraction,
            kelvins(330.0),
            &[0.6, 0.4],
        )
        .unwrap_or_else(|e| panic!("a fraction of {fraction} should have a temperature: {e}"));
        assert!(
            (result.beta - fraction).abs() < 1e-7,
            "asked for {fraction}, the answer's own fraction is {}",
            result.beta
        );
        assert_eq!(result.phase, azoth_eos::Phase::TwoPhase);
    }
}
