//! Spec-driven tests for the `characterization.tbp_density` model.
//!
//! Every expected value is a row of `validation/neqsim/captures/tbp_closure_probe.tsv`.

use azoth_characterization::{TbpClosureKind, TbpDensityResult, model_gen, tbp_density};
use azoth_core::units::{kelvins, kilograms_per_mole};
use azoth_test_support as common;

const MODEL_ID: &str = "characterization.tbp_density";

fn call(case: &azoth_core::spec::TestCase) -> TbpDensityResult {
    let closure: TbpClosureKind = case
        .string("closure")
        .unwrap_or("riazi_daubert_1980")
        .parse()
        .expect("the case names a closure");
    tbp_density(
        closure,
        kelvins(common::input(case, "boiling_point")),
        kilograms_per_mole(common::input(case, "molar_mass")),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        common::assert_close(
            result.density.value,
            common::expected(case, "density"),
            case.tolerance,
            &format!("{}::{} (density)", spec.id, case.id),
        );
        common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
    }
}

#[test]
fn the_three_unsupported_closures_are_refused_by_name() {
    // Each refusal is an `[[unported]]` row, and each names the class that would close it. A
    // port that quietly answered with the 1980 pair for all four would return a number here.
    for (closure, key) in [
        (
            TbpClosureKind::RiaziDaubert1987,
            "closure=riazi_daubert_1987",
        ),
        (TbpClosureKind::Soreide, "closure=soreide"),
        (TbpClosureKind::TbpModel, "closure=tbp_model"),
    ] {
        let refused = tbp_density(closure, kelvins(500.0), kilograms_per_mole(0.2));
        let error = refused.expect_err(&format!("{key} must be refused"));
        assert!(
            error.to_string().contains(key),
            "{key}: the refusal should name its own key, got {error}"
        );
    }
}

#[test]
fn the_default_closure_is_the_one_that_inverts() {
    // The 1980 pair is `Default`, so an absent `closure` cannot silently pick an unsupported one.
    assert_eq!(TbpClosureKind::default(), TbpClosureKind::RiaziDaubert1980);
}
