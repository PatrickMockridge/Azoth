//! Spec-driven tests for the `eos.viscosity` model.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::databank;
use azoth_eos::model_gen;
use azoth_eos::viscosity;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.viscosity";

fn from_case(case: &azoth_core::spec::TestCase) -> azoth_eos::mixture::Mixture {
    let names = case
        .list("components")
        .expect("the case declares components");
    databank::mixture_of(names, None)
        .expect("the case's components resolve")
        .0
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let mixture = from_case(case);
        let result = viscosity(
            &mixture,
            kelvins(common::input(case, "T")),
            pascals(common::input(case, "P")),
            case.vector("z").expect("z"),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        common::assert_close(
            result.mu.value,
            common::expected(case, "mu"),
            case.tolerance,
            &format!("{}::{} (mu)", spec.id, case.id),
        );
        common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
    }
}
