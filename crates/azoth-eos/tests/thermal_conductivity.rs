//! Spec-driven tests for the `eos.thermal_conductivity` model.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::databank;
use azoth_eos::model_gen;
use azoth_eos::thermal_conductivity;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.thermal_conductivity";

fn from_case(
    case: &azoth_core::spec::TestCase,
) -> (azoth_eos::mixture::Mixture, azoth_eos::IdealGasModel) {
    let names = case
        .list("components")
        .expect("the case declares components");
    databank::mixture_of(names, None).expect("the case's components resolve")
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let (mixture, ideal_gas) = from_case(case);
        let result = thermal_conductivity(
            &mixture,
            &ideal_gas,
            kelvins(common::input(case, "T")),
            pascals(common::input(case, "P")),
            case.vector("z").expect("z"),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        common::assert_close(
            result.k.value,
            common::expected(case, "k"),
            case.tolerance,
            &format!("{}::{} (k)", spec.id, case.id),
        );
        common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
    }
}
