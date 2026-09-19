//! Spec-driven tests for the `eos.vh_flash` model.

use azoth_core::units::{cubic_meters_per_mole, joules_per_mole};
use azoth_eos::Cubic;
use azoth_eos::databank;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.vh_flash";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::VhFlashResult {
    let (mixture, ideal_gas) = databank::mixture_of(
        case.list("components").expect("components"),
        Cubic::Pr,
        None,
    )
    .unwrap_or_else(|e| panic!("case `{}` should resolve but failed: {e}", case.id));
    azoth_eos::vh_flash::vh_flash(
        &mixture,
        &ideal_gas,
        cubic_meters_per_mole(common::input(case, "V")),
        joules_per_mole(common::input(case, "H")),
        case.vector("z").expect("z"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = azoth_eos::model_gen::model(MODEL_ID).expect("the model");
    assert!(!spec.cases.is_empty(), "the model should have cases");
    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        for name in ["P", "T"] {
            let got = if name == "P" {
                result.pressure.value
            } else {
                result.temperature.value
            };
            common::assert_close(
                got,
                common::expected(case, name),
                case.tolerance,
                &format!("{context} ({name})"),
            );
        }
        common::assert_consistent(&result, context);
    }
}
