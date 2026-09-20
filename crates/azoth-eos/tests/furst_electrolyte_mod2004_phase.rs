//! Spec-driven tests for the `eos.furst_electrolyte_mod2004_phase` model.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::furst_electrolyte_mod2004_phase::furst_electrolyte_mod2004_phase;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.furst_electrolyte_mod2004_phase";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::FurstElectrolyteMod2004PhaseResult {
    let components: Vec<String> = case
        .list("components")
        .expect("components")
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    furst_electrolyte_mod2004_phase(
        &components,
        kelvins(common::input(case, "T")),
        pascals(common::input(case, "P")),
        case.vector("x").expect("x"),
        common::input_str(case, "compressed_phase"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = azoth_eos::model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        common::assert_close(
            result.z_factor,
            case.expected_value("z_factor").expect("z_factor"),
            case.tolerance,
            &format!("{context} (z_factor)"),
        );
        let expected = case.expected_vector("ln_phi").expect("ln_phi");
        for (i, (&got, &want)) in result.ln_phi.iter().zip(expected).enumerate() {
            common::assert_close(
                got,
                want,
                case.tolerance,
                &format!("{context} (ln_phi[{i}])"),
            );
        }
        common::assert_consistent(&result, context);
    }
}
