//! Spec-driven tests for the `eos.pcsaft_rahmat_phase` model.
//!
//! The spec's case pins the two implementations to each other. These tests are for what a
//! case cannot say: the checks that refuse rather than answer, and the difference between
//! the model and the solve it is built on.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::pcsaft_rahmat_phase::pcsaft_rahmat_phase;

/// The spec's case, through the model rather than through the solve: names resolved from
/// the databank, the checks run, the result the model's own.
#[test]
fn every_case_in_the_spec() {
    let spec =
        azoth_eos::model_gen::model("eos.pcsaft_rahmat_phase").expect("the model is registered");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let names: Vec<String> = azoth_test_support::list_input(case, "components")
            .iter()
            .map(|name| (*name).to_string())
            .collect();
        let result = pcsaft_rahmat_phase(
            &names,
            kelvins(azoth_test_support::input(case, "T")),
            pascals(azoth_test_support::input(case, "P")),
            case.vector("z").expect("the case declares z"),
            azoth_test_support::input_str(case, "compressed_phase"),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        let label = |field: &str| format!("{}::{} ({field})", spec.id, case.id);
        azoth_test_support::assert_close(
            result.z_factor,
            azoth_test_support::expected(case, "z_factor"),
            case.tolerance,
            &label("z_factor"),
        );
        azoth_test_support::assert_close(
            result.v.value,
            azoth_test_support::expected(case, "v"),
            case.tolerance,
            &label("v"),
        );
        let want = case
            .expected_vector("ln_phi")
            .expect("the case declares ln_phi");
        assert_eq!(
            result.ln_phi.len(),
            want.len(),
            "one coefficient per component"
        );
        for (i, (got, wanted)) in result.ln_phi.iter().zip(want).enumerate() {
            azoth_test_support::assert_close(
                *got,
                *wanted,
                case.tolerance,
                &label(&format!("ln_phi[{i}]")),
            );
        }
    }
}

/// The composition is checked rather than renormalised, and the branch is checked rather
/// than defaulted: both are the shape of failure that answers for a different fluid.
#[test]
fn a_bad_composition_or_branch_is_refused() {
    let names = ["methane".to_string()];
    let error = pcsaft_rahmat_phase(&names, kelvins(300.0), pascals(5.0e6), &[0.9], "vapour")
        .expect_err("a composition that does not sum to one");
    assert!(error.to_string().contains("sum"), "{error}");

    let error = pcsaft_rahmat_phase(&names, kelvins(300.0), pascals(5.0e6), &[1.0], "gas")
        .expect_err("a branch that is not one of the two");
    assert!(error.to_string().contains("not a side"), "{error}");
}
