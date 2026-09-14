//! Spec-driven tests for `eos.rk_alpha_ab`.

use azoth_eos::rk_alpha_ab;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.rk_alpha_ab";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::RkAlphaAbResult {
    rk_alpha_ab(common::input(case, "Tr"), common::input(case, "Pr"))
        .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
}

fn field(result: &azoth_eos::RkAlphaAbResult, name: &str) -> f64 {
    match name {
        "alpha" => result.alpha,
        "a_reduced" => result.a_reduced,
        _ => result.b_reduced,
    }
}

#[test]
fn every_case_in_the_spec() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    common::assert_skips_are_explained(spec);

    let mut executed = 0;
    for case in spec.all_tests() {
        if !case.is_active() {
            continue;
        }
        match case.kind {
            "worked_example" | "reference" => {
                let result = call(case);
                for name in ["alpha", "a_reduced", "b_reduced"] {
                    common::assert_close(
                        field(&result, name),
                        common::expected(case, name),
                        case.tolerance,
                        &format!("{}::{} ({name})", spec.id, case.id),
                    );
                }
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "Tr" => Some(common::input(case, "Tr")),
                        "Pr" => Some(common::input(case, "Pr")),
                        "alpha" => Some(result.alpha),
                        "a_reduced" => Some(result.a_reduced),
                        "b_reduced" => Some(result.b_reduced),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
            }
            "property" => panic!("{}::{}: unhandled property", spec.id, case.id),
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 2,
        "expected several active cases, ran {executed}"
    );
}
