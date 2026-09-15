//! Spec-driven tests for `eos.twu_kappa`.

use azoth_eos::spec_gen;
use azoth_eos::twu_kappa;
use azoth_test_support as common;

const CALC_ID: &str = "eos.twu_kappa";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::TwuKappaResult {
    twu_kappa(common::input(case, "omega"))
        .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
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
                common::assert_close(
                    result.kappa,
                    common::expected(case, "kappa"),
                    case.tolerance,
                    &format!("{}::{} (kappa)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "omega" => Some(common::input(case, "omega")),
                        "kappa" => Some(result.kappa),
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
