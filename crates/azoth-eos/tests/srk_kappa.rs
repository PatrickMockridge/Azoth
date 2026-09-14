//! Spec-driven tests for `eos.srk_kappa`.

use azoth_eos::spec_gen;
use azoth_eos::srk_kappa;
use azoth_test_support as common;

const CALC_ID: &str = "eos.srk_kappa";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::SrkKappaResult {
    srk_kappa(common::input(case, "omega"))
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

/// The worked example claims the result is the double nearest the exact decimal, which
/// is what justifies a 1e-12 tolerance on decimal-input arithmetic.
#[test]
fn the_worked_example_is_exact_in_binary() {
    let r = srk_kappa(0.152).unwrap();
    assert_eq!(
        r.kappa.to_bits(),
        0.715181696_f64.to_bits(),
        "the worked example is no longer bit-exact; the spec's derivation needs revising"
    );
}
