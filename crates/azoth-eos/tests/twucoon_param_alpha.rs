//! Spec-driven tests for `eos.twucoon_param_alpha`.

use azoth_eos::spec_gen;
use azoth_eos::twucoon_param_alpha;
use azoth_test_support as common;

const CALC_ID: &str = "eos.twucoon_param_alpha";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::TwucoonParamAlphaResult {
    twucoon_param_alpha(
        common::input(case, "a"),
        common::input(case, "b"),
        common::input(case, "c"),
        common::input(case, "Tr"),
    )
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
                    result.alpha,
                    common::expected(case, "alpha"),
                    case.tolerance,
                    &format!("{}::{} (alpha)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "a" => Some(common::input(case, "a")),
                        "b" => Some(common::input(case, "b")),
                        "c" => Some(common::input(case, "c")),
                        "Tr" => Some(common::input(case, "Tr")),
                        "alpha" => Some(result.alpha),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
            }
            "property" => match case.property {
                Some("unit_round_trip") => {}
                other => panic!("{}::{}: unhandled property {other:?}", spec.id, case.id),
            },
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 2,
        "expected several active cases, ran {executed}"
    );
}
