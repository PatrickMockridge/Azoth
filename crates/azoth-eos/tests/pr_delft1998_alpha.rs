//! Spec-driven tests for `eos.pr_delft1998_alpha`.

use azoth_eos::pr_delft1998_alpha;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.pr_delft1998_alpha";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::PrDelft1998AlphaResult {
    pr_delft1998_alpha(common::input(case, "omega"), common::input(case, "Tr"))
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
                        "omega" => Some(common::input(case, "omega")),
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
        executed >= 3,
        "expected the worked example, the reference case and the round trip, ran {executed}"
    );
}
