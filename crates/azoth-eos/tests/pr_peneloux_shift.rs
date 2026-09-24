//! Spec-driven tests for `eos.pr_peneloux_shift`.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::pr_peneloux_shift;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.pr_peneloux_shift";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::PrPenelouxShiftResult {
    // `z_ra` is optional and a case that states none is asking for the correlation, which
    // is NeqSim's own fallback for a zero column.
    let z_ra = case.input("z_ra");
    pr_peneloux_shift(
        common::input(case, "omega"),
        kelvins(common::input(case, "Tc")),
        pascals(common::input(case, "Pc")),
        z_ra,
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
                    result.c.value,
                    common::expected(case, "c"),
                    case.tolerance,
                    &format!("{}::{} (c)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "omega" => Some(common::input(case, "omega")),
                        "Tc" => Some(common::input(case, "Tc")),
                        "Pc" => Some(common::input(case, "Pc")),
                        "z_ra" => case.input("z_ra"),
                        "c" => Some(result.c.value),
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
