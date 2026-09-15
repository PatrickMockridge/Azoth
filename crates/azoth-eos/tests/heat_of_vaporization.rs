//! Spec-driven tests for `eos.heat_of_vaporization`.

use azoth_core::units::kelvins;
use azoth_eos::heat_of_vaporization;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.heat_of_vaporization";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::HeatOfVaporizationResult {
    heat_of_vaporization(
        common::input(case, "c0"),
        common::input(case, "c1"),
        common::input(case, "c2"),
        common::input(case, "c3"),
        kelvins(common::input(case, "Tc")),
        kelvins(common::input(case, "T")),
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
                    result.hov.value,
                    common::expected(case, "hov"),
                    case.tolerance,
                    &format!("{}::{} (hov)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "c0" => Some(common::input(case, "c0")),
                        "c1" => Some(common::input(case, "c1")),
                        "c2" => Some(common::input(case, "c2")),
                        "c3" => Some(common::input(case, "c3")),
                        "Tc" => Some(common::input(case, "Tc")),
                        "T" => Some(common::input(case, "T")),
                        "hov" => Some(result.hov.value),
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
