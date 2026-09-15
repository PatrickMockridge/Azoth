//! Spec-driven tests for `eos.antoine_vapor_pressure`.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::antoine_vapor_pressure;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.antoine_vapor_pressure";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::AntoineVaporPressureResult {
    antoine_vapor_pressure(
        common::input(case, "A"),
        common::input(case, "B"),
        common::input(case, "C"),
        common::input(case, "D"),
        common::input(case, "E"),
        common::input_str(case, "form").parse().unwrap(),
        kelvins(common::input(case, "Tc")),
        pascals(common::input(case, "Pc")),
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
                    result.p_sat.value,
                    common::expected(case, "p_sat"),
                    case.tolerance,
                    &format!("{}::{} (p_sat)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "A" => Some(common::input(case, "A")),
                        "B" => Some(common::input(case, "B")),
                        "C" => Some(common::input(case, "C")),
                        "D" => Some(common::input(case, "D")),
                        "E" => Some(common::input(case, "E")),
                        "Tc" => Some(common::input(case, "Tc")),
                        "Pc" => Some(common::input(case, "Pc")),
                        "T" => Some(common::input(case, "T")),
                        "p_sat" => Some(result.p_sat.value),
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
    assert!(executed >= 2, "expected several active cases, ran {executed}");
}
