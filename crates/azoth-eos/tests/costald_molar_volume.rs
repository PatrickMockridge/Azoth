//! Spec-driven tests for `eos.costald_molar_volume`.

use azoth_core::units::{
    cubic_meters_per_mole, kelvins, kilograms_per_cubic_meter, kilograms_per_mole,
};
use azoth_eos::costald_molar_volume;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.costald_molar_volume";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::CostaldMolarVolumeResult {
    costald_molar_volume(
        common::input(case, "omega"),
        kelvins(common::input(case, "Tc")),
        cubic_meters_per_mole(common::input(case, "Vc")),
        kilograms_per_mole(common::input(case, "M")),
        kilograms_per_cubic_meter(common::input(case, "rho_normal")),
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
                    result.v.value,
                    common::expected(case, "v"),
                    case.tolerance,
                    &format!("{}::{} (v)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "omega" => Some(common::input(case, "omega")),
                        "Tc" => Some(common::input(case, "Tc")),
                        "Vc" => Some(common::input(case, "Vc")),
                        "M" => Some(common::input(case, "M")),
                        "rho_normal" => Some(common::input(case, "rho_normal")),
                        "T" => Some(common::input(case, "T")),
                        "v" => Some(result.v.value),
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
