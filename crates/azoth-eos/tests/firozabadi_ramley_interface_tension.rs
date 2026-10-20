//! Spec-driven tests for `eos.firozabadi_ramley_interface_tension`.

use azoth_core::units::{kelvins, kilograms_per_cubic_meter};
use azoth_eos::firozabadi_ramley_interface_tension;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.firozabadi_ramley_interface_tension";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::FirozabadiRamleyInterfaceTensionResult {
    firozabadi_ramley_interface_tension(
        kelvins(common::input(case, "T")),
        kelvins(common::input(case, "pseudo_critical_temperature")),
        kilograms_per_cubic_meter(common::input(case, "rho_1")),
        kilograms_per_cubic_meter(common::input(case, "rho_2")),
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
                    result.sigma.value,
                    common::expected(case, "sigma"),
                    case.tolerance,
                    &format!("{}::{} (sigma)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "T" => Some(common::input(case, "T")),
                        "pseudo_critical_temperature" => {
                            Some(common::input(case, "pseudo_critical_temperature"))
                        }
                        "rho_1" => Some(common::input(case, "rho_1")),
                        "rho_2" => Some(common::input(case, "rho_2")),
                        "sigma" => Some(result.sigma.value),
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
