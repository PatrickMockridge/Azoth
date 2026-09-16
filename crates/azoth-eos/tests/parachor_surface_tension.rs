//! Spec-driven tests for `eos.parachor_surface_tension`.

use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_eos::parachor_surface_tension;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.parachor_surface_tension";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::ParachorSurfaceTensionResult {
    parachor_surface_tension(
        common::input(case, "parachor"),
        kilograms_per_cubic_meter(common::input(case, "rho_l")),
        kilograms_per_cubic_meter(common::input(case, "rho_v")),
        kilograms_per_mole(common::input(case, "M")),
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
                        "parachor" => Some(common::input(case, "parachor")),
                        "rho_l" => Some(common::input(case, "rho_l")),
                        "rho_v" => Some(common::input(case, "rho_v")),
                        "M" => Some(common::input(case, "M")),
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
