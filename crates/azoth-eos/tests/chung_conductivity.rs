//! Spec-driven tests for `eos.chung_conductivity`.

use azoth_core::units::{
    cubic_meters_per_mole, joules_per_mole_kelvin, kelvins, kilograms_per_mole,
};
use azoth_eos::chung_conductivity;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.chung_conductivity";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::ChungConductivityResult {
    chung_conductivity(
        joules_per_mole_kelvin(common::input(case, "Cv0")),
        kilograms_per_mole(common::input(case, "M")),
        common::input(case, "omega"),
        kelvins(common::input(case, "Tc")),
        cubic_meters_per_mole(common::input(case, "Vc")),
        common::input(case, "dipole"),
        common::input(case, "kappa"),
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
                    result.k.value,
                    common::expected(case, "k"),
                    case.tolerance,
                    &format!("{}::{} (k)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "Cv0" => Some(common::input(case, "Cv0")),
                        "M" => Some(common::input(case, "M")),
                        "omega" => Some(common::input(case, "omega")),
                        "Tc" => Some(common::input(case, "Tc")),
                        "Vc" => Some(common::input(case, "Vc")),
                        "dipole" => Some(common::input(case, "dipole")),
                        "kappa" => Some(common::input(case, "kappa")),
                        "T" => Some(common::input(case, "T")),
                        "k" => Some(result.k.value),
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
