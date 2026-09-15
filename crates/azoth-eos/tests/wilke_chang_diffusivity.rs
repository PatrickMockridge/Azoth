//! Spec-driven tests for `eos.wilke_chang_diffusivity`.

use azoth_core::units::{cubic_meters_per_mole, kelvins, kilograms_per_mole, pascal_seconds};
use azoth_eos::spec_gen;
use azoth_eos::wilke_chang_diffusivity;
use azoth_test_support as common;

const CALC_ID: &str = "eos.wilke_chang_diffusivity";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::WilkeChangDiffusivityResult {
    wilke_chang_diffusivity(
        common::input(case, "phi"),
        kilograms_per_mole(common::input(case, "M")),
        kelvins(common::input(case, "T")),
        pascal_seconds(common::input(case, "eta")),
        cubic_meters_per_mole(common::input(case, "VA")),
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
                    result.d.value,
                    common::expected(case, "d"),
                    case.tolerance,
                    &format!("{}::{} (d)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "phi" => Some(common::input(case, "phi")),
                        "M" => Some(common::input(case, "M")),
                        "T" => Some(common::input(case, "T")),
                        "eta" => Some(common::input(case, "eta")),
                        "VA" => Some(common::input(case, "VA")),
                        "d" => Some(result.d.value),
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
