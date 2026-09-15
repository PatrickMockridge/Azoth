//! Spec-driven tests for `eos.hayduk_minhas_diffusivity`.

use azoth_core::units::{cubic_meters_per_mole, kelvins, pascal_seconds};
use azoth_eos::{HaydukMinhasForm, hayduk_minhas_diffusivity, spec_gen};
use azoth_test_support as common;

const CALC_ID: &str = "eos.hayduk_minhas_diffusivity";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::HaydukMinhasDiffusivityResult {
    let form: HaydukMinhasForm = common::input_str(case, "form").parse().unwrap();
    hayduk_minhas_diffusivity(
        form,
        cubic_meters_per_mole(common::input(case, "VA")),
        kelvins(common::input(case, "T")),
        pascal_seconds(common::input(case, "eta")),
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
                        "VA" => Some(common::input(case, "VA")),
                        "T" => Some(common::input(case, "T")),
                        "eta" => Some(common::input(case, "eta")),
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
