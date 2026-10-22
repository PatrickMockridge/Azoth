//! Spec-driven tests for `eos.scale_mass`.

use azoth_core::units::moles_per_cubic_meter;
use azoth_eos::scale_mass::scale_mass;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.scale_mass";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::results::ScaleMassResult {
    scale_mass(
        common::input_str(case, "mineral"),
        moles_per_cubic_meter(common::input(case, "cation_concentration")),
        moles_per_cubic_meter(common::input(case, "anion_concentration")),
        common::input(case, "saturation_index"),
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
                let context = &format!("{}::{}", spec.id, case.id);
                for (field, value) in [
                    ("precipitated_mass", result.precipitated_mass.value),
                    ("precipitated_moles", result.precipitated_moles.value),
                    ("fraction_excess", result.fraction_excess),
                    (
                        "limiting_concentration",
                        result.limiting_concentration.value,
                    ),
                ] {
                    common::assert_close(
                        value,
                        common::expected(case, field),
                        case.tolerance,
                        &format!("{context} ({field})"),
                    );
                }
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "cation_concentration" => Some(common::input(case, "cation_concentration")),
                        "anion_concentration" => Some(common::input(case, "anion_concentration")),
                        _ => None,
                    },
                    context,
                );
                common::assert_consistent(&result, context);
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

/// **The 0.99 cap binds above SI 2.** A state twice as supersaturated gives the same mass.
#[test]
fn the_cap_makes_every_high_index_equal() {
    let run = |si: f64| {
        scale_mass(
            "baso4",
            moles_per_cubic_meter(10.0),
            moles_per_cubic_meter(20.0),
            si,
        )
        .expect("the state runs")
    };
    assert_eq!(run(2.5).precipitated_mass, run(5.0).precipitated_mass);
    assert_eq!(run(2.5).fraction_excess, 0.99);
}

#[test]
fn an_unknown_mineral_is_refused() {
    let error = scale_mass(
        "gypsum",
        moles_per_cubic_meter(10.0),
        moles_per_cubic_meter(20.0),
        2.5,
    )
    .expect_err("gypsum is not one of the class's five arms");
    assert!(error.to_string().contains("caco3"), "{error}");
}
