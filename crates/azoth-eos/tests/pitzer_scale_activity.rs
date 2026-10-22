//! Spec-driven tests for `eos.pitzer_scale_activity`.

use azoth_core::units::kelvins;
use azoth_eos::pitzer_scale_activity::pitzer_scale_activity;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.pitzer_scale_activity";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::results::PitzerScaleActivityResult {
    pitzer_scale_activity(
        kelvins(common::input(case, "temperature")),
        common::input(case, "ionic_strength_molal"),
        common::input(case, "debye_huckel_a"),
        common::input_str(case, "salt"),
        common::input(case, "salt_molality"),
        common::input_str(case, "ion"),
        common::input_str(case, "cation"),
        common::input_str(case, "anion"),
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
                    (
                        "mean_activity_coefficient",
                        result.mean_activity_coefficient,
                    ),
                    ("ion_activity_coefficient", result.ion_activity_coefficient),
                    (
                        "activity_coefficient_product",
                        result.activity_coefficient_product,
                    ),
                    ("gamma_sodium", result.gamma_sodium),
                    ("gamma_chloride", result.gamma_chloride),
                    ("gamma_calcium", result.gamma_calcium),
                    ("gamma_barium", result.gamma_barium),
                    ("gamma_strontium", result.gamma_strontium),
                    ("gamma_ferrous", result.gamma_ferrous),
                    ("gamma_sulphate", result.gamma_sulphate),
                    ("gamma_carbonate", result.gamma_carbonate),
                    ("gamma_bicarbonate", result.gamma_bicarbonate),
                ] {
                    if let Some(expected) = case.expected_value(field) {
                        common::assert_close(
                            value,
                            expected,
                            case.tolerance,
                            &format!("{context} ({field})"),
                        );
                    }
                }
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "temperature" => Some(common::input(case, "temperature")),
                        "ionic_strength_molal" => Some(common::input(case, "ionic_strength_molal")),
                        "salt_molality" => Some(common::input(case, "salt_molality")),
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

/// **The mean coefficient ignores the configured strength.** It computes its own ionic
/// strength from the molality, so two calls with the same molality and different strengths
/// agree on the mean and differ on the mapped ions.
#[test]
fn the_mean_uses_its_own_strength() {
    let run = |strength: f64| {
        pitzer_scale_activity(
            kelvins(298.15),
            strength,
            0.5107140231833494,
            "nacl",
            1.0,
            "sodium",
            "sodium",
            "chloride",
        )
        .expect("the state runs")
    };
    let dilute = run(0.01);
    let concentrated = run(1.0);
    assert_eq!(
        dilute.mean_activity_coefficient,
        concentrated.mean_activity_coefficient
    );
    assert_ne!(dilute.gamma_calcium, concentrated.gamma_calcium);
}

#[test]
fn an_unknown_salt_or_ion_is_refused() {
    let salt_error = pitzer_scale_activity(
        kelvins(298.15),
        1.0,
        0.5107140231833494,
        "kcl",
        1.0,
        "sodium",
        "sodium",
        "chloride",
    )
    .expect_err("kcl is not one of the eight binary salts");
    assert!(salt_error.to_string().contains("eight binary salts"));

    let ion_error = pitzer_scale_activity(
        kelvins(298.15),
        1.0,
        0.5107140231833494,
        "nacl",
        1.0,
        "lithium",
        "sodium",
        "chloride",
    )
    .expect_err("lithium is not one of the nine ions");
    assert!(ion_error.to_string().contains("nine ions"));
}
