//! Spec-driven tests for `eos.scale_inhibitor_performance`.

use azoth_core::units::{kelvins, kilograms_per_cubic_meter};
use azoth_eos::results::ScaleInhibitorPerformanceResult;
use azoth_eos::scale_inhibitor_performance::scale_inhibitor_performance;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.scale_inhibitor_performance";

fn call(case: &azoth_core::spec::TestCase) -> ScaleInhibitorPerformanceResult {
    scale_inhibitor_performance(
        common::input_str(case, "scale_type"),
        common::input_str(case, "inhibitor_chemistry"),
        kelvins(common::input(case, "temperature")),
        common::input(case, "saturation_ratio"),
        kilograms_per_cubic_meter(common::input(case, "total_dissolved_solids")),
        kilograms_per_cubic_meter(common::input(case, "calcium")),
        kilograms_per_cubic_meter(common::input(case, "available_dose")),
    )
    .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
}

fn run(
    scale_type: &str,
    chemistry: &str,
    temperature_c: f64,
    saturation_ratio: f64,
    tds_mg_l: f64,
    calcium_mg_l: f64,
    dose_mg_l: f64,
) -> ScaleInhibitorPerformanceResult {
    scale_inhibitor_performance(
        scale_type,
        chemistry,
        kelvins(temperature_c + 273.15),
        saturation_ratio,
        kilograms_per_cubic_meter(tds_mg_l * 1e-3),
        kilograms_per_cubic_meter(calcium_mg_l * 1e-3),
        kilograms_per_cubic_meter(dose_mg_l * 1e-3),
    )
    .expect("the state runs")
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
                executed += 1;
                let result = call(case);
                let context = &format!("{}::{}", spec.id, case.id);
                for (field, value) in [
                    (
                        "minimum_inhibitor_concentration",
                        result.minimum_inhibitor_concentration.value,
                    ),
                    ("recommended_dose", result.recommended_dose.value),
                    ("efficiency", result.efficiency),
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
                        "saturation_ratio" => Some(common::input(case, "saturation_ratio")),
                        "total_dissolved_solids" => {
                            Some(common::input(case, "total_dissolved_solids"))
                        }
                        "calcium" => Some(common::input(case, "calcium")),
                        "available_dose" => Some(common::input(case, "available_dose")),
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
            other => panic!("{}::{} has an unhandled kind `{other}`", spec.id, case.id),
        }
    }
    assert!(executed >= 2, "the spec should exercise more than one row");
}

/// **The four warnings are the class's predicates**, each on its own threshold, and the
/// calcium penalty belongs to phosphonate alone.
#[test]
fn the_warning_predicates() {
    let hot = run("caco3", "phosphonate", 180.0, 5.0, 50000.0, 1000.0, 0.0);
    assert!(hot.thermal_degradation_risk);

    let calcic = run("caco3", "phosphonate", 60.0, 5.0, 50000.0, 5000.0, 0.0);
    assert!(calcic.calcium_phosphonate_precipitation_risk);

    let acrylic = run("caco3", "polyacrylate", 60.0, 5.0, 50000.0, 5000.0, 0.0);
    assert!(!acrylic.calcium_phosphonate_precipitation_risk);

    let extreme = run("caco3", "phosphonate", 60.0, 200.0, 50000.0, 1000.0, 0.0);
    assert!(extreme.extreme_supersaturation_risk);

    let cold = run("caco3", "phosphate_ester", 30.0, 5.0, 50000.0, 1000.0, 0.0);
    assert!(cold.low_temperature_solubility_risk);
}

/// **The efficiency follows the dose against the MIC**: none below it, and the class's
/// 0.85-at-the-MIC and approaching-one curve above it.
#[test]
fn the_efficiency_follows_the_dose() {
    let none = run("caco3", "phosphonate", 60.0, 10.0, 50000.0, 1000.0, 0.0);
    assert_eq!(none.efficiency, 0.0);
    assert!(!none.adequate);

    let half = run(
        "caco3",
        "phosphonate",
        60.0,
        10.0,
        50000.0,
        1000.0,
        2.951544993495972,
    );
    assert!((half.efficiency - 0.425).abs() < 1e-12);
    assert!(!half.adequate);

    let huge = run("caco3", "phosphonate", 60.0, 10.0, 50000.0, 1000.0, 1.0e6);
    assert!(huge.adequate);
    assert!((huge.efficiency - 1.0).abs() < 1e-12);
}

/// **An unknown scale type or chemistry is refused** rather than defaulted to the class's
/// switch fall-through, which would silently dose the wrong chemistry.
#[test]
fn an_unknown_enum_is_refused() {
    let unknown_scale = scale_inhibitor_performance(
        "cacl2",
        "phosphonate",
        kelvins(333.15),
        5.0,
        kilograms_per_cubic_meter(50.0),
        kilograms_per_cubic_meter(1.0),
        kilograms_per_cubic_meter(0.0),
    )
    .expect_err("`cacl2` is not a scale type");
    assert!(unknown_scale.to_string().contains("scale_type"));

    let unknown_chemistry = scale_inhibitor_performance(
        "caco3",
        "silicate",
        kelvins(333.15),
        5.0,
        kilograms_per_cubic_meter(50.0),
        kilograms_per_cubic_meter(1.0),
        kilograms_per_cubic_meter(0.0),
    )
    .expect_err("`silicate` is not an inhibitor chemistry");
    assert!(
        unknown_chemistry
            .to_string()
            .contains("inhibitor_chemistry")
    );
}
