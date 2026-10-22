//! Spec-driven tests for `eos.scale_prediction`.

use azoth_core::units::{kelvins, kilograms_per_cubic_meter, pascals};
use azoth_eos::scale_prediction::scale_prediction;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.scale_prediction";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::results::ScalePredictionResult {
    let mass = |name: &str| kilograms_per_cubic_meter(common::input(case, name));
    scale_prediction(
        mass("calcium"),
        mass("barium"),
        mass("strontium"),
        mass("iron"),
        mass("magnesium"),
        mass("sodium"),
        mass("bicarbonate"),
        mass("sulphate"),
        mass("total_dissolved_solids"),
        kelvins(common::input(case, "temperature")),
        pascals(common::input(case, "pressure")),
        pascals(common::input(case, "co2_partial_pressure")),
        common::input(case, "ph"),
        case.flag("auto_ph").expect("the case states auto_ph"),
        case.flag("second_order_pressure_correction")
            .expect("the case states second_order_pressure_correction"),
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
                executed += 1;
                let result = call(case);
                let context = &format!("{}::{}", spec.id, case.id);
                for (field, value) in [
                    ("ionic_strength", result.ionic_strength.value),
                    ("debye_huckel_a", result.debye_huckel_a),
                    (
                        "divalent_activity_coefficient",
                        result.divalent_activity_coefficient,
                    ),
                    (
                        "monovalent_activity_coefficient",
                        result.monovalent_activity_coefficient,
                    ),
                    ("ph_used", result.ph_used),
                    ("free_calcium", result.free_calcium.value),
                    ("free_sulphate", result.free_sulphate.value),
                    ("free_carbonate", result.free_carbonate.value),
                    ("free_bicarbonate", result.free_bicarbonate.value),
                    ("total_barium", result.total_barium.value),
                    ("total_strontium", result.total_strontium.value),
                    ("total_iron", result.total_iron.value),
                    ("total_sodium", result.total_sodium.value),
                    ("ksp_calcite", result.ksp_calcite),
                    ("ksp_barite", result.ksp_barite),
                    ("ksp_celestite", result.ksp_celestite),
                    ("ksp_anhydrite", result.ksp_anhydrite),
                    ("ksp_siderite", result.ksp_siderite),
                    ("log10_ksp_calcite", result.log10_ksp_calcite),
                    ("log10_ksp_barite", result.log10_ksp_barite),
                    ("log10_ksp_celestite", result.log10_ksp_celestite),
                    ("log10_ksp_anhydrite", result.log10_ksp_anhydrite),
                    ("log10_ksp_siderite", result.log10_ksp_siderite),
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
                for (field, value) in [
                    ("si_caco3", result.si_caco3),
                    ("si_baso4", result.si_baso4),
                    ("si_srso4", result.si_srso4),
                    ("si_caso4", result.si_caso4),
                    ("si_feco3", result.si_feco3),
                ] {
                    // An unlisted index is not pinned; `an_absent_ion_gives_an_absent_index`
                    // holds the absence rule.
                    if let Some(expected) = case.expected_value(field) {
                        assert_eq!(
                            value,
                            Some(expected),
                            "{context}: `{field}` should be {expected}"
                        );
                    }
                }
                if let Some(expected) = case.expected_value("scale_risk_count") {
                    assert_eq!(
                        result.scale_risk_count, expected as u32,
                        "{context}: the risk count"
                    );
                    assert_eq!(
                        result.has_scaling_risk,
                        result.scale_risk_count > 0,
                        "{context}: the flag follows the count"
                    );
                }
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "temperature" => Some(common::input(case, "temperature")),
                        "pressure" => Some(common::input(case, "pressure")),
                        "co2_partial_pressure" => Some(common::input(case, "co2_partial_pressure")),
                        "calcium" => Some(common::input(case, "calcium")),
                        "barium" => Some(common::input(case, "barium")),
                        "strontium" => Some(common::input(case, "strontium")),
                        "iron" => Some(common::input(case, "iron")),
                        "magnesium" => Some(common::input(case, "magnesium")),
                        "sodium" => Some(common::input(case, "sodium")),
                        "bicarbonate" => Some(common::input(case, "bicarbonate")),
                        "sulphate" => Some(common::input(case, "sulphate")),
                        "total_dissolved_solids" => {
                            Some(common::input(case, "total_dissolved_solids"))
                        }
                        "ph" => Some(common::input(case, "ph")),
                        "scale_risk_count" => Some(f64::from(result.scale_risk_count)),
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

/// **The auto-pH path is the class's, defect included.** A manual pH is the physical route, and
/// the port's estimate is about -45 at these conditions; this test says so rather than letting
/// a future "fix" pass quietly.
#[test]
fn auto_ph_is_the_class_defect() {
    let manual = scale_prediction(
        kilograms_per_cubic_meter(0.4),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.15),
        kilograms_per_cubic_meter(0.01),
        kilograms_per_cubic_meter(35.0),
        kelvins(353.15),
        pascals(1.0e7),
        pascals(2.0e5),
        6.5,
        false,
        false,
    )
    .expect("the state runs");
    let automatic = scale_prediction(
        kilograms_per_cubic_meter(0.4),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.15),
        kilograms_per_cubic_meter(0.01),
        kilograms_per_cubic_meter(35.0),
        kelvins(353.15),
        pascals(1.0e7),
        pascals(2.0e5),
        6.5,
        true,
        false,
    )
    .expect("the state runs");
    assert!((automatic.ph_used - -45.30001384171815).abs() < 1e-9);
    assert!(automatic.free_carbonate.value < 1e-50);
    assert!(automatic.si_caco3.expect("calcite has both ions") < -50.0);
    assert!(manual.si_caco3.expect("calcite has both ions") > -1.0);
}

/// **An absent ion gives an absent index.** The class returns NaN where the mineral's ion is
/// missing, and the result carries no value there.
#[test]
fn an_absent_ion_gives_an_absent_index() {
    let result = scale_prediction(
        kilograms_per_cubic_meter(0.4),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.15),
        kilograms_per_cubic_meter(0.01),
        kilograms_per_cubic_meter(35.0),
        kelvins(353.15),
        pascals(1.0e7),
        pascals(2.0e5),
        6.5,
        false,
        false,
    )
    .expect("the state runs");
    assert!(result.si_caco3.is_some());
    assert!(result.si_caso4.is_some());
    assert!(result.si_baso4.is_none());
    assert!(result.si_srso4.is_none());
    assert!(result.si_feco3.is_none());
}

/// **The pressure correction begins just above one atmosphere.** At 1.013 bara the class
/// bypasses it; at 1.014 bara it applies, which is a discontinuity of its own making.
#[test]
fn the_pressure_bypass_is_at_1_013_bara() {
    let run = |pressure: f64| {
        scale_prediction(
            kilograms_per_cubic_meter(0.1),
            kilograms_per_cubic_meter(0.2),
            kilograms_per_cubic_meter(0.0),
            kilograms_per_cubic_meter(0.0),
            kilograms_per_cubic_meter(0.0),
            kilograms_per_cubic_meter(0.0),
            kilograms_per_cubic_meter(0.1),
            kilograms_per_cubic_meter(2.0),
            kilograms_per_cubic_meter(40.0),
            kelvins(333.15),
            pascals(pressure),
            pascals(0.5e5),
            6.5,
            false,
            false,
        )
        .expect("the state runs")
    };
    let at_base = run(1.013e5);
    let above = run(1.014e5);
    assert!((at_base.ksp_barite - 2.4740496629926317e-10).abs() < 1e-20);
    assert!((above.ksp_barite - 2.474052771229189e-10).abs() < 1e-20);
    assert!(above.ksp_barite > at_base.ksp_barite);
}
