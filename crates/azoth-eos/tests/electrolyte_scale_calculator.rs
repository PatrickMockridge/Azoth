//! Spec-driven tests for `eos.electrolyte_scale_calculator`.

use azoth_core::units::{kelvins, kilograms_per_cubic_meter, pascals};
use azoth_eos::electrolyte_scale_calculator::electrolyte_scale_calculator;
use azoth_eos::results::ElectrolyteScaleCalculatorResult;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.electrolyte_scale_calculator";

fn call(case: &azoth_core::spec::TestCase) -> ElectrolyteScaleCalculatorResult {
    let mass = |name: &str| kilograms_per_cubic_meter(common::input(case, name));
    electrolyte_scale_calculator(
        mass("calcium"),
        mass("barium"),
        mass("strontium"),
        mass("magnesium"),
        mass("sodium"),
        mass("potassium"),
        mass("iron"),
        mass("chloride"),
        mass("sulphate"),
        mass("bicarbonate"),
        mass("carbonate"),
        kelvins(common::input(case, "temperature")),
        pascals(common::input(case, "pressure")),
        common::input(case, "ph"),
    )
    .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
}

fn seawater(temperature_c: f64, ph: f64) -> ElectrolyteScaleCalculatorResult {
    let mas = |value: f64| kilograms_per_cubic_meter(value * 1e-3);
    electrolyte_scale_calculator(
        mas(412.0),
        mas(0.02),
        mas(7.9),
        mas(1290.0),
        mas(10780.0),
        mas(399.0),
        mas(0.002),
        mas(19350.0),
        mas(2710.0),
        mas(142.0),
        mas(0.0),
        kelvins(temperature_c + 273.15),
        pascals(1.013e5),
        ph,
    )
    .expect("the seawater state runs")
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
                    ("debye_hueckel_a", result.debye_hueckel_a),
                    (
                        "divalent_activity_coefficient",
                        result.divalent_activity_coefficient,
                    ),
                    (
                        "monovalent_activity_coefficient",
                        result.monovalent_activity_coefficient,
                    ),
                    ("si_caco3", result.si_caco3),
                    ("si_baso4", result.si_baso4),
                    ("si_caso4", result.si_caso4),
                    ("si_srso4", result.si_srso4),
                    ("si_brucite", result.si_brucite),
                    ("hydroxide_activity", result.hydroxide_activity),
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
                        "pressure" => Some(common::input(case, "pressure")),
                        "calcium" => Some(common::input(case, "calcium")),
                        "barium" => Some(common::input(case, "barium")),
                        "strontium" => Some(common::input(case, "strontium")),
                        "magnesium" => Some(common::input(case, "magnesium")),
                        "sodium" => Some(common::input(case, "sodium")),
                        "potassium" => Some(common::input(case, "potassium")),
                        "iron" => Some(common::input(case, "iron")),
                        "chloride" => Some(common::input(case, "chloride")),
                        "sulphate" => Some(common::input(case, "sulphate")),
                        "bicarbonate" => Some(common::input(case, "bicarbonate")),
                        "carbonate" => Some(common::input(case, "carbonate")),
                        "ph" => Some(common::input(case, "ph")),
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

/// **Brucite's index rises two decades per pH unit**, because two hydroxides enter the
/// activity product; the upstream test's own property.
#[test]
fn brucite_rises_two_decades_per_ph_unit() {
    let low = seawater(10.0, 9.0).si_brucite;
    let high = seawater(10.0, 10.0).si_brucite;
    assert!((high - low - 2.0).abs() < 0.05, "{low} -> {high}");
}

/// **The hydroxide activity is the water ion product**, so pH 14 at 25 C is 1.0 within the
/// class's own `pKw` fit.
#[test]
fn hydroxide_activity_matches_water_ion_product() {
    let result = seawater(25.0, 14.0);
    assert!((result.hydroxide_activity - 1.0).abs() < 0.05);
}

/// **Brucite solubility is retrograde**: the index rises with temperature at fixed pH.
#[test]
fn brucite_solubility_is_retrograde() {
    let cold = seawater(5.0, 10.0).si_brucite;
    let warm = seawater(40.0, 10.0).si_brucite;
    assert!(warm > cold, "{cold} -> {warm}");
}

/// **The `-10` sentinel is unreachable.** An absent ion takes the `1e-15` floor rather than
/// the sentinel branch, so its index is a large negative number and not exactly `-10`.
#[test]
fn the_sentinel_is_unreachable() {
    let mas = |value: f64| kilograms_per_cubic_meter(value * 1e-3);
    let result = electrolyte_scale_calculator(
        mas(400.0),
        mas(0.0),
        mas(0.0),
        mas(0.0),
        mas(10000.0),
        mas(0.0),
        mas(0.0),
        mas(19000.0),
        mas(0.0),
        mas(0.0),
        mas(0.0),
        kelvins(298.15),
        pascals(1.013e5),
        7.0,
    )
    .expect("the state runs");
    for (label, value) in [
        ("si_caco3", result.si_caco3),
        ("si_baso4", result.si_baso4),
        ("si_caso4", result.si_caso4),
        ("si_srso4", result.si_srso4),
    ] {
        assert_ne!(value, -10.0, "{label} took the unreachable sentinel");
        assert!(value < -5.0, "{label} = {value}");
    }
}

/// **The pressure correction begins just above one atmosphere**, where the class's guard
/// stops bypassing it; a larger Ksp means a smaller index.
#[test]
fn the_pressure_correction_runs_above_1_013_bara() {
    let run = |pressure: f64| {
        let mas = |value: f64| kilograms_per_cubic_meter(value * 1e-3);
        electrolyte_scale_calculator(
            mas(420.0),
            mas(0.05),
            mas(8.0),
            mas(1280.0),
            mas(10800.0),
            mas(400.0),
            mas(0.0),
            mas(19400.0),
            mas(2700.0),
            mas(142.0),
            mas(0.0),
            kelvins(333.15),
            pascals(pressure),
            6.5,
        )
        .expect("the state runs")
    };
    let at_base = run(1.013e5);
    let above = run(1.014e5);
    assert!(above.si_caco3 < at_base.si_caco3);
}
