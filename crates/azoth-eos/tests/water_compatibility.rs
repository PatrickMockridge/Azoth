//! Spec-driven tests for `eos.water_compatibility`.

use azoth_core::spec::{RecordField, TestCase};
use azoth_core::units::{kelvins, kilograms_per_cubic_meter, pascals};
use azoth_eos::model_gen;
use azoth_eos::results::{WaterCompatibilityResult, WaterMixingResult};
use azoth_eos::water_compatibility::{WaterAnalysis, water_compatibility};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.water_compatibility";

fn analysis(row: &[RecordField]) -> WaterAnalysis {
    let number = |name: &str| common::record_number(row, name).expect("every field is stated");
    WaterAnalysis {
        calcium: kilograms_per_cubic_meter(number("calcium")),
        barium: kilograms_per_cubic_meter(number("barium")),
        strontium: kilograms_per_cubic_meter(number("strontium")),
        iron: kilograms_per_cubic_meter(number("iron")),
        magnesium: kilograms_per_cubic_meter(number("magnesium")),
        sodium: kilograms_per_cubic_meter(number("sodium")),
        bicarbonate: kilograms_per_cubic_meter(number("bicarbonate")),
        sulphate: kilograms_per_cubic_meter(number("sulphate")),
        total_dissolved_solids: kilograms_per_cubic_meter(number("total_dissolved_solids")),
        temperature: kelvins(number("temperature")),
        pressure: pascals(number("pressure")),
        co2_partial_pressure: pascals(number("co2_partial_pressure")),
    }
}

fn call(case: &TestCase) -> WaterCompatibilityResult {
    let formation: Vec<WaterAnalysis> = common::record_input(case, "formation_water")
        .iter()
        .map(|row| analysis(row))
        .collect();
    let injection: Vec<WaterAnalysis> = common::record_input(case, "injection_water")
        .iter()
        .map(|row| analysis(row))
        .collect();
    let percents = case
        .vector("injection_water_percent")
        .expect("every case states the percentages");
    water_compatibility(&formation, &injection, percents)
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

fn assert_row(expected: &[RecordField], row: &WaterMixingResult, at: &str, tol: f64) {
    common::assert_close(
        row.injection_water_percent,
        common::record_number(expected, "injection_water_percent")
            .expect("every row states its percentage"),
        tol,
        &format!("{at}::injection_water_percent"),
    );
    for (label, got) in [
        ("si_caco3", row.si_caco3),
        ("si_baso4", row.si_baso4),
        ("si_srso4", row.si_srso4),
        ("si_caso4", row.si_caso4),
        ("si_feco3", row.si_feco3),
    ] {
        match (got, common::record_number(expected, label)) {
            (Some(got), Some(want)) => {
                common::assert_close(got, want, tol, &format!("{at}::{label}"))
            }
            (None, None) => {}
            (got, want) => panic!("{at}::{label}: {got:?} against {want:?}"),
        }
    }
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        for field in ["worst_case_percent", "worst_case_saturation_index"] {
            if let Some(expected) = case.expected_value(field) {
                let actual = match field {
                    "worst_case_percent" => result.worst_case_percent,
                    _ => result.worst_case_saturation_index,
                };
                common::assert_close(
                    actual.unwrap_or_else(|| panic!("{context}: {field} is absent")),
                    expected,
                    case.tolerance,
                    &format!("{context} ({field})"),
                );
            }
        }
        if let Some(expected) = case.expected_string("worst_case_scale") {
            assert_eq!(
                result.worst_case_scale.as_deref(),
                Some(expected),
                "{context}: the worst-case mineral"
            );
        }
        let expected = common::expected_record(case, "mixtures");
        assert_eq!(
            result.mixtures.len(),
            expected.len(),
            "{context}: the mixture count"
        );
        for (index, (row, declared)) in result.mixtures.iter().zip(expected).enumerate() {
            assert_row(
                declared,
                row,
                &format!("{context}::mixtures[{index}]"),
                case.tolerance,
            );
        }
        common::assert_consistent(&result, context);
    }
}

/// **A water is exactly one row.** A list of two is refused rather than read for its first
/// entry, because the choice of which one to read would be this port's and not the class's.
#[test]
fn a_two_row_water_is_refused() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    let case = &spec.cases[0];
    let formation: Vec<WaterAnalysis> = common::record_input(case, "formation_water")
        .iter()
        .map(|row| analysis(row))
        .collect();
    let doubled = [formation[0].clone(), formation[0].clone()];
    let error = water_compatibility(&doubled, &formation, &[0.0])
        .expect_err("two formation waters are not one");
    assert!(error.to_string().contains("formation_water"));
}

/// **A percentage outside the class's range is mixed as stated and warned about.** The class
/// carries no clamp, so the model reproduces the number and says it is outside the range.
#[test]
fn a_percentage_outside_the_range_warns_and_still_mixes() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    let case = &spec.cases[0];
    let formation: Vec<WaterAnalysis> = common::record_input(case, "formation_water")
        .iter()
        .map(|row| analysis(row))
        .collect();
    let injection: Vec<WaterAnalysis> = common::record_input(case, "injection_water")
        .iter()
        .map(|row| analysis(row))
        .collect();
    let result = water_compatibility(&formation, &injection, &[150.0]).expect("the state runs");
    assert_eq!(result.mixtures.len(), 1);
    assert!(
        result
            .warnings
            .iter()
            .any(|warning| warning.field.is_none())
    );
}
