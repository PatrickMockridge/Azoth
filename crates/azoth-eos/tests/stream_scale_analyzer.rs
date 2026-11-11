//! Spec-driven tests for `eos.stream_scale_analyzer`.

use azoth_core::spec::{RecordField, TestCase};
use azoth_core::units::{cubic_meters_per_second, kelvins, kilograms_per_cubic_meter, pascals};
use azoth_eos::model_gen;
use azoth_eos::results::{StreamScaleAnalyzerResult, StreamScaleMineralRow};
use azoth_eos::stream_scale_analyzer::stream_scale_analyzer;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.stream_scale_analyzer";

fn call(case: &TestCase) -> StreamScaleAnalyzerResult {
    let mass = |name: &str| kilograms_per_cubic_meter(common::input(case, name));
    stream_scale_analyzer(
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
        common::input_str(case, "activity_model"),
        cubic_meters_per_second(common::input(case, "water_flow")),
        case.input("bdot_ion_size"),
        case.input("bdot_b"),
        case.input("pitzer_ionic_strength_molal"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

fn assert_row(expected: &[RecordField], row: &StreamScaleMineralRow, at: &str, tol: f64) {
    if let Some(mineral) = common::record_name(expected, "mineral") {
        assert_eq!(row.mineral, mineral, "{at}: the mineral");
    }
    for (field, got) in [
        ("initial_si", row.initial_si),
        ("final_si", row.final_si),
        ("precipitated_mass", row.precipitated_mass.value),
        ("scale_rate", row.scale_rate.value),
    ] {
        if let Some(want) = common::record_number(expected, field) {
            common::assert_close(got, want, tol, &format!("{at}::{field}"));
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
        for (field, got) in [
            ("total_scale_mass", result.total_scale_mass.value),
            ("total_scale_rate", result.total_scale_rate.value),
        ] {
            if let Some(want) = case.expected_value(field) {
                common::assert_close(got, want, case.tolerance, &format!("{context} ({field})"));
            }
        }
        if let Some(expected) = case.expected_string("dominant_scale") {
            assert_eq!(
                result.dominant_scale, expected,
                "{context}: the dominant mineral"
            );
        }
        assert_eq!(
            result.has_scaling,
            result.total_scale_mass.value > 0.0,
            "{context}: the flag follows the total"
        );

        let expected = case.record("minerals").unwrap_or(&[]);
        if case.record("minerals").is_some() {
            assert_eq!(
                result.minerals.len(),
                expected.len(),
                "{context}: the mineral count"
            );
        }
        for (index, (row, declared)) in result.minerals.iter().zip(expected).enumerate() {
            assert_row(
                declared,
                row,
                &format!("{context}::minerals[{index}]"),
                case.tolerance,
            );
        }
        common::assert_consistent(&result, context);
    }
}

/// **Every rate is the precipitated mass times the throughput**, per mineral and in total,
/// so the class's mg/L and L/day cancel to the same number.
#[test]
fn the_rate_is_the_mass_times_the_flow() {
    let result = stream_scale_analyzer(
        kilograms_per_cubic_meter(1.0),
        kilograms_per_cubic_meter(1.0),
        kilograms_per_cubic_meter(3.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.2),
        kilograms_per_cubic_meter(1.5),
        kilograms_per_cubic_meter(60.0),
        kelvins(343.15),
        pascals(1.0e7),
        pascals(1.0e5),
        0.0,
        false,
        false,
        "davies",
        cubic_meters_per_second(0.0011574074074074073),
        None,
        None,
        None,
    )
    .expect("the state runs");
    let flow = 0.0011574074074074073;
    for row in &result.minerals {
        common::assert_close(
            row.scale_rate.value,
            row.precipitated_mass.value * flow,
            1e-15,
            &row.mineral,
        );
    }
    common::assert_close(
        result.total_scale_rate.value,
        result.total_scale_mass.value * flow,
        1e-15,
        "the total",
    );
}

/// **Nothing precipitated reports `none`**, not the first registered mineral.
#[test]
fn a_clean_water_has_no_dominant_mineral() {
    let result = stream_scale_analyzer(
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.001),
        kelvins(298.15),
        pascals(1.013e5),
        pascals(1.0e5),
        7.0,
        false,
        false,
        "davies",
        cubic_meters_per_second(1.1574074074074073e-5),
        None,
        None,
        None,
    )
    .expect("the state runs");
    assert!(result.minerals.is_empty());
    assert_eq!(result.dominant_scale, "none");
    assert!(!result.has_scaling);
    assert_eq!(result.total_scale_rate.value, 0.0);
}
