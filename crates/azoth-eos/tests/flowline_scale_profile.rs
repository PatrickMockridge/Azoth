//! Spec-driven tests for `eos.flowline_scale_profile`.

use azoth_core::spec::{RecordField, TestCase};
use azoth_core::units::{kelvins, kilograms_per_cubic_meter, pascals};
use azoth_eos::flowline_scale_profile::flowline_scale_profile;
use azoth_eos::model_gen;
use azoth_eos::results::{FlowlineScaleProfileResult, FlowlineScaleSegment};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.flowline_scale_profile";

fn call(case: &TestCase) -> FlowlineScaleProfileResult {
    let input = |name: &str| common::input(case, name);
    let mass = |name: &str| kilograms_per_cubic_meter(input(name));
    flowline_scale_profile(
        mass("calcium"),
        mass("barium"),
        mass("strontium"),
        mass("iron"),
        mass("magnesium"),
        mass("sodium"),
        mass("bicarbonate"),
        mass("sulphate"),
        mass("total_dissolved_solids"),
        pascals(input("co2_partial_pressure")),
        case.flag("auto_ph").unwrap_or(false),
        kelvins(input("inlet_temperature")),
        pascals(input("inlet_pressure")),
        kelvins(input("outlet_temperature")),
        pascals(input("outlet_pressure")),
        input("segments") as usize,
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

fn assert_row(expected: &[RecordField], row: &FlowlineScaleSegment, at: &str, tol: f64) {
    common::assert_close(
        row.distance_fraction,
        common::record_number(expected, "distance_fraction")
            .expect("every row states its fraction"),
        tol,
        &format!("{at}::distance_fraction"),
    );
    common::assert_close(
        row.temperature.value,
        common::record_number(expected, "temperature").expect("every row states its temperature"),
        tol,
        &format!("{at}::temperature"),
    );
    common::assert_close(
        row.pressure.value,
        common::record_number(expected, "pressure").expect("every row states its pressure"),
        tol,
        &format!("{at}::pressure"),
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
        for field in [
            "max_caco3",
            "max_baso4",
            "max_srso4",
            "max_caso4",
            "max_feco3",
        ] {
            if let Some(expected) = case.expected_value(field) {
                let actual = match field {
                    "max_caco3" => result.max_caco3,
                    "max_baso4" => result.max_baso4,
                    "max_srso4" => result.max_srso4,
                    "max_caso4" => result.max_caso4,
                    _ => result.max_feco3,
                };
                common::assert_close(
                    actual.unwrap_or_else(|| panic!("{context}: {field} is absent")),
                    expected,
                    case.tolerance,
                    &format!("{context} ({field})"),
                );
            }
        }
        let expected = common::expected_record(case, "segments");
        assert_eq!(
            result.segments.len(),
            expected.len(),
            "{context}: the point count"
        );
        for (index, (row, declared)) in result.segments.iter().zip(expected).enumerate() {
            assert_row(
                declared,
                row,
                &format!("{context}::segments[{index}]"),
                case.tolerance,
            );
        }
        common::assert_consistent(&result, context);
    }
}

/// **One segment is clamped to two**, so the grid is three points rather than two: the inlet,
/// the midpoint and the outlet.
#[test]
fn a_single_segment_is_clamped_to_two() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    let case = spec
        .cases
        .iter()
        .find(|case| case.id == "one_segment_is_clamped_to_two")
        .expect("the clamp case");
    let result = call(case);
    assert_eq!(result.segments.len(), 3);
    assert!((result.segments[1].distance_fraction - 0.5).abs() < 1e-12);
}

/// **The indices follow the predictor's absent-ion rule**: a water with no barium or strontium
/// has no barite or celestite index at any point, and a mineral no point defines has no
/// maximum rather than the class's `-infinity`.
#[test]
fn absent_ions_leave_absent_indices_and_no_maximum() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    let case = spec
        .cases
        .iter()
        .find(|case| case.id == "high_sulphate_120c")
        .expect("the high-sulphate case");
    let result = call(case);
    assert!(result.segments.iter().all(|row| row.si_baso4.is_none()));
    assert!(result.max_baso4.is_none());
    assert!(result.max_feco3.is_none());
    assert!(result.max_caso4.is_some_and(|value| value > 0.0));
}
