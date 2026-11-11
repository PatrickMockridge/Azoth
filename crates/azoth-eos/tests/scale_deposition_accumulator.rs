//! Spec-driven tests for `eos.scale_deposition_accumulator`.

use azoth_core::units::{
    kilograms_per_cubic_meter, kilograms_per_square_meter_second, meters, seconds,
};
use azoth_eos::model_gen;
use azoth_eos::results::ScaleDepositionAccumulatorResult;
use azoth_eos::scale_deposition_accumulator::scale_deposition_accumulator;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.scale_deposition_accumulator";

fn call(case: &azoth_core::spec::TestCase) -> ScaleDepositionAccumulatorResult {
    let mass = |name: &str| kilograms_per_cubic_meter(common::input(case, name));
    scale_deposition_accumulator(
        mass("calcium"),
        mass("bicarbonate"),
        mass("sulphate"),
        mass("barium"),
        mass("sodium"),
        mass("total_dissolved_solids"),
        common::input(case, "ph"),
        meters(common::input(case, "diameter")),
        meters(common::input(case, "length")),
        case.vector("segment_temperatures")
            .expect("the case states its temperatures"),
        case.vector("segment_pressures")
            .expect("the case states its pressures"),
        seconds(common::input(case, "service_time")),
        kilograms_per_cubic_meter(common::input(case, "scale_density")),
        kilograms_per_square_meter_second(common::input(case, "base_rate")),
        common::input(case, "inhibitor_efficiency"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        for (field, got) in [
            ("total_deposition_mass", result.total_deposition_mass.value),
            ("max_thickness", result.max_thickness.value),
            ("time_to_blockage", result.time_to_blockage.value),
        ] {
            if let Some(want) = case.expected_value(field) {
                common::assert_close(got, want, case.tolerance, &format!("{context} ({field})"));
            }
        }
        let expected = case.record("segments").unwrap_or(&[]);
        if case.record("segments").is_some() {
            assert_eq!(
                result.segments.len(),
                expected.len(),
                "{context}: the segment count"
            );
        }
        for (index, (row, declared)) in result.segments.iter().zip(expected).enumerate() {
            for (field, got) in [
                ("saturation_index", row.saturation_index),
                ("deposition_mass", row.deposition_mass.value),
                ("thickness", row.thickness.value),
            ] {
                if let Some(want) = common::record_number(declared, field) {
                    common::assert_close(
                        got,
                        want,
                        case.tolerance,
                        &format!("{context}: segments[{index}] {field}"),
                    );
                }
            }
        }
        common::assert_consistent(&result, context);
    }
}

/// **The segment count is the shorter profile**, so a one-value temperature list with an
/// eleven-value pressure list evaluates one segment over the whole length - the class's own
/// `min` of the two lists.
#[test]
fn the_shorter_profile_governs() {
    let result = scale_deposition_accumulator(
        kilograms_per_cubic_meter(1.0),
        kilograms_per_cubic_meter(0.5),
        kilograms_per_cubic_meter(0.1),
        kilograms_per_cubic_meter(0.01),
        kilograms_per_cubic_meter(20.0),
        kilograms_per_cubic_meter(50.0),
        6.5,
        meters(0.1),
        meters(1000.0),
        &[333.15],
        &[1.0e7; 11],
        seconds(31_557_600.0),
        kilograms_per_cubic_meter(2700.0),
        kilograms_per_square_meter_second(1.5844043908606944e-8),
        0.0,
    )
    .expect("the state runs");
    assert_eq!(result.segments.len(), 1);
}

/// **A pipe with no wall area gives an empty result**, and an empty profile does the same.
#[test]
fn a_pipe_without_area_has_no_segments() {
    for (diameter, temperatures) in [(0.0, &[333.15f64][..]), (0.1, &[][..])] {
        let result = scale_deposition_accumulator(
            kilograms_per_cubic_meter(1.0),
            kilograms_per_cubic_meter(0.5),
            kilograms_per_cubic_meter(0.1),
            kilograms_per_cubic_meter(0.01),
            kilograms_per_cubic_meter(20.0),
            kilograms_per_cubic_meter(50.0),
            6.5,
            meters(diameter),
            meters(1000.0),
            temperatures,
            &[1.0e7],
            seconds(31_557_600.0),
            kilograms_per_cubic_meter(2700.0),
            kilograms_per_square_meter_second(1.5844043908606944e-8),
            0.0,
        )
        .expect("the state runs");
        assert!(result.segments.is_empty());
        assert_eq!(result.total_deposition_mass.value, 0.0);
        assert!(result.time_to_blockage.value.is_infinite());
    }
}

/// **Nothing deposited gives an infinite blockage time**, not a division by zero.
#[test]
fn an_undersaturated_segment_never_blocks() {
    let result = scale_deposition_accumulator(
        kilograms_per_cubic_meter(0.01),
        kilograms_per_cubic_meter(0.01),
        kilograms_per_cubic_meter(0.1),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(20.0),
        kilograms_per_cubic_meter(50.0),
        4.0,
        meters(0.1),
        meters(1000.0),
        &[333.15],
        &[1.0e7],
        seconds(31_557_600.0),
        kilograms_per_cubic_meter(2700.0),
        kilograms_per_square_meter_second(1.5844043908606944e-8),
        0.0,
    )
    .expect("the state runs");
    assert_eq!(result.max_thickness.value, 0.0);
    assert!(result.time_to_blockage.value.is_infinite());
}
