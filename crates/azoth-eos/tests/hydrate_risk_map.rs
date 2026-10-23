//! Spec-driven tests for `eos.hydrate_risk_map`.

use azoth_core::spec::{FieldValue, RecordField};
use azoth_core::units::kelvin_intervals;
use azoth_eos::hydrate_risk_map::hydrate_risk_map;
use azoth_eos::model_gen;
use azoth_eos::results::{HydrateRiskMapResult, HydrateRiskPoint};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.hydrate_risk_map";

fn call(case: &azoth_core::spec::TestCase) -> HydrateRiskMapResult {
    let components: Vec<String> = common::list_input(case, "components")
        .iter()
        .map(|name| (*name).to_owned())
        .collect();
    hydrate_risk_map(
        &components,
        case.vector("z").expect("the case states its composition"),
        case.vector("distances").expect("the case states distances"),
        case.vector("pressures").expect("the case states pressures"),
        case.vector("temperatures")
            .expect("the case states temperatures"),
        case.string("eos"),
        case.string("hydrate_model"),
        case.input("high_risk_subcooling").map(kelvin_intervals),
        case.input("medium_risk_subcooling").map(kelvin_intervals),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

fn assert_row(expected: &[RecordField], row: &HydrateRiskPoint, at: &str, tol: f64) {
    for (label, got) in [
        ("distance", row.distance.value),
        ("pressure", row.pressure.value),
        ("actual_temperature", row.actual_temperature.value),
    ] {
        common::assert_close(
            got,
            common::record_number(expected, label)
                .unwrap_or_else(|| panic!("{at}::{label} is absent")),
            tol,
            &format!("{at}::{label}"),
        );
    }
    for (label, got) in [
        (
            "hydrate_temperature",
            row.hydrate_temperature.map(|value| value.value),
        ),
        ("subcooling", row.subcooling.map(|value| value.value)),
    ] {
        match (got, common::record_number(expected, label)) {
            (Some(got), Some(want)) => {
                common::assert_close(got, want, tol, &format!("{at}::{label}"))
            }
            (None, None) => {}
            (got, want) => panic!("{at}::{label}: {got:?} against {want:?}"),
        }
    }
    let expected_risk = expected
        .iter()
        .find(|field| field.name == "risk_level")
        .and_then(|field| match field.value {
            FieldValue::Name(name) => Some(name),
            FieldValue::Number(_) | FieldValue::Flag(_) => None,
        })
        .unwrap_or_else(|| panic!("{at}::risk_level is absent"));
    assert_eq!(row.risk_level, expected_risk, "{at}::risk_level");
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        for (field, value) in [
            (
                "critical_point_count",
                f64::from(result.critical_point_count),
            ),
            ("unknown_point_count", f64::from(result.unknown_point_count)),
            ("total_points", f64::from(result.total_points)),
            (
                "minimum_subcooling",
                result
                    .minimum_subcooling
                    .map_or(f64::NAN, |value| value.value),
            ),
        ] {
            if let Some(want) = case.expected_value(field) {
                common::assert_close(value, want, case.tolerance, &format!("{context} ({field})"));
            }
        }
        if let Some(expected) = case.expected_string("overall_risk") {
            assert_eq!(
                result.overall_risk.as_str(),
                expected,
                "{context}: the overall risk"
            );
        }
        let expected = common::expected_record(case, "points");
        assert_eq!(
            result.points.len(),
            expected.len(),
            "{context}: the point count"
        );
        for (index, (row, declared)) in result.points.iter().zip(expected).enumerate() {
            assert_row(
                declared,
                row,
                &format!("{context}::points[{index}]"),
                case.tolerance,
            );
        }
        // The flag is self-asserted: `TestCase` carries no boolean expectations.
        let evaluated = !matches!(case.id, "dry" | "water_zero");
        assert_eq!(
            result.all_points_evaluated, evaluated,
            "{context}: all points evaluated"
        );
        common::assert_consistent(&result, context);
    }
}

/// **Three vectors of different lengths are refused** rather than read to the shortest.
#[test]
fn mismatched_profiles_are_refused() {
    let components: Vec<String> = ["methane", "ethane", "propane", "water"]
        .iter()
        .map(|name| (*name).to_owned())
        .collect();
    let error = hydrate_risk_map(
        &components,
        &[0.82, 0.06, 0.03, 0.09],
        &[0.0, 1000.0],
        &[1.0e7],
        &[300.0],
        Some("srk"),
        Some("pvtsim"),
        None,
        None,
    )
    .expect_err("the profile lengths disagree");
    assert!(error.to_string().contains("distances"));
}
