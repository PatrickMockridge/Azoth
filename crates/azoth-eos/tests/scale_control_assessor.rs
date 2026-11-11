//! Spec-driven tests for `eos.scale_control_assessor`.

use azoth_core::spec::{RecordField, TestCase};
use azoth_core::units::{kelvins, kilograms_per_cubic_meter, pascals};
use azoth_eos::model_gen;
use azoth_eos::results::ScaleControlAssessorResult;
use azoth_eos::scale_control_assessor::{ScaleInhibitor, scale_control_assessor};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.scale_control_assessor";

fn inhibitor(row: &[RecordField]) -> ScaleInhibitor {
    let number = |name: &str| common::record_number(row, name).expect("every field is stated");
    ScaleInhibitor {
        scale_type: common::record_name(row, "scale_type")
            .expect("every field is stated")
            .to_owned(),
        inhibitor_chemistry: common::record_name(row, "inhibitor_chemistry")
            .expect("every field is stated")
            .to_owned(),
        temperature: kelvins(number("temperature")),
        saturation_ratio: number("saturation_ratio"),
        total_dissolved_solids: kilograms_per_cubic_meter(number("total_dissolved_solids")),
        calcium: kilograms_per_cubic_meter(number("calcium")),
        available_dose: kilograms_per_cubic_meter(number("available_dose")),
    }
}

fn call(case: &TestCase) -> ScaleControlAssessorResult {
    let mass = |name: &str| kilograms_per_cubic_meter(common::input(case, name));
    let inhibitors: Vec<ScaleInhibitor> = common::record_input(case, "inhibitors")
        .iter()
        .map(|row| inhibitor(row))
        .collect();
    scale_control_assessor(
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
        &inhibitors,
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

/// An expected field the case omits is not pinned; the absence rule is held by
/// `absent_scales_and_uninhibited_scales`, not by every row that leaves an index out.
fn assert_index(expected: &[RecordField], row: &Option<f64>, at: &str, tol: f64) {
    if let Some(want) = common::record_number(expected, "saturation_index") {
        common::assert_close(
            row.unwrap_or_else(|| panic!("{at}: the case pins an index")),
            want,
            tol,
            at,
        );
    }
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);

        let thermo = case.record("thermodynamic_si").unwrap_or(&[]);
        if !thermo.is_empty() {
            assert_eq!(
                result.thermodynamic_si.len(),
                thermo.len(),
                "{context}: thermo rows"
            );
        }
        for (index, (row, declared)) in result.thermodynamic_si.iter().zip(thermo).enumerate() {
            if let Some(scale) = common::record_name(declared, "scale") {
                assert_eq!(row.scale, scale, "{context}: thermo[{index}] scale");
            }
            assert_index(
                declared,
                &row.saturation_index,
                &format!("{context}: thermo[{index}] index"),
                case.tolerance,
            );
        }

        let risk = case.record("kinetic_risk_index").unwrap_or(&[]);
        if !risk.is_empty() {
            assert_eq!(
                result.kinetic_risk_index.len(),
                risk.len(),
                "{context}: risk rows"
            );
        }
        for (index, (row, declared)) in result.kinetic_risk_index.iter().zip(risk).enumerate() {
            if let Some(scale) = common::record_name(declared, "scale") {
                assert_eq!(row.scale, scale, "{context}: risk[{index}] scale");
            }
            if let Some(want) = common::record_number(declared, "risk_index") {
                common::assert_close(
                    row.risk_index
                        .unwrap_or_else(|| panic!("{context}: risk[{index}] is pinned")),
                    want,
                    case.tolerance,
                    &format!("{context}: risk[{index}] index"),
                );
            }
        }

        if let Some(expected) = case.expected_value("worst_kinetic_risk_index") {
            common::assert_close(
                result
                    .worst_kinetic_risk_index
                    .expect("the case pins a worst"),
                expected,
                case.tolerance,
                &format!("{context}: the worst index"),
            );
        }
        assert_eq!(
            result.kinetically_controlled,
            result
                .worst_kinetic_risk_index
                .is_some_and(|value| value < 0.5),
            "{context}: the verdict follows the 0.5 threshold"
        );
        common::assert_consistent(&result, context);
    }
}

/// **An absent ion stays absent in both maps**, and a scale with no inhibitor carries its
/// index unchanged.
#[test]
fn absent_scales_and_uninhibited_scales() {
    let result = scale_control_assessor(
        kilograms_per_cubic_meter(2.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(20.0),
        kilograms_per_cubic_meter(0.8),
        kilograms_per_cubic_meter(0.5),
        kilograms_per_cubic_meter(80.0),
        kelvins(353.15),
        pascals(1.0e7),
        pascals(2.0e5),
        7.0,
        &[],
    )
    .expect("the state runs");
    assert!(result.thermodynamic_si[1].saturation_index.is_none());
    assert!(result.kinetic_risk_index[1].risk_index.is_none());
    assert_eq!(
        result.thermodynamic_si[0].saturation_index,
        result.kinetic_risk_index[0].risk_index
    );
}

/// **A duplicate inhibitor row keeps the last**, as the class's map put does.
#[test]
fn the_last_inhibitor_row_wins() {
    let row = |dose: f64| ScaleInhibitor {
        scale_type: "caco3".to_owned(),
        inhibitor_chemistry: "phosphonate".to_owned(),
        temperature: kelvins(353.15),
        saturation_ratio: 1.6851624454329173,
        total_dissolved_solids: kilograms_per_cubic_meter(80.0),
        calcium: kilograms_per_cubic_meter(2.0),
        available_dose: kilograms_per_cubic_meter(dose),
    };
    let run = |inhibitors: &[ScaleInhibitor]| {
        scale_control_assessor(
            kilograms_per_cubic_meter(2.0),
            kilograms_per_cubic_meter(0.05),
            kilograms_per_cubic_meter(0.0),
            kilograms_per_cubic_meter(0.0),
            kilograms_per_cubic_meter(0.0),
            kilograms_per_cubic_meter(20.0),
            kilograms_per_cubic_meter(0.8),
            kilograms_per_cubic_meter(0.5),
            kilograms_per_cubic_meter(80.0),
            kelvins(353.15),
            pascals(1.0e7),
            pascals(2.0e5),
            7.0,
            inhibitors,
        )
        .expect("the state runs")
    };
    let zero_then_dosed = run(&[row(0.0), row(0.02)]);
    let dosed = run(&[row(0.02)]);
    assert_eq!(
        zero_then_dosed.kinetic_risk_index[0].risk_index,
        dosed.kinetic_risk_index[0].risk_index
    );
}

/// **The 0.999 efficiency cap is what bounds the residual**: a dose past the MIC cannot
/// push the risk index more than three decades below the thermodynamic index.
#[test]
fn the_efficiency_cap_bounds_the_residual() {
    let huge = ScaleInhibitor {
        scale_type: "caco3".to_owned(),
        inhibitor_chemistry: "phosphonate".to_owned(),
        temperature: kelvins(353.15),
        saturation_ratio: 1.6851624454329173,
        total_dissolved_solids: kilograms_per_cubic_meter(80.0),
        calcium: kilograms_per_cubic_meter(2.0),
        available_dose: kilograms_per_cubic_meter(1000.0),
    };
    let result = scale_control_assessor(
        kilograms_per_cubic_meter(2.0),
        kilograms_per_cubic_meter(0.05),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(20.0),
        kilograms_per_cubic_meter(0.8),
        kilograms_per_cubic_meter(0.5),
        kilograms_per_cubic_meter(80.0),
        kelvins(353.15),
        pascals(1.0e7),
        pascals(2.0e5),
        7.0,
        &[huge],
    )
    .expect("the state runs");
    let thermodynamic = result.thermodynamic_si[0]
        .saturation_index
        .expect("calcite has ions");
    let risk = result.kinetic_risk_index[0]
        .risk_index
        .expect("calcite has ions");
    assert!((thermodynamic - risk - 3.0).abs() < 1e-9);
}
