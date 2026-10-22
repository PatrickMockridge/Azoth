//! Spec-driven tests for `eos.multi_mineral_scale_equilibrium`.

use azoth_core::WarningCode;
use azoth_core::spec::{FieldValue, RecordField, TestCase};
use azoth_core::units::{cubic_meters, kelvins, kilograms_per_cubic_meter, pascals};
use azoth_eos::model_gen;
use azoth_eos::multi_mineral_scale_equilibrium::multi_mineral_scale_equilibrium;
use azoth_eos::results::{
    ActivityCoefficientRow, MineralEquilibriumRow, MultiMineralScaleEquilibriumResult,
};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.multi_mineral_scale_equilibrium";

fn call(case: &TestCase) -> MultiMineralScaleEquilibriumResult {
    let input = |name: &str| common::input(case, name);
    let mass = |name: &str| kilograms_per_cubic_meter(input(name));
    multi_mineral_scale_equilibrium(
        mass("calcium"),
        mass("barium"),
        mass("strontium"),
        mass("iron"),
        mass("magnesium"),
        mass("sodium"),
        mass("bicarbonate"),
        mass("sulphate"),
        mass("total_dissolved_solids"),
        kelvins(input("temperature")),
        pascals(input("pressure")),
        pascals(input("co2_partial_pressure")),
        input("ph"),
        case.flag("auto_ph").unwrap_or(false),
        case.flag("second_order_pressure_correction")
            .unwrap_or(false),
        common::input_str(case, "activity_model"),
        case.input("bdot_ion_size"),
        case.input("bdot_b"),
        case.input("pitzer_ionic_strength_molal"),
        case.input("water_volume")
            .map(|litres| cubic_meters(litres * 1.0e-3)),
        case.input("max_iterations").map(|value| value as usize),
        case.input("relaxation"),
        case.input("convergence_threshold")
            .map(azoth_core::units::moles_per_cubic_meter),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

fn record_flag(row: &[RecordField], name: &str) -> Option<bool> {
    row.iter()
        .find(|field| field.name == name)
        .and_then(|field| match field.value {
            FieldValue::Flag(value) => Some(value),
            FieldValue::Number(_) | FieldValue::Name(_) => None,
        })
}

fn assert_mineral_row(expected: &[RecordField], row: &MineralEquilibriumRow, at: &str, tol: f64) {
    assert_eq!(
        row.name,
        common::record_name(expected, "name").expect("every mineral row is named"),
        "{at}: the mineral"
    );
    assert_eq!(
        row.common_name,
        common::record_name(expected, "common_name")
            .expect("every mineral row carries its common name"),
        "{at}: the common name"
    );
    for (label, got) in [
        ("initial_saturation_index", row.initial_saturation_index),
        ("final_saturation_index", row.final_saturation_index),
        (
            "activity_coefficient_product",
            row.activity_coefficient_product,
        ),
        ("precipitated_amount", row.precipitated_amount.value),
        ("precipitated_mass", row.precipitated_mass.value),
        ("total_precipitated_mass", row.total_precipitated_mass.value),
    ] {
        common::assert_close(
            got,
            common::record_number(expected, label)
                .unwrap_or_else(|| panic!("{at}: the row states {label}")),
            tol,
            &format!("{at}::{label}"),
        );
    }
    assert_eq!(
        row.limiting_ion_exhausted,
        record_flag(expected, "limiting_ion_exhausted").expect("the row states the flag"),
        "{at}: the limiting flag"
    );
}

fn assert_coefficient_row(
    expected: &[RecordField],
    row: &ActivityCoefficientRow,
    at: &str,
    tol: f64,
) {
    assert_eq!(
        row.ion,
        common::record_name(expected, "ion").expect("every coefficient row names its ion"),
        "{at}: the ion"
    );
    common::assert_close(
        row.activity_coefficient,
        common::record_number(expected, "activity_coefficient").expect("the row states its value"),
        tol,
        &format!("{at}::activity_coefficient"),
    );
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        for field in [
            "total_mass",
            "divalent_activity_coefficient",
            "maximum_ion_balance_residual",
            "residual_calcium",
            "residual_barium",
            "residual_strontium",
            "residual_ferrous",
            "residual_sulphate",
            "residual_carbonate",
        ] {
            if let Some(expected) = case.expected_value(field) {
                let actual = match field {
                    "total_mass" => result.total_mass.value,
                    "divalent_activity_coefficient" => result.divalent_activity_coefficient,
                    "maximum_ion_balance_residual" => result.maximum_ion_balance_residual.value,
                    "residual_calcium" => result.residual_calcium.value,
                    "residual_barium" => result.residual_barium.value,
                    "residual_strontium" => result.residual_strontium.value,
                    "residual_ferrous" => result.residual_ferrous.value,
                    "residual_sulphate" => result.residual_sulphate.value,
                    _ => result.residual_carbonate.value,
                };
                common::assert_close(
                    actual,
                    expected,
                    case.tolerance,
                    &format!("{context} ({field})"),
                );
            }
        }
        if let Some(expected) = case.expected_record("activity_coefficients") {
            assert_eq!(
                result.activity_coefficients.len(),
                expected.len(),
                "{context}: the coefficient count"
            );
            for (index, (row, declared)) in result
                .activity_coefficients
                .iter()
                .zip(expected)
                .enumerate()
            {
                assert_coefficient_row(
                    declared,
                    row,
                    &format!("{context}::activity_coefficients[{index}]"),
                    case.tolerance,
                );
            }
        }
        if let Some(expected) = case.expected_record("minerals") {
            assert_eq!(
                result.minerals.len(),
                expected.len(),
                "{context}: the mineral count"
            );
            for (index, (row, declared)) in result.minerals.iter().zip(expected).enumerate() {
                assert_mineral_row(
                    declared,
                    row,
                    &format!("{context}::minerals[{index}]"),
                    case.tolerance,
                );
            }
        }
        common::assert_consistent(&result, context);
    }
}

/// **A mineral whose free ion is below the class's `1e-18` floor has no row.** On the
/// high-salinity brine the free carbonate is `4.9e-23` mol/L, so calcite is never registered and
/// the assemblage is three rows rather than four; its absence is the registration, not a zero
/// mass.
#[test]
fn a_mineral_below_the_registration_floor_has_no_row() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    let case = spec
        .cases
        .iter()
        .find(|case| case.id == "high_salinity_davies_precipitates_nothing")
        .expect("the high-salinity case");
    let result = call(case);
    assert_eq!(result.minerals.len(), 3, "calcite is not registered");
    assert!(
        result.minerals.iter().all(|row| row.name != "CaCO3"),
        "the calcite row is absent rather than zero"
    );
    assert!(result.maximum_ion_balance_residual.value < 1e-12);
}

/// **A capped solve is returned with its flag and a warning**, not refused: the class reports
/// the flag beside the last iterate, and so does this.
#[test]
fn a_capped_solve_says_so() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    let case = spec
        .cases
        .iter()
        .find(|case| case.id == "iteration_limit_one_step")
        .expect("the iteration-cap case");
    let result = call(case);
    assert_eq!(result.iterations, 1);
    assert!(result.iteration_limit_reached);
    assert!(
        result
            .warnings
            .iter()
            .any(|warning| warning.code == WarningCode::SolverNotConverged),
        "a capped solve carries its warning"
    );
    assert!(
        result.maximum_complementarity_violation > 1.0,
        "one step leaves a large violation"
    );
}

/// **The water volume scales the grams and nothing else.** The per-litre mass and the
/// composition are the same state at one litre and at a thousand.
#[test]
fn the_water_volume_scales_only_the_grams() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    let case = spec
        .cases
        .iter()
        .find(|case| case.id == "water_volume_1000l")
        .expect("the volume case");
    let thousand_litres = call(case);
    let ten_litres = multi_mineral_scale_equilibrium(
        kilograms_per_cubic_meter(common::input(case, "calcium")),
        kilograms_per_cubic_meter(common::input(case, "barium")),
        kilograms_per_cubic_meter(common::input(case, "strontium")),
        kilograms_per_cubic_meter(common::input(case, "iron")),
        kilograms_per_cubic_meter(common::input(case, "magnesium")),
        kilograms_per_cubic_meter(common::input(case, "sodium")),
        kilograms_per_cubic_meter(common::input(case, "bicarbonate")),
        kilograms_per_cubic_meter(common::input(case, "sulphate")),
        kilograms_per_cubic_meter(common::input(case, "total_dissolved_solids")),
        kelvins(common::input(case, "temperature")),
        pascals(common::input(case, "pressure")),
        pascals(common::input(case, "co2_partial_pressure")),
        common::input(case, "ph"),
        true,
        false,
        "davies",
        None,
        None,
        None,
        Some(cubic_meters(0.01)),
        None,
        None,
        None,
    )
    .expect("the state runs");
    let barite_thousand = &thousand_litres.minerals[0];
    let barite_ten = &ten_litres.minerals[0];
    assert!(
        (barite_thousand.precipitated_mass.value - barite_ten.precipitated_mass.value).abs()
            < 1e-12
    );
    common::assert_close(
        barite_ten.total_precipitated_mass.value,
        barite_ten.precipitated_mass.value * 0.01,
        1e-12,
        "ten litres of the per-litre mass",
    );
}

#[test]
fn an_unknown_activity_model_is_refused() {
    let error = multi_mineral_scale_equilibrium(
        kilograms_per_cubic_meter(1.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(0.0),
        kilograms_per_cubic_meter(60.0),
        kelvins(343.15),
        pascals(1.0e7),
        pascals(0.0),
        7.0,
        false,
        false,
        "debye_huckel",
        None,
        None,
        None,
        Some(cubic_meters(1.0e-3)),
        None,
        None,
        None,
    )
    .expect_err("`debye_huckel` is not one of the three");
    assert!(error.to_string().contains("pitzer_binary"));
}
