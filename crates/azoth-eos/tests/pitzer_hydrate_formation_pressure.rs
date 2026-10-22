//! Spec-driven tests for `eos.pitzer_hydrate_formation_pressure`.

use azoth_core::AzothError;
use azoth_core::spec::TestCase;
use azoth_core::units::kelvins;
use azoth_eos::model_gen;
use azoth_eos::pitzer_hydrate_formation_pressure::pitzer_hydrate_formation_pressure;
use azoth_eos::results::PitzerHydrateFormationPressureResult;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.pitzer_hydrate_formation_pressure";

fn call(case: &TestCase) -> PitzerHydrateFormationPressureResult {
    let names = case
        .list("components")
        .expect("every case names its components");
    let moles = case.vector("moles").expect("every case states its moles");
    let cubic = common::input_str(case, "cubic")
        .parse()
        .expect("the case's cubic parses");
    pitzer_hydrate_formation_pressure(
        names,
        cubic,
        moles,
        kelvins(common::input(case, "T")),
        case.string("structure"),
        case.input("minimum_pressure"),
        case.input("maximum_pressure"),
        case.string("pitzer_dataset"),
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
        common::assert_close(
            result.pressure.value,
            common::expected(case, "pressure"),
            case.tolerance,
            &format!("{context} (pressure)"),
        );
        common::assert_close(
            result.water_activity,
            common::expected(case, "water_activity"),
            case.tolerance,
            &format!("{context} (water_activity)"),
        );
        assert_eq!(
            Some(result.structure.as_str()),
            case.expected_string("structure"),
            "{context}: the stable structure"
        );
        assert!(
            result.residual.abs() < 1e-7,
            "{context}: the residual at the answer is {}",
            result.residual
        );
        common::assert_consistent(&result, context);
    }
}

/// **A fixed temperature below the guests' Henry floor is refused**, as `validateFeed` refuses
/// it: CO2's fitted range starts at 274.19 K.
#[test]
fn a_temperature_below_the_henry_floor_is_refused() {
    let error = pitzer_hydrate_formation_pressure(
        &["CO2", "water"],
        azoth_eos::Cubic::Srk,
        &[10.0, 55.508435061791985],
        kelvins(273.5),
        None,
        None,
        None,
        None,
    )
    .expect_err("273.5 K is below CO2's 274.19 K floor");
    assert!(matches!(error, AzothError::OutOfRange { .. }));
    assert!(error.to_string().contains("274.19"));
}

/// **A fixed temperature outside the class's 273.15-323.15 K window is refused**, and so is a
/// pressure bound that is not positive.
#[test]
fn the_windows_are_enforced() {
    let hot = pitzer_hydrate_formation_pressure(
        &["CO2", "water"],
        azoth_eos::Cubic::Srk,
        &[10.0, 55.508435061791985],
        kelvins(330.0),
        None,
        None,
        None,
        None,
    )
    .expect_err("330 K is outside the class's window");
    assert!(matches!(hot, AzothError::OutOfRange { .. }));

    let zero_lower = pitzer_hydrate_formation_pressure(
        &["CO2", "water"],
        azoth_eos::Cubic::Srk,
        &[10.0, 55.508435061791985],
        kelvins(278.0),
        None,
        Some(0.0),
        None,
        None,
    )
    .expect_err("zero pressure is not a bound");
    assert!(matches!(zero_lower, AzothError::OutOfRange { .. }));
}
