//! Spec-driven tests for `eos.pitzer_hydrate_formation_temperature`.

use azoth_core::AzothError;
use azoth_core::spec::TestCase;
use azoth_core::units::pascals;
use azoth_eos::model_gen;
use azoth_eos::pitzer_hydrate_formation_temperature::pitzer_hydrate_formation_temperature;
use azoth_eos::results::PitzerHydrateFormationTemperatureResult;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.pitzer_hydrate_formation_temperature";

fn call(case: &TestCase) -> PitzerHydrateFormationTemperatureResult {
    let names = case
        .list("components")
        .expect("every case names its components");
    let moles = case.vector("moles").expect("every case states its moles");
    let cubic = common::input_str(case, "cubic")
        .parse()
        .expect("the case's cubic parses");
    pitzer_hydrate_formation_temperature(
        names,
        cubic,
        moles,
        pascals(common::input(case, "P")),
        case.string("structure"),
        case.input("minimum_temperature"),
        case.input("maximum_temperature"),
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
            result.temperature.value,
            common::expected(case, "temperature"),
            case.tolerance,
            &format!("{context} (temperature)"),
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

/// **A brine beyond the Henry-limited window has no bracket**, and the model says so with the
/// solver's own error rather than a clamped temperature: 20 wt% NaCl at 30 bara is the capture's
/// `co2_nacl_20wt_30bara` refusal.
#[test]
fn a_salt_beyond_the_henry_window_is_solver_not_converged() {
    let error = pitzer_hydrate_formation_temperature(
        &["CO2", "water", "Na+", "K+", "Ca++", "Mg++", "Cl-"],
        azoth_eos::Cubic::Srk,
        &[10.0, 55.508435061791985, 4.0, 0.0, 0.0, 0.0, 4.0],
        pascals(3.0e6),
        None,
        None,
        None,
        None,
    )
    .expect_err("20 wt% NaCl has no bracket in 274.19-323.15 K");
    assert!(matches!(error, AzothError::SolverNotConverged { .. }));
}

/// **A feed without a guest is refused by the resolver**, and one without water too: the class's
/// `validateFeed` says the same in its own words.
#[test]
fn a_feed_without_a_guest_or_without_water_is_refused() {
    let no_guest = pitzer_hydrate_formation_temperature(
        &["water"],
        azoth_eos::Cubic::Srk,
        &[55.508435061791985],
        pascals(3.0e6),
        None,
        None,
        None,
        None,
    )
    .expect_err("water alone has no guest");
    assert!(no_guest.to_string().contains("hydrate former"));

    let no_water = pitzer_hydrate_formation_temperature(
        &["CO2"],
        azoth_eos::Cubic::Srk,
        &[10.0],
        pascals(3.0e6),
        None,
        None,
        None,
        None,
    )
    .expect_err("CO2 alone has no water");
    assert!(no_water.to_string().contains("no water"));
}

/// **Pinning a structure evaluates it alone**, so the pinned answer differs from the stable one
/// where the comparison would have chosen the other structure.
#[test]
fn pinning_a_structure_changes_the_answer() {
    let run = |structure: Option<&str>| {
        pitzer_hydrate_formation_temperature(
            &["CO2", "water"],
            azoth_eos::Cubic::Srk,
            &[10.0, 55.508435061791985],
            pascals(3.0e6),
            structure,
            None,
            None,
            None,
        )
        .expect("the state runs")
    };
    let auto = run(None);
    let pinned = run(Some("sii"));
    assert_eq!(auto.structure.as_str(), "structure_i");
    assert_eq!(pinned.structure.as_str(), "structure_ii");
    assert!(pinned.temperature.value < auto.temperature.value);
}

/// **The temperature and pressure solves are the same boundary**: solving T at 30 bara and then
/// P at that temperature returns 30 bara, which is the upstream test's round trip.
#[test]
fn the_two_solves_are_one_boundary() {
    use azoth_eos::pitzer_hydrate_formation_pressure::pitzer_hydrate_formation_pressure;

    let names = ["CO2", "water", "Na+", "K+", "Ca++", "Mg++", "Cl-"];
    let moles = [
        10.0,
        55.508435061791985,
        0.9006088115566122,
        0.0,
        0.0,
        0.0,
        0.9006088115566122,
    ];
    let temperature = pitzer_hydrate_formation_temperature(
        &names,
        azoth_eos::Cubic::Srk,
        &moles,
        pascals(3.0e6),
        None,
        None,
        None,
        None,
    )
    .expect("the temperature solves");
    let pressure = pitzer_hydrate_formation_pressure(
        &names,
        azoth_eos::Cubic::Srk,
        &moles,
        azoth_core::units::kelvins(temperature.temperature.value),
        None,
        None,
        None,
        None,
    )
    .expect("the pressure solves back");
    common::assert_close(
        pressure.pressure.value,
        3.0e6,
        3.0e-3,
        "the round trip back to 30 bara",
    );
}
