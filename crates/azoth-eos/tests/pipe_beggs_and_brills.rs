//! Spec-driven tests for `eos.pipe_beggs_and_brills`.

use azoth_core::units::{kelvins, meters, pascals};
use azoth_eos::model_gen;
use azoth_eos::pipe_beggs_and_brills::pipe_beggs_and_brills;
use azoth_eos::results::PipeBeggsAndBrillsResult;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.pipe_beggs_and_brills";

fn call(case: &azoth_core::spec::TestCase) -> PipeBeggsAndBrillsResult {
    let components: Vec<String> = common::list_input(case, "components")
        .iter()
        .map(|name| (*name).to_owned())
        .collect();
    pipe_beggs_and_brills(
        &components,
        common::input(case, "inlet_n"),
        case.vector("inlet_z")
            .expect("the case states its composition"),
        pascals(common::input(case, "inlet_p")),
        kelvins(common::input(case, "inlet_t")),
        meters(common::input(case, "length")),
        meters(common::input(case, "diameter")),
        meters(common::input(case, "elevation")),
        common::input(case, "number_of_increments") as usize,
        meters(common::input(case, "roughness")),
        Some(common::input_str(case, "eos")),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

fn assert_vector(case: &azoth_core::spec::TestCase, context: &str, name: &str, got: &[f64]) {
    if let Some(expected) = case.expected_vector(name) {
        assert_eq!(got.len(), expected.len(), "{context}: {name} length");
        for (index, (value, want)) in got.iter().zip(expected).enumerate() {
            common::assert_close(
                *value,
                *want,
                case.tolerance,
                &format!("{context}: {name}[{index}]"),
            );
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
        assert_vector(case, context, "pressure_profile", &result.pressure_profile);
        assert_vector(
            case,
            context,
            "pressure_drop_profile",
            &result.pressure_drop_profile,
        );
        assert_vector(
            case,
            context,
            "temperature_profile",
            &result.temperature_profile,
        );
        assert_vector(
            case,
            context,
            "liquid_holdup_profile",
            &result.liquid_holdup_profile,
        );
        assert_vector(
            case,
            context,
            "mixture_superficial_velocity_profile",
            &result.mixture_superficial_velocity_profile,
        );
        for (field, value) in [
            ("outlet_pressure", result.outlet_pressure.value),
            ("total_pressure_drop", result.total_pressure_drop.value),
        ] {
            if let Some(want) = case.expected_value(field) {
                common::assert_close(value, want, case.tolerance, &format!("{context} ({field})"));
            }
        }
        common::assert_consistent(&result, context);
    }
}

fn run(names: &[&str], z: &[f64], flow: f64) -> PipeBeggsAndBrillsResult {
    let components: Vec<String> = names.iter().map(|name| (*name).to_owned()).collect();
    pipe_beggs_and_brills(
        &components,
        flow,
        z,
        pascals(1.0e7),
        kelvins(333.15),
        meters(1000.0),
        meters(0.1),
        meters(0.0),
        10,
        meters(1.0e-5),
        Some("srk"),
    )
    .expect("the state runs")
}

/// **The three liquid-rich branches the wet-gas line never reaches** are reached by the
/// flows the capture uses, and each reports one regime over its whole length.
#[test]
fn the_branches_are_reached() {
    let segregated = run(&["methane", "CO2", "water"], &[0.85, 0.05, 0.10], 80.0);
    assert!(
        segregated
            .flow_regime_profile
            .iter()
            .all(|r| r == "segregated")
    );

    let intermittent = run(&["methane", "CO2", "water"], &[0.2, 0.05, 0.75], 320.0);
    assert!(
        intermittent
            .flow_regime_profile
            .iter()
            .all(|r| r == "intermittent")
    );

    let distributed = run(&["methane", "CO2", "water"], &[0.5, 0.05, 0.45], 800.0);
    assert!(
        distributed
            .flow_regime_profile
            .iter()
            .all(|r| r == "distributed")
    );
}

/// **A line that cannot deliver refuses**, as the class's negative outlet pressure does.
#[test]
fn a_line_that_cannot_deliver_is_refused() {
    let components: Vec<String> = ["methane", "CO2", "water"]
        .iter()
        .map(|n| (*n).to_owned())
        .collect();
    let error = pipe_beggs_and_brills(
        &components,
        1600.0,
        &[0.85, 0.05, 0.10],
        pascals(1.0e7),
        kelvins(333.15),
        meters(1000.0),
        meters(0.1),
        meters(0.0),
        10,
        meters(1.0e-5),
        Some("srk"),
    )
    .expect_err("the line does not deliver at this flow");
    assert!(error.to_string().contains("outlet pressure"));
}
