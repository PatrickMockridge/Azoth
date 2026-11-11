//! Spec-driven tests for `eos.closed_loop_deposition_solver`.

use azoth_core::units::{
    kelvins, kilograms_per_cubic_meter, kilograms_per_square_meter_second, meters, pascals, seconds,
};
use azoth_eos::closed_loop_deposition_solver::closed_loop_deposition_solver;
use azoth_eos::model_gen;
use azoth_eos::results::ClosedLoopDepositionSolverResult;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.closed_loop_deposition_solver";

fn call(case: &azoth_core::spec::TestCase) -> ClosedLoopDepositionSolverResult {
    let components: Vec<String> = common::list_input(case, "components")
        .iter()
        .map(|name| (*name).to_owned())
        .collect();
    let density = |name: &str| kilograms_per_cubic_meter(common::input(case, name));
    closed_loop_deposition_solver(
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
        density("calcium"),
        density("bicarbonate"),
        density("sulphate"),
        density("barium"),
        density("sodium"),
        density("total_dissolved_solids"),
        common::input(case, "ph"),
        seconds(common::input(case, "service_time")),
        density("scale_density"),
        kilograms_per_square_meter_second(common::input(case, "base_rate")),
        common::input(case, "inhibitor_efficiency"),
        meters(common::input(case, "tolerance")),
        common::input(case, "max_iterations") as usize,
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
        assert_vector(case, context, "diameter_history", &result.diameter_history);
        assert_vector(
            case,
            context,
            "max_thickness_history",
            &result.max_thickness_history,
        );
        assert_vector(case, context, "velocity_history", &result.velocity_history);
        for (field, value) in [
            ("iterations_taken", result.iterations_taken as f64),
            (
                "final_effective_diameter",
                result.final_effective_diameter.value,
            ),
        ] {
            if let Some(want) = case.expected_value(field) {
                common::assert_close(value, want, case.tolerance, &format!("{context} ({field})"));
            }
        }
        // The boolean is self-asserted: `TestCase` carries no flag expectations, so the one
        // case that stops at its cap is named here.
        if case.id == "max_iterations_2" {
            assert!(
                !result.converged,
                "{context}: the cap was reached unconverged"
            );
        } else {
            assert!(result.converged, "{context}: the loop converged");
        }
        common::assert_consistent(&result, context);
    }
}

/// **A pipe that refuses stops the solve**: the class swallows the failure, and this model
/// propagates it rather than evaluating a partial profile.
#[test]
fn a_pipe_that_refuses_stops_the_solve() {
    let components: Vec<String> = ["methane", "CO2", "water"]
        .iter()
        .map(|name| (*name).to_owned())
        .collect();
    let density = kilograms_per_cubic_meter;
    let error = closed_loop_deposition_solver(
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
        density(1.0),
        density(0.5),
        density(0.1),
        density(0.01),
        density(20.0),
        density(50.0),
        6.5,
        seconds(31_557_600.0),
        density(2700.0),
        kilograms_per_square_meter_second(1.5844043908606944e-8),
        0.0,
        meters(1.0e-4),
        10,
        Some("srk"),
    )
    .expect_err("the line does not deliver at this flow");
    assert!(error.to_string().contains("outlet pressure"));
}
