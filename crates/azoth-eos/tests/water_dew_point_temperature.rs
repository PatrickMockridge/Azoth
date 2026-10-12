//! Spec-driven tests for the `eos.water_dew_point_temperature` model.

// The oracle values are NeqSim's ground truth, so their full digits are the point.
#![allow(clippy::excessive_precision)]

use azoth_core::units::{kelvins, pascals};
use azoth_eos::mixture::RootSide;
use azoth_eos::water_dew_point_temperature::water_dew_point_temperature;
use azoth_eos::{Cubic, model_gen};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.water_dew_point_temperature";

fn names(case: &azoth_core::spec::TestCase) -> Vec<String> {
    case.list("components")
        .expect("the case declares components")
        .iter()
        .map(|name| (*name).to_string())
        .collect()
}

fn solve(case: &azoth_core::spec::TestCase) -> azoth_eos::results::WaterDewPointTemperatureResult {
    water_dew_point_temperature(
        &names(case),
        case.vector("z").expect("the case declares z"),
        pascals(common::input(case, "P")),
        Cubic::Srk,
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty());

    for case in spec.cases {
        let result = solve(case);
        common::assert_close(
            result.temperature.value,
            common::expected(case, "temperature"),
            case.tolerance,
            &format!("{}::{} (T)", spec.id, case.id),
        );
    }
}

/// **The reported residual is the equation, re-derived from outside.**
///
/// This is the property that makes the case values mean something. The fixed point's gain is
/// a fixed `100 K` per unit residual, so it approaches the root slowly and stops on its step
/// cap with a residual of about `1.5e-6` rather than on the `1e-7` it tests for - which is
/// why the cases pin the *temperature* and not the residual. What they cannot catch is a port
/// that returned a temperature from somewhere else; this can, because it evaluates
/// `z_w - phi_aq / phi_vap` with the crate's own phase-state call and compares it against what
/// the model says it was.
///
/// **The evaluation is one step behind the answer, and that is the class's own arrangement**:
/// `run()` computes the residual, steps `T += 100 r`, and exits, so the temperature a caller
/// reads is *past* the state the residual belongs to. Re-deriving at the reported temperature
/// instead is a relative `2.9e-4` out, which is what this test found before it was written
/// this way.
#[test]
fn the_reported_residual_is_the_equation_one_step_behind_the_answer() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");

    for case in spec.cases {
        let names = names(case);
        let z = case.vector("z").expect("the case declares z");
        let p = pascals(common::input(case, "P"));
        let result = solve(case);

        let borrowed: Vec<&str> = names.iter().map(String::as_str).collect();
        let (mixture, _) = azoth_eos::databank::mixture_of(&borrowed, Cubic::Srk, None)
            .expect("the case's fluid should resolve");
        let index = borrowed
            .iter()
            .position(|name| *name == "water")
            .expect("every case names water");

        let evaluated_at = kelvins(result.temperature.value - 100.0 * result.residual);
        let reduced = mixture
            .reduced_parameters(evaluated_at, p)
            .expect("the state the residual belongs to");
        let vapour = mixture
            .phase_state(&reduced, z, RootSide::Vapour)
            .expect("the feed's vapour root");
        let mut x = vec![1.0e-30_f64; names.len()];
        x[index] = 1.0;
        let aqueous = mixture
            .phase_state(&reduced, &x, RootSide::Liquid)
            .expect("the incipient phase's liquid root");

        let re_derived = z[index] - (aqueous.ln_phi[index] - vapour.ln_phi[index]).exp();
        common::assert_close(
            result.residual,
            re_derived,
            1e-12,
            &format!("{}::{} (residual against the equation)", spec.id, case.id),
        );
    }
}

/// **The answer is not the seed, and the walk says it is not.**
///
/// The fixed point starts at the water component's `melting_point_temperature` (`273.15` K)
/// and returns a colder temperature; a port that returned its seed, or that stopped on the
/// first step, would pass a bare `temperature` comparison against nothing but fail this.
#[test]
fn the_answer_is_reached_from_the_seed_rather_than_being_it() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");

    for case in spec.cases {
        let result = solve(case);
        assert!(
            result.temperature.value < 273.15 - 5.0,
            "case `{}`: the answer {} K is not a walk away from the 273.15 K seed",
            case.id,
            result.temperature.value
        );
        // The gain converges too slowly for `1e-7` inside the class's own `10000`-step cap,
        // so the walk runs out of steps - which is NeqSim's behaviour too and not an error,
        // because the temperature has settled to the oracle's digits regardless.
        assert_eq!(result.iterations, 10000, "case `{}`: cap", case.id);
        assert!(!result.converged, "case `{}`: converged", case.id);
    }
}

/// **Colder for less water, and colder for less pressure** - the two directions the cases
/// differ along, so a port that swapped them would fail rather than pass three numbers.
#[test]
fn the_dew_point_moves_the_way_the_state_does() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    let by_id = |id: &str| {
        spec.cases
            .iter()
            .find(|case| case.id == id)
            .unwrap_or_else(|| panic!("case `{id}` should be in the spec"))
    };

    let wet = solve(by_id("neqsims_own_test_gas_at_35_ppm_and_70_bara"));
    let dry = solve(by_id("the_same_gas_at_22_ppm_and_70_bara"));
    let low = solve(by_id("the_same_gas_at_35_ppm_and_20_bara"));

    assert!(
        dry.temperature.value < wet.temperature.value,
        "less water should condense later: {} K at 22 ppm against {} K at 35 ppm",
        dry.temperature.value,
        wet.temperature.value
    );
    assert!(
        low.temperature.value < wet.temperature.value,
        "less pressure should condense later: {} K at 20 bara against {} K at 70 bara",
        low.temperature.value,
        wet.temperature.value
    );
}

/// **A fluid with neither `water` nor `MEG` is refused**, where NeqSim writes a zero
/// temperature; and a zero-mole target is refused rather than walked.
#[test]
fn the_two_fluids_the_class_walks_are_the_two_this_refuses() {
    let names = ["methane".to_string(), "ethane".to_string()];
    let z = [0.9, 0.1];
    let refused = water_dew_point_temperature(&names, &z, pascals(7.0e6), Cubic::Srk);
    assert!(refused.is_err(), "a fluid with no water should be refused");

    let names = ["methane".to_string(), "water".to_string()];
    let z = [1.0, 0.0];
    let refused = water_dew_point_temperature(&names, &z, pascals(7.0e6), Cubic::Srk);
    assert!(
        refused.is_err(),
        "a zero-mole water is the phase this seeds from and should be refused"
    );
}
