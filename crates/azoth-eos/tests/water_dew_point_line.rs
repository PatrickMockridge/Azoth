//! Spec-driven tests for the `eos.water_dew_point_line` model.

// The oracle values are NeqSim's ground truth, so their full digits are the point.
#![allow(clippy::excessive_precision)]

use azoth_core::units::pascals;
use azoth_eos::Cubic;
use azoth_eos::model_gen;
use azoth_eos::results::WaterDewPointLineResult;
use azoth_eos::water_dew_point_line::water_dew_point_line;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.water_dew_point_line";

fn names(case: &azoth_core::spec::TestCase) -> Vec<String> {
    case.list("components")
        .expect("the case declares components")
        .iter()
        .map(|name| (*name).to_string())
        .collect()
}

fn solve(case: &azoth_core::spec::TestCase) -> WaterDewPointLineResult {
    water_dew_point_line(
        &names(case),
        case.vector("z").expect("the case declares z"),
        pascals(common::input(case, "minimum_pressure")),
        pascals(common::input(case, "maximum_pressure")),
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
        let context = &format!("{}::{}", spec.id, case.id);

        let want_t = case
            .expected_vector("temperature")
            .expect("the case declares temperature");
        assert_eq!(
            result.temperature.len(),
            want_t.len(),
            "{context}: the line has {} points, not {}",
            result.temperature.len(),
            want_t.len()
        );
        for (i, (&got, &want)) in result.temperature.iter().zip(want_t).enumerate() {
            common::assert_close(
                got,
                want,
                case.tolerance,
                &format!("{context} (temperature[{i}])"),
            );
        }

        let want_p = case
            .expected_vector("pressure")
            .expect("the case declares pressure");
        for (i, (&got, &want)) in result.pressure.iter().zip(want_p).enumerate() {
            common::assert_close(
                got,
                want,
                case.tolerance,
                &format!("{context} (pressure[{i}])"),
            );
        }

        assert_eq!(
            result.failed_points,
            common::expected(case, "failed_points") as u32,
            "{context}: failed points"
        );
    }
}

/// **The line is a curve, and a port that solved one point and multiplied it would fail here.**
///
/// The ten temperatures rise from the cold end to a maximum near the middle and fall back at the
/// top, so the spread is a fact about the sweep and not about the arithmetic at any one point.
#[test]
fn the_line_is_a_curve_and_not_one_point_repeated() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    for case in spec.cases {
        let result = solve(case);
        let low = result
            .temperature
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let high = result
            .temperature
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            high - low > 1.0,
            "{}: the line spans {:e} K, so it is not a curve",
            case.id,
            high - low
        );
    }
}

/// **The cold end is warm by a declared margin, and the other nine are exact.**
///
/// At `1` bara this library answers `224.4312` K where NeqSim answers `224.2818` - a `0.149` K,
/// relative `6.7e-4`, that sits inside the model family's `1e-3` capture tolerance, so the case
/// matches and this test is what stops the divergence hiding inside the tolerance. It is a
/// flash-level difference in a marginal state: at NeqSim's own answer this library's flash
/// reports two phases, gas and aqueous, and it loses the aqueous phase only `0.149` K warmer.
/// The nine remaining points are exact to `1e-9`, and either number moving has to move this.
#[test]
fn the_cold_end_is_warm_by_a_declared_margin_and_the_rest_are_exact() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    for case in spec.cases {
        let result = solve(case);
        let want = case
            .expected_vector("temperature")
            .expect("the case declares temperature");

        let warm = result.temperature[0] - want[0];
        assert!(
            (warm - 0.1494).abs() < 1.0e-3,
            "{}: the cold end is warm by {warm} K against the declared 0.1494",
            case.id
        );
        for (i, (&got, &oracle)) in result.temperature.iter().zip(want).enumerate().skip(1) {
            let relative = (got - oracle).abs() / oracle.abs();
            assert!(
                relative <= 1.0e-9,
                "{}: point {i} is {relative:e} from NeqSim, over the 1e-9 the other nine hold to",
                case.id
            );
        }
    }
}

/// **A sweep that runs backwards is refused rather than swept.**
///
/// The class steps negatively and answers a reversed curve; the port refuses, because a minimum
/// above its maximum is a caller's error and not a line.
#[test]
fn a_sweep_that_runs_backwards_is_refused() {
    let names: Vec<String> = ["methane", "water"]
        .iter()
        .map(|s| (*s).to_string())
        .collect();
    let z = [0.99, 0.01];
    let refused = water_dew_point_line(&names, &z, pascals(2.0e7), pascals(1.0e5), Cubic::Srk);
    assert!(refused.is_err(), "a minimum above its maximum is refused");

    let refused = water_dew_point_line(&names, &z, pascals(0.0), pascals(1.0e5), Cubic::Srk);
    assert!(refused.is_err(), "a zero pressure is not a state");
}
