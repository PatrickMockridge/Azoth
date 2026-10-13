//! Spec-driven tests for the `eos.water_dew_point_multiphase` model.

// The oracle values are NeqSim's ground truth, so their full digits are the point.
#![allow(clippy::excessive_precision)]

use azoth_core::units::pascals;
use azoth_eos::Cubic;
use azoth_eos::model_gen;
use azoth_eos::results::WaterDewPointMultiphaseResult;
use azoth_eos::water_dew_point_multiphase::water_dew_point_multiphase;
use azoth_eos::water_dew_point_temperature::water_dew_point_temperature;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.water_dew_point_multiphase";

fn names(case: &azoth_core::spec::TestCase) -> Vec<String> {
    case.list("components")
        .expect("the case declares components")
        .iter()
        .map(|name| (*name).to_string())
        .collect()
}

fn solve(case: &azoth_core::spec::TestCase) -> WaterDewPointMultiphaseResult {
    water_dew_point_multiphase(
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
        let context = &format!("{}::{}", spec.id, case.id);
        common::assert_close(
            result.temperature.value,
            common::expected(case, "temperature"),
            case.tolerance,
            &format!("{context} (T)"),
        );
        // The bracket has to have closed, and the fluid at the answer has to be the two
        // phases the search was looking for - a dew point with one phase is a temperature
        // nothing appeared at, whatever the number says.
        assert!(
            result.bracket_width <= 1.0e-4,
            "{context}: the bracket closed to {} K, over the 1e-4 it bisects to",
            result.bracket_width
        );
        assert_eq!(result.phase_count, 2, "{context}: phases at the answer");
    }
}

/// **The two routes are two ids because they disagree, and the disagreement has a sign.**
///
/// `eos.water_dew_point_temperature` prescribes its incipient phase and walks to the fugacity
/// equality; this one asks a multiphase flash where the phase appears. Upstream keeps them as
/// separate entry points and `WaterDewPointAnalyser` selects between them by name, so a port
/// that collapsed them - or that made one call the other - would have to pick one of two
/// answers and lose the other.
///
/// **They cross, and the crossing is measured here rather than assumed.** At 22 and 35 ppm the
/// prescribed route answers the *warmer* temperature, by `1.68` and `0.51` K; at 200 ppm it
/// answers the *colder*, by `0.008` K. So "which route is conservative" is a function of the
/// water content and not a property of either - which is why this test asserts a difference
/// and states the signs rather than pinning one of them.
#[test]
fn the_prescribed_route_and_this_one_are_different_answers() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");

    for case in spec.cases {
        let names = names(case);
        let z = case.vector("z").expect("the case declares z");
        let p = pascals(common::input(case, "P"));
        let bracketed = solve(case);
        let prescribed = water_dew_point_temperature(&names, z, p, Cubic::Srk)
            .expect("the prescribed route solves the same feed");

        let difference = prescribed.temperature.value - bracketed.temperature.value;
        assert!(
            difference.abs() > 1.0e-3,
            "{}: the two routes agree to {:e} K, so one of them is not carrying its own rule",
            case.id,
            difference
        );
    }
}

/// **A fluid with no water is refused rather than bracketed.**
///
/// Upstream's coarse search would step to its `200` K floor and return that, which reads as a
/// dew point and is not one - so the refusal is the port's and the spec states it.
#[test]
fn a_feed_with_no_water_is_refused() {
    let names: Vec<String> = ["methane", "ethane"].iter().map(|s| (*s).to_string()).collect();
    let refused = water_dew_point_multiphase(&names, &[0.9, 0.1], pascals(7.0e6), Cubic::Srk);
    assert!(refused.is_err(), "a feed with no water has no aqueous phase");

    let names: Vec<String> = ["methane", "water"].iter().map(|s| (*s).to_string()).collect();
    let refused = water_dew_point_multiphase(&names, &[1.0, 0.0], pascals(7.0e6), Cubic::Srk);
    assert!(
        refused.is_err(),
        "water at a zero mole fraction is not a phase the bracket can find"
    );
}
