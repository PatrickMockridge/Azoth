//! Spec-driven tests for the `eos.multiflash_stability_test` model.

// The oracle values are NeqSim's ground truth, so their full digits are the point.
#![allow(clippy::excessive_precision)]

use azoth_core::units::{kelvins, pascals};
use azoth_eos::Cubic;
use azoth_eos::databank::mixture_of;
use azoth_eos::model_gen;
use azoth_eos::multiflash_stability_test::multiflash_stability_test;
use azoth_eos::results::{MultiflashSeedPath, MultiflashStabilityTestResult};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.multiflash_stability_test";

fn names(case: &azoth_core::spec::TestCase) -> Vec<String> {
    case.list("components")
        .expect("the case declares components")
        .iter()
        .map(|name| (*name).to_string())
        .collect()
}

fn solve(case: &azoth_core::spec::TestCase) -> MultiflashStabilityTestResult {
    let names = names(case);
    let borrowed: Vec<&str> = names.iter().map(String::as_str).collect();
    let (mixture, _) =
        mixture_of(&borrowed, Cubic::Srk, None).expect("the case's fluid should resolve");
    multiflash_stability_test(
        &mixture,
        kelvins(common::input(case, "T")),
        pascals(common::input(case, "P")),
        case.vector("z").expect("the case declares z"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty());

    for case in spec.cases {
        let result = solve(case);
        let context = format!("{}::{}", spec.id, case.id);

        let expected = case
            .expected_string("path")
            .unwrap_or_else(|| panic!("{context}: the case should declare a path"));
        let actual = match result.path {
            MultiflashSeedPath::Stable => "stable",
            MultiflashSeedPath::WilsonLiquid => "wilson_liquid",
            MultiflashSeedPath::WilsonVapour => "wilson_vapour",
            MultiflashSeedPath::PureComponent => "pure_component",
        };
        assert_eq!(actual, expected, "{context}: path");

        if let Some(index) = case.expected_value("seed_component_index") {
            common::assert_close(
                result.seed_component_index as f64,
                index,
                case.tolerance,
                &format!("{context} (seed_component_index)"),
            );
        }
        if let Some(fraction) = case.expected_value("seed_fraction") {
            common::assert_close(result.seed_fraction, fraction, case.tolerance, &context);
        }
    }
}

/// **The pure-water trial reproduces NeqSim's own stationary point.**
///
/// `WaterDewPointSeedProbe` reads NeqSim's own `TPmultiflash.stabilityAnalysis` adding an
/// aqueous phase at `x_water = 0.999993878160488` on this gas at its dew point. That number
/// is the whole oracle for this model - it is what the missing seeding has to reach - so it
/// is pinned here rather than only through the case's coarse fields.
#[test]
fn the_pure_water_trial_lands_on_neqsims_own_composition() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    let case = spec
        .cases
        .iter()
        .find(|case| case.id == "the_water_dew_point_gas_at_its_own_dew_point")
        .expect("the dew-point case should be in the spec");

    let result = solve(case);
    assert_eq!(result.path, MultiflashSeedPath::PureComponent);
    assert_eq!(result.seed_component_index, 7);
    let water = *result
        .seed_composition
        .last()
        .expect("the seed composition is one entry per component");
    common::assert_close(
        water,
        0.999993878160488,
        1.0e-06,
        "the aqueous composition NeqSim's own seeding reaches",
    );
}

/// **A stable feed says so, and hands back nothing to add.**
///
/// The warm case is ten kelvin from the dew point, so it is the threshold doing the work
/// rather than a model that always finds a phase: if this returned a composition the case
/// above would prove nothing about the boundary.
#[test]
fn a_stable_feed_offers_no_phase() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    let case = spec
        .cases
        .iter()
        .find(|case| case.id == "the_same_gas_ten_kelvin_warm_of_its_dew_point")
        .expect("the warm case should be in the spec");

    let result = solve(case);
    assert_eq!(result.path, MultiflashSeedPath::Stable);
    assert!(result.seed_composition.is_empty());
    assert_eq!(result.seed_fraction, 0.0);
    assert_eq!(result.seed_component_index, -1);
    // Every trial was attempted and every one was at or above the plane.
    assert!(!result.tm.is_empty());
    assert!(result.tm.iter().all(|distance| *distance >= -1.0e-08), "{:?}", result.tm);
}
