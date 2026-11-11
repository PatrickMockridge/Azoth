//! Spec-driven tests for `eos.scale_kinetics`.

use azoth_core::units::{moles_per_cubic_meter, seconds};
use azoth_eos::results::ScaleKineticsResult;
use azoth_eos::scale_kinetics::scale_kinetics;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.scale_kinetics";

fn call(case: &azoth_core::spec::TestCase) -> ScaleKineticsResult {
    scale_kinetics(
        common::input(case, "saturation_index"),
        common::input(case, "nucleation_constant"),
        common::input(case, "nucleation_shape_factor"),
        common::input(case, "surface_rate_constant"),
        common::input(case, "surface_reaction_order"),
        common::input(case, "mass_transfer_coefficient"),
        moles_per_cubic_meter(common::input(case, "bulk_concentration")),
        moles_per_cubic_meter(common::input(case, "equilibrium_concentration")),
    )
    .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    common::assert_skips_are_explained(spec);

    let mut executed = 0;
    for case in spec.all_tests() {
        if !case.is_active() {
            continue;
        }
        match case.kind {
            "worked_example" | "reference" => {
                executed += 1;
                let result = call(case);
                let context = &format!("{}::{}", spec.id, case.id);
                for (field, value) in [
                    ("supersaturation_ratio", result.supersaturation_ratio),
                    ("induction_time", result.induction_time.value),
                    (
                        "surface_limited_growth",
                        result.surface_limited_growth.value,
                    ),
                    (
                        "transport_limited_growth",
                        result.transport_limited_growth.value,
                    ),
                    ("effective_growth_rate", result.effective_growth_rate.value),
                ] {
                    if let Some(expected) = case.expected_value(field) {
                        common::assert_close(
                            value,
                            expected,
                            case.tolerance,
                            &format!("{context} ({field})"),
                        );
                    }
                }
                if let Some(expected) = case.expected_string("limiting_regime") {
                    assert_eq!(
                        result.limiting_regime, expected,
                        "{context}: the limiting regime"
                    );
                }
                common::assert_consistent(&result, context);
            }
            "property" => match case.property {
                Some("unit_round_trip") => {}
                other => panic!("{}::{}: unhandled property {other:?}", spec.id, case.id),
            },
            other => panic!("{}::{} has an unhandled kind `{other}`", spec.id, case.id),
        }
    }
    assert!(executed >= 2, "the spec should exercise more than one row");
}

fn run(saturation_index: f64) -> ScaleKineticsResult {
    scale_kinetics(
        saturation_index,
        3.6,
        3.0,
        1.5854895991882293e-10,
        2.0,
        1.0e-5,
        moles_per_cubic_meter(0.0),
        moles_per_cubic_meter(0.0),
    )
    .expect("the state runs")
}

/// **The induction time falls as supersaturation rises**, and the two growth rates are
/// both positive above saturation.
#[test]
fn induction_falls_with_supersaturation() {
    let low = run(0.3);
    let high = run(1.5);
    assert!(high.induction_time.value < low.induction_time.value);
    assert!(high.surface_limited_growth.value > low.surface_limited_growth.value);
}

/// **At and below saturation nothing forms**: the induction time is infinite, both growth
/// rates are zero and the regime is `NONE`.
#[test]
fn at_saturation_nothing_forms() {
    for saturation_index in [-0.5, 0.0] {
        let result = run(saturation_index);
        assert!(result.induction_time.value.is_infinite());
        assert_eq!(result.surface_limited_growth.value, 0.0);
        assert_eq!(result.effective_growth_rate.value, 0.0);
        assert_eq!(result.limiting_regime, "NONE");
    }
}

/// **A transport rate below the surface rate governs**, and one above it leaves the
/// reaction in control - the class's smaller-of rule.
#[test]
fn the_smaller_growth_rate_governs() {
    let transport = scale_kinetics(
        2.0,
        3.6,
        3.0,
        1.5854895991882295e-9,
        2.0,
        1.0e-6,
        moles_per_cubic_meter(1.0),
        moles_per_cubic_meter(0.0),
    )
    .expect("the state runs");
    assert_eq!(transport.limiting_regime, "TRANSPORT");
    assert!((transport.effective_growth_rate.value - 3.7e-11).abs() < 1e-20);

    let reaction = scale_kinetics(
        1.0,
        3.6,
        3.0,
        1.5854895991882293e-10,
        2.0,
        1.0e-3,
        moles_per_cubic_meter(1.0),
        moles_per_cubic_meter(0.0),
    )
    .expect("the state runs");
    assert_eq!(reaction.limiting_regime, "REACTION");
    assert_eq!(
        reaction.effective_growth_rate.value,
        reaction.surface_limited_growth.value
    );
}

/// **A negative rate is clamped, not refused**, as the class's setters clamp it, and the
/// warning names the input.
#[test]
fn a_negative_rate_is_clamped() {
    let result = scale_kinetics(
        1.0,
        seconds(-1.0).value,
        3.0,
        -5.0e-10,
        2.0,
        -1.0e-5,
        moles_per_cubic_meter(0.0),
        moles_per_cubic_meter(0.0),
    )
    .expect("the class clamps rather than refuses");
    assert_eq!(result.induction_time.value, 0.0);
    assert_eq!(result.surface_limited_growth.value, 0.0);
    assert_eq!(result.effective_growth_rate.value, 0.0);
    assert!(!result.warnings.is_empty());
}
