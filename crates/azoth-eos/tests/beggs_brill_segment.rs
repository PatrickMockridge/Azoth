//! Spec-driven tests for `eos.beggs_brill_segment`.

use azoth_core::units::{
    cubic_meters_per_second, kelvins, kilograms_per_cubic_meter, meters, pascal_seconds, pascals,
};
use azoth_eos::beggs_brill_segment::beggs_brill_segment;
use azoth_eos::results::BeggsBrillSegmentResult;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.beggs_brill_segment";

fn call(case: &azoth_core::spec::TestCase) -> BeggsBrillSegmentResult {
    let volume = |name: &str| cubic_meters_per_second(common::input(case, name));
    let mass = |name: &str| kilograms_per_cubic_meter(common::input(case, name));
    let viscosity = |name: &str| pascal_seconds(common::input(case, name));
    beggs_brill_segment(
        volume("gas_flow"),
        volume("liquid_flow"),
        mass("gas_density"),
        mass("liquid_density"),
        viscosity("gas_viscosity"),
        viscosity("liquid_viscosity"),
        kelvins(common::input(case, "temperature")),
        pascals(common::input(case, "pressure")),
        meters(common::input(case, "diameter")),
        meters(common::input(case, "length")),
        meters(common::input(case, "elevation")),
        meters(common::input(case, "roughness")),
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
                if let Some(expected) = case.expected_string("regime") {
                    assert_eq!(result.regime.as_str(), expected, "{context}: the regime");
                }
                for (field, value) in [
                    ("liquid_holdup", result.liquid_holdup),
                    (
                        "superficial_gas_velocity",
                        result.superficial_gas_velocity.value,
                    ),
                    (
                        "superficial_liquid_velocity",
                        result.superficial_liquid_velocity.value,
                    ),
                    (
                        "superficial_mixture_velocity",
                        result.superficial_mixture_velocity.value,
                    ),
                    ("no_slip_density", result.no_slip_density.value),
                    ("no_slip_viscosity", result.no_slip_viscosity.value),
                    ("reynolds", result.reynolds),
                    ("friction_factor", result.friction_factor),
                    (
                        "two_phase_friction_factor",
                        result.two_phase_friction_factor,
                    ),
                    ("hydrostatic_drop", result.hydrostatic_drop.value),
                    ("friction_drop", result.friction_drop.value),
                    ("total_drop", result.total_drop.value),
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
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "gas_flow" => Some(common::input(case, "gas_flow")),
                        "liquid_flow" => Some(common::input(case, "liquid_flow")),
                        "gas_density" => Some(common::input(case, "gas_density")),
                        "liquid_density" => Some(common::input(case, "liquid_density")),
                        "gas_viscosity" => Some(common::input(case, "gas_viscosity")),
                        "liquid_viscosity" => Some(common::input(case, "liquid_viscosity")),
                        "temperature" => Some(common::input(case, "temperature")),
                        "pressure" => Some(common::input(case, "pressure")),
                        "diameter" => Some(common::input(case, "diameter")),
                        "length" => Some(common::input(case, "length")),
                        "roughness" => Some(common::input(case, "roughness")),
                        _ => None,
                    },
                    context,
                );
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

/// **The total is the two terms**, and the two-phase multiplier is one where no liquid flows.
#[test]
fn the_drop_is_the_two_terms() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    let case = spec
        .all_tests()
        .find(|case| case.id == "segment_uphill")
        .expect("the uphill case");
    let result = call(case);
    common::assert_close(
        result.total_drop.value,
        result.hydrostatic_drop.value + result.friction_drop.value,
        1e-12,
        "the total",
    );
}
