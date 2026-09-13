//! Spec-driven tests for `thermal.conduction_plane_wall`.

use azoth_core::units::{
    Area, Length, TemperatureInterval, kelvin_intervals, meters, square_meters,
    watts_per_meter_kelvin,
};
use azoth_core::{AzothError, CalcResult};
use azoth_test_support as common;
use azoth_thermal::conduction_plane_wall;
use azoth_thermal::spec_gen;

const CALC_ID: &str = "thermal.conduction_plane_wall";

fn call(case: &azoth_core::spec::TestCase) -> azoth_thermal::ConductionPlaneWallResult {
    conduction_plane_wall(
        watts_per_meter_kelvin(common::input(case, "k")),
        square_meters(common::input(case, "A")),
        kelvin_intervals(common::input(case, "dT")),
        meters(common::input(case, "L")),
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
                let result = call(case);
                common::assert_close(
                    result.q.value,
                    common::expected(case, "q"),
                    case.tolerance,
                    &format!("{}::{} (q)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "k" => Some(common::input(case, "k")),
                        "A" => Some(common::input(case, "A")),
                        "dT" => Some(common::input(case, "dT")),
                        "L" => Some(common::input(case, "L")),
                        "q" => Some(result.q.value),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
            }
            "property" => match case.property {
                Some("monotonic") => monotonic(),
                Some("unit_round_trip") => unit_round_trip(),
                other => panic!("{}::{}: unhandled property {other:?}", spec.id, case.id),
            },
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 4,
        "expected several active cases, ran {executed}"
    );
}

/// `q` rises with each of `k`, `A` and `dT`, and falls with `L`.
///
/// These are the four partial derivatives of the equation. A point-value test can
/// pass with a wrong formula that happens to be right at one set of inputs; an
/// inverted ratio or a sign error cannot survive all four directions at once.
#[test]
fn monotonic() {
    let at = |k: f64, a: f64, dt: f64, l: f64| {
        conduction_plane_wall(
            watts_per_meter_kelvin(k),
            square_meters(a),
            kelvin_intervals(dt),
            meters(l),
        )
        .unwrap()
        .q
        .value
    };
    let base = at(45.0, 2.0, 30.0, 0.05);

    assert!(at(90.0, 2.0, 30.0, 0.05) > base, "q must increase with k");
    assert!(at(45.0, 4.0, 30.0, 0.05) > base, "q must increase with A");
    assert!(at(45.0, 2.0, 60.0, 0.05) > base, "q must increase with dT");
    assert!(at(45.0, 2.0, 30.0, 0.10) < base, "q must decrease with L");
}

/// The same physical wall described in SI and in US customary units.
///
/// Three of the four inputs are converted: the area to square feet, the thickness
/// to feet, and the temperature difference to a Fahrenheit interval - which is the
/// interesting one, because a *difference* of 30 K is 54 degF, while an absolute
/// 30 K is -243.15 degC. The conductivity stays in SI because `uom` carries no US
/// customary thermal conductivity, and the Python half of this test is kept
/// identical rather than converting one more quantity in one language only.
#[test]
fn unit_round_trip() {
    let si = conduction_plane_wall(
        watts_per_meter_kelvin(45.0),
        square_meters(2.0),
        kelvin_intervals(30.0),
        meters(0.05),
    )
    .unwrap()
    .q
    .value;

    let us = conduction_plane_wall(
        watts_per_meter_kelvin(45.0),
        Area::new::<uom::si::area::square_foot>(2.0 / 0.3048 / 0.3048),
        TemperatureInterval::new::<uom::si::temperature_interval::degree_fahrenheit>(
            30.0 * 9.0 / 5.0,
        ),
        Length::new::<uom::si::length::foot>(0.05 / 0.3048),
    )
    .unwrap()
    .q
    .value;

    common::assert_same_state(si, us, 1e-12, "q from SI vs US customary state");
}

#[test]
fn a_non_positive_conductivity_is_an_error() {
    let err = conduction_plane_wall(
        watts_per_meter_kelvin(0.0),
        square_meters(2.0),
        kelvin_intervals(30.0),
        meters(0.05),
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("k"));
}

#[test]
fn a_zero_thickness_is_an_error_rather_than_infinite_heat_flow() {
    // `L` is a divisor, so zero is the singular limit. Returning infinity would be
    // arithmetic; refusing is the honest answer, because a zero-thickness wall is
    // not a state the plane-wall model describes.
    let err = conduction_plane_wall(
        watts_per_meter_kelvin(45.0),
        square_meters(2.0),
        kelvin_intervals(30.0),
        meters(0.0),
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("L"));
}

#[test]
fn a_negative_temperature_difference_reverses_the_flow_instead_of_failing() {
    // The one modelling decision this calc makes, and the reason `dT` has no lower
    // bound: it models a difference across a slab, not a named hot face, so the
    // signed answer is more useful than a refusal. A result carrying no warnings is
    // the point - a negative heat flow is a direction, not a problem.
    let r = conduction_plane_wall(
        watts_per_meter_kelvin(45.0),
        square_meters(2.0),
        kelvin_intervals(-30.0),
        meters(0.05),
    )
    .unwrap();
    assert!((r.q.value + 54000.0).abs() < 1e-9, "got {}", r.q.value);
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}
