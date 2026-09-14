//! Spec-driven tests for `hydraulics.reynolds_number`.
//!
//! Every case here comes from `specs/calcs/hydraulics/reynolds_number.toml`. Add
//! a test to that file and it runs here and in Python with no new test code.

use azoth_test_support as common;

use azoth_core::spec::TestCase;
use azoth_core::units::{kilograms_per_cubic_meter, meters, meters_per_second, pascal_seconds};
use azoth_core::{AzothError, CalcResult, FlowRegime, WarningCode};
use azoth_hydraulics::{regime_warning, reynolds_number, spec_gen};

const CALC_ID: &str = "hydraulics.reynolds_number";

fn call(case: &TestCase) -> azoth_hydraulics::ReynoldsNumberResult {
    reynolds_number(
        kilograms_per_cubic_meter(common::input(case, "rho")),
        meters_per_second(common::input(case, "v")),
        meters(common::input(case, "D")),
        pascal_seconds(common::input(case, "mu")),
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
                    result.re,
                    common::expected(case, "re"),
                    case.tolerance,
                    &format!("{}::{} (re)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |q| match q {
                        "rho" => Some(common::input(case, "rho")),
                        "v" => Some(common::input(case, "v")),
                        "D" => Some(common::input(case, "D")),
                        "mu" => Some(common::input(case, "mu")),
                        "re" => Some(result.re),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
                if case.kind == "worked_example" {
                    assert_eq!(
                        result.regime,
                        FlowRegime::from_reynolds_number(result.re),
                        "regime must follow from the Reynolds number"
                    );
                }
            }
            "property" => match case.property {
                Some("unit_round_trip") => unit_round_trip(),
                other => panic!("{}::{}: unhandled property {other:?}", spec.id, case.id),
            },
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 3,
        "expected several active cases, ran {executed}"
    );
}

/// The same physical state described in SI and in US customary units must give
/// the same Reynolds number.
///
/// Reynolds number is dimensionless, so this is the sharpest check available
/// that unit handling is right: any stray conversion factor that survived into
/// the ratio would show up as a discrepancy rather than cancelling.
#[test]
fn unit_round_trip() {
    let si = reynolds_number(
        kilograms_per_cubic_meter(998.0),
        meters_per_second(1.5),
        meters(0.1),
        pascal_seconds(1.002e-3),
    )
    .expect("SI state should compute");

    // The identical state, with length and velocity expressed in feet.
    // 0.1 m = 0.32808399... ft, 1.5 m/s = 4.92125984... ft/s.
    let us = reynolds_number(
        kilograms_per_cubic_meter(998.0),
        uom::si::f64::Velocity::new::<uom::si::velocity::foot_per_second>(1.5 / 0.3048),
        uom::si::f64::Length::new::<uom::si::length::foot>(0.1 / 0.3048),
        pascal_seconds(1.002e-3),
    )
    .expect("US customary state should compute");

    common::assert_same_state(si.re, us.re, 1e-12, "Re from SI vs US customary");
    assert_eq!(si.regime, us.regime);
}

#[test]
fn hard_bounds_raise_rather_than_warn() {
    // Each of these makes the Reynolds number undefined or meaningless, so they
    // are errors. A warning here would return inf dressed as a result.
    //
    // Written out one at a time rather than tabulated: the table version needed
    // a `fn() -> AzothError` per row, which is more machinery than four
    // assertions are worth.
    let rho = |value| kilograms_per_cubic_meter(value);
    let vel = |value| meters_per_second(value);
    let dia = |value| meters(value);
    let vis = |value| pascal_seconds(value);

    let expect = |field: &str, err: AzothError| {
        assert_eq!(
            err.field(),
            Some(field),
            "error for `{field}` named the wrong field: {err}"
        );
    };

    expect(
        "rho",
        reynolds_number(rho(0.0), vel(1.5), dia(0.1), vis(1e-3)).unwrap_err(),
    );
    expect(
        "v",
        reynolds_number(rho(998.0), vel(-1.0), dia(0.1), vis(1e-3)).unwrap_err(),
    );
    expect(
        "D",
        reynolds_number(rho(998.0), vel(1.5), dia(0.0), vis(1e-3)).unwrap_err(),
    );
    expect(
        "mu",
        reynolds_number(rho(998.0), vel(1.5), dia(0.1), vis(0.0)).unwrap_err(),
    );
}

#[test]
fn laminar_flow_gets_no_regime_warning() {
    // The transitional band is the only regime that warns. Laminar and turbulent
    // are both determinate, so warning about them would be noise.
    let laminar = reynolds_number(
        kilograms_per_cubic_meter(1000.0),
        meters_per_second(0.01),
        meters(0.05),
        pascal_seconds(1e-2),
    )
    .unwrap();
    assert_eq!(laminar.regime, FlowRegime::Laminar);
    assert!(!laminar.has_warning(WarningCode::TransitionalFlow));
    assert!(laminar.is_clean(), "{:?}", laminar.warnings);
    assert!(regime_warning(laminar.regime).is_none());
}

#[test]
fn transitional_flow_warns_with_the_specific_code() {
    // The spec names TRANSITIONAL_FLOW for this band, so a caller can react to
    // indeterminate friction without string-matching a message.
    let t = reynolds_number(
        kilograms_per_cubic_meter(1000.0),
        meters_per_second(0.1),
        meters(0.05),
        pascal_seconds(1.5e-3),
    )
    .unwrap();
    assert_eq!(t.regime, FlowRegime::Transitional);
    assert!(t.has_warning(WarningCode::TransitionalFlow));
    assert!(
        !t.has_warning(WarningCode::OutOfValidRange),
        "the band check names its own code; the generic one must not also fire"
    );
    assert!(regime_warning(t.regime).is_some());
}

#[test]
fn nan_input_is_rejected() {
    // NaN comparisons are all false, so without the explicit check in
    // RangeCheck::violated a NaN would pass every bound.
    let err = reynolds_number(
        kilograms_per_cubic_meter(f64::NAN),
        meters_per_second(1.5),
        meters(0.1),
        pascal_seconds(1e-3),
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
}
