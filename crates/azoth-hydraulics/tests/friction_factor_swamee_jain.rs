//! Spec-driven tests for `hydraulics.friction_factor_swamee_jain`.

mod common;

use azoth_core::units::{kilograms_per_cubic_meter, meters, meters_per_second, pascal_seconds};
use azoth_core::{AzothError, CalcResult, WarningCode};
use azoth_hydraulics::spec_gen::TestCase;
use azoth_hydraulics::{friction_factor_swamee_jain, reynolds_number};

const CALC_ID: &str = "hydraulics.friction_factor_swamee_jain";

fn call(case: &TestCase) -> azoth_hydraulics::SwameeJainResult {
    friction_factor_swamee_jain(
        common::input(case, "re"),
        common::input(case, "relative_roughness"),
    )
    .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = common::spec(CALC_ID);
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
                    result.f,
                    common::expected(case, "f"),
                    case.tolerance,
                    &format!("{}::{} (f)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |q| match q {
                        "re" => Some(common::input(case, "re")),
                        "relative_roughness" => Some(common::input(case, "relative_roughness")),
                        "f" => Some(result.f),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
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
        executed >= 4,
        "expected several active cases, ran {executed}"
    );
}

#[test]
fn unit_round_trip() {
    let relative_roughness = 4.6e-4;

    let re_si = reynolds_number(
        kilograms_per_cubic_meter(998.0),
        meters_per_second(1.5),
        meters(0.1),
        pascal_seconds(1.002e-3),
    )
    .unwrap()
    .re;

    let re_us = reynolds_number(
        kilograms_per_cubic_meter(998.0),
        uom::si::f64::Velocity::new::<uom::si::velocity::foot_per_second>(1.5 / 0.3048),
        uom::si::f64::Length::new::<uom::si::length::foot>(0.1 / 0.3048),
        pascal_seconds(1.002e-3),
    )
    .unwrap()
    .re;

    let f_si = friction_factor_swamee_jain(re_si, relative_roughness)
        .unwrap()
        .f;
    let f_us = friction_factor_swamee_jain(re_us, relative_roughness)
        .unwrap()
        .f;

    common::assert_same_state(f_si, f_us, 1e-12, "f from SI vs US customary state");
}

#[test]
fn outside_the_published_range_warns_but_computes() {
    // Below Re = 5000 the equation is outside the range the paper fitted it to.
    // It is still well defined, so this is a warning, not an error - but the
    // accuracy claim no longer covers the answer.
    let low = friction_factor_swamee_jain(1000.0, 4.6e-4).unwrap();
    assert!(low.has_warning(WarningCode::OutOfValidRange));
    assert!(low.f.is_finite());

    // Above the upper Reynolds bound.
    let high = friction_factor_swamee_jain(1e9, 4.6e-4).unwrap();
    assert!(high.has_warning(WarningCode::OutOfValidRange));
}

#[test]
fn inside_the_published_range_is_clean() {
    // A mid-range case must produce no warnings at all. If this ever fires, the
    // bounds or the inclusive flags have drifted.
    let r = friction_factor_swamee_jain(1e5, 4.6e-4).unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

#[test]
fn zero_reynolds_number_is_an_error() {
    // Re**0.9 is zero at the origin, so 5.74/Re**0.9 is singular.
    let err = friction_factor_swamee_jain(0.0, 4.6e-4).unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("re"));
}

#[test]
fn the_explicit_form_needs_no_iteration() {
    // Deliberate contrast with Colebrook: the same inputs produce a result with
    // no iteration count, because there is no iteration. The two result types
    // differ for this reason.
    let a = friction_factor_swamee_jain(1e5, 4.6e-4).unwrap();
    let b = friction_factor_swamee_jain(1e5, 4.6e-4).unwrap();
    assert_eq!(a.f.to_bits(), b.f.to_bits());
}
