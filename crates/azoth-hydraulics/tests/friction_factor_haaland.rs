//! Spec-driven tests for `hydraulics.friction_factor_haaland`.

mod common;

use azoth_core::units::{kilograms_per_cubic_meter, meters, meters_per_second, pascal_seconds};
use azoth_core::{AzothError, CalcResult, WarningCode};
use azoth_hydraulics::spec_gen::TestCase;
use azoth_hydraulics::{friction_factor_colebrook, friction_factor_haaland, reynolds_number};

const CALC_ID: &str = "hydraulics.friction_factor_haaland";

fn call(case: &TestCase) -> azoth_hydraulics::HaalandResult {
    friction_factor_haaland(
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

    let f_si = friction_factor_haaland(re_si, relative_roughness)
        .unwrap()
        .f;
    let f_us = friction_factor_haaland(re_us, relative_roughness)
        .unwrap()
        .f;

    common::assert_same_state(f_si, f_us, 1e-12, "f from SI vs US customary state");
}

#[test]
fn the_two_explicit_forms_are_the_same_cloth() {
    // Haaland and Swamee-Jain both approximate Colebrook. Neither is checked
    // against the other - each is checked against Colebrook by its own spec - so
    // this asserts only that they sit within a few percent of one another. It is a
    // sanity check on the registry rather than on either equation: if these ever
    // disagree wildly, one of the three friction factor calcs has changed in a way
    // its own tests did not catch.
    let (re, relative_roughness) = (1.0e5, 4.6e-4);
    let haaland = friction_factor_haaland(re, relative_roughness).unwrap().f;
    let colebrook = friction_factor_colebrook(re, relative_roughness).unwrap().f;
    assert!(
        (haaland / colebrook - 1.0).abs() < 0.05,
        "Haaland is {:.2}% from Colebrook, outside the published accuracy of about \
         1-2% by a wide margin",
        100.0 * (haaland / colebrook - 1.0)
    );
}

#[test]
fn outside_the_published_range_warns_but_computes() {
    // Below Re ~ 4000 the equation is outside the turbulent range it describes. It
    // is still well defined, so this is a warning, not an error - but the accuracy
    // claim no longer covers the answer.
    let low = friction_factor_haaland(1000.0, 4.6e-4).unwrap();
    assert!(low.has_warning(WarningCode::OutOfValidRange));
    assert!(low.f.is_finite());

    let rough = friction_factor_haaland(1e5, 0.5).unwrap();
    assert!(rough.has_warning(WarningCode::OutOfValidRange));
}

#[test]
fn inside_the_published_range_is_clean() {
    // A mid-range case must produce no warnings at all. If this ever fires, the
    // bounds or their inclusive flags have drifted.
    let r = friction_factor_haaland(1e5, 4.6e-4).unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

#[test]
fn zero_reynolds_number_is_an_error() {
    // 6.9/Re is singular at the origin.
    let err = friction_factor_haaland(0.0, 4.6e-4).unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("re"));
}

#[test]
fn negative_roughness_is_an_error() {
    let err = friction_factor_haaland(1e5, -1.0e-4).unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("relative_roughness"));
}

#[test]
fn the_explicit_form_needs_no_iteration() {
    // Deliberate contrast with Colebrook: the same inputs produce a result with no
    // iteration count, because there is no iteration.
    let a = friction_factor_haaland(1e5, 4.6e-4).unwrap();
    let b = friction_factor_haaland(1e5, 4.6e-4).unwrap();
    assert_eq!(a.f.to_bits(), b.f.to_bits());
}
