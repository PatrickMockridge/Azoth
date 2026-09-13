//! Spec-driven tests for `hydraulics.friction_factor_colebrook`.

use azoth_test_support as common;

use azoth_core::spec::TestCase;
use azoth_core::units::{kilograms_per_cubic_meter, meters, meters_per_second, pascal_seconds};
use azoth_core::{AzothError, CalcResult, WarningCode};
use azoth_hydraulics::{
    friction_factor_colebrook, friction_factor_swamee_jain, fully_rough_limit, reynolds_number,
    spec_gen,
};

const CALC_ID: &str = "hydraulics.friction_factor_colebrook";

fn call(case: &TestCase) -> azoth_hydraulics::ColebrookResult {
    friction_factor_colebrook(
        common::input(case, "re"),
        common::input(case, "relative_roughness"),
    )
    .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
}

fn resolve(case: &TestCase, f: f64) -> impl Fn(&str) -> Option<f64> + '_ {
    move |q| match q {
        "re" => Some(common::input(case, "re")),
        "relative_roughness" => Some(common::input(case, "relative_roughness")),
        "f" => Some(f),
        _ => None,
    }
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
                assert!(
                    result.converged,
                    "{}::{}: the solver reported non-convergence, which should have been \
                     returned as an error",
                    spec.id, case.id
                );
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
                    resolve(case, result.f),
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

/// The same physical state described in SI and US customary units must give the
/// same friction factor.
///
/// Both inputs here are dimensionless, so the round trip runs the whole chain -
/// velocity and diameter in feet, through Reynolds number, into the friction
/// factor - rather than testing a conversion on a ratio that has no units to
/// convert.
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

    let f_si = friction_factor_colebrook(re_si, relative_roughness)
        .unwrap()
        .f;
    let f_us = friction_factor_colebrook(re_us, relative_roughness)
        .unwrap()
        .f;

    common::assert_same_state(f_si, f_us, 1e-12, "f from SI vs US customary state");
}

#[test]
fn converges_to_the_fully_rough_asymptote() {
    // As Re grows the roughness term dominates and the answer approaches an
    // explicit limit. This confirms the solver lands on the right fixed point
    // rather than on some nearby one that also satisfies its own stopping rule.
    let rr = 0.01;
    let solved = friction_factor_colebrook(1e12, rr).unwrap().f;
    let asymptote = fully_rough_limit(rr);
    common::assert_close(solved, asymptote, 1e-8, "fully rough limit at Re = 1e12");
}

#[test]
fn smooth_pipe_agrees_with_the_textbook_value() {
    // At zero roughness Colebrook reduces to the smooth-pipe correlation, for
    // which the accepted value at Re = 1e5 is about 0.0180. This pins the
    // implementation to an independent value rather than to itself.
    let f = friction_factor_colebrook(1e5, 0.0).unwrap().f;
    assert!(
        (0.0178..0.0182).contains(&f),
        "smooth-pipe f at Re=1e5 was {f}, outside the textbook range"
    );
}

#[test]
fn iteration_count_is_deterministic() {
    // Deterministic iteration count is what makes the Python and Rust results
    // comparable rather than merely close.
    let a = friction_factor_colebrook(1e5, 4.6e-4).unwrap();
    let b = friction_factor_colebrook(1e5, 4.6e-4).unwrap();
    assert_eq!(a.iterations, b.iterations);
    assert_eq!(a.f.to_bits(), b.f.to_bits(), "f must be bit-identical");
}

#[test]
fn below_transition_warns_but_still_returns_a_number() {
    // Out of the fitted range, not undefined: a warning, not an error.
    let r = friction_factor_colebrook(3000.0, 4.6e-4).unwrap();
    assert!(r.has_warning(WarningCode::OutOfValidRange));
    assert!(r.f.is_finite());
    assert!(r.converged);
}

#[test]
fn zero_reynolds_number_is_an_error() {
    // 2.51/(Re*sqrt(f)) is singular at Re = 0, so this is undefined rather than
    // out of range.
    let err = friction_factor_colebrook(0.0, 4.6e-4).unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("re"));
}

#[test]
fn negative_roughness_is_an_error() {
    let err = friction_factor_colebrook(1e5, -1e-4).unwrap_err();
    assert_eq!(err.field(), Some("relative_roughness"));
}

#[test]
fn comparably_accurate_to_swamee_jain_within_the_published_claim() {
    // The two correlations in this slice have to be consistent with each other.
    // Swamee and Jain claim their explicit form is within 1% of Colebrook over
    // 5000 < Re < 1e8, so that published figure is the tolerance here rather
    // than a number chosen to make the test pass.
    //
    // The points are strictly inside the range, not on its boundary - see
    // `swamee_jain_exceeds_its_claim_at_the_lower_bound` for why.
    for (re, rr) in [(1e4, 4.6e-4), (1e5, 4.6e-4), (1e6, 1e-3), (1e7, 5e-3)] {
        let f_c = friction_factor_colebrook(re, rr).unwrap().f;
        let f_sj = friction_factor_swamee_jain(re, rr).unwrap().f;
        common::assert_close(
            f_sj,
            f_c,
            1e-2,
            &format!("swamee-jain vs colebrook at Re={re}, rr={rr}"),
        );
    }
}

#[test]
fn swamee_jain_exceeds_its_claim_at_the_lower_bound() {
    // Measured, not assumed: the explicit approximation degrades as Re
    // approaches the lower end of its stated range. At Re = 5000 with
    // epsilon/D = 4.6e-4 it is 1.39% from Colebrook, against a headline claim of
    // 1%; by Re = 1e4 it is back to 0.59%.
    //
    // This is recorded rather than hidden, for two reasons. It documents that
    // the published 1% is not uniform across the claimed range, and it means a
    // future change that silently made the low-Re behaviour *better or worse*
    // would be noticed rather than absorbed by a loosened tolerance.
    let (re, rr) = (5000.0, 4.6e-4);
    let f_c = friction_factor_colebrook(re, rr).unwrap().f;
    let f_sj = friction_factor_swamee_jain(re, rr).unwrap().f;
    let rel = (f_sj - f_c).abs() / f_c;
    assert!(
        (0.013..0.015).contains(&rel),
        "swamee-jain is {rel:.4} from colebrook at Re=5000; expected about 0.0139. \
         If this moved, the low-Re behaviour of the approximation changed - check \
         whether that was intended before adjusting this range."
    );
}
