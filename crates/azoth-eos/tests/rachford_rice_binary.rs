//! Spec-driven tests for `eos.rachford_rice_binary`.

use azoth_core::{AzothError, CalcResult};
use azoth_eos::rachford_rice_binary;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.rachford_rice_binary";

/// The Rachford-Rice sum, written out here rather than imported.
///
/// The calc returns the closed form; this is the equation the closed form solves.
/// Keeping them separate is the point - the residual test below compares the answer
/// against the *equation*, not against the algebra that produced it, so an error in
/// the derivation cannot hide by being present in both.
fn rachford_rice_sum(z1: f64, k1: f64, k2: f64, beta: f64) -> f64 {
    z1 * (k1 - 1.0) / (1.0 + beta * (k1 - 1.0))
        + (1.0 - z1) * (k2 - 1.0) / (1.0 + beta * (k2 - 1.0))
}

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::RachfordRiceBinaryResult {
    rachford_rice_binary(
        common::input(case, "z1"),
        common::input(case, "K1"),
        common::input(case, "K2"),
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
                    result.beta,
                    common::expected(case, "beta"),
                    case.tolerance,
                    &format!("{}::{} (beta)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "z1" => Some(common::input(case, "z1")),
                        "K1" => Some(common::input(case, "K1")),
                        "K2" => Some(common::input(case, "K2")),
                        "beta" => Some(result.beta),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
            }
            "property" => match case.property {
                Some("consistency_with") => {}
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

/// The returned `beta` solves the Rachford-Rice equation, over a broad sweep.
///
/// The check that does not depend on knowing the right answer, which is what matters
/// for a closed form with five roundings in it. Swept across compositions and
/// K-value pairs, **including pairs whose `beta` falls outside `[0, 1]`** - the
/// equation still has a solution there and the residual must still vanish. A solver
/// that clamped `beta` to the physical interval, or that bailed out when the root
/// left it, fails here while passing every point-value case that happens to be
/// two-phase.
#[test]
fn the_returned_beta_solves_the_equation() {
    let mut checked = 0_u32;
    let mut outside = 0_u32;
    for z1_step in 1..10 {
        for k1 in [1.5, 2.0, 4.0, 8.0, 30.0] {
            for k2 in [0.02, 0.2, 0.5, 0.8, 0.95] {
                let z1 = f64::from(z1_step) / 10.0;
                let r = rachford_rice_binary(z1, k1, k2).unwrap();
                let residual = rachford_rice_sum(z1, k1, k2, r.beta);
                assert!(
                    residual.abs() < 1e-12,
                    "z1 = {z1}, K1 = {k1}, K2 = {k2}: beta = {} leaves a residual of {residual:e}",
                    r.beta
                );
                if !(0.0..=1.0).contains(&r.beta) {
                    outside += 1;
                }
                checked += 1;
            }
        }
    }
    assert!(
        checked > 100,
        "the sweep should be broad, checked {checked}"
    );
    assert!(
        outside > 0,
        "the sweep should reach single-phase feeds; none produced a beta outside [0, 1]"
    );
}

/// A feed whose K-values all exceed 1 does not split, and says so.
///
/// The warning is on the *output* rather than on the inputs because whether the feed
/// splits depends on all three, not on any one of them. Clamping `beta` to `[0, 1]`
/// would pass a value check and fail here, and would also discard the information
/// that the feed is subcooled liquid rather than merely not-split.
#[test]
fn a_feed_that_does_not_split_carries_a_warning() {
    let r = rachford_rice_binary(0.5, 2.0, 1.5).unwrap();
    assert!(
        r.beta < 0.0,
        "both K > 1 should give a negative beta, got {}",
        r.beta
    );
    assert!(!r.is_clean(), "a beta outside [0, 1] must be visible");

    // And the mirror image: both K < 1 gives beta > 1, the superheated-vapour case.
    let vapour = rachford_rice_binary(0.5, 0.8, 0.5).unwrap();
    assert!(
        vapour.beta > 1.0,
        "both K < 1 should give beta > 1, got {}",
        vapour.beta
    );
    assert!(!vapour.is_clean());
}

/// A two-phase feed is silent.
#[test]
fn a_two_phase_feed_carries_no_warning() {
    for (z1, k1, k2) in [(0.6, 4.0, 0.25), (0.3, 5.0, 0.2), (0.5, 3.0, 0.3)] {
        let r = rachford_rice_binary(z1, k1, k2).unwrap();
        assert!(
            (0.0..=1.0).contains(&r.beta),
            "expected a two-phase feed for z1 = {z1}, K1 = {k1}, K2 = {k2}"
        );
        assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
    }
}

#[test]
fn a_k_value_of_one_is_an_error() {
    // `K - 1` is a divisor, and the component degenerates besides: with K1 = 1 it
    // distributes equally and drops out of the sum, leaving a one-component problem.
    for (k1, k2) in [(1.0, 0.5), (2.0, 1.0)] {
        let err = rachford_rice_binary(0.5, k1, k2).unwrap_err();
        assert!(matches!(err, AzothError::OutOfRange { .. }));
    }
}

#[test]
fn a_non_positive_k_value_is_an_error() {
    for (k1, k2) in [(0.0, 0.5), (-1.0, 0.5), (2.0, 0.0)] {
        let err = rachford_rice_binary(0.5, k1, k2).unwrap_err();
        assert!(matches!(err, AzothError::OutOfRange { .. }));
    }
}

#[test]
fn a_mole_fraction_outside_zero_to_one_is_an_error() {
    for z1 in [-0.01, 1.01] {
        let err = rachford_rice_binary(z1, 4.0, 0.25).unwrap_err();
        assert!(matches!(err, AzothError::OutOfRange { .. }));
        assert_eq!(err.field(), Some("z1"));
    }
}
