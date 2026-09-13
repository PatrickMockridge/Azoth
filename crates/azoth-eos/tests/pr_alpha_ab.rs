//! Spec-driven tests for `eos.pr_alpha_ab`.

use azoth_core::AzothError;
use azoth_eos::spec_gen;
use azoth_eos::{OMEGA_A, OMEGA_B, pr_alpha_ab};
use azoth_test_support as common;

const CALC_ID: &str = "eos.pr_alpha_ab";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::PrAlphaAbResult {
    pr_alpha_ab(
        common::input(case, "kappa"),
        common::input(case, "Tr"),
        common::input(case, "Pr"),
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
                for field in ["alpha", "a_reduced", "b_reduced"] {
                    let actual = match field {
                        "alpha" => result.alpha,
                        "a_reduced" => result.a_reduced,
                        _ => result.b_reduced,
                    };
                    common::assert_close(
                        actual,
                        common::expected(case, field),
                        case.tolerance,
                        &format!("{}::{} ({field})", spec.id, case.id),
                    );
                }
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "kappa" => Some(common::input(case, "kappa")),
                        "Tr" => Some(common::input(case, "Tr")),
                        "Pr" => Some(common::input(case, "Pr")),
                        "alpha" => Some(result.alpha),
                        "a_reduced" => Some(result.a_reduced),
                        "b_reduced" => Some(result.b_reduced),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
            }
            "property" => match case.property {
                Some("monotonic") => monotonic(),
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

/// The reduced cubic at a given `A` and `B`, as `(p, q)` of its depressed form.
///
/// `z**3 - (1-B)z**2 + (A - 3B**2 - 2B)z - (A*B - B**2 - B**3) = 0` written as
/// `z**3 + c2*z**2 + c1*z + c0` and depressed by `z = w - c2/3`. A triple root
/// leaves `p` and `q` both zero, and nothing else does - which is what makes this
/// a sharper check on the constants than comparing a root would be.
fn depressed(a: f64, b: f64) -> (f64, f64) {
    let (c2, c1, c0) = (
        -(1.0 - b),
        a - 3.0 * b * b - 2.0 * b,
        -(a * b - b * b - b * b * b),
    );
    let p = c1 - c2 * c2 / 3.0;
    let q = 2.0 * c2 * c2 * c2 / 27.0 - c2 * c1 / 3.0 + c0;
    (p, q)
}

/// The Omega constants are the values that make the critical point a triple root.
///
/// This is the claim the whole spec rests on, and it is a *derivable* one: the two
/// constants must satisfy `Omega_a = 3u**2 + 2u + (1-u)**2/3` and
/// `Omega_a*u - u**2 - u**3 = (1-u)**3/27` simultaneously, which fixes them. The
/// spec states the pair in full so the worked example is retraceable without the
/// paper, and this test is what makes that claim checkable rather than asserted.
///
/// Measured, not assumed: with the constants this crate ships, `p` comes out
/// **exactly zero** in f64 and `q` is 6.9e-18 - one ulp of the coefficient scale.
#[test]
fn the_omegas_make_the_critical_point_a_triple_root() {
    let (p, q) = depressed(OMEGA_A, OMEGA_B);
    assert!(
        p.abs() < 1e-12 && q.abs() < 1e-12,
        "the critical point is not a triple root: p = {p:e}, q = {q:e}"
    );
}

/// ...and the values the paper *prints* do not, which is why this crate does not use them.
///
/// This is the trap the spec's `notes` exist for. `0.45724` and
/// `0.07780` are roundings of the pair above, and rounding them breaks the
/// condition by three orders of magnitude more than the tolerance - because a
/// triple root is cubically ill-conditioned, so a 5e-6 error in the coefficient
/// moves the root by about its cube root, 1.7e-2.
///
/// Asserting the failure rather than only the success is the point: it pins *how*
/// wrong the printed pair is, so an implementer who "corrects" the constants to
/// match the paper makes this test fail and is told why, instead of silently
/// shifting every downstream number by 4.55% at the critical point.
#[test]
fn the_printed_omegas_do_not_and_that_is_why_they_are_not_used() {
    let (p, q) = depressed(0.45724, 0.07780);
    assert!(
        p.abs() > 1e-9,
        "expected the printed Omegas to break the triple-root condition; p = {p:e}"
    );
    assert!(
        q.abs() > 1e-9,
        "expected the printed Omegas to break the triple-root condition; q = {q:e}"
    );
}

/// `alpha` falls with `Tr`, `a_reduced` falls faster, `b_reduced` falls, and all
/// three rise with `Pr`.
///
/// Four relations, and each is what a different error breaks: a sign error in
/// `1 - sqrt(Tr)`, a missing `Tr**2` divisor, a missing `Tr` divisor, and a
/// swapped `Tr`/`Pr`. A point-value case can pass with any of those present.
#[test]
fn monotonic() {
    let at = |tr: f64, pr: f64| pr_alpha_ab(0.60282728832, tr, pr).unwrap();

    let base = at(0.8, 0.25);
    assert!(
        at(0.7, 0.25).alpha > base.alpha,
        "alpha must fall as Tr rises"
    );
    assert!(
        at(0.7, 0.25).a_reduced > base.a_reduced,
        "a_reduced must fall as Tr rises"
    );
    assert!(
        at(0.7, 0.25).b_reduced > base.b_reduced,
        "b_reduced must fall as Tr rises"
    );
    assert!(
        at(0.8, 0.5).alpha == base.alpha,
        "alpha does not depend on Pr"
    );
    assert!(
        at(0.8, 0.5).a_reduced > base.a_reduced,
        "a_reduced must rise with Pr"
    );
    assert!(
        at(0.8, 0.5).b_reduced > base.b_reduced,
        "b_reduced must rise with Pr"
    );

    // And the two reduced parameters are linear in Pr, which is what `* Pr` in one
    // place and not an exponent means. A quadratic would pass every inequality
    // above and fail this.
    let double = at(0.8, 0.5);
    assert!((double.a_reduced - 2.0 * base.a_reduced).abs() < 1e-15);
    assert!((double.b_reduced - 2.0 * base.b_reduced).abs() < 1e-15);
}

#[test]
fn a_non_positive_reduced_temperature_is_an_error() {
    // Tr appears as a square root and as a squared divisor, so zero is singular
    // and a negative value makes the root undefined. Neither is a limiting case.
    for tr in [0.0, -0.5] {
        let err = pr_alpha_ab(0.60282728832, tr, 0.25).unwrap_err();
        assert!(matches!(err, AzothError::OutOfRange { .. }));
        assert_eq!(err.field(), Some("Tr"), "for Tr = {tr}");
    }
}

#[test]
fn a_non_positive_reduced_pressure_is_an_error() {
    // Pr = 0 drives B to zero, and B is a divisor in the fugacity expression this
    // feeds. Refusing here means the failure is named where the reason is visible.
    for pr in [0.0, -0.1] {
        let err = pr_alpha_ab(0.60282728832, 0.8, pr).unwrap_err();
        assert!(matches!(err, AzothError::OutOfRange { .. }));
        assert_eq!(err.field(), Some("Pr"), "for Pr = {pr}");
    }
}

#[test]
fn alpha_is_exactly_one_at_the_critical_temperature_whatever_kappa_is() {
    // The structural fact the critical-point case rests on, and the reason that
    // case can pin the constants with no arithmetic in the way: `1 - sqrt(1)` is
    // exactly zero, so the square is exactly 1 for any finite kappa.
    for kappa in [-0.5, 0.0, 0.60282728832, 3.0] {
        let r = pr_alpha_ab(kappa, 1.0, 1.0).unwrap();
        assert_eq!(r.alpha.to_bits(), 1.0_f64.to_bits(), "for kappa = {kappa}");
    }
}
