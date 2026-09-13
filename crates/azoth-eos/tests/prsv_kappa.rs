//! Spec-driven tests for `eos.prsv_kappa`.

use azoth_core::AzothError;
use azoth_eos::spec_gen;
use azoth_eos::{pr_alpha_ab, prsv_kappa};
use azoth_test_support as common;

const CALC_ID: &str = "eos.prsv_kappa";

/// A deliberately arbitrary, non-zero `kappa1`.
///
/// This crate ships no fitted values for the parameter - a table of them is the
/// databank it does not have - so tests pick a round number and say so rather than
/// borrowing one that looks like data.
const SOME_KAPPA1: f64 = 0.05;

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::PrsvKappaResult {
    prsv_kappa(
        common::input(case, "omega"),
        common::input(case, "Tr"),
        common::input(case, "kappa1"),
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
                    result.kappa,
                    common::expected(case, "kappa"),
                    case.tolerance,
                    &format!("{}::{} (kappa)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "omega" => Some(common::input(case, "omega")),
                        "Tr" => Some(common::input(case, "Tr")),
                        "kappa1" => Some(common::input(case, "kappa1")),
                        "kappa" => Some(result.kappa),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
            }
            "property" => match case.property {
                Some("consistency_with") | Some("monotonic") => {}
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

/// The PRSV coefficient feeds Peng-Robinson's alpha function unchanged.
///
/// The whole argument for splitting the coefficients out of the equation of state
/// rests on this: PRSV modifies the coefficient's temperature dependence and
/// nothing else, so `eos.pr_alpha_ab` served it without one line changing. If that
/// were wrong, `pr_alpha_ab` would have needed a variant and the decomposition
/// would have been a fork after all - which is what this test is here to catch.
///
/// Cross-calc rather than single-calc, which is why it is declared in the spec and
/// implemented by hand in both test files instead of generated.
#[test]
fn the_prsv_coefficient_feeds_the_pr_alpha_function() {
    for (omega, tr, kappa1) in [(0.152, 0.8, SOME_KAPPA1), (0.01142, 1.2, -0.03)] {
        let kappa = prsv_kappa(omega, tr, kappa1).unwrap().kappa;
        let ab = pr_alpha_ab(kappa, tr, 0.25).unwrap();

        // The alpha function is the published relation, restated here so the test
        // checks it rather than trusting the calc that computes it.
        let expected_alpha = (1.0 + kappa * (1.0 - tr.sqrt())).powi(2);
        assert!(
            (ab.alpha - expected_alpha).abs() < 1e-12,
            "alpha does not follow from the PRSV coefficient: {} vs {expected_alpha}",
            ab.alpha
        );
        // And the reduced parameters must be built from *that* alpha, not a stale
        // one - a calc that took the coefficient but ignored it would pass the line
        // above and fail here.
        assert!(
            (ab.a_reduced - azoth_eos::OMEGA_A * expected_alpha * 0.25 / (tr * tr)).abs() < 1e-12
        );
    }
}

/// The coefficient rises with `omega`, for any `kappa1`.
///
/// The derivative of the acentric-only polynomial has a negative discriminant and
/// so is positive everywhere, and the `kappa1` term does not involve `omega` at all.
/// Monotonicity therefore holds whatever the fitted parameter is, which is what
/// makes it safe to assert without one.
#[test]
fn monotonic() {
    for kappa1 in [0.0, SOME_KAPPA1, -0.5] {
        for pair in [(0.0, 0.1), (0.1, 0.2), (0.2, 0.35), (0.35, 0.6), (0.6, 1.0)] {
            let (low, high) = pair;
            let at_low = prsv_kappa(low, 0.8, kappa1).unwrap().kappa;
            let at_high = prsv_kappa(high, 0.8, kappa1).unwrap().kappa;
            assert!(
                at_high > at_low,
                "kappa must rise with omega at kappa1 = {kappa1}: {low} -> {high} \
                 gave {at_low} -> {at_high}"
            );
        }
    }
}

/// The `kappa1` term vanishes at `Tr = 0.7`, whatever `kappa1` is.
///
/// `(1 + sqrt(Tr))` is never zero, so the whole term is zero exactly when
/// `0.7 - Tr` is. That makes the coefficient collapse to its acentric-only part at
/// the temperature where the acentric factor is defined - the anchor PRSV is built
/// around - and it is the sharpest cheap test of the temperature term.
#[test]
fn the_kappa1_term_vanishes_at_the_anchor() {
    let reference = prsv_kappa(0.152, 0.7, 0.0).unwrap().kappa;
    for kappa1 in [0.0, SOME_KAPPA1, -3.2, 100.0, 1e6] {
        let at_anchor = prsv_kappa(0.152, 0.7, kappa1).unwrap().kappa;
        assert_eq!(
            at_anchor.to_bits(),
            reference.to_bits(),
            "at Tr = 0.7 the coefficient should not depend on kappa1, but kappa1 = \
             {kappa1} gave {at_anchor} against {reference} for kappa1 = 0"
        );
    }
}

/// With `kappa1 = 0` the coefficient does not depend on `Tr` at all.
#[test]
fn a_zero_kappa1_is_temperature_independent() {
    let reference = prsv_kappa(0.152, 0.7, 0.0).unwrap().kappa;
    for tr in [0.3, 0.5, 0.7, 0.9, 1.0, 1.5] {
        let kappa = prsv_kappa(0.152, tr, 0.0).unwrap().kappa;
        assert_eq!(
            kappa.to_bits(),
            reference.to_bits(),
            "with kappa1 = 0 the coefficient should not move with Tr, but Tr = {tr} gave {kappa}"
        );
    }
}

#[test]
fn a_non_positive_reduced_temperature_is_an_error() {
    for tr in [0.0, -0.5] {
        let err = prsv_kappa(0.152, tr, SOME_KAPPA1).unwrap_err();
        assert!(matches!(err, AzothError::OutOfRange { .. }));
        assert_eq!(err.field(), Some("Tr"), "for Tr = {tr}");
    }
}

#[test]
fn kappa1_is_unconstrained_in_sign() {
    // Both signs occur in published fits, so neither may be refused. Asserting it
    // keeps a future bound from being added on the assumption that it is like a
    // coefficient that must be positive.
    for kappa1 in [-1.0, -0.03, 0.0, 0.05, 12.0] {
        assert!(prsv_kappa(0.152, 0.8, kappa1).is_ok(), "kappa1 = {kappa1}");
    }
}
