//! Spec-driven tests for `eos.vdw1f_mix_binary`.

use azoth_core::AzothError;
use azoth_eos::spec_gen;
use azoth_eos::vdw1f_mix_binary;
use azoth_test_support as common;

const CALC_ID: &str = "eos.vdw1f_mix_binary";

/// A deliberately arbitrary `k12`.
///
/// This library ships no fitted values for it - a table of binary parameters is the
/// databank it does not have - so tests pick a round number and say so rather than
/// borrowing one that looks like data.
const SOME_K12: f64 = 0.05;

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::Vdw1fMixBinaryResult {
    vdw1f_mix_binary(
        common::input(case, "z1"),
        common::input(case, "a1"),
        common::input(case, "a2"),
        common::input(case, "b1"),
        common::input(case, "b2"),
        common::input(case, "k12"),
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
                for field in ["a_mix", "b_mix"] {
                    let actual = if field == "a_mix" {
                        result.a_mix
                    } else {
                        result.b_mix
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
                        "z1" => Some(common::input(case, "z1")),
                        "a1" => Some(common::input(case, "a1")),
                        "a2" => Some(common::input(case, "a2")),
                        "b1" => Some(common::input(case, "b1")),
                        "b2" => Some(common::input(case, "b2")),
                        "k12" => Some(common::input(case, "k12")),
                        "a_mix" => Some(result.a_mix),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
            }
            "property" => match case.property {
                Some("monotonic") => {}
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

/// A pure component reproduces its own parameters, for arbitrary inputs.
///
/// The check that a mixing rule is a mixing rule, and the sharpest cheap one: at
/// `z1 = 1` every term carrying `z2` vanishes and the mixture parameters must be
/// component 1's own; at `z1 = 0`, component 2's.
///
/// Swept over parameters rather than fixed at the spec's case, because it is the
/// *structure* being tested rather than a value. An implementation with a
/// transposed weight, or with `(1 - z1)` written as `z1` in the cross term, passes
/// a mid-range composition and fails here.
#[test]
fn a_pure_component_reproduces_its_own_parameters() {
    for (a1, a2, b1, b2, k12) in [
        (
            0.20206500174625697,
            0.08448417260831159,
            0.02431127309496514,
            0.025932024634629486,
            SOME_K12,
        ),
        (0.115, 0.31, 0.018, 0.041, 0.4),
        (1.0, 2.0, 3.0, 4.0, 0.0),
    ] {
        let pure1 = vdw1f_mix_binary(1.0, a1, a2, b1, b2, k12).unwrap();
        assert_eq!(pure1.a_mix.to_bits(), a1.to_bits(), "z1 = 1 should give a1");
        assert_eq!(pure1.b_mix.to_bits(), b1.to_bits(), "z1 = 1 should give b1");

        let pure2 = vdw1f_mix_binary(0.0, a1, a2, b1, b2, k12).unwrap();
        assert_eq!(pure2.a_mix.to_bits(), a2.to_bits(), "z1 = 0 should give a2");
        assert_eq!(pure2.b_mix.to_bits(), b2.to_bits(), "z1 = 0 should give b2");
    }
}

/// `b_mix` is the mole-fraction-weighted mean, exactly.
///
/// Worth asserting bit-for-bit rather than approximately, because `b_mix` is linear
/// in `z1` and has no square root in it: `z1*b1 + (1 - z1)*b2` is either computed
/// that way or it is not, and a `k12` that leaked into the covolume would show up
/// here and nowhere else. vdW1f has no `l12`.
#[test]
fn b_mix_is_linear_in_composition() {
    let (b1, b2) = (0.02431127309496514, 0.025932024634629486);
    for z1 in [0.0, 0.25, 0.5, 0.6, 0.75, 1.0] {
        let r = vdw1f_mix_binary(z1, 0.2, 0.1, b1, b2, SOME_K12).unwrap();
        assert_eq!(r.b_mix.to_bits(), (z1 * b1 + (1.0 - z1) * b2).to_bits());
    }
}

/// Each mixture parameter moves monotonically with `z1`, in the direction its own
/// pure parameters order.
///
/// The direction is the substance of the test rather than an incidental. A weighted
/// average moves toward whichever component's fraction is rising, so `a_mix` rises
/// with `z1` when `a1 > a2` and falls when it does not - and the two need not agree,
/// because nothing relates one component's attraction to its covolume. The propane-
/// like and methane-like pair below has `a1 > a2` but `b1 < b2`, so the two
/// parameters move in *opposite* directions. An implementation with a transposed
/// weight would move both the same way and fail.
#[test]
fn each_mixture_parameter_moves_with_its_own_ordering() {
    let cases = [
        // (a1, a2, b1, b2) - a1 > a2 but b1 < b2
        (
            0.20206500174625697,
            0.08448417260831159,
            0.02431127309496514,
            0.025932024634629486,
        ),
        // both larger in component 1
        (0.31, 0.115, 0.041, 0.018),
        // both smaller in component 1
        (0.05, 0.4, 0.01, 0.09),
    ];
    for (a1, a2, b1, b2) in cases {
        let mut last: Option<(f64, f64)> = None;
        for step in 0..=10 {
            let z1 = f64::from(step) / 10.0;
            let r = vdw1f_mix_binary(z1, a1, a2, b1, b2, SOME_K12).unwrap();
            if let Some((prev_a, prev_b)) = last {
                if a1 > a2 {
                    assert!(r.a_mix > prev_a, "a_mix must rise with z1, at {z1}");
                } else {
                    assert!(r.a_mix < prev_a, "a_mix must fall with z1, at {z1}");
                }
                if b1 > b2 {
                    assert!(r.b_mix > prev_b, "b_mix must rise with z1, at {z1}");
                } else {
                    assert!(r.b_mix < prev_b, "b_mix must fall with z1, at {z1}");
                }
            }
            last = Some((r.a_mix, r.b_mix));
        }
    }
}

/// A `k12` outside `[0, 2]` can drive `a_mix` negative, and that is refused.
///
/// The bound is on `a_mix` rather than on `k12` because a `k12` outside the band is
/// not wrong in itself - most compositions still give a physical `a_mix` - and
/// refusing the parameter would reject calls that are fine.
#[test]
fn an_unphysical_mixture_attraction_is_refused() {
    // a1 = a2 = 1 with k12 = 3 gives a12 = -2, and the quadratic form is negative
    // at the midpoint: 0.25 - 1 + 0.25 = -0.5.
    let err = vdw1f_mix_binary(0.5, 1.0, 1.0, 0.1, 0.1, 3.0).unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("a_mix"));

    // But the same k12 at a composition that stays physical is fine - which is the
    // reason the bound is not on k12.
    assert!(vdw1f_mix_binary(0.99, 1.0, 1.0, 0.1, 0.1, 3.0).is_ok());
}

#[test]
fn a_mole_fraction_outside_zero_to_one_is_an_error() {
    for z1 in [-0.01, 1.01, 2.0] {
        let err = vdw1f_mix_binary(z1, 0.2, 0.1, 0.02, 0.03, SOME_K12).unwrap_err();
        assert!(matches!(err, AzothError::OutOfRange { .. }));
        assert_eq!(err.field(), Some("z1"), "for z1 = {z1}");
    }
}
