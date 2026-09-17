//! Spec-driven tests for the `eos.rachford_rice` model.

use azoth_core::CalcResult;
use azoth_eos::rachford_rice::rachford_rice;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.rachford_rice";

/// A K-value below this is an ion and takes no part in the split.
const ION_THRESHOLD: f64 = 1.0e-30;

/// The bracket `g` has its root in, or `None` when its K-values do not straddle one.
///
/// Written here rather than taken from `azoth_eos::flash_iteration`: the point of the
/// comparison below is that two procedures written apart find the same root, and a
/// shared helper would hide the agreement it is meant to show.
fn bracket(k: &[f64]) -> Option<(f64, f64)> {
    let (mut lo, mut hi) = (f64::NEG_INFINITY, f64::INFINITY);
    for &value in k {
        if value < ION_THRESHOLD {
            continue;
        }
        if value > 1.0 {
            lo = lo.max(1.0 / (1.0 - value));
        } else if value < 1.0 {
            hi = hi.min(1.0 / (1.0 - value));
        } else {
            return None;
        }
    }
    (lo.is_finite() && hi.is_finite()).then_some((lo, hi))
}

/// `g` bisected on its bracket until the interval closes to `1e-14`.
fn bisect(z: &[f64], k: &[f64]) -> Option<f64> {
    let (mut lo, mut hi) = bracket(k)?;
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if hi - lo <= 1.0e-14 {
            break;
        }
        let g: f64 = z
            .iter()
            .zip(k)
            .filter(|&(_, &ki)| ki >= ION_THRESHOLD)
            .map(|(&zi, &ki)| zi * (ki - 1.0) / (1.0 + mid * (ki - 1.0)))
            .sum();
        if g > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

#[test]
fn every_case_in_the_spec() {
    let spec = azoth_eos::model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let context = &format!("{}::{}", spec.id, case.id);
        let result = rachford_rice(
            case.vector("z").expect("the case declares z"),
            case.vector("K").expect("the case declares K"),
        )
        .unwrap_or_else(|e| panic!("{context} should compute but failed: {e}"));
        common::assert_close(
            result.beta,
            common::expected(case, "beta"),
            case.tolerance,
            &format!("{context} (beta)"),
        );
        common::assert_consistent(&result, context);
    }
}

/// The two procedures agree wherever both are defined.
///
/// The model runs Nielsen & Lia's reformulation; `eos.pt_flash` bisects the same equation
/// on the same bracket. They are separate code written from separate sources, so a
/// bracket or a sign error in either shows up here as a disagreement about which root the
/// equation has - and the tolerance is the model's `1e-10` stopping rule, not a tighter
/// one it does not claim.
#[test]
fn the_root_agrees_with_the_bisection_wherever_both_are_defined() {
    let spec = azoth_eos::model_gen::model(MODEL_ID).expect("the model");
    let mut compared = 0;
    for (z, k) in spec
        .cases
        .iter()
        .map(|case| {
            (
                case.vector("z").unwrap().to_vec(),
                case.vector("K").unwrap().to_vec(),
            )
        })
        .chain([
            // The two states NeqSim answers with a clamp and the bracket is defined for.
            (vec![0.1, 0.9], vec![5.799172708809655, 0.14913889826410245]),
            (vec![0.5, 0.5], vec![1.5, 0.9]),
        ])
    {
        let Some(expected) = bisect(&z, &k).filter(|_| bracket(&k).is_some()) else {
            continue;
        };
        let got = rachford_rice(&z, &k)
            .expect("a bracketed feed computes")
            .beta;
        assert!(
            (got - expected).abs() <= 1.0e-10,
            "z={z:?}, K={k:?}: the model gives {got} and the bisection {expected}"
        );
        compared += 1;
    }
    assert!(compared >= 5, "only {compared} states were compared");
}

/// Every state `eos.pt_flash` settles on is a state both procedures agree at.
///
/// The flash's own cases record the K-values it converges to, so this is the model
/// checked against the solver the flash actually runs, at the iterates a real flash
/// produces - which is the comparison the model's registration is for.
#[test]
fn the_root_matches_the_bisection_on_every_flash_case() {
    let flash = azoth_eos::model_gen::model("eos.pt_flash").expect("the flash");
    let mut compared = 0;
    for case in flash.cases {
        let context = &format!("eos.pt_flash::{}", case.id);
        let z = case.vector("z").expect("the case declares z");
        let k = case
            .expected_vector("k")
            .unwrap_or_else(|| panic!("{context} records the K-values it converged to"));
        let Some(expected) = bisect(z, k) else {
            continue;
        };
        let got = rachford_rice(z, k)
            .unwrap_or_else(|e| panic!("{context} should compute but failed: {e}"))
            .beta;
        assert!(
            (got - expected).abs() <= 1.0e-10,
            "{context}: the model gives {got} and the flash's bisection {expected}"
        );
        common::assert_close(
            got,
            common::expected(case, "beta"),
            case.tolerance,
            &format!("{context} (beta)"),
        );
        compared += 1;
    }
    assert!(compared >= 3, "only {compared} flash cases were compared");
}

/// A feed with no root is NeqSim's clamp, not an error.
///
/// The two cases are the ends: every K below one is a liquid, every K above one is a
/// vapour, and the equation has nothing to solve in either.
#[test]
fn a_feed_that_cannot_split_comes_back_clamped() {
    for (z, k, expected) in [
        (vec![0.5, 0.5], vec![0.2, 0.3], 1.0e-12),
        (vec![0.5, 0.5], vec![3.0, 5.0], 1.0 - 1.0e-12),
    ] {
        let got = rachford_rice(&z, &k).expect("a single-phase feed is not an error");
        assert!(
            (got.beta - expected).abs() < 1.0e-15,
            "z={z:?}, K={k:?}: expected the clamp {expected}, got {}",
            got.beta
        );
        assert!(
            bracket(&k).is_none(),
            "z={z:?}, K={k:?}: the clamp is the answer only where no bracket exists"
        );
    }
}

/// A root outside `[0, 1]` is reported rather than clamped away.
///
/// This is the difference from NeqSim, and it is the difference the whole model turns
/// on: `RachfordRice.calcBeta` answers `g(0) < 0` and `g(1) > 0` with
/// `phaseFractionMinimumLimit`, so the vapour fraction below comes back as `1e-12` from
/// both of its solvers - measured, `validation/neqsim/RachfordRiceProbe.java`. A negative
/// vapour fraction is the negative flash, which `eos.pt_flash` reports and explains, so
/// restoring either of NeqSim's two single-phase returns fails this test.
#[test]
fn a_root_outside_the_unit_interval_is_reported_not_clamped() {
    // Subcooled: the root is below zero.
    let z = [0.1, 0.9];
    let k = [5.799172708809655, 0.14913889826410245];
    let below = rachford_rice(&z, &k)
        .expect("a bracketed feed computes")
        .beta;
    assert!(
        below < -0.06 && below > -0.08,
        "the negative flash should be near -0.07, got {below}"
    );
    assert!(
        (below - bisect(&z, &k).expect("the bracket exists")).abs() <= 1.0e-10,
        "the reported root is not the equation's"
    );

    // Superheated: the root is above one. `K = [1.5, 0.9]` puts it at exactly four -
    // `0.5(K1-1)/(1+3b) + 0.5(K2-1)/(1-0.1b) = 0` rearranges to `b = 4`.
    let above = rachford_rice(&[0.5, 0.5], &[1.5, 0.9])
        .expect("a bracketed feed computes")
        .beta;
    assert!(
        (above - 4.0).abs() <= 1.0e-10,
        "the superheated root should be four, got {above}"
    );
}

/// An ion stays in the liquid and does not enter the equation.
///
/// NeqSim's threshold and NeqSim's reason: a component with `K < 1e-30` is skipped by
/// the sum, by the bracket and by the rescaling alike. The value is
/// `validation/neqsim/RachfordRiceProbe.java`'s, and it is the same number the
/// methane/butane pair alone would give once `z` is read over the participating
/// components.
#[test]
fn an_ion_is_skipped() {
    let with_ion = rachford_rice(
        &[0.1, 0.45, 0.45],
        &[1.0e-40, 7.304244305324782, 0.33749596785762953],
    )
    .expect("an ion does not stop the split")
    .beta;
    assert!(
        (with_ion - 0.6754007405823642).abs() <= 1.0e-12,
        "the ion case gives {with_ion}"
    );

    // Every component an ion: there is no split to find, and the model says so with
    // the clamp rather than by dividing by zero.
    let none = rachford_rice(&[0.5, 0.5], &[1.0e-40, 1.0e-50])
        .expect("a feed of ions is not an error")
        .beta;
    assert!((none - 1.0e-12).abs() < 1.0e-15, "got {none}");
}

#[test]
fn a_k_value_that_is_not_positive_is_refused() {
    for bad in [0.0, -0.5] {
        let err = rachford_rice(&[0.5, 0.5], &[bad, 2.0]).expect_err("K must be positive");
        assert!(
            matches!(err, azoth_core::AzothError::OutOfRange { .. }),
            "K={bad} gave {err:?}"
        );
    }
}

#[test]
fn a_feed_and_a_k_vector_of_different_lengths_are_refused() {
    let err = rachford_rice(&[0.5, 0.5], &[2.0]).expect_err("the lengths must agree");
    assert!(matches!(err, azoth_core::AzothError::InvalidInput { .. }));
    let empty = rachford_rice(&[], &[]).expect_err("a feed of nothing has no split");
    assert!(matches!(empty, azoth_core::AzothError::InvalidInput { .. }));
}

/// The result is the model's, and the model's id is the one the spec registers.
#[test]
fn the_result_carries_the_models_own_id() {
    let result = rachford_rice(&[0.6, 0.4], &[7.304244305324782, 0.33749596785762953]).unwrap();
    assert_eq!(
        <azoth_eos::RachfordRiceResult as CalcResult>::CALC_ID,
        MODEL_ID
    );
    assert!(
        result.warnings.is_empty(),
        "an ordinary state warns about nothing"
    );
}
