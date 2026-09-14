//! Spec-driven tests for the `eos.bubble_pressure` model.

use azoth_core::units::{kelvins, pascals};
use azoth_core::{AzothError, CalcResult};
use azoth_eos::mixture::Mixture;
use azoth_eos::{bubble_pressure, databank, model_gen, pt_flash};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.bubble_pressure";

/// The methane/n-butane pair the spec's cases use, resolved through the databank so
/// the pair the sweeps run and the pair the cases run are the same fluid.
fn methane_butane() -> Mixture {
    databank::mixture_of(&["methane", "n-butane"])
        .expect("the pair resolves")
        .0
}

fn ternary() -> Mixture {
    databank::mixture_of(&["methane", "propane", "n-butane"])
        .expect("the trio resolves")
        .0
}

fn mixture_from_case(case: &azoth_core::spec::TestCase) -> Mixture {
    let names = case
        .list("components")
        .expect("the case declares components");
    databank::mixture_of(names)
        .expect("the case's components resolve")
        .0
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty());

    for case in spec.cases {
        let mixture = mixture_from_case(case);
        let result = bubble_pressure(
            &mixture,
            kelvins(common::input(case, "T")),
            case.vector("x").expect("x"),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        let context = &format!("{}::{}", spec.id, case.id);
        common::assert_close(
            result.pressure.value,
            common::expected(case, "pressure"),
            case.tolerance,
            &format!("{context} (pressure)"),
        );
        for name in ["z_liquid", "z_vapour", "min_t_over_tc"] {
            common::assert_close(
                match name {
                    "z_liquid" => result.z_liquid,
                    "z_vapour" => result.z_vapour,
                    _ => result.min_t_over_tc,
                },
                common::expected(case, name),
                case.tolerance,
                &format!("{context} ({name})"),
            );
        }
        assert_eq!(
            result.iterations as f64,
            common::expected(case, "iterations")
        );
        assert!(
            result.residual <= spec.algorithm.expect("a procedure").tolerance,
            "{context}: residual {:e} exceeds the declared tolerance",
            result.residual
        );

        for name in ["incipient", "k"] {
            let actual: &[f64] = if name == "incipient" {
                &result.incipient
            } else {
                &result.k
            };
            let expected = case.expected_vector(name).expect("declared");
            assert_eq!(actual.len(), expected.len(), "{name}: length");
            for (i, (a, e)) in actual.iter().zip(expected).enumerate() {
                common::assert_close(*a, *e, case.tolerance, &format!("{context} ({name}[{i}])"));
            }
        }
    }
}

/// The bubble point is where a flash's vapour fraction vanishes.
///
/// The strongest check available without external data, and the two models are
/// computed by entirely different procedures: one searches for a state where a
/// phase split exists, the other for the pressure at which it stops existing. A
/// flash of the *held* composition at the returned pressure must therefore come
/// back with `beta` at zero and the same two compositions.
#[test]
fn the_flash_agrees_at_the_returned_pressure() {
    for (mixture, t, held) in [
        (methane_butane(), 300.0, vec![0.2, 0.8]),
        (methane_butane(), 330.0, vec![0.2, 0.8]),
        (ternary(), 320.0, vec![0.2, 0.3, 0.5]),
    ] {
        let boundary = bubble_pressure(&mixture, kelvins(t), &held).unwrap();
        let flash = pt_flash(
            &mixture,
            kelvins(t),
            pascals(boundary.pressure.value),
            &held,
        )
        .unwrap();

        let beta = flash
            .beta
            .unwrap_or_else(|| panic!("T={t}: the flash found no split at the bubble point"));
        assert!(
            beta.abs() < 1e-9,
            "T={t}: beta is {beta} at the bubble pressure, not zero"
        );
        for (i, &xi) in held.iter().enumerate() {
            assert!(
                (flash.x[i] - xi).abs() < 1e-9,
                "T={t}: the flash's liquid {i} is {} not {xi}",
                flash.x[i]
            );
            assert!(
                (flash.y[i] - boundary.incipient[i]).abs() < 1e-9,
                "T={t}: the flash's vapour {i} is {} not {}",
                flash.y[i],
                boundary.incipient[i]
            );
        }
    }
}

/// A mixture above its critical condition has no bubble point, and says so.
///
/// The trap: `sum_i x_i K_i = 1` is satisfied by `K_i = 1` at *every* pressure, and
/// the residual `S - 1` **cannot** detect it because it is a weighted sum whose
/// terms cancel - measured at 1.5e-13 while the individual K-values were 1e-7 from
/// one. So the guard is on the K-values, and this is the check that it fires.
#[test]
fn a_mixture_with_no_bubble_point_is_refused() {
    for (t, held) in [
        (350.0, vec![0.5, 0.3, 0.2]),
        (340.0, vec![0.5, 0.3, 0.2]),
        (330.0, vec![0.5, 0.5]),
    ] {
        let mixture = if held.len() == 3 {
            ternary()
        } else {
            methane_butane()
        };
        let err = bubble_pressure(&mixture, kelvins(t), &held)
            .expect_err(&format!("T={t} {held:?} has no bubble point"));
        assert!(matches!(err, AzothError::OutOfRange { .. }), "T={t}");
        assert_eq!(err.field(), Some("min_t_over_tc"), "T={t}");
    }
}

/// The guard must not fire on a genuine bubble point, however slow the iteration.
///
/// The threshold's whole justification is the gap between these states and the ones
/// above: genuine boundaries converge with `max |ln K|` at 1.10 or more against the
/// degenerate cluster's `9.1e-04` at most.
#[test]
fn a_genuine_bubble_point_is_never_refused() {
    for (t, held) in [
        (300.0, vec![0.2, 0.8]),
        (280.0, vec![0.2, 0.8]),
        (320.0, vec![0.2, 0.3, 0.5]),
        (300.0, vec![0.2, 0.3, 0.5]),
    ] {
        let mixture = if held.len() == 3 {
            ternary()
        } else {
            methane_butane()
        };
        let result = bubble_pressure(&mixture, kelvins(t), &held)
            .unwrap_or_else(|e| panic!("T={t} {held:?} should have a bubble point: {e}"));
        let max_ln_k = result
            .k
            .iter()
            .map(|value| value.ln().abs())
            .fold(0.0_f64, f64::max);
        assert!(
            max_ln_k > 1e-3,
            "T={t} {held:?}: max |ln K| is {max_ln_k:e}, which the guard would reject - \
             the threshold and the measurement behind it disagree"
        );
    }
}

#[test]
fn a_single_component_is_refused_and_points_at_the_right_calc() {
    let propane = databank::mixture_of(&["propane"])
        .expect("propane resolves")
        .0;
    let err = bubble_pressure(&propane, kelvins(300.0), &[1.0]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }));
    assert!(
        err.to_string().contains("pure_saturation"),
        "the refusal should name the calc that does answer it: {err}"
    );
}

#[test]
fn a_malformed_composition_is_refused() {
    let mixture = methane_butane();
    for x in [
        vec![0.2, 0.9],
        vec![0.2, 0.7],
        vec![0.6, -0.6],
        vec![0.2, 0.8, 0.0],
    ] {
        assert!(
            bubble_pressure(&mixture, kelvins(300.0), &x).is_err(),
            "x = {x:?}"
        );
    }
}

#[test]
fn a_non_positive_temperature_is_refused() {
    let mixture = methane_butane();
    for t in [0.0, -1.0] {
        assert!(
            bubble_pressure(&mixture, kelvins(t), &[0.2, 0.8]).is_err(),
            "T={t}"
        );
    }
}

/// The three compositions the iteration produces are consistent with each other.
///
/// `sum_i x_i K_i = 1` is the equation being solved, so it holds at the answer by
/// construction - but `K_i = y_i / x_i` and `sum_i y_i = 1` are properties the
/// result has to satisfy separately, and a K-vector that was right while the
/// composition was not would pass a point-value case.
#[test]
fn the_result_is_self_consistent() {
    let result = bubble_pressure(&methane_butane(), kelvins(300.0), &[0.2, 0.8]).unwrap();
    let sum_y: f64 = result.incipient.iter().sum();
    assert!(
        (sum_y - 1.0).abs() < 1e-12,
        "the incipient phase sums to {sum_y}"
    );
    // `K_i = y_i / x_i` up to the normalisation, which is enforced to the
    // residual: the iteration sets `y_i = x_i K_i / s` with `s - 1` at the
    // tolerance, so the identity holds relatively rather than absolutely. An
    // absolute bound here would be a bound on how large `K` is allowed to be.
    for i in 0..2 {
        let held = [0.2, 0.8][i];
        let ratio = result.incipient[i] / held;
        assert!(
            (ratio / result.k[i] - 1.0).abs() < 1e-11,
            "K[{i}] is {} but y/x is {ratio}",
            result.k[i]
        );
    }
    assert!(
        result.is_clean(),
        "unexpected warnings: {:?}",
        result.warnings
    );
}
