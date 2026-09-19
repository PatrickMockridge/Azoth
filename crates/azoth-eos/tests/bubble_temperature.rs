//! Spec-driven tests for the `eos.bubble_temperature` model.

use azoth_core::units::{kelvins, pascals};
use azoth_core::{AzothError, CalcResult};
use azoth_eos::Cubic;
use azoth_eos::mixture::Mixture;
use azoth_eos::{bubble_temperature, databank, model_gen, pt_flash};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.bubble_temperature";

/// The methane/n-butane pair the spec's cases use, resolved through the databank so
/// the pair the sweeps run and the pair the cases run are the same fluid.
fn methane_butane() -> Mixture {
    databank::mixture_of(&["methane", "n-butane"], Cubic::Pr, None)
        .expect("the pair resolves")
        .0
}

fn ternary() -> Mixture {
    databank::mixture_of(&["methane", "propane", "n-butane"], Cubic::Pr, None)
        .expect("the trio resolves")
        .0
}

fn mixture_from_case(case: &azoth_core::spec::TestCase) -> Mixture {
    let names = case
        .list("components")
        .expect("the case declares components");
    databank::mixture_of(names, Cubic::Pr, None)
        .expect("the case's components resolve")
        .0
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty());

    for case in spec.cases {
        let mixture = mixture_from_case(case);
        let result = bubble_temperature(
            &mixture,
            pascals(common::input(case, "P")),
            case.vector("x").expect("x"),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        let context = &format!("{}::{}", spec.id, case.id);
        common::assert_close(
            result.temperature.value,
            common::expected(case, "temperature"),
            case.tolerance,
            &format!("{context} (T)"),
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
/// computed by entirely different procedures: one searches for a state where a phase
/// split exists, the other for the temperature at which it stops existing.
#[test]
fn the_flash_agrees_at_the_returned_temperature() {
    for (mixture, p, held) in [
        (methane_butane(), 3950960.4937437344, vec![0.2, 0.8]),
        (ternary(), 4664798.1363369655, vec![0.2, 0.3, 0.5]),
    ] {
        let boundary = bubble_temperature(&mixture, pascals(p), &held).unwrap();
        let flash = pt_flash(
            &mixture,
            kelvins(boundary.temperature.value),
            pascals(p),
            &held,
        )
        .unwrap();

        let beta = flash
            .beta
            .unwrap_or_else(|| panic!("P={p}: the flash found no split at the bubble point"));
        assert!(
            beta.abs() < 1e-8,
            "P={p}: beta is {beta} at the bubble temperature, not zero"
        );
        for (i, &xi) in held.iter().enumerate() {
            assert!(
                (flash.x[i] - xi).abs() < 1e-8,
                "P={p}: the flash's liquid {i} is {} not {xi}",
                flash.x[i]
            );
            // The incipient composition is a successive-substitution estimate from
            // `y_i = x_i K_i / S`, accurate to ~1e-8 rather than to the flash's exact
            // split: the Newton step moves temperature so `S` is small before the
            // composition fully settles, where the pressure search's fixed-point
            // update settles both together.
            assert!(
                (flash.y[i] - boundary.incipient[i]).abs() < 1e-7,
                "P={p}: the flash's vapour {i} is {} not {}",
                flash.y[i],
                boundary.incipient[i]
            );
        }
    }
}

/// A mixture above its cricondenbar has no bubble point, and says so.
///
/// The trap is the same one the pressure search guards against: `sum_i x_i K_i = 1`
/// is satisfied by `K_i = 1` at every temperature, and the residual cannot detect it
/// because its terms cancel. So the guard is on the K-values.
#[test]
fn a_mixture_with_no_bubble_point_is_refused() {
    for p in [6.0e6, 8.0e6] {
        let err = bubble_temperature(&methane_butane(), pascals(p), &[0.2, 0.8])
            .expect_err(&format!("P={p} has no bubble point"));
        assert!(matches!(err, AzothError::OutOfRange { .. }), "P={p}");
        assert_eq!(err.field(), Some("min_t_over_tc"), "P={p}");
    }
}

/// The guard must not fire on a genuine bubble point, however slow the iteration.
#[test]
fn a_genuine_bubble_point_is_never_refused() {
    for (p, held) in [
        (3950960.4937437344, vec![0.2, 0.8]),
        (4664798.1363369655, vec![0.2, 0.3, 0.5]),
    ] {
        let mixture = if held.len() == 3 {
            ternary()
        } else {
            methane_butane()
        };
        let result = bubble_temperature(&mixture, pascals(p), &held)
            .unwrap_or_else(|e| panic!("P={p} {held:?} should have a bubble point: {e}"));
        let max_ln_k = result
            .k
            .iter()
            .map(|value| value.ln().abs())
            .fold(0.0_f64, f64::max);
        assert!(
            max_ln_k > 1e-3,
            "P={p} {held:?}: max |ln K| is {max_ln_k:e}, which the guard would reject - \
             the threshold and the measurement behind it disagree"
        );
    }
}

#[test]
fn a_single_component_is_refused_and_points_at_the_right_calc() {
    let propane = databank::mixture_of(&["propane"], Cubic::Pr, None)
        .expect("propane resolves")
        .0;
    let err = bubble_temperature(&propane, pascals(1.0e6), &[1.0]).unwrap_err();
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
            bubble_temperature(&mixture, pascals(3.0e6), &x).is_err(),
            "x = {x:?}"
        );
    }
}

#[test]
fn a_non_positive_pressure_is_refused() {
    let mixture = methane_butane();
    for p in [0.0, -1.0] {
        assert!(
            bubble_temperature(&mixture, pascals(p), &[0.2, 0.8]).is_err(),
            "P={p}"
        );
    }
}

/// The three compositions the iteration produces are consistent with each other.
#[test]
fn the_result_is_self_consistent() {
    let result =
        bubble_temperature(&methane_butane(), pascals(3950960.4937437344), &[0.2, 0.8]).unwrap();
    let sum_y: f64 = result.incipient.iter().sum();
    assert!(
        (sum_y - 1.0).abs() < 1e-12,
        "the incipient phase sums to {sum_y}"
    );
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
