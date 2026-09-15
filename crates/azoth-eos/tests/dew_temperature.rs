//! Spec-driven tests for the `eos.dew_temperature` model.

use azoth_core::units::{kelvins, pascals};
use azoth_core::{AzothError, CalcResult};
use azoth_eos::mixture::Mixture;
use azoth_eos::{databank, dew_temperature, model_gen, pt_flash};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.dew_temperature";

/// The methane/n-butane pair the spec's cases use, resolved through the databank so
/// the pair the sweeps run and the pair the cases run are the same fluid.
fn methane_butane() -> Mixture {
    databank::mixture_of(&["methane", "n-butane"], None)
        .expect("the pair resolves")
        .0
}

fn ternary() -> Mixture {
    databank::mixture_of(&["methane", "propane", "n-butane"], None)
        .expect("the trio resolves")
        .0
}

fn mixture_from_case(case: &azoth_core::spec::TestCase) -> Mixture {
    let names = case
        .list("components")
        .expect("the case declares components");
    databank::mixture_of(names, None)
        .expect("the case's components resolve")
        .0
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty());

    for case in spec.cases {
        let mixture = mixture_from_case(case);
        let result = dew_temperature(
            &mixture,
            pascals(common::input(case, "P")),
            case.vector("y").expect("y"),
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

/// The dew point is where a flash's vapour fraction reaches one.
///
/// The strongest check available without external data, and the two models are
/// computed by entirely different procedures: one searches for a state where a phase
/// split exists, the other for the temperature at which it stops existing.
#[test]
fn the_flash_agrees_at_the_returned_temperature() {
    for (mixture, p, held) in [
        (methane_butane(), 1567347.3649394324, vec![0.8, 0.2]),
        (ternary(), 1036545.9819627925, vec![0.5, 0.3, 0.2]),
    ] {
        let boundary = dew_temperature(&mixture, pascals(p), &held).unwrap();
        let flash = pt_flash(
            &mixture,
            kelvins(boundary.temperature.value),
            pascals(p),
            &held,
        )
        .unwrap();

        let beta = flash
            .beta
            .unwrap_or_else(|| panic!("P={p}: the flash found no split at the dew point"));
        assert!(
            (beta - 1.0).abs() < 1e-8,
            "P={p}: beta is {beta} at the dew temperature, not one"
        );
        for (i, &yi) in held.iter().enumerate() {
            assert!(
                (flash.y[i] - yi).abs() < 1e-8,
                "P={p}: the flash's vapour {i} is {} not {yi}",
                flash.y[i]
            );
            // The incipient (liquid) composition is a successive-substitution estimate
            // from `x_i = y_i / (K_i S)`, accurate to ~1e-9 rather than to the flash's
            // exact split, for the same reason the bubble estimate is ~1e-8.
            assert!(
                (flash.x[i] - boundary.incipient[i]).abs() < 1e-7,
                "P={p}: the flash's liquid {i} is {} not {}",
                flash.x[i],
                boundary.incipient[i]
            );
        }
    }
}

/// The guard must not fire on a genuine dew point, however slow the iteration.
#[test]
fn a_genuine_dew_point_is_never_refused() {
    for (p, held) in [
        (1567347.3649394324, vec![0.8, 0.2]),
        (1036545.9819627925, vec![0.5, 0.3, 0.2]),
    ] {
        let mixture = if held.len() == 3 {
            ternary()
        } else {
            methane_butane()
        };
        let result = dew_temperature(&mixture, pascals(p), &held)
            .unwrap_or_else(|e| panic!("P={p} {held:?} should have a dew point: {e}"));
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
    let propane = databank::mixture_of(&["propane"], None)
        .expect("propane resolves")
        .0;
    let err = dew_temperature(&propane, pascals(1.0e6), &[1.0]).unwrap_err();
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
        vec![0.8, 0.1],
        vec![0.8, 0.3],
        vec![1.6, -0.6],
        vec![0.8, 0.2, 0.0],
    ] {
        assert!(
            dew_temperature(&mixture, pascals(1.5e6), &x).is_err(),
            "x = {x:?}"
        );
    }
}

#[test]
fn a_non_positive_pressure_is_refused() {
    let mixture = methane_butane();
    for p in [0.0, -1.0] {
        assert!(
            dew_temperature(&mixture, pascals(p), &[0.8, 0.2]).is_err(),
            "P={p}"
        );
    }
}

/// The three compositions the iteration produces are consistent with each other.
#[test]
fn the_result_is_self_consistent() {
    let result =
        dew_temperature(&methane_butane(), pascals(1567347.3649394324), &[0.8, 0.2]).unwrap();
    let sum_x: f64 = result.incipient.iter().sum();
    assert!(
        (sum_x - 1.0).abs() < 1e-12,
        "the incipient phase sums to {sum_x}"
    );
    for i in 0..2 {
        let held = [0.8, 0.2][i];
        // `K_i = y_i / x_i`, so `x_i = y_i / K_i` and `x_i * K_i / y_i = 1` up to the
        // normalisation residual.
        let ratio = result.incipient[i] * result.k[i] / held;
        assert!(
            (ratio - 1.0).abs() < 1e-11,
            "K[{i}] is {} but y/x is {}",
            result.k[i],
            held / result.incipient[i]
        );
    }
    assert!(
        result.is_clean(),
        "unexpected warnings: {:?}",
        result.warnings
    );
}
