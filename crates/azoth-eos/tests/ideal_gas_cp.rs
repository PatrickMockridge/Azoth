//! Spec-driven tests for `eos.ideal_gas_cp`.
//!
//! The calc is a port of NeqSim's `Component.getCp0`, and it is dimensional: each
//! coefficient carries a power of temperature, which is what lets the `CPA`-`CPE`
//! columns of NeqSim's `COMP.csv` be read as they are stored.

use azoth_core::units::kelvins;
use azoth_core::{AzothError, CalcResult, WarningCode};
use azoth_eos::{ideal_gas_cp, spec_gen};
use azoth_test_support as common;

const CALC_ID: &str = "eos.ideal_gas_cp";

/// Methane's coefficients as NeqSim ships them, and the set most of these tests use.
const METHANE: (f64, f64, f64, f64, f64) =
    (37.978352, -0.07461815, 0.000301881, -2.83e-07, 9.070574e-11);

fn cp_of(
    cp_a: f64,
    cp_b: f64,
    cp_c: f64,
    cp_d: f64,
    cp_e: f64,
    t: f64,
) -> azoth_eos::IdealGasCpResult {
    ideal_gas_cp(cp_a, cp_b, cp_c, cp_d, cp_e, kelvins(t))
        .unwrap_or_else(|e| panic!("t = {t} should compute but failed: {e}"))
}

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::IdealGasCpResult {
    cp_of(
        common::input(case, "cp_a"),
        common::input(case, "cp_b"),
        common::input(case, "cp_c"),
        common::input(case, "cp_d"),
        common::input(case, "cp_e"),
        common::input(case, "T"),
    )
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
        if case.kind != "worked_example" && case.kind != "reference" {
            continue;
        }
        executed += 1;
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);

        common::assert_close(
            result.cp.value,
            common::expected(case, "cp"),
            case.tolerance,
            &format!("{context} (cp)"),
        );
        common::assert_consistent(&result, context);
        common::assert_warnings_agree_with_spec(
            spec,
            &result.warnings,
            |quantity| match quantity {
                "cp_a" => Some(common::input(case, "cp_a")),
                "cp_b" => Some(common::input(case, "cp_b")),
                "cp_c" => Some(common::input(case, "cp_c")),
                "cp_d" => Some(common::input(case, "cp_d")),
                "cp_e" => Some(common::input(case, "cp_e")),
                "T" => Some(common::input(case, "T")),
                "cp" => Some(result.cp.value),
                _ => None,
            },
            context,
        );
    }
    assert!(executed >= 3, "expected several cases, ran {executed}");
}

/// The polynomial is homogeneous of degree one in its coefficients.
///
/// `cp` is linear in all five by construction, so doubling every coefficient doubles
/// the answer and changing the sign of all five negates it. That is not a property of
/// any particular coefficient set - it holds for all of them - and it is what catches
/// an implementation which, say, multiplies a coefficient by the temperature twice or
/// drops a term, both of which pass at the values one worked example happens to use.
#[test]
fn the_polynomial_is_linear_in_its_coefficients() {
    let t = 500.0;
    let base = cp_of(METHANE.0, METHANE.1, METHANE.2, METHANE.3, METHANE.4, t);
    let doubled = cp_of(
        2.0 * METHANE.0,
        2.0 * METHANE.1,
        2.0 * METHANE.2,
        2.0 * METHANE.3,
        2.0 * METHANE.4,
        t,
    );
    let negated = cp_of(
        -METHANE.0, -METHANE.1, -METHANE.2, -METHANE.3, -METHANE.4, t,
    );

    assert!((doubled.cp.value - 2.0 * base.cp.value).abs() < 1e-9);
    assert!((negated.cp.value + base.cp.value).abs() < 1e-9);
}

/// At 1 K every power of the temperature is one.
///
/// The one temperature at which the five coefficients can be checked against the
/// answer by addition alone. An implementation using a reduced temperature - a
/// different scale, or a divisor - agrees elsewhere after rescaling its coefficients,
/// and fails here.
#[test]
fn at_one_kelvin_the_polynomial_is_the_sum_of_its_coefficients() {
    let (a, b, c, d, e) = (4.0, -20.0, 10.0, 1.0, 0.5);
    let result = cp_of(a, b, c, d, e, 1.0);
    assert!(
        (result.cp.value - (a + b + c + d + e)).abs() < 1e-12,
        "at T = 1 K the polynomial should be the sum of its coefficients, got {}",
        result.cp.value
    );
}

/// A polynomial that has turned over reports it rather than passing in silence.
///
/// The failure this catches is the one that matters: a negative heat capacity fed
/// into `integral Cp dT` produces an enthalpy wrong by an amount nobody can see.
#[test]
fn a_non_positive_heat_capacity_warns() {
    let result = cp_of(4.0, -20.0, 10.0, 1.0, 0.5, 1.0);
    assert!(result.cp.value < 0.0);
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code == WarningCode::OutOfValidRange),
        "a negative heat capacity should say so: {:?}",
        result.warnings
    );

    // And methane's own coefficients at 300 K, which are ordinary, are clean.
    let ordinary = cp_of(METHANE.0, METHANE.1, METHANE.2, METHANE.3, METHANE.4, 300.0);
    assert!(
        ordinary.is_clean(),
        "unexpected warnings: {:?}",
        ordinary.warnings
    );
}

#[test]
fn a_non_positive_temperature_is_refused() {
    for t in [0.0, -1.0] {
        let err = ideal_gas_cp(
            METHANE.0,
            METHANE.1,
            METHANE.2,
            METHANE.3,
            METHANE.4,
            kelvins(t),
        )
        .unwrap_err();
        assert!(matches!(err, AzothError::OutOfRange { .. }), "T = {t}");
        assert_eq!(err.field(), Some("T"), "T = {t}");
    }
}

/// All five coefficients matter, individually.
///
/// Five separate single-coefficient perturbations, each of which must move the answer.
/// An implementation that silently dropped one - the commonest way a polynomial goes
/// wrong, and `cp_e` is what a four-term implementation drops - would be
/// indistinguishable at most coefficient sets, because the term it drops is a small
/// correction.
#[test]
fn every_coefficient_changes_the_answer() {
    let t = 700.0;
    let base = cp_of(METHANE.0, METHANE.1, METHANE.2, METHANE.3, METHANE.4, t)
        .cp
        .value;
    let perturbations = [
        cp_of(
            METHANE.0 + 0.5,
            METHANE.1,
            METHANE.2,
            METHANE.3,
            METHANE.4,
            t,
        )
        .cp
        .value,
        cp_of(
            METHANE.0,
            METHANE.1 + 0.5,
            METHANE.2,
            METHANE.3,
            METHANE.4,
            t,
        )
        .cp
        .value,
        cp_of(
            METHANE.0,
            METHANE.1,
            METHANE.2 - 0.4,
            METHANE.3,
            METHANE.4,
            t,
        )
        .cp
        .value,
        cp_of(
            METHANE.0,
            METHANE.1,
            METHANE.2,
            METHANE.3 + 0.4,
            METHANE.4,
            t,
        )
        .cp
        .value,
        cp_of(
            METHANE.0,
            METHANE.1,
            METHANE.2,
            METHANE.3,
            METHANE.4 + 1.0e-11,
            t,
        )
        .cp
        .value,
    ];
    for (index, perturbed) in perturbations.into_iter().enumerate() {
        assert!(
            (perturbed - base).abs() > 1e-9,
            "coefficient {index} has no effect on the answer"
        );
    }
}
