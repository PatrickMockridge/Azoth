//! Spec-driven tests for `eos.ideal_gas_cp`.

use azoth_core::units::kelvins;
use azoth_core::{AzothError, CalcResult, WarningCode};
use azoth_eos::{ideal_gas_cp, spec_gen};
use azoth_test_support as common;

const CALC_ID: &str = "eos.ideal_gas_cp";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::IdealGasCpResult {
    ideal_gas_cp(
        common::input(case, "a"),
        common::input(case, "b"),
        common::input(case, "c"),
        common::input(case, "d"),
        kelvins(common::input(case, "T")),
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
        if case.kind != "worked_example" && case.kind != "reference" {
            continue;
        }
        executed += 1;
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);

        common::assert_close(
            result.cp_over_r,
            common::expected(case, "cp_over_r"),
            case.tolerance,
            &format!("{context} (cp_over_r)"),
        );
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
                "a" => Some(common::input(case, "a")),
                "b" => Some(common::input(case, "b")),
                "c" => Some(common::input(case, "c")),
                "d" => Some(common::input(case, "d")),
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
/// `cp_over_r` is linear in `(a, b, c, d)` by construction, so doubling every
/// coefficient doubles the answer and changing the sign of all four negates it. That
/// is not a property of any particular coefficient set - it holds for all of them -
/// and it is what catches an implementation which, say, multiplies `d` by `theta` twice
/// or drops a term, both of which pass at the values one worked example happens to use.
#[test]
fn the_polynomial_is_linear_in_its_coefficients() {
    let (a, b, c, d, t) = (4.0, 1.0, -0.5, 0.1, 500.0);
    let base = ideal_gas_cp(a, b, c, d, kelvins(t)).unwrap();
    let doubled = ideal_gas_cp(2.0 * a, 2.0 * b, 2.0 * c, 2.0 * d, kelvins(t)).unwrap();
    let negated = ideal_gas_cp(-a, -b, -c, -d, kelvins(t)).unwrap();

    assert!((doubled.cp_over_r - 2.0 * base.cp_over_r).abs() < 1e-15);
    assert!((negated.cp_over_r + base.cp_over_r).abs() < 1e-15);
    assert!((doubled.cp.value - 2.0 * base.cp.value).abs() < 1e-12);
}

/// The reference temperature is what makes a table's printed numbers dimensionless.
///
/// At `T = 1000 K` the reduced temperature is exactly one, so the polynomial collapses
/// to `a + b + c + d`. That is the one temperature at which the four coefficients can
/// be checked against the answer by addition alone, and it is the check that catches an
/// implementation using `T/100` or `T` rather than `T/1000` - all of which agree at any
/// other temperature after rescaling the coefficients, and none of which agree here.
#[test]
fn the_reference_temperature_is_a_thousand_kelvin() {
    let result = ideal_gas_cp(4.0, -20.0, 10.0, 1.0, kelvins(1000.0)).unwrap();
    assert!(
        (result.cp_over_r - (4.0 - 20.0 + 10.0 + 1.0)).abs() < 1e-15,
        "at T = 1000 K the polynomial should be a + b + c + d, got {}",
        result.cp_over_r
    );
    assert!((result.cp.value - result.cp_over_r * azoth_eos::MOLAR_GAS_CONSTANT).abs() < 1e-12);
}

/// A polynomial that has turned over reports it rather than passing in silence.
///
/// The failure this catches is the one that matters: a negative heat capacity fed
/// into `integral Cp dT` produces an enthalpy wrong by an amount nobody can see.
#[test]
fn a_non_positive_heat_capacity_warns() {
    let result = ideal_gas_cp(4.0, -20.0, 10.0, 1.0, kelvins(1000.0)).unwrap();
    assert!(result.cp.value < 0.0);
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code == WarningCode::OutOfValidRange),
        "a negative heat capacity should say so: {:?}",
        result.warnings
    );

    // And a positive one at the same coefficients, inside the range, is clean.
    let ordinary = ideal_gas_cp(4.0, 1.0, -0.5, 0.1, kelvins(500.0)).unwrap();
    assert!(
        ordinary.is_clean(),
        "unexpected warnings: {:?}",
        ordinary.warnings
    );
}

#[test]
fn a_non_positive_temperature_is_refused() {
    for t in [0.0, -1.0] {
        let err = ideal_gas_cp(4.0, 1.0, -0.5, 0.1, kelvins(t)).unwrap_err();
        assert!(matches!(err, AzothError::OutOfRange { .. }), "T = {t}");
        assert_eq!(err.field(), Some("T"), "T = {t}");
    }
}

/// All four coefficients matter, individually.
///
/// Four separate single-coefficient perturbations, each of which must move the answer.
/// An implementation that silently dropped one - the commonest way a four-term
/// polynomial goes wrong - would be indistinguishable at most coefficient sets,
/// because the term it drops is a small correction.
#[test]
fn every_coefficient_changes_the_answer() {
    let base = ideal_gas_cp(4.0, 1.0, -0.5, 0.1, kelvins(700.0)).unwrap();
    for (index, perturbed) in [
        ideal_gas_cp(4.5, 1.0, -0.5, 0.1, kelvins(700.0)).unwrap(),
        ideal_gas_cp(4.0, 1.5, -0.5, 0.1, kelvins(700.0)).unwrap(),
        ideal_gas_cp(4.0, 1.0, -0.9, 0.1, kelvins(700.0)).unwrap(),
        ideal_gas_cp(4.0, 1.0, -0.5, 0.5, kelvins(700.0)).unwrap(),
    ]
    .into_iter()
    .enumerate()
    {
        assert!(
            (perturbed.cp_over_r - base.cp_over_r).abs() > 1e-6,
            "coefficient {index} has no effect on the answer"
        );
    }
}
