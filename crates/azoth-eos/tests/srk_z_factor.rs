//! Spec-driven tests for `eos.srk_z_factor`.

use azoth_core::AzothError;
use azoth_eos::spec_gen;
use azoth_eos::{Cubic, RootStructure, srk_z_factor};
use azoth_test_support as common;

const CALC_ID: &str = "eos.srk_z_factor";

/// The SRK cubic as `f(z)`, from the crate's geometry - so a wrong coefficient in the
/// crate is wrong here too, which is the point: this checks the returned root against
/// the published polynomial.
fn cubic(a: f64, b: f64, z: f64) -> f64 {
    let (c2, c1, c0) = Cubic::Srk.z_coefficients(a, b);
    ((z + c2) * z + c1) * z + c0
}

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::SrkZFactorResult {
    srk_z_factor(
        common::input(case, "a_reduced"),
        common::input(case, "b_reduced"),
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
                let a = common::input(case, "a_reduced");
                let b = common::input(case, "b_reduced");
                for field in ["z_min", "z_max"] {
                    let actual = if field == "z_min" {
                        result.z_min
                    } else {
                        result.z_max
                    };
                    common::assert_close(
                        actual,
                        common::expected(case, field),
                        case.tolerance,
                        &format!("{}::{} ({field})", spec.id, case.id),
                    );
                    let residual = cubic(a, b, actual);
                    assert!(
                        residual.abs() < 1e-12,
                        "{}::{}: {field} = {actual} does not satisfy the cubic: f = {residual:e}",
                        spec.id,
                        case.id
                    );
                }
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
            }
            "property" => panic!("{}::{}: unhandled property", spec.id, case.id),
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 2,
        "expected several active cases, ran {executed}"
    );
}

/// The admissible root count is one or three, and never two.
///
/// The twin of `pr_z_factor`'s sweep. `f(b_reduced)` is exactly `-2*b_reduced**2`, so `B`
/// lies below all three roots or between the middle and the largest one, and there is no
/// in-between case; a sign error in the constant term would make a count of two
/// reachable.
#[test]
fn the_admissible_count_is_never_two() {
    let mut seen_one = 0_u32;
    let mut seen_three = 0_u32;
    for a_step in 0..=120 {
        for b_step in 1..=98 {
            let a = 0.6 * f64::from(a_step) / 120.0;
            let b = 0.002 + 0.098 * f64::from(b_step) / 98.0;
            let r = srk_z_factor(a, b).unwrap();
            match r.root_structure {
                RootStructure::One => seen_one += 1,
                RootStructure::Three => seen_three += 1,
            }
            assert!(
                cubic(a, b, b) < 0.0,
                "f(B) should be -2*B**2 < 0 for A = {a}, B = {b}"
            );
        }
    }
    assert!(
        seen_one > 0 && seen_three > 0,
        "the sweep should reach both cases"
    );
}

#[test]
fn a_non_positive_b_reduced_is_an_error() {
    for b in [0.0, -0.01] {
        let err = srk_z_factor(0.2, b).unwrap_err();
        assert!(matches!(err, AzothError::OutOfRange { .. }));
        assert_eq!(err.field(), Some("b_reduced"), "for B = {b}");
    }
}

#[test]
fn a_negative_a_reduced_is_an_error() {
    let err = srk_z_factor(-0.01, 0.024).unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("a_reduced"));
}
