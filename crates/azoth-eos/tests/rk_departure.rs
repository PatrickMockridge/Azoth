//! Spec-driven tests for `eos.rk_departure`.

use azoth_core::AzothError;
use azoth_eos::rk_departure;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.rk_departure";

const FIELDS: [&str; 4] = ["ln_phi", "h_dep_rt", "s_dep_r", "cp_dep_r"];

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::RkDepartureResult {
    rk_departure(
        common::input(case, "a_reduced"),
        common::input(case, "b_reduced"),
        common::input(case, "z"),
    )
    .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
}

fn field(result: &azoth_eos::RkDepartureResult, name: &str) -> f64 {
    match name {
        "ln_phi" => result.ln_phi,
        "h_dep_rt" => result.h_dep_rt,
        "s_dep_r" => result.s_dep_r,
        _ => result.cp_dep_r,
    }
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
                for name in FIELDS {
                    common::assert_close(
                        field(&result, name),
                        common::expected(case, name),
                        case.tolerance,
                        &format!("{}::{} ({name})", spec.id, case.id),
                    );
                }
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "a_reduced" => Some(common::input(case, "a_reduced")),
                        "b_reduced" => Some(common::input(case, "b_reduced")),
                        "z" => Some(common::input(case, "z")),
                        "z_minus_b_reduced" => {
                            Some(common::input(case, "z") - common::input(case, "b_reduced"))
                        }
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
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

/// `h_dep_rt - s_dep_r` equals `ln_phi` - the Gibbs identity, exact in real arithmetic.
#[test]
fn the_gibbs_identity_holds() {
    for a in [0.05, 0.1866943088347048, 0.4] {
        for b in [0.005, 0.02707510936404929, 0.08] {
            for z in [0.05, 0.3, 0.8, 0.95] {
                if z <= b {
                    continue;
                }
                let r = rk_departure(a, b, z).unwrap();
                let identity = r.h_dep_rt - r.s_dep_r;
                assert!(
                    (identity - r.ln_phi).abs() < 1e-12,
                    "A = {a}, B = {b}, z = {z}: h - s = {identity} against ln_phi = {}",
                    r.ln_phi
                );
            }
        }
    }
}

#[test]
fn a_z_at_or_below_b_reduced_is_an_error_naming_the_difference() {
    for z in [0.02707510936404929, 0.01, 0.0, -1.0] {
        let err = rk_departure(0.2, 0.02707510936404929, z).unwrap_err();
        assert!(matches!(err, AzothError::OutOfRange { .. }));
        assert_eq!(err.field(), Some("z_minus_b_reduced"), "for z = {z}");
    }
}
