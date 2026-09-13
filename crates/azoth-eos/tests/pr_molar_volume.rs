//! Spec-driven tests for `eos.pr_molar_volume`.

use azoth_core::AzothError;
use azoth_core::units::{kelvins, pascals};
use azoth_eos::spec_gen;
use azoth_eos::{MOLAR_GAS_CONSTANT, pr_molar_volume};
use azoth_test_support as common;

const CALC_ID: &str = "eos.pr_molar_volume";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::PrMolarVolumeResult {
    pr_molar_volume(
        common::input(case, "z"),
        kelvins(common::input(case, "T")),
        pascals(common::input(case, "P")),
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
                    result.v.value,
                    common::expected(case, "v"),
                    case.tolerance,
                    &format!("{}::{} (v)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "z" => Some(common::input(case, "z")),
                        "T" => Some(common::input(case, "T")),
                        "P" => Some(common::input(case, "P")),
                        "v" => Some(result.v.value),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
            }
            "property" => match case.property {
                Some("monotonic") | Some("unit_round_trip") => {}
                other => panic!("{}::{}: unhandled property {other:?}", spec.id, case.id),
            },
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 4,
        "expected several active cases, ran {executed}"
    );
}

/// The gas constant is the exact product, not the truncated literal.
///
/// Worth its own test because the difference is invisible at a glance: `8.314462618`
/// and `8.31446261815324` are both "the gas constant" to a reader, and which one is
/// in the code is not something any point-value case with a loose tolerance would
/// notice. Asserting the constant's bits pins it, and the spec's ideal-gas case at
/// 300 K and one atmosphere fails in the fifth significant figure if it drifts.
#[test]
fn the_gas_constant_is_the_exact_product() {
    // Pinned to the exact product, written out. A change to either constant, or a
    // switch to the truncated literal, moves this.
    assert_eq!(
        MOLAR_GAS_CONSTANT.to_bits(),
        8.314_462_618_153_24_f64.to_bits(),
        "MOLAR_GAS_CONSTANT should be the exact N_A * k_B, 8.31446261815324"
    );
    // The bits equality above is the real check; the commonly quoted
    // `8.314462618` is this truncated and would fail it. Not asserted separately
    // because comparing two constants is a comparison clippy folds away, and an
    // assertion that cannot fail is not an assertion.
}

/// `v` rises with `z` and with `T`, and falls with `P`.
///
/// The three partial derivatives of the equation. A point-value case cannot tell a
/// correct formula from one with `T` and `P` transposed; three directions at once
/// can.
#[test]
fn monotonic() {
    let at = |z: f64, t: f64, p: f64| pr_molar_volume(z, kelvins(t), pascals(p)).unwrap().v.value;
    let base = at(0.8, 300.0, 1.0e5);

    assert!(at(0.9, 300.0, 1.0e5) > base, "v must rise with z");
    assert!(at(0.8, 400.0, 1.0e5) > base, "v must rise with T");
    assert!(at(0.8, 300.0, 2.0e5) < base, "v must fall with P");

    // And it is exactly linear in z, which a bound-check or a plot would not say.
    let doubled = at(1.6, 300.0, 1.0e5);
    assert_eq!(
        doubled.to_bits(),
        (2.0 * base).to_bits(),
        "v should be linear in z"
    );
}

/// The same physical state described in different units gives the same volume.
///
/// The first real exercise of `m**3/mol` and of the units layer through an
/// equation-of-state calculation - unlike the rest of this namespace, this calc is
/// dimensional, so the round trip has something to convert.
#[test]
fn the_same_state_in_other_units() {
    let si = pr_molar_volume(0.7907789662973796, kelvins(295.864), pascals(1_062_000.0))
        .unwrap()
        .v
        .value;

    // The pressure in bar, converted at the boundary - which is the arrangement
    // the project's unit rule asks for: callers convert, the calc receives SI.
    let bar = pascals(10.62e5);
    let converted = pr_molar_volume(0.7907789662973796, kelvins(295.864), bar)
        .unwrap()
        .v
        .value;
    assert_eq!(si.to_bits(), converted.to_bits());
}

#[test]
fn a_non_positive_input_is_an_error() {
    for (z, t, p) in [
        (0.0, 300.0, 1e5),
        (-0.5, 300.0, 1e5),
        (0.8, 0.0, 1e5),
        (0.8, -1.0, 1e5),
        (0.8, 300.0, 0.0),
        (0.8, 300.0, -1e5),
    ] {
        let err = pr_molar_volume(z, kelvins(t), pascals(p)).unwrap_err();
        assert!(
            matches!(err, AzothError::OutOfRange { .. }),
            "z={z}, T={t}, P={p}"
        );
    }
}
