//! Spec-driven tests for `eos.scale_saturation_ratio`.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::scale_saturation_ratio;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.scale_saturation_ratio";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::results::ScaleSaturationRatioResult {
    scale_saturation_ratio(
        common::input_str(case, "salt"),
        common::input(case, "x1"),
        common::input(case, "x2"),
        common::input(case, "x_water"),
        common::input(case, "gamma1"),
        common::input(case, "gamma2"),
        common::input(case, "water_activity"),
        case.input("h3o_molality"),
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
                let context = format!("{}::{}", spec.id, case.id);
                common::assert_close(
                    result.saturation_ratio,
                    common::expected(case, "saturation_ratio"),
                    case.tolerance,
                    &format!("{context} (saturation ratio)"),
                );
                common::assert_consistent(&result, &context);
            }
            "property" => match case.property {
                Some("unit_round_trip") => {}
                other => panic!("{}::{}: unhandled property {other:?}", spec.id, case.id),
            },
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 2,
        "expected several active cases, ran {executed}"
    );
}

/// **The ratio is one at saturation and the sign of `ln SR` is the whole verdict.**
///
/// The two captured states bracket it: the brine's NaCl is `1.28e-07` - far under, so it will
/// not precipitate - and its FeCO3 is `2.58e+05`, far over. A port that had the ratio inverted
/// would report both the other way round and the numbers would still be plausible.
#[test]
fn one_is_saturation() {
    let at = |x: f64| {
        scale_saturation_ratio(
            "NaCl",
            x,
            x,
            0.994547937931656,
            1.0,
            1.0,
            1.0,
            None,
            kelvins(298.15),
            pascals(1.0e6),
        )
        .expect("computes")
    };
    // Bisect the ion mole fraction until the ratio is one, which is what saturation means.
    // The saturation mole fraction is `0.1107...` at this `x_water` - above `0.1`, which is
    // where the first bracket put it and why the bisection returned the boundary's `0.815`.
    let (mut low, mut high) = (1.0e-9, 0.5);
    for _ in 0..200 {
        let mid = 0.5 * (low + high);
        if at(mid).saturation_ratio < 1.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    let saturated = 0.5 * (low + high);
    common::assert_close(
        at(saturated).saturation_ratio,
        1.0,
        1.0e-12,
        "at saturation",
    );
    assert!(
        at(saturated * 0.5).saturation_ratio < 1.0,
        "under-saturated"
    );
    assert!(at(saturated * 2.0).saturation_ratio > 1.0, "over-saturated");
}

/// A salt the table does not carry is refused, and `FeS` without a hydrogen molality is too.
#[test]
fn the_refusals() {
    let error = scale_saturation_ratio(
        "unobtainium",
        1.0e-4,
        1.0e-4,
        0.99,
        1.0,
        1.0,
        1.0,
        None,
        kelvins(298.15),
        pascals(1.0e6),
    )
    .expect_err("not a row");
    assert!(
        matches!(error, azoth_core::AzothError::InvalidInput { .. }),
        "{error:?}"
    );

    let error = scale_saturation_ratio(
        "FeS",
        1.0e-4,
        1.0e-4,
        0.99,
        1.0,
        1.0,
        1.0,
        None,
        kelvins(298.15),
        pascals(1.0e6),
    )
    .expect_err("FeS needs the hydrogen molality");
    assert!(
        matches!(error, azoth_core::AzothError::InvalidInput { .. }),
        "{error:?}"
    );
}
