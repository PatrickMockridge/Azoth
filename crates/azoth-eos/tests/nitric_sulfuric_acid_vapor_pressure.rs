//! Spec-driven tests for `eos.nitric_sulfuric_acid_vapor_pressure`.

use azoth_core::units::kelvins;
use azoth_core::{AzothError, WarningCode};
use azoth_eos::nitric_sulfuric_acid_vapor_pressure;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.nitric_sulfuric_acid_vapor_pressure";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::NitricSulfuricAcidVaporPressureResult {
    nitric_sulfuric_acid_vapor_pressure(kelvins(common::input(case, "T")))
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
                let context = &format!("{}::{}", spec.id, case.id);
                for name in ["p_water", "p_nitric_acid", "p_sulfuric_acid"] {
                    let got = match name {
                        "p_water" => result.p_water.value,
                        "p_nitric_acid" => result.p_nitric_acid.value,
                        _ => result.p_sulfuric_acid.value,
                    };
                    common::assert_close(
                        got,
                        common::expected(case, name),
                        case.tolerance,
                        &format!("{context} ({name})"),
                    );
                }
                common::assert_consistent(&result, context);
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "T" => Some(common::input(case, "T")),
                        _ => None,
                    },
                    context,
                );
            }
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 2,
        "expected several active cases, ran {executed}"
    );
}

/// The stated range is a warning, and outside it the answer is still produced.
///
/// The arithmetic is well defined past 298 K, and inspecting the extrapolation is a
/// legitimate thing for a caller to be doing - which is why the bound is a warning rather
/// than a refusal, and why this asserts both halves of that.
#[test]
fn outside_the_stated_range_is_a_warning_and_not_a_failure() {
    let inside = nitric_sulfuric_acid_vapor_pressure(kelvins(273.15)).unwrap();
    assert!(
        !inside
            .warnings
            .iter()
            .any(|w| w.code == WarningCode::OutOfValidRange),
        "273.15 K is inside the stated range"
    );

    let outside = nitric_sulfuric_acid_vapor_pressure(kelvins(340.0)).unwrap();
    assert!(
        outside
            .warnings
            .iter()
            .any(|w| w.code == WarningCode::OutOfValidRange),
        "340 K is outside it, and the caller should be told: {:?}",
        outside.warnings
    );
}

/// The nitric-acid form is a pole at 43 K, and it is refused rather than returned.
///
/// The three correlations differ in where they break, and only this one breaks *inside*
/// a range a caller might ask for: `B/(T - 43.0)` is infinite at 43 K and negative below,
/// so below it the "vapour pressure" would be a negative number rather than a failure.
#[test]
fn the_nitric_acid_pole_is_refused() {
    for bad in [43.0, 40.0, 1.0] {
        let err = nitric_sulfuric_acid_vapor_pressure(kelvins(bad)).unwrap_err();
        assert!(
            matches!(err, AzothError::OutOfRange { .. }),
            "{bad}: {err:?}"
        );
    }
}

/// Water at the ice point is 611 Pa, and that is what fixes the unit factor.
///
/// The water correlation is written in millibars while every other vapour pressure in
/// this library is in pascals, and the factor is the one thing in the port that a reader
/// cannot check by looking at the equation. A measured value can.
#[test]
fn water_at_the_ice_point_matches_its_measured_vapour_pressure() {
    let r = nitric_sulfuric_acid_vapor_pressure(kelvins(273.15)).unwrap();
    assert!(
        (r.p_water.value - 611.0).abs() < 2.0,
        "water at 273.15 K is 611 Pa, not {}",
        r.p_water.value
    );
}
