//! Spec-driven tests for `hydraulics.control_valve_cv`.

use azoth_core::units::pascals;
use azoth_core::{AzothError, CalcResult};
use azoth_hydraulics::spec_gen;
use azoth_hydraulics::{CV_TO_SI, control_valve_cv};
use azoth_test_support as common;

const CALC_ID: &str = "hydraulics.control_valve_cv";

fn call(case: &azoth_core::spec::TestCase) -> azoth_hydraulics::ControlValveCvResult {
    control_valve_cv(
        common::input(case, "Cv"),
        pascals(common::input(case, "dP")),
        common::input(case, "SG"),
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
                    result.q.value,
                    common::expected(case, "q"),
                    case.tolerance,
                    &format!("{}::{} (q)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "Cv" => Some(common::input(case, "Cv")),
                        "dP" => Some(common::input(case, "dP")),
                        "SG" => Some(common::input(case, "SG")),
                        "q" => Some(result.q.value),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
            }
            "property" => match case.property {
                Some("monotonic") => monotonic(),
                Some("unit_round_trip") => unit_round_trip(),
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

/// Flow rises with the coefficient and the drop, falls with gravity, and is zero
/// when the drop is.
///
/// The directions place each input in the expression: a wrong power, a misplaced
/// square root or an inverted ratio all show up as a direction changing, and the
/// zero case is the one a form that divided by `dP` would get wrong.
#[test]
fn monotonic() {
    let at = |cv: f64, dp: f64, sg: f64| control_valve_cv(cv, pascals(dp), sg).unwrap().q.value;
    let base = at(10.0, 100000.0, 1.0);

    assert!(at(20.0, 100000.0, 1.0) > base, "q must rise with Cv");
    assert!(at(10.0, 200000.0, 1.0) > base, "q must rise with dP");
    assert!(at(10.0, 100000.0, 1.5) < base, "q must fall as SG rises");
    assert_eq!(at(10.0, 0.0, 1.0), 0.0, "no drop, no flow");
}

/// The same valve described with the drop in psi rather than pascals.
///
/// `Cv` and `SG` are dimensionless and cannot round-trip, so `dP` is the only input
/// here that carries a unit. Converting it is a real check rather than a restatement
/// of the coefficient convention: it runs the pressure through `uom`'s conversion
/// and back into the relation, which a calc that had baked `psi` into its constant
/// would fail.
#[test]
fn unit_round_trip() {
    use uom::si::pressure::psi;

    let si_dp = pascals(100000.0);
    let si = control_valve_cv(10.0, si_dp, 1.0).unwrap().q.value;
    let us = control_valve_cv(
        10.0,
        azoth_core::units::Pressure::new::<psi>(si_dp.get::<psi>()),
        1.0,
    )
    .unwrap()
    .q
    .value;

    common::assert_same_state(si, us, 1e-12, "q from Pa vs psi");
}

/// The definition of `Cv`, stated directly: a coefficient of 1 passes 1 gpm of
/// water at 1 psi.
///
/// This is the check on the conversion constant rather than on the arithmetic. If
/// `CV_TO_SI` were wrong by any factor, the flow at unit coefficient and unit drop
/// would not be one gallon per minute, and the spec's derivation would not close.
#[test]
fn a_coefficient_of_one_passes_one_gallon_per_minute() {
    let one_psi = azoth_hydraulics::control_valve_cv::psi_in_pascals();
    let r = control_valve_cv(1.0, pascals(one_psi), 1.0).unwrap();
    let gallon_m3 = azoth_hydraulics::control_valve_cv::GALLON_M3;
    assert!(
        (r.q.value - gallon_m3 / 60.0).abs() < 1e-18,
        "one Cv at one psi should pass one gpm, got {} m^3/s",
        r.q.value
    );
}

/// The constant the two implementations use must be the same number.
///
/// Python computes `CV_TO_SI` from the definitions; Rust stores it as a literal
/// because `const` cannot call `sqrt`. The cross-language agreement tests compare
/// *results*, which would hide a constant that differed by less than the tolerance
/// on one case and not on another - so the constant itself is pinned here.
#[test]
fn the_conversion_constant_is_the_expected_value() {
    assert!((CV_TO_SI - 7.598_054_212_083_37e-7).abs() < 1e-21);
}

#[test]
fn a_non_positive_coefficient_is_an_error() {
    let err = control_valve_cv(0.0, pascals(100000.0), 1.0).unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("Cv"));
}

#[test]
fn a_negative_drop_is_an_error_rather_than_reversed_flow() {
    let err = control_valve_cv(10.0, pascals(-1000.0), 1.0).unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("dP"));
}

#[test]
fn a_non_positive_specific_gravity_is_an_error() {
    let err = control_valve_cv(10.0, pascals(100000.0), 0.0).unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("SG"));
}

#[test]
fn a_zero_drop_gives_zero_flow_and_no_warning() {
    // Allowed, and correct rather than merely tolerable: no pressure difference
    // drives no flow.
    let r = control_valve_cv(10.0, pascals(0.0), 1.0).unwrap();
    assert_eq!(r.q.value, 0.0);
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}
