//! Spec-driven tests for `hydraulics.darcy_weisbach`.

use azoth_test_support as common;

use azoth_core::spec::TestCase;
use azoth_core::units::{kilograms_per_cubic_meter, meters, meters_per_second, pascal_seconds};
use azoth_core::{AzothError, CalcResult, FlowRegime, WarningCode};
use azoth_hydraulics::{add_fitting_loss, darcy_weisbach, spec_gen};

const CALC_ID: &str = "hydraulics.darcy_weisbach";

/// The spec declares `mu` optional, so its absence is a legitimate case rather
/// than a missing argument.
fn call(case: &TestCase) -> azoth_hydraulics::DarcyWeisbachResult {
    darcy_weisbach(
        common::input(case, "f"),
        meters(common::input(case, "L")),
        meters(common::input(case, "D")),
        kilograms_per_cubic_meter(common::input(case, "rho")),
        meters_per_second(common::input(case, "v")),
        case.input("mu").map(pascal_seconds),
    )
    .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    common::assert_skips_are_explained(spec);

    let mut executed = 0;
    let mut skipped = 0;
    for case in spec.all_tests() {
        if !case.is_active() {
            skipped += 1;
            continue;
        }
        match case.kind {
            "worked_example" | "reference" => {
                let result = call(case);
                common::assert_close(
                    result.dp.value,
                    common::expected(case, "dp"),
                    case.tolerance,
                    &format!("{}::{} (dp)", spec.id, case.id),
                );
                if let Some(expected_re) = case.expected_value("re") {
                    let re = result.re.unwrap_or_else(|| {
                        panic!("{}::{}: re expected but missing", spec.id, case.id)
                    });
                    common::assert_close(
                        re,
                        expected_re,
                        case.tolerance,
                        &format!("{}::{} (re)", spec.id, case.id),
                    );
                }
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |q| match q {
                        "f" => Some(common::input(case, "f")),
                        "L" => Some(common::input(case, "L")),
                        "D" => Some(common::input(case, "D")),
                        "rho" => Some(common::input(case, "rho")),
                        "v" => Some(common::input(case, "v")),
                        // None when mu was omitted, which is exactly the case
                        // the spec's range check has to report as unchecked.
                        "re" => result.re,
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
            }
            "property" => match case.property {
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
    // The two skipped cases are the unverifiable Crane and Perry's references.
    // They are skipped with a reason, and checked as such by
    // assert_skips_are_explained, but they must not have vanished silently.
    assert_eq!(
        skipped, 2,
        "expected the two source_needed references to be skipped"
    );
}

#[test]
fn unit_round_trip() {
    // The same physical state with length and velocity in feet must give the
    // same pressure drop.
    let si = darcy_weisbach(
        0.02,
        meters(100.0),
        meters(0.1),
        kilograms_per_cubic_meter(998.0),
        meters_per_second(1.5),
        None,
    )
    .unwrap();

    let us = darcy_weisbach(
        0.02,
        uom::si::f64::Length::new::<uom::si::length::foot>(100.0 / 0.3048),
        uom::si::f64::Length::new::<uom::si::length::foot>(0.1 / 0.3048),
        kilograms_per_cubic_meter(998.0),
        uom::si::f64::Velocity::new::<uom::si::velocity::foot_per_second>(1.5 / 0.3048),
        None,
    )
    .unwrap();

    common::assert_same_state(
        si.dp.value,
        us.dp.value,
        1e-12,
        "dp from SI vs US customary",
    );
}

#[test]
fn omitting_viscosity_marks_the_regime_unchecked() {
    // This is the heart of the optional-input design: "checked and fine" and
    // "never checked" must not look the same.
    let unchecked = darcy_weisbach(
        0.02,
        meters(100.0),
        meters(0.1),
        kilograms_per_cubic_meter(998.0),
        meters_per_second(1.5),
        None,
    )
    .unwrap();
    assert!(unchecked.re.is_none());
    assert!(unchecked.regime.is_none());
    assert!(unchecked.has_warning(WarningCode::RangeCheckSkipped));

    // Supplying mu must not change the pressure drop at all - it only enables
    // the check.
    let checked = darcy_weisbach(
        0.02,
        meters(100.0),
        meters(0.1),
        kilograms_per_cubic_meter(998.0),
        meters_per_second(1.5),
        Some(pascal_seconds(1.002e-3)),
    )
    .unwrap();
    assert_eq!(
        unchecked.dp.value.to_bits(),
        checked.dp.value.to_bits(),
        "supplying mu must not change dp"
    );
    assert!(checked.re.is_some());
    assert_eq!(checked.regime, Some(FlowRegime::Turbulent));
    assert!(
        !checked.has_warning(WarningCode::RangeCheckSkipped),
        "with mu supplied the regime check did run, so it must not report as skipped"
    );
}

#[test]
fn transitional_flow_warns_through_this_calc_too() {
    // Selected so Re lands in the 2000-4000 band: q = v*A with a 0.05 m pipe at
    // a low velocity.
    let r = darcy_weisbach(
        0.032,
        meters(10.0),
        meters(0.05),
        kilograms_per_cubic_meter(1000.0),
        meters_per_second(0.1),
        Some(pascal_seconds(1.5e-3)),
    )
    .unwrap();
    assert_eq!(r.regime, Some(FlowRegime::Transitional));
    assert!(r.has_warning(WarningCode::TransitionalFlow));
}

#[test]
fn zero_diameter_is_an_error() {
    // L/D is singular at D = 0.
    let err = darcy_weisbach(
        0.02,
        meters(100.0),
        meters(0.0),
        kilograms_per_cubic_meter(998.0),
        meters_per_second(1.5),
        None,
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("D"));
}

#[test]
fn zero_length_is_allowed_and_gives_zero_drop() {
    // Zero length is physical - a degenerate pipe - so it is not a hard bound.
    let r = darcy_weisbach(
        0.02,
        meters(0.0),
        meters(0.1),
        kilograms_per_cubic_meter(998.0),
        meters_per_second(1.5),
        None,
    )
    .unwrap();
    assert_eq!(r.dp.value, 0.0);
}

#[test]
fn pressure_drop_scales_with_density_and_velocity_squared() {
    // The equation is linear in rho and quadratic in v. Checking the scaling
    // catches a transcription error that a single point would not.
    let base = darcy_weisbach(
        0.02,
        meters(100.0),
        meters(0.1),
        kilograms_per_cubic_meter(998.0),
        meters_per_second(1.5),
        None,
    )
    .unwrap()
    .dp
    .value;

    let denser = darcy_weisbach(
        0.02,
        meters(100.0),
        meters(0.1),
        kilograms_per_cubic_meter(1996.0),
        meters_per_second(1.5),
        None,
    )
    .unwrap()
    .dp
    .value;
    common::assert_close(denser, 2.0 * base, 1e-12, "dp is linear in rho");

    let faster = darcy_weisbach(
        0.02,
        meters(100.0),
        meters(0.1),
        kilograms_per_cubic_meter(998.0),
        meters_per_second(3.0),
        None,
    )
    .unwrap()
    .dp
    .value;
    common::assert_close(faster, 4.0 * base, 1e-12, "dp is quadratic in v");
}

#[test]
fn fitting_loss_adds_the_velocity_head_term() {
    // The CLI's composition, exercised here so the CLI inherits a tested path.
    // dP_total = f*(L/D)*(rho v^2/2) + K*(rho v^2/2).
    let (f, l, d, rho, v) = (0.018, 100.0, 0.1, 998.0, 1.5);
    let straight = darcy_weisbach(
        f,
        meters(l),
        meters(d),
        kilograms_per_cubic_meter(rho),
        meters_per_second(v),
        None,
    )
    .unwrap()
    .dp
    .value;

    let k_total = 0.684;
    let total = add_fitting_loss(straight, k_total, rho, v);
    let velocity_head = rho * v * v / 2.0;

    common::assert_close(
        total,
        straight + k_total * velocity_head,
        1e-15,
        "fitting loss is K times the velocity head",
    );

    // Pin the actual magnitude rather than asserting a threshold chosen by feel.
    // For these inputs: velocity head = 998*1.5**2/2 = 1122.75 Pa, the fitting
    // term is 0.684 * 1122.75 = 767.961 Pa, and the straight-pipe term is
    // 0.018 * 1000 * 1122.75 = 20209.5 Pa, so the fittings add 3.8%.
    let fitting_term = k_total * velocity_head;
    common::assert_close(fitting_term, 767.961, 1e-9, "fitting term in Pa");
    let ratio = fitting_term / straight;
    assert!(
        (0.037..0.039).contains(&ratio),
        "fitting loss is {:.4} of the straight-pipe loss, expected about 0.038",
        ratio
    );
}
