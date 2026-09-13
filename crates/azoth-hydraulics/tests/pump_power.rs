//! Spec-driven tests for `hydraulics.pump_power`.

use azoth_core::units::{
    Length, MassDensity, VolumeRate, cubic_meters_per_second, kilograms_per_cubic_meter, meters,
};
use azoth_core::{AzothError, CalcResult};
use azoth_hydraulics::spec_gen;
use azoth_hydraulics::{STANDARD_GRAVITY_M_S2, pump_power};
use azoth_test_support as common;

const CALC_ID: &str = "hydraulics.pump_power";

fn call(case: &azoth_core::spec::TestCase) -> azoth_hydraulics::PumpPowerResult {
    pump_power(
        kilograms_per_cubic_meter(common::input(case, "rho")),
        cubic_meters_per_second(common::input(case, "q")),
        meters(common::input(case, "H")),
        common::input(case, "eta"),
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
                    result.power.value,
                    common::expected(case, "power"),
                    case.tolerance,
                    &format!("{}::{} (power)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "rho" => Some(common::input(case, "rho")),
                        "q" => Some(common::input(case, "q")),
                        "H" => Some(common::input(case, "H")),
                        "eta" => Some(common::input(case, "eta")),
                        "power" => Some(result.power.value),
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

/// Power rises with density, flow and head, and falls with efficiency.
///
/// The four partial derivatives. A point-value test can pass with a wrong formula
/// that happens to be right at one operating point; an inverted ratio or a sign
/// error cannot survive all four directions at once.
#[test]
fn monotonic() {
    let at = |rho: f64, q: f64, h: f64, eta: f64| {
        pump_power(
            kilograms_per_cubic_meter(rho),
            cubic_meters_per_second(q),
            meters(h),
            eta,
        )
        .unwrap()
        .power
        .value
    };
    let base = at(998.0, 0.01, 30.0, 0.75);

    assert!(
        at(1200.0, 0.01, 30.0, 0.75) > base,
        "power must rise with rho"
    );
    assert!(at(998.0, 0.02, 30.0, 0.75) > base, "power must rise with q");
    assert!(at(998.0, 0.01, 60.0, 0.75) > base, "power must rise with H");
    assert!(
        at(998.0, 0.01, 30.0, 0.50) > base,
        "power must rise as eta falls"
    );
}

/// The same operating point described in SI and in US customary units.
///
/// All three dimensioned inputs are converted, which is the point: density, flow
/// and head each carry a unit and each is converted separately, so an accidental
/// factor between them shows up as a different power rather than cancelling out.
///
/// The US values are obtained through `uom`'s own `get::<U>()` rather than by
/// writing the conversion factors out here. That is deliberate: the first version
/// of this test restated them, and the factors were only accurate to about 7e-8,
/// which failed at the 1e-12 tolerance the rest of the suite uses. Restating what
/// the units library already knows is exactly the mistake that produced the `mm`
/// bug - see `azoth_core::units` - and it hid here as a test that looked like it
/// was checking the calc and was really checking the test author's arithmetic.
#[test]
fn unit_round_trip() {
    use uom::si::length::foot;
    use uom::si::mass_density::pound_per_cubic_foot;
    use uom::si::volume_rate::cubic_foot_per_second;

    let si_rho = kilograms_per_cubic_meter(998.0);
    let si_q = cubic_meters_per_second(0.01);
    let si_h = meters(30.0);

    let si = pump_power(si_rho, si_q, si_h, 0.75).unwrap().power.value;

    let us = pump_power(
        MassDensity::new::<pound_per_cubic_foot>(si_rho.get::<pound_per_cubic_foot>()),
        VolumeRate::new::<cubic_foot_per_second>(si_q.get::<cubic_foot_per_second>()),
        Length::new::<foot>(si_h.get::<foot>()),
        0.75,
    )
    .unwrap()
    .power
    .value;

    common::assert_same_state(si, us, 1e-12, "power from SI vs US customary state");
}

#[test]
fn standard_gravity_is_the_defined_value() {
    // Pinned in both languages, so a change to one of them fails a test rather than
    // putting a 0.5% error between the two implementations that neither reports.
    assert_eq!(STANDARD_GRAVITY_M_S2, 9.80665);
}

#[test]
fn an_inefficiency_of_one_reduces_to_the_hydraulic_power() {
    // The limiting case: a lossless pump delivers exactly rho * g * q * H. This
    // separates the two halves of the equation, so a failure says which half broke.
    let hydraulic = 998.0 * STANDARD_GRAVITY_M_S2 * 0.01 * 30.0;
    let r = pump_power(
        kilograms_per_cubic_meter(998.0),
        cubic_meters_per_second(0.01),
        meters(30.0),
        1.0,
    )
    .unwrap();
    assert!((r.power.value - hydraulic).abs() < 1e-9);
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

#[test]
fn a_negative_flow_is_an_error_rather_than_reversed_flow() {
    let err = pump_power(
        kilograms_per_cubic_meter(998.0),
        cubic_meters_per_second(-0.01),
        meters(30.0),
        0.75,
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("q"));
}

#[test]
fn an_efficiency_above_one_is_an_error() {
    // Not an engineering uncertainty to warn about: more hydraulic power out than
    // shaft power in is a violation of the first law.
    let err = pump_power(
        kilograms_per_cubic_meter(998.0),
        cubic_meters_per_second(0.01),
        meters(30.0),
        1.01,
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("eta"));
}

#[test]
fn a_zero_efficiency_is_an_error_rather_than_infinite_power() {
    let err = pump_power(
        kilograms_per_cubic_meter(998.0),
        cubic_meters_per_second(0.01),
        meters(30.0),
        0.0,
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("eta"));
}

#[test]
fn zero_flow_needs_zero_power_and_is_not_a_warning() {
    // Allowed, and correct rather than merely tolerable: no flow moves no fluid.
    let r = pump_power(
        kilograms_per_cubic_meter(998.0),
        cubic_meters_per_second(0.0),
        meters(30.0),
        0.75,
    )
    .unwrap();
    assert_eq!(r.power.value, 0.0);
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}
