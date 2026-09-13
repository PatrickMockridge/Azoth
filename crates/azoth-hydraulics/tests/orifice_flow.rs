//! Spec-driven tests for `hydraulics.orifice_flow`.

use azoth_core::units::{
    Length, MassDensity, Pressure, VolumeRate, kilograms_per_cubic_meter, millimeters, pascals,
};
use azoth_core::{AzothError, CalcResult};
use azoth_hydraulics::orifice_flow;
use azoth_hydraulics::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "hydraulics.orifice_flow";

fn call(case: &azoth_core::spec::TestCase) -> azoth_hydraulics::OrificeFlowResult {
    orifice_flow(
        millimeters(common::input(case, "d")),
        pascals(common::input(case, "dP")),
        kilograms_per_cubic_meter(common::input(case, "rho")),
        common::input(case, "Cd"),
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
                        "d" => Some(common::input(case, "d")),
                        "dP" => Some(common::input(case, "dP")),
                        "rho" => Some(common::input(case, "rho")),
                        "Cd" => Some(common::input(case, "Cd")),
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

/// Flow rises with the bore, the differential and the coefficient, and is zero when
/// the differential is.
///
/// The directions locate each input in the expression: a wrong power on `d`, a
/// missing square root or an inverted ratio all show up as a direction changing.
/// The zero case is the one a form that divided by `dP` would get wrong.
#[test]
fn monotonic() {
    let at = |d: f64, dp: f64, rho: f64, cd: f64| {
        orifice_flow(
            millimeters(d),
            pascals(dp),
            kilograms_per_cubic_meter(rho),
            cd,
        )
        .unwrap()
        .q
        .value
    };
    let base = at(50.0, 25000.0, 998.0, 0.62);

    assert!(at(100.0, 25000.0, 998.0, 0.62) > base, "q must rise with d");
    assert!(at(50.0, 50000.0, 998.0, 0.62) > base, "q must rise with dP");
    assert!(
        at(50.0, 25000.0, 500.0, 0.62) > base,
        "q must rise as rho falls"
    );
    assert!(at(50.0, 25000.0, 998.0, 0.80) > base, "q must rise with Cd");

    assert_eq!(at(50.0, 0.0, 998.0, 0.62), 0.0, "no difference, no flow");
}

/// The same orifice and operating point described in SI and in US customary units.
///
/// This is the test the `mm` fix was waiting for. `d` is the one input in the whole
/// registry whose spec unit is not its own SI base unit, so this is the first calc
/// through which the millimetre conversion runs end to end - the case that was wrong
/// by a factor of 1000 until M1, and that no calculation had exercised until now.
///
/// The US values are obtained through `uom`'s own `get::<U>()`. Writing the factors
/// out by hand is how the pump_power round trip first failed: they were good to 7e-8
/// and the suite asserts to 1e-12.
#[test]
fn unit_round_trip() {
    use uom::si::length::inch;
    use uom::si::mass_density::pound_per_cubic_foot;
    use uom::si::pressure::psi;

    let si_d = millimeters(50.0);
    let si_dp = pascals(25000.0);
    let si_rho = kilograms_per_cubic_meter(998.0);

    let si = orifice_flow(si_d, si_dp, si_rho, 0.62).unwrap().q.value;

    let us = orifice_flow(
        Length::new::<inch>(si_d.get::<inch>()),
        Pressure::new::<psi>(si_dp.get::<psi>()),
        MassDensity::new::<pound_per_cubic_foot>(si_rho.get::<pound_per_cubic_foot>()),
        0.62,
    )
    .unwrap()
    .q
    .value;

    common::assert_same_state(si, us, 1e-12, "q from SI vs US customary state");
}

#[test]
fn an_ideal_coefficient_gives_the_ideal_flow() {
    // The limiting case, which separates the coefficient from the rest of the
    // expression: if the worked example and this both fail, the fault is in the area
    // or the root; if only the worked example fails, it is in the coefficient.
    let ideal = std::f64::consts::PI * 0.05 * 0.05 / 4.0 * (2.0_f64 * 25000.0 / 998.0).sqrt();
    let r = orifice_flow(
        millimeters(50.0),
        pascals(25000.0),
        kilograms_per_cubic_meter(998.0),
        1.0,
    )
    .unwrap();
    assert!((r.q.value - ideal).abs() < 1e-15, "got {}", r.q.value);
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

#[test]
fn a_zero_bore_is_an_error() {
    let err = orifice_flow(
        millimeters(0.0),
        pascals(25000.0),
        kilograms_per_cubic_meter(998.0),
        0.62,
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("d"));
}

#[test]
fn a_negative_differential_is_an_error_rather_than_reversed_flow() {
    let err = orifice_flow(
        millimeters(50.0),
        pascals(-25000.0),
        kilograms_per_cubic_meter(998.0),
        0.62,
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("dP"));
}

#[test]
fn a_coefficient_above_one_is_an_error() {
    // More flow than an obstruction of that area can pass at that difference is not
    // an uncertainty to warn about but a violation of the relation.
    let err = orifice_flow(
        millimeters(50.0),
        pascals(25000.0),
        kilograms_per_cubic_meter(998.0),
        1.01,
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("Cd"));
}

#[test]
fn a_millimetre_bore_matches_the_same_bore_in_metres() {
    // The `mm` conversion stated directly, so a regression in `to_si` shows up here
    // as well as in the units module's own tests. A 50 mm bore and a 0.05 m bore are
    // the same orifice and must give the same flow; if the millimetre path returned
    // millimetres rather than metres, this is the assertion that reports it.
    let in_mm = orifice_flow(
        millimeters(50.0),
        pascals(25000.0),
        kilograms_per_cubic_meter(998.0),
        0.62,
    )
    .unwrap();
    let in_m = orifice_flow(
        Length::new::<uom::si::length::meter>(0.05),
        pascals(25000.0),
        kilograms_per_cubic_meter(998.0),
        0.62,
    )
    .unwrap();
    assert_eq!(in_mm.q.value.to_bits(), in_m.q.value.to_bits());
}

#[test]
fn volume_rate_converts_back_to_the_declared_unit() {
    // The output is a `VolumeRate`, and the Python side rebuilds it from the SI
    // magnitude plus the unit string. Pinned here so a change to either half is
    // visible from the Rust side too.
    let r = orifice_flow(
        millimeters(50.0),
        pascals(25000.0),
        kilograms_per_cubic_meter(998.0),
        0.62,
    )
    .unwrap();
    let _: VolumeRate = r.q;
    assert!((r.q.value - 0.008616706712060997).abs() < 1e-15);
}
