//! Spec-driven tests for `hydraulics.choked_flow_area`.

use azoth_core::units::{
    MassDensity, MassRate, Pressure, kilograms_per_cubic_meter, kilograms_per_second, pascals,
};
use azoth_core::{AzothError, CalcResult, WarningCode};
use azoth_hydraulics::choked_flow_area;
use azoth_hydraulics::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "hydraulics.choked_flow_area";

fn call(case: &azoth_core::spec::TestCase) -> azoth_hydraulics::ChokedFlowAreaResult {
    choked_flow_area(
        kilograms_per_second(common::input(case, "m_dot")),
        pascals(common::input(case, "P0")),
        kilograms_per_cubic_meter(common::input(case, "rho0")),
        common::input(case, "k"),
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
                    result.a.value,
                    common::expected(case, "a"),
                    case.tolerance,
                    &format!("{}::{} (a)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "m_dot" => Some(common::input(case, "m_dot")),
                        "P0" => Some(common::input(case, "P0")),
                        "rho0" => Some(common::input(case, "rho0")),
                        "k" => Some(common::input(case, "k")),
                        "a" => Some(result.a.value),
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

/// The area rises with mass flow and falls with pressure, density and `k`.
///
/// The last of those is the least obvious and the most worth having: a higher
/// isentropic exponent means a larger geometric factor and so a higher critical flux,
/// so the same mass flow needs *less* area. A form that got the exponent's sign or
/// its numerator wrong would move that direction.
#[test]
fn monotonic() {
    let at = |m: f64, p: f64, r: f64, k: f64| {
        choked_flow_area(
            kilograms_per_second(m),
            pascals(p),
            kilograms_per_cubic_meter(r),
            k,
        )
        .unwrap()
        .a
        .value
    };
    let base = at(1.0, 1.0e6, 10.0, 1.4);

    assert!(at(2.0, 1.0e6, 10.0, 1.4) > base, "a must rise with m_dot");
    assert!(at(1.0, 5.0e5, 10.0, 1.4) > base, "a must rise as P0 falls");
    assert!(at(1.0, 1.0e6, 5.0, 1.4) > base, "a must rise as rho0 falls");
    assert!(
        at(1.0, 1.0e6, 10.0, 1.6666666666666667) < base,
        "a must fall as k rises"
    );
    assert_eq!(at(0.0, 1.0e6, 10.0, 1.4), 0.0, "no flow, no area");
}

/// The same throat described in SI and in US customary units.
///
/// All three dimensioned inputs are converted, through `uom`'s own `get::<U>()` rather
/// than through factors written out here - writing them out is how the pump_power
/// round trip first failed, at 7e-8 against a 1e-12 tolerance.
#[test]
fn unit_round_trip() {
    use uom::si::mass_density::pound_per_cubic_foot;
    use uom::si::mass_rate::pound_per_second;
    use uom::si::pressure::psi;

    let si_m = kilograms_per_second(1.0);
    let si_p = pascals(1.0e6);
    let si_r = kilograms_per_cubic_meter(10.0);
    let si = choked_flow_area(si_m, si_p, si_r, 1.4).unwrap().a.value;

    let us = choked_flow_area(
        MassRate::new::<pound_per_second>(si_m.get::<pound_per_second>()),
        Pressure::new::<psi>(si_p.get::<psi>()),
        MassDensity::new::<pound_per_cubic_foot>(si_r.get::<pound_per_cubic_foot>()),
        1.4,
    )
    .unwrap()
    .a
    .value;

    common::assert_same_state(si, us, 1e-12, "a from SI vs US customary units");
}

/// An exponent above the monatomic limit warns but still computes.
///
/// This is the first warning bound in five calcs, and the reason is that `k` is a
/// property of a substance with a known range rather than a coefficient in a
/// convention. The value is still returned: the relation evaluates perfectly well at
/// `k = 1.8`, and refusing would be less useful than saying so.
#[test]
fn an_exponent_above_the_monatomic_limit_warns_but_computes() {
    let r = choked_flow_area(
        kilograms_per_second(1.0),
        pascals(1.0e6),
        kilograms_per_cubic_meter(10.0),
        1.8,
    )
    .unwrap();
    assert!(r.has_warning(WarningCode::OutOfValidRange));
    assert!(r.a.value.is_finite() && r.a.value > 0.0);
}

#[test]
fn the_monatomic_limit_itself_is_clean() {
    // 5/3 is the boundary and is inclusive: it is argon and helium, not a mistake.
    let r = choked_flow_area(
        kilograms_per_second(1.0),
        pascals(1.0e6),
        kilograms_per_cubic_meter(10.0),
        5.0 / 3.0,
    )
    .unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

#[test]
fn an_exponent_of_one_is_an_error_rather_than_a_division_by_zero() {
    // `k - 1` is a denominator in the exponent, so k = 1 is singular - and physically
    // it is a substance whose specific heats are equal, which no gas is.
    let err = choked_flow_area(
        kilograms_per_second(1.0),
        pascals(1.0e6),
        kilograms_per_cubic_meter(10.0),
        1.0,
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("k"));
}

#[test]
fn non_positive_pressure_or_density_is_an_error() {
    for (p, r) in [(0.0, 10.0), (1.0e6, 0.0)] {
        let err = choked_flow_area(
            kilograms_per_second(1.0),
            pascals(p),
            kilograms_per_cubic_meter(r),
            1.4,
        )
        .unwrap_err();
        assert!(matches!(err, AzothError::OutOfRange { .. }));
    }
}
