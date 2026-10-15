//! Spec-driven tests for `eos.iapws_saturation_pressure`.
//!
//! The oracle is `validation/neqsim/captures/iapws_if97_probe.tsv`, which carries both
//! saturation directions: `psat_t(T)` at five temperatures and `tsat_p(p)` at six pressures.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::iapws_if97;
use azoth_eos::iapws_saturation_pressure;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.iapws_saturation_pressure";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::IapwsSaturationPressureResult {
    iapws_saturation_pressure(kelvins(common::input(case, "T")))
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
                    result.p_sat.value,
                    common::expected(case, "p_sat"),
                    case.tolerance,
                    &format!("{}::{} (p_sat)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "T" => Some(common::input(case, "T")),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
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
        executed >= 4,
        "expected several active cases, ran {executed}"
    );
}

/// **The curve's own domain, refused at both ends and not extrapolated to.**
///
/// Below the triple point there is no liquid water to saturate; at the critical point the two
/// phases meet and the curve ends. Both are `error` in the spec, because outside them the
/// question is a different one rather than a bad number - so the two ends are a test and not
/// a note. 647.096 K *is* the end, so it is refused too.
#[test]
fn the_domain_is_refused_at_both_ends() {
    for refused in [0.0, 200.0, 273.0, 647.096, 700.0, 1000.0] {
        assert!(
            iapws_saturation_pressure(kelvins(refused)).is_err(),
            "{refused} K is outside the Region 4 curve's [273.15, 647.096) K and must be refused"
        );
    }
    for accepted in [273.15, 300.0, 647.0] {
        assert!(
            iapws_saturation_pressure(kelvins(accepted)).is_ok(),
            "{accepted} K is inside the domain and must be answered"
        );
    }
}

/// The pressure rises with the temperature, strictly, over the whole curve.
///
/// Three captured points cannot show this: a coefficient wrong in its curvature can be right
/// at three temperatures and wrong between them. The Clausius-Clapeyron slope is positive
/// everywhere on a vapour-pressure curve, so this is the property the equation has to have.
#[test]
fn the_pressure_rises_with_the_temperature() {
    let at = |t: f64| iapws_saturation_pressure(kelvins(t)).unwrap().p_sat.value;

    let mut previous = at(273.15);
    for step in 1..=200 {
        let t = 273.15 + (647.0 - 273.15) * f64::from(step) / 200.0;
        let current = at(t);
        assert!(
            current > previous,
            "the saturation pressure must rise with T: at {t} K it is {current} Pa, \
             not above {previous} Pa at the step below"
        );
        previous = current;
    }
}

/// **The two directions agree, which is what pins the inverse the spec does not expose.**
///
/// `p4_T` and `T4_p` are separate fits in the release and not one equation inverted, so their
/// composition being the identity is a fact about the fits rather than a tautology - and it is
/// the only check this tree has on `tsat_p`, whose values are in the same capture and whose id
/// is owed. Measured, they agree to `9e-14` relative at the worst of the five pressures both
/// are defined at, which is why the bound is `1e-12`.
///
/// **The capture's sixth pressure is not round-tripped, and that is a finding.** `0.000611`
/// MPa is 611 Pa, which is **below the backward equation's own lower limit of 611.213 Pa** -
/// the release states 273.15 K as that fit's floor - so `tsat_p(611 Pa)` is an extrapolation
/// of it and answers 273.145211 K, five thousandths of a kelvin *under* the forward fit's own
/// floor of 273.15 K. The forward equation then refuses it, correctly. The two fits are
/// consistent wherever both are defined; it is the capture's grid that reaches past one of
/// them, and the value is kept there because it is what the grid asked for.
#[test]
fn the_two_directions_round_trip() {
    for pressure in [0.000611213, 0.000611657, 0.101325, 1.0, 5.0, 10.0, 22.064] {
        let temperature = iapws_if97::saturation_temperature(pressure);
        let back = iapws_saturation_pressure(kelvins(temperature))
            .expect("the backward fit lands inside the forward fit's domain")
            .p_sat
            .value;
        let relative = (back / (pressure * 1.0e6) - 1.0).abs();
        assert!(
            relative <= 1.0e-12,
            "psat(tsat({pressure} MPa)) = {} MPa is {relative:e} from it, over the 1e-12 the \
             two fits hold to",
            back / 1.0e6
        );
    }
    // And the pair disagrees with the module's own MPa scale by exactly the factor the model
    // applies, so the conversion cannot drift away from the module it wraps.
    let direct = iapws_if97::saturation_pressure(373.15);
    let through = iapws_saturation_pressure(kelvins(373.15))
        .unwrap()
        .p_sat
        .value;
    assert!(
        (through - direct * 1.0e6).abs() <= pascals(direct * 1.0e6).value * 1.0e-15,
        "the model's Pa and the module's MPa are not the same number: {through} against {}",
        direct * 1.0e6
    );
}
