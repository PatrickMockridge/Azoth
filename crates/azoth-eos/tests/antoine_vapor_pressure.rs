//! Spec-driven tests for `eos.antoine_vapor_pressure`.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::antoine_vapor_pressure;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.antoine_vapor_pressure";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::AntoineVaporPressureResult {
    antoine_vapor_pressure(
        common::input(case, "A"),
        common::input(case, "B"),
        common::input(case, "C"),
        common::input(case, "D"),
        common::input(case, "E"),
        common::input_str(case, "form").parse().unwrap(),
        kelvins(common::input(case, "Tc")),
        pascals(common::input(case, "Pc")),
        kelvins(common::input(case, "T")),
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
                        "A" => Some(common::input(case, "A")),
                        "B" => Some(common::input(case, "B")),
                        "C" => Some(common::input(case, "C")),
                        "D" => Some(common::input(case, "D")),
                        "E" => Some(common::input(case, "E")),
                        "Tc" => Some(common::input(case, "Tc")),
                        "Pc" => Some(common::input(case, "Pc")),
                        "T" => Some(common::input(case, "T")),
                        "p_sat" => Some(result.p_sat.value),
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
        executed >= 2,
        "expected several active cases, ran {executed}"
    );
}

/// `E` decides the form, and the label alone does not.
///
/// Twenty rows in NeqSim's `COMP.csv` carry DIPPR-101 coefficients under the label
/// `log`, which names the two-term exponential instead. Reading the label alone returns
/// `i-pentane` at 8.3e38 bar where the correlation gives 0.918, so the selection takes
/// the fifth coefficient as well - NeqSim's own rule,
/// `Component.usesDipprVaporPressureCorrelation`.
#[test]
fn a_non_zero_exponent_selects_the_dippr_form() {
    use azoth_eos::AntoineForm::*;
    use azoth_eos::form_from_type;

    // The fix: the same label, two different forms, decided by `E`.
    assert_eq!(form_from_type("log", 2.0), Some(Dippr101));
    assert_eq!(form_from_type("log", 0.0), Some(Exp));
    assert_eq!(form_from_type("exp", 2.0), Some(Dippr101));

    // `pow10` and `pow10KPa` keep precedence: their coefficients are log10-based and
    // would not survive the exponential form, so a non-zero `E` does not outrank them.
    assert_eq!(form_from_type("pow10", 2.0), Some(Pow10));
    assert_eq!(form_from_type("pow10KPa", 2.0), Some(Pow10Kpa));

    // `loglog`/`log10` still fall through to Wagner, which NeqSim has not changed.
    assert_eq!(form_from_type("loglog", 0.0), Some(Wagner));
    assert_eq!(form_from_type("log10", 0.0), Some(Wagner));
}

/// The DIPPR form is `exp(A + B/T + C ln T + D T**E)` in pascals.
///
/// Worth a test of its own because the unit factor is the one thing a reader cannot see
/// from the equation: NeqSim returns this one in pascals already, `exp(...) / 100000` in
/// bar, so unlike the `exp` and `pow10` forms there is no `1e5` factor on either side.
#[test]
fn the_dippr_form_returns_pascals_without_a_factor() {
    use azoth_eos::AntoineForm;
    use azoth_eos::antoine_vapor_pressure;

    // i-pentane, at the state validation/eos/i_pentane_antoine_dippr101_against_neqsim.json
    // states: NeqSim at commit 8922111 returns 0.10329497716461689 bar.
    let r = antoine_vapor_pressure(
        72.35,
        -5010.9,
        -7.883,
        8.979e-06,
        2.0,
        AntoineForm::Dippr101,
        kelvins(460.43),
        pascals(3_381_200.0),
        kelvins(248.15),
    )
    .unwrap();
    assert!(
        (r.p_sat.value - 10329.497716461689_f64).abs() < 1e-9,
        "got {}",
        r.p_sat.value
    );
}

/// **The Wagner form above `Tc` is refused, and it used to be `NaN`.**
///
/// `x = 1 - T/Tc` is negative above the critical temperature, and `x.powf(1.5)` is not
/// real there. Rust's `powf` returns `NaN`, which flowed through the `exp` into `p_sat`
/// and came back as a *result* - so a caller checking `p_sat > 0` saw `false` with nothing
/// to attribute it to. The Python reference was worse in kind rather than degree: `**`
/// promotes a negative base to a **complex**, so it raised `TypeError` where Rust returned
/// a number, and the two kernels disagreed about what the function even is.
///
/// The fix is the third of the policy's options rather than the first: this is a refusal
/// and not a clamp, because there is no saturation pressure above `Tc` to report.
#[test]
fn the_wagner_form_is_refused_above_the_critical_temperature() {
    let tc = 190.56;
    // Below `Tc` it is an ordinary number, so the refusal is about the domain and not
    // about the coefficients.
    let below = antoine_vapor_pressure(
        -7.0,
        1.5,
        -2.0,
        0.5,
        0.0,
        "wagner".parse().unwrap(),
        kelvins(tc),
        pascals(4.5992e6),
        kelvins(150.0),
    )
    .expect("150 K is below the critical temperature");
    assert!(below.p_sat.value.is_finite() && below.p_sat.value > 0.0);

    // And above it, a refusal naming the state.
    for t in [190.57, 200.0, 300.0] {
        let error = antoine_vapor_pressure(
            -7.0,
            1.5,
            -2.0,
            0.5,
            0.0,
            "wagner".parse().unwrap(),
            kelvins(tc),
            pascals(4.5992e6),
            kelvins(t),
        )
        .expect_err("above Tc there is no saturation pressure");
        assert!(
            error.to_string().contains("critical temperature"),
            "{error}"
        );
    }
}
