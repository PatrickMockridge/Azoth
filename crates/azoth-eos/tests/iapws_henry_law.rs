//! Spec-driven tests for `eos.iapws_henry_law`.

use azoth_core::units::kelvins;
use azoth_eos::spec_gen;
use azoth_eos::{Gas, HenryStatus, gas_from_name, iapws_henry_law};
use azoth_test_support as common;

const CALC_ID: &str = "eos.iapws_henry_law";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::IapwsHenryLawResult {
    let gas: Gas = common::input_str(case, "gas")
        .parse()
        .expect("the case names a row");
    iapws_henry_law(gas, kelvins(common::input(case, "T")))
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
                    result.henry.value,
                    common::expected(case, "henry"),
                    case.tolerance,
                    &format!("{}::{} (henry)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "gas" => None,
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

/// The reported derivative is the equation's, checked against the constant it describes.
///
/// A dropped or doubled factor in the logarithmic derivative would not show in any of the
/// captured constants: they pin `kH`, and this pins `d(ln kH)/dT` by differencing `kH`
/// itself. The two expressions share the saturation series and differ in every term after
/// it, so this is the only test here that reads the second half of the method.
#[test]
fn the_derivative_is_the_logarithmic_derivative_of_the_constant() {
    let h = 1.0e-3;
    for (gas, temperature) in [
        (Gas::Ch4, 298.15),
        (Gas::N2, 298.15),
        (Gas::Co2, 300.0),
        (Gas::H2s, 293.15),
    ] {
        let middle = iapws_henry_law(gas, kelvins(temperature)).expect("in domain");
        let below = iapws_henry_law(gas, kelvins(temperature - h)).expect("in domain");
        let above = iapws_henry_law(gas, kelvins(temperature + h)).expect("in domain");
        let difference = (above.henry.value.ln() - below.henry.value.ln()) / (2.0 * h);
        assert!(
            (difference - middle.d_ln_henry_d_t).abs() < 1.0e-09,
            "{gas:?} at {temperature} K: the equation gives {}, differencing gives {difference}",
            middle.d_ln_henry_d_t,
        );
        // And `ln_henry` is the logarithm of the constant it is reported beside.
        assert!((middle.ln_henry - middle.henry.value.ln()).abs() < 1.0e-12);
    }
}

/// **The fitted window is returned, and the domain is refused.** The two are different
/// questions and NeqSim answers them in different places.
///
/// `ComponentGE.getEffectiveHenryCoefficient` asks `isUsable` and turns a `false` into the
/// insoluble limit - so a methane phase at 275 K never uses the 275.46 K row's own number.
/// That decision is the phase's, and this model returns the number and says which range it
/// came from; what it refuses is a temperature where the guideline is not water at all,
/// which is `getHenryCoefficientBar`'s domain test.
#[test]
fn the_fitted_window_is_reported_and_the_domain_is_refused() {
    // Methane's row is fitted over 275.46-633.11 K. Below the window, inside liquid water:
    let cold = iapws_henry_law(Gas::Ch4, kelvins(273.20)).expect("liquid water");
    assert_eq!(cold.status, HenryStatus::GuidelineExtrapolation);
    assert!(
        cold.henry.value > 0.0,
        "the number is returned, not withheld"
    );
    assert!(
        (cold.rms_log_residual - 0.038_6).abs() < 1e-12,
        "the row's fit residual travels with the row"
    );

    // Nitrogen's starts at 278.12 K, so the same temperature is inside methane's window
    // and outside nitrogen's - the window is a property of the row.
    let warm = iapws_henry_law(Gas::N2, kelvins(298.15)).expect("liquid water");
    assert_eq!(warm.status, HenryStatus::WithinFittedRange);
    assert_eq!(
        iapws_henry_law(Gas::N2, kelvins(277.0))
            .expect("liquid water")
            .status,
        HenryStatus::GuidelineExtrapolation
    );

    // The domain's two edges: the triple point is in, the critical temperature is out.
    assert!(iapws_henry_law(Gas::Ch4, kelvins(273.15)).is_ok());
    assert!(iapws_henry_law(Gas::Ch4, kelvins(647.096)).is_err());
    assert!(iapws_henry_law(Gas::Ch4, kelvins(273.14)).is_err());
}

/// The row a name means, which is NeqSim's `findGas`.
///
/// A row answers to its formula and to its name; **the table's vocabulary is not this
/// library's component list**, and the three rows that are not components here are the
/// reason the resolution cannot be a databank lookup.
#[test]
fn a_row_answers_to_its_formula_and_to_its_name() {
    assert_eq!(gas_from_name("CH4"), Some(Gas::Ch4));
    assert_eq!(gas_from_name("methane"), Some(Gas::Ch4));
    assert_eq!(gas_from_name("  Methane "), Some(Gas::Ch4));
    assert_eq!(gas_from_name("h2s"), Some(Gas::H2s));
    assert_eq!(gas_from_name("hydrogen sulphide"), Some(Gas::H2s));
    assert_eq!(gas_from_name("hydrogen sulfide"), Some(Gas::H2s));
    assert_eq!(gas_from_name("sf6"), Some(Gas::Sf6));

    // Rows with no component behind them in this library's data.
    assert_eq!(gas_from_name("krypton"), Some(Gas::Kr));
    assert_eq!(gas_from_name("xenon"), Some(Gas::Xe));
    assert_eq!(gas_from_name("carbon monoxide"), Some(Gas::Co));

    // And a substance the guideline does not carry.
    assert_eq!(gas_from_name("propane"), None);
    assert_eq!(gas_from_name("water"), None);

    // Every row is reachable by the spelling the spec names it with.
    for gas in [
        Gas::He,
        Gas::Ne,
        Gas::Ar,
        Gas::Kr,
        Gas::Xe,
        Gas::H2,
        Gas::N2,
        Gas::O2,
        Gas::Co,
        Gas::Co2,
        Gas::H2s,
        Gas::Ch4,
        Gas::C2h6,
        Gas::Sf6,
    ] {
        assert_eq!(gas_from_name(gas.as_str()), Some(gas));
        assert_eq!(gas.as_str().parse::<Gas>(), Ok(gas));
    }
}
