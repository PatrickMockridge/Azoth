//! Spec-driven tests for `eos.tbp_fraction_properties`.

use azoth_core::{AzothError, CalcResult};
use azoth_eos::spec_gen;
use azoth_eos::tbp_fraction_properties;
use azoth_test_support as common;

const CALC_ID: &str = "eos.tbp_fraction_properties";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::results::TbpFractionPropertiesResult {
    tbp_fraction_properties(
        common::input(case, "molar_mass"),
        common::input(case, "density"),
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
                for (field, actual, expected) in [
                    ("tc", result.tc.value, common::expected(case, "tc")),
                    ("pc", result.pc.value, common::expected(case, "pc")),
                    (
                        "boiling_temperature",
                        result.boiling_temperature.value,
                        common::expected(case, "boiling_temperature"),
                    ),
                    (
                        "acentric_factor",
                        result.acentric_factor,
                        common::expected(case, "acentric_factor"),
                    ),
                    (
                        "attraction_exponent",
                        result.attraction_exponent,
                        common::expected(case, "attraction_exponent"),
                    ),
                ] {
                    common::assert_close(
                        actual,
                        expected,
                        case.tolerance,
                        &format!("{context} ({field})"),
                    );
                }
                common::assert_consistent(&result, &context);
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "molar_mass" => Some(common::input(case, "molar_mass")),
                        "density" => Some(common::input(case, "density")),
                        "acentric_factor" => Some(result.acentric_factor),
                        _ => None,
                    },
                    &context,
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

/// **Both branches, at the gram-per-mole thresholds NeqSim states them in.**
///
/// The cut's molar mass is stated here in kg/mol and the correlations branch in g/mol, so
/// reading either threshold in the caller's unit is a port that agrees at every light cut and
/// diverges only at heavy ones. The tests either side of each switch are the ones that catch
/// it, and this names what they are catching: `539` and `541` g/mol straddle the boiling
/// point's cubic-to-power-law switch, `1119` and `1121` the coefficient set's.
#[test]
fn the_two_switches_are_at_540_and_1120_grams_per_mole() {
    // The boiling point: the cubic at 539 g/mol, the power law at 541, and they differ.
    let light = tbp_fraction_properties(0.539, 890.0).expect("computes");
    let heavy = tbp_fraction_properties(0.541, 900.0).expect("computes");
    let cubic = 2.0e-6 * 539.0_f64.powi(3) - 0.0035 * 539.0_f64.powi(2) + 2.4003 * 539.0 + 171.74;
    let power_law = 97.58 * 541.0_f64.powf(0.3323) * 0.9_f64.powf(0.04609);
    assert!((light.boiling_temperature.value - cubic).abs() < 1e-9);
    assert!((heavy.boiling_temperature.value - power_law).abs() < 1e-9);

    // The coefficient set: at 1119 g/mol the oil set's critical temperature is 1250.4 K and
    // at 1121 the heavy set's is 994.7 - a step, not a trend.
    let oil = tbp_fraction_properties(1.119, 990.0).expect("computes");
    let set = tbp_fraction_properties(1.121, 1000.0).expect("computes");
    assert!((oil.tc.value - 1250.39775284238).abs() < 1e-6);
    assert!((set.tc.value - 994.658850370848).abs() < 1e-6);
    assert!(set.tc.value < oil.tc.value);

    // And the switch is on the **molar mass alone**, not on which set gave a sane answer:
    // the heavy set's acentric factor here is `-44.6`, and it is returned.
    assert!(set.acentric_factor < -1.0);
}

/// **The unphysical heavy branch is reported rather than passed on.**
///
/// Above `mw = 1120` g/mol the heavy set gives a boiling point above the critical
/// temperature, so `tc/tb - 1` is negative and the acentric factor comes out large and
/// negative - measured, `-44.6358323245519` at `mw = 1121`, `d = 1.0`. The correlation is
/// reproduced rather than repaired, and the output's range is what says so.
#[test]
fn an_unphysical_acentric_factor_is_warned_about() {
    let result = tbp_fraction_properties(1.121, 1000.0).expect("computes");
    assert!(result.boiling_temperature.value > result.tc.value);
    assert!(
        result.has_warning(azoth_core::WarningCode::OutOfValidRange),
        "the warning list is {:?}",
        result.warnings
    );

    // The physical region is clean, so the warning is about the branch and not the calc.
    let light = tbp_fraction_properties(0.5, 880.0).expect("computes");
    assert!(light.is_clean(), "{:?}", light.warnings);
}

/// A cut with no molar mass or no density is refused.
#[test]
fn a_cut_with_no_mass_or_no_density_is_refused() {
    for (molar_mass, density) in [(0.0, 880.0), (0.5, 0.0), (-0.5, 880.0), (0.5, -1.0)] {
        let error = tbp_fraction_properties(molar_mass, density)
            .expect_err("a cut needs both numbers positive");
        assert!(
            matches!(
                error,
                AzothError::OutOfRange { .. } | AzothError::InvalidInput { .. }
            ),
            "mw {molar_mass}, d {density}: {error:?}"
        );
    }
}
