//! Spec-driven tests for `eos.de_waard_milliams_corrosion`.

use azoth_core::units::{
    kelvins, kilograms_per_square_meter_second, meters_per_second, pascals, seconds,
};
use azoth_eos::de_waard_milliams_corrosion::de_waard_milliams_corrosion;
use azoth_eos::results::DeWaardMilliamsCorrosionResult;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.de_waard_milliams_corrosion";

fn call(case: &azoth_core::spec::TestCase) -> DeWaardMilliamsCorrosionResult {
    let pressure = |name: &str| pascals(common::input(case, name));
    de_waard_milliams_corrosion(
        kelvins(common::input(case, "temperature")),
        pressure("co2_partial_pressure"),
        common::input(case, "ph"),
        pressure("h2s_partial_pressure"),
        meters_per_second(common::input(case, "flow_velocity")),
        common::input(case, "inhibitor_efficiency"),
        common::input(case, "glycol_fraction"),
        kilograms_per_square_meter_second(common::input(case, "sulfur_deposition_rate")),
        seconds(common::input(case, "design_life")),
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
                executed += 1;
                let result = call(case);
                let context = &format!("{}::{}", spec.id, case.id);
                if let Some(expected) = case.expected_string("corrosion_severity") {
                    assert_eq!(
                        result.corrosion_severity.as_str(),
                        expected,
                        "{context}: the severity"
                    );
                }
                if let Some(expected) = case.expected_string("fes_morphology") {
                    assert_eq!(
                        result.fes_morphology.as_str(),
                        expected,
                        "{context}: the morphology"
                    );
                }
                for (field, value) in [
                    ("baseline_rate", result.baseline_rate.value),
                    ("ph_correction_factor", result.ph_correction_factor),
                    ("scale_correction_factor", result.scale_correction_factor),
                    ("glycol_correction_factor", result.glycol_correction_factor),
                    ("flow_correction_factor", result.flow_correction_factor),
                    ("inhibitor_factor", result.inhibitor_factor),
                    ("corrected_rate", result.corrected_rate.value),
                    ("sulfur_rate", result.sulfur_rate.value),
                    ("fes_rate", result.fes_rate.value),
                    ("co2_corrected_rate", result.co2_corrected_rate.value),
                    ("total_corrosion_rate", result.total_corrosion_rate.value),
                    ("feco3_saturation_index", result.feco3_saturation_index),
                    ("corrosion_allowance", result.corrosion_allowance.value),
                ] {
                    if let Some(expected) = case.expected_value(field) {
                        common::assert_close(
                            value,
                            expected,
                            case.tolerance,
                            &format!("{context} ({field})"),
                        );
                    }
                }
                // The two flags are self-asserted: `TestCase` carries no boolean expectations.
                let sour = matches!(
                    case.id,
                    "sour_threshold_above" | "fes_cap" | "morph_50" | "morph_150"
                );
                assert_eq!(result.is_sour_service, sour, "{context}: sour service");
                let scale_forming = matches!(
                    case.id,
                    "ph_above_sat" | "scale_160_ph7" | "fes_cap" | "morph_150" | "feco3_high_ph"
                );
                assert_eq!(
                    result.feco3_scale_forming, scale_forming,
                    "{context}: FeCO3 scale forming"
                );
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "temperature" => Some(common::input(case, "temperature")),
                        "co2_partial_pressure" => Some(common::input(case, "co2_partial_pressure")),
                        "h2s_partial_pressure" => Some(common::input(case, "h2s_partial_pressure")),
                        "flow_velocity" => Some(common::input(case, "flow_velocity")),
                        "design_life" => Some(common::input(case, "design_life")),
                        _ => None,
                    },
                    context,
                );
                common::assert_consistent(&result, context);
            }
            "property" => match case.property {
                Some("unit_round_trip") => {}
                other => panic!("{}::{}: unhandled property {other:?}", spec.id, case.id),
            },
            other => panic!("{}::{} has an unhandled kind `{other}`", spec.id, case.id),
        }
    }
    assert!(executed >= 2, "the spec should exercise more than one row");
}

/// **A temperature below absolute zero is refused** rather than fed to the correlation.
#[test]
fn a_non_positive_temperature_is_refused() {
    let error = de_waard_milliams_corrosion(
        kelvins(0.0),
        pascals(2.0e5),
        4.5,
        pascals(0.0),
        meters_per_second(2.0),
        0.0,
        0.0,
        kilograms_per_square_meter_second(0.0),
        seconds(31_557_600.0),
    )
    .expect_err("zero kelvin is not a state");
    assert!(error.to_string().contains("temperature"));
}
