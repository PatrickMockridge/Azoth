//! Spec-driven tests for `eos.surf_cooldown`.

use azoth_core::units::{
    kelvin_intervals, kelvins, meters, pascals, seconds, watts_per_meter_kelvin,
    watts_per_square_meter_kelvin,
};
use azoth_eos::model_gen;
use azoth_eos::results::SurfCooldownResult;
use azoth_eos::surf_cooldown::surf_cooldown;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.surf_cooldown";

fn call(case: &azoth_core::spec::TestCase) -> SurfCooldownResult {
    let components: Vec<String> = common::list_input(case, "components")
        .iter()
        .map(|name| (*name).to_owned())
        .collect();
    surf_cooldown(
        &components,
        case.vector("z").expect("the case states its composition"),
        pascals(common::input(case, "operating_pressure")),
        kelvins(common::input(case, "flash_temperature")),
        meters(common::input(case, "internal_diameter")),
        meters(common::input(case, "wall_thickness")),
        meters(common::input(case, "insulation_thickness")),
        meters(common::input(case, "coating_thickness")),
        watts_per_meter_kelvin(common::input(case, "insulation_conductivity")),
        kelvins(common::input(case, "seabed_temperature")),
        case.input("operating_temperature").map(kelvins),
        case.input("overall_u_value")
            .map(watts_per_square_meter_kelvin),
        case.input("external_htc")
            .map(watts_per_square_meter_kelvin),
        case.input("hydrate_margin").map(kelvin_intervals),
        case.input("required_no_touch_time").map(seconds),
        case.input("total_time").map(seconds),
        case.input("time_step").map(seconds),
        case.flag("require_water"),
        case.string("eos"),
        case.string("hydrate_model"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        if let Some(expected) = case.expected_string("verdict") {
            assert_eq!(result.verdict.as_str(), expected, "{context}: the verdict");
        }
        for (field, got) in [
            ("water_mole_fraction", Some(result.water_mole_fraction)),
            (
                "no_touch_time",
                result.no_touch_time.map(|value| value.value),
            ),
            (
                "initial_fluid_temperature",
                Some(result.initial_fluid_temperature.value),
            ),
            (
                "hydrate_equilibrium_temperature",
                result
                    .hydrate_equilibrium_temperature
                    .map(|value| value.value),
            ),
            ("hydrate_margin", Some(result.hydrate_margin.value)),
            ("seabed_temperature", Some(result.seabed_temperature.value)),
            ("overall_u_value", Some(result.overall_u_value.value)),
            ("time_constant", Some(result.time_constant.value)),
            ("fluid_density", Some(result.fluid_density.value)),
            (
                "fluid_specific_heat",
                Some(result.fluid_specific_heat.value),
            ),
            ("internal_diameter", Some(result.internal_diameter.value)),
            (
                "insulation_thickness",
                Some(result.insulation_thickness.value),
            ),
        ] {
            match (got, case.expected_value(field)) {
                (Some(got), Some(want)) => {
                    common::assert_close(got, want, case.tolerance, &format!("{context} ({field})"))
                }
                (None, None) => {}
                (got, want) => panic!("{context} ({field}): {got:?} against {want:?}"),
            }
        }
        // The flag is self-asserted: `TestCase` carries no boolean expectations.
        assert_eq!(
            result.water_present,
            case.id != "dry",
            "{context}: water present"
        );
        common::assert_consistent(&result, context);
    }
}

/// **A fluid with no water is refused when the caller requires one.**
#[test]
fn an_absent_water_is_refused_when_required() {
    let components: Vec<String> = ["methane", "ethane", "propane"]
        .iter()
        .map(|name| (*name).to_owned())
        .collect();
    let error = surf_cooldown(
        &components,
        &[0.90, 0.07, 0.03],
        pascals(1.2e7),
        kelvins(338.15),
        meters(0.254),
        meters(0.0159),
        meters(0.06),
        meters(0.006),
        watts_per_meter_kelvin(0.22),
        kelvins(277.15),
        None,
        Some(watts_per_square_meter_kelvin(2.5)),
        None,
        None,
        None,
        None,
        None,
        Some(true),
        Some("srk"),
        Some("pvtsim"),
    )
    .expect_err("a dry fluid is refused when water is required");
    assert!(error.to_string().contains("water"));
}
