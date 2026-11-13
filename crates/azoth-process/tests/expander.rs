//! `process.expander`'s stepped branch against `validation/neqsim/captures/process_expander_steps.tsv`.

use azoth_core::units::{kelvins, pascals};
use azoth_process::models::expander::expander;

fn expander_with(steps: usize) -> azoth_process::models::ExpanderResult {
    expander(
        &[
            "methane".to_string(),
            "ethane".to_string(),
            "propane".to_string(),
            "n-butane".to_string(),
        ],
        690.6184184064726,
        &[0.82, 0.10, 0.05, 0.03],
        pascals(9_000_000.0),
        kelvins(298.15),
        pascals(3_000_000.0),
        1.0,
        Some(true),
        Some(0.80),
        Some(steps),
        None,
        None,
        None,
        None,
    )
    .expect("the captured state expands")
}

fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= tolerance,
        "{what}: {actual} against {expected}, relative error {error:.3e} > {tolerance:.3e}"
    );
}

#[test]
fn the_five_step_default_matches_the_capture() {
    let result = expander_with(5);
    close(
        result.outlet_t.value,
        246.26461237189093,
        1.0e-4,
        "outlet_t",
    );
    close(
        result.outlet_h.value,
        -3093.134500487529,
        1.0e-4,
        "outlet_h",
    );
    close(result.power.value, -999915.6921234869, 1.0e-4, "power");
}

#[test]
fn a_zero_step_count_is_one_step() {
    // The class's `Math.max(1, getNumberOfCompressorCalcSteps())`.
    let zero = expander_with(0);
    let one = expander_with(1);
    assert_eq!(zero.outlet_t.value, one.outlet_t.value);
    assert_eq!(zero.outlet_h.value, one.outlet_h.value);
    assert_eq!(zero.power.value, one.power.value);
}

#[test]
fn a_stated_power_reports_its_efficiency() {
    let result = expander(
        &[
            "methane".to_string(),
            "ethane".to_string(),
            "propane".to_string(),
            "n-butane".to_string(),
        ],
        690.618_418_406_472_6,
        &[0.82, 0.10, 0.05, 0.03],
        pascals(9_000_000.0),
        kelvins(298.15),
        pascals(3_000_000.0),
        0.8,
        None,
        None,
        None,
        None,
        Some(azoth_core::units::watts(-985_780.729_484_916_9)),
        None,
        None,
    )
    .expect("the captured state expands");
    close(
        result.outlet_t.value,
        246.508_123_701_218_98,
        5.0e-4,
        "outlet_t",
    );
    close(
        result
            .reported_isentropic_efficiency
            .expect("an efficiency"),
        0.800_000_006_064_420_7,
        5.0e-4,
        "efficiency",
    );
}

#[test]
fn the_output_temperature_solve_finds_two_hundred_and_fifty() {
    let result = expander(
        &[
            "methane".to_string(),
            "ethane".to_string(),
            "propane".to_string(),
            "n-butane".to_string(),
        ],
        690.618_418_406_472_6,
        &[0.82, 0.10, 0.05, 0.03],
        pascals(9_000_000.0),
        kelvins(298.15),
        pascals(3_000_000.0),
        0.8,
        None,
        None,
        None,
        None,
        None,
        Some(kelvins(250.0)),
        None,
    )
    .expect("the captured state expands");
    close(
        result.outlet_t.value,
        250.000_002_665_501_09,
        5.0e-4,
        "outlet_t",
    );
    close(
        result
            .reported_isentropic_efficiency
            .expect("an efficiency"),
        0.635_490_578_922_437_7,
        5.0e-4,
        "efficiency",
    );
}

#[test]
fn both_modes_are_refused() {
    let result = expander(
        &["methane".to_string()],
        1.0,
        &[1.0],
        pascals(6_000_000.0),
        kelvins(320.0),
        pascals(3_000_000.0),
        0.75,
        None,
        None,
        None,
        None,
        Some(azoth_core::units::watts(-1000.0)),
        Some(kelvins(280.0)),
        None,
    );
    assert!(result.is_err());
}

#[test]
fn the_modes_are_the_isentropic_branchs() {
    let result = expander(
        &["methane".to_string()],
        1.0,
        &[1.0],
        pascals(6_000_000.0),
        kelvins(320.0),
        pascals(3_000_000.0),
        0.75,
        Some(true),
        Some(0.75),
        None,
        None,
        Some(azoth_core::units::watts(-1000.0)),
        None,
        None,
    );
    assert!(result.is_err());
}

#[test]
fn the_captured_radial_design_is_reproduced() {
    let result = expander(
        &["methane".to_string(), "n-butane".to_string()],
        1.0,
        &[0.9, 0.1],
        pascals(6_000_000.0),
        kelvins(320.0),
        pascals(3_000_000.0),
        0.75,
        None,
        None,
        None,
        None,
        None,
        None,
        Some(true),
    )
    .expect("the design computes");

    assert_eq!(
        result.design_expander_type.as_deref(),
        Some("radial_inflow")
    );
    assert_eq!(result.design_load_type.as_deref(), Some("brake"));
    assert_eq!(result.design_number_of_stages, Some(1.0));
    let design = result;
    close(
        design.design_wheel_diameter.unwrap().value,
        0.475_642_771_656_904_36,
        1e-6,
        "wheel",
    );
    close(
        design.design_tip_speed.unwrap().value,
        231.299_961_133_479_8,
        1e-6,
        "tip",
    );
    close(
        design.design_shaft_diameter.unwrap().value,
        0.030,
        1e-9,
        "shaft",
    );
    close(
        design.design_rated_speed.unwrap(),
        9287.440_325_802_469,
        1e-6,
        "rated",
    );
    close(
        design.design_max_continuous_speed.unwrap(),
        9751.812_342_092_593,
        1e-6,
        "mcs",
    );
    close(
        design.design_trip_speed.unwrap(),
        10_216.184_358_382_716,
        1e-6,
        "trip",
    );
    close(
        design.design_first_critical_speed.unwrap(),
        4643.720_162_901_234,
        1e-6,
        "critical",
    );
    close(
        design.design_casing_pressure.unwrap().value,
        6.6e6,
        1e-9,
        "casing p",
    );
    close(
        design.design_casing_temperature.unwrap().value,
        350.000_000_000_000_06,
        1e-12,
        "casing t",
    );
    close(
        design.design_casing_wall_thickness.unwrap().value,
        0.020_707_358_792_104_134,
        1e-6,
        "wall",
    );
    close(
        design.design_bearing_span.unwrap().value,
        1.189_106_929_142_261,
        1e-6,
        "span",
    );
    close(
        design.design_wheel_weight.unwrap().value,
        154.802_876_690_851_1,
        1e-6,
        "wheel kg",
    );
    close(
        design.design_shaft_weight.unwrap().value,
        9.897_242_615_258_442,
        1e-6,
        "shaft kg",
    );
    close(
        design.design_generator_weight.unwrap().value,
        0.0,
        1e-12,
        "generator",
    );
    assert_eq!(design.design_seal_type.as_deref(), Some("dry_gas"));
    assert_eq!(design.design_bearing_type.as_deref(), Some("sleeve"));
    close(
        design.design_inner_diameter.unwrap().value,
        0.618_335_603_153_975_7,
        1e-6,
        "id",
    );
    close(
        design.design_outer_diameter.unwrap().value,
        0.659_750_320_738_184,
        1e-6,
        "od",
    );
    close(
        design.design_tantan_length.unwrap().value,
        2.378_213_858_284_522,
        1e-6,
        "tantan",
    );
    close(
        design.design_weight_vessel.unwrap().value,
        1164.170_147_573_361_4,
        1e-6,
        "vessel kg",
    );
    close(
        design.design_weight_internals.unwrap().value,
        211.140_982_313_364_87,
        1e-6,
        "internals kg",
    );
    close(
        design.design_weight_piping.unwrap().value,
        906.517_775_984_658_8,
        1e-6,
        "piping kg",
    );
    close(
        design.design_weight_electro_instrument.unwrap().value,
        50.0,
        1e-12,
        "electro kg",
    );
    close(
        design.design_weight_structural_steel.unwrap().value,
        412.593_338_966_017_94,
        1e-6,
        "structural kg",
    );
    close(
        design.design_weight_total.unwrap().value,
        2744.422_244_837_403,
        1e-6,
        "total kg",
    );
    close(
        design.design_module_length.unwrap().value,
        5.472_767_322_855_653,
        1e-6,
        "module l",
    );
    close(
        design.design_module_width.unwrap().value,
        2.319_500_641_476_367_7,
        1e-6,
        "module w",
    );
    close(
        design.design_module_height.unwrap().value,
        1.659_750_320_738_183_8,
        1e-6,
        "module h",
    );
}

#[test]
fn the_captured_axial_design_selects_the_generator() {
    let result = expander(
        &["methane".to_string(), "n-butane".to_string()],
        250_000.0,
        &[0.9, 0.1],
        pascals(6_000_000.0),
        kelvins(320.0),
        pascals(3_000_000.0),
        0.75,
        None,
        None,
        None,
        None,
        None,
        None,
        Some(true),
    )
    .expect("the design computes");

    assert_eq!(result.design_expander_type.as_deref(), Some("axial"));
    assert_eq!(result.design_load_type.as_deref(), Some("generator"));
    assert_eq!(result.design_number_of_stages, Some(2.0));
    close(
        result.design_wheel_diameter.unwrap().value,
        1.0,
        1e-9,
        "wheel",
    );
    close(
        result.design_shaft_diameter.unwrap().value,
        0.595,
        1e-9,
        "shaft",
    );
    close(
        result.design_generator_weight.unwrap().value,
        278_783.127_062_265_7,
        1e-6,
        "generator kg",
    );
    close(
        result.design_weight_total.unwrap().value,
        572_166.940_002_067_7,
        1e-6,
        "total kg",
    );
}
