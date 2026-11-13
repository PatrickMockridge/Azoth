//! `process.gas_dispersion` against `validation/neqsim/captures/process_gas_dispersion.tsv`.

use azoth_core::units::{kelvins, kilograms_per_second, meters, meters_per_second, pascals};
use azoth_process::models::gas_dispersion::{GasDispersionResult, gas_dispersion};

fn names(list: &[&str]) -> Vec<String> {
    list.iter().map(|name| (*name).to_string()).collect()
}

fn close(actual: f64, expected: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= 1.0e-6,
        "{what}: {actual} against {expected}, relative error {error:.3e}"
    );
}

fn standard() -> GasDispersionResult {
    gas_dispersion(
        &names(&["methane"]),
        &[1.0],
        kilograms_per_second(1.0),
        Some(kelvins(288.15)),
        Some(pascals(101_325.0)),
        Some(meters_per_second(5.0)),
        Some("D"),
        Some(meters(0.1)),
        Some(false),
        Some(meters(0.0)),
        None,
        None,
        None,
        None,
        None,
    )
    .expect("the captured state is valid")
}

#[test]
fn the_methane_plume_reproduces_the_capture() {
    let result = standard();
    assert_eq!(result.selected_model, "gaussian_plume");
    close(
        result.source_density.value,
        0.680_058_687_200_945_9,
        "source density",
    );
    close(
        result.air_density.value,
        1.225_012_265_990_694_6,
        "air density",
    );
    close(result.fuel_mole_fraction, 1.0, "fuel mole fraction");
    close(
        result.lower_flammable_limit.expect("an LFL"),
        0.044,
        "the LFL",
    );
    close(
        result
            .distance_to_flammable_endpoint
            .expect("an endpoint")
            .value,
        30.074_944_990_742_615,
        "the endpoint distance",
    );
    close(
        result.distance_to_lfl.expect("an LFL distance").value,
        21.291_430_521_781_088,
        "the LFL distance",
    );
    close(
        result.flammable_cloud_volume.value,
        31.641_537_822_943_3,
        "the cloud volume",
    );
    assert!(result.has_flammable_cloud);
    assert!(!result.has_toxic_endpoint);
    assert_eq!(result.screening_basis, "gaussian_plume_briggs");
}

#[test]
fn carbon_dioxide_takes_the_dense_model_and_says_it_has_no_fuel() {
    let result = gas_dispersion(
        &names(&["CO2"]),
        &[1.0],
        kilograms_per_second(5.0),
        Some(kelvins(288.15)),
        Some(pascals(101_325.0)),
        Some(meters_per_second(5.0)),
        Some("D"),
        Some(meters(0.1)),
        Some(false),
        Some(meters(0.0)),
        None,
        None,
        None,
        None,
        None,
    )
    .expect("the captured state is valid");
    assert_eq!(result.selected_model, "heavy_gas");
    assert!(result.lower_flammable_limit.is_none());
    assert!(result.distance_to_lfl.is_none());
    assert!(result.distance_to_half_lfl.is_none());
    assert!(result.distance_to_flammable_endpoint.is_none());
    assert_eq!(result.flammable_cloud_volume.value, 0.0);
    assert!(!result.has_flammable_cloud);
    assert_eq!(result.screening_basis, "heavy_gas_britter_mcquaid");
}

#[test]
fn the_sour_gas_toxic_endpoint_is_its_own_scan() {
    let result = gas_dispersion(
        &names(&["methane", "H2S"]),
        &[0.98, 0.02],
        kilograms_per_second(2.0),
        Some(kelvins(288.15)),
        Some(pascals(101_325.0)),
        Some(meters_per_second(5.0)),
        Some("D"),
        Some(meters(0.1)),
        Some(false),
        Some(meters(0.0)),
        None,
        None,
        None,
        Some("H2S"),
        Some(100.0e-6),
    )
    .expect("the captured state is valid");
    close(
        result.lower_flammable_limit.expect("an LFL"),
        0.043_912_175_648_702_596,
        "Le Chatelier's mixture LFL",
    );
    close(
        result.toxic_distance.expect("a toxic distance").value,
        90.610_373_323_904_55,
        "the toxic distance",
    );
    assert_eq!(result.toxic_component, "H2S");
    assert_eq!(result.toxic_threshold, Some(100.0e-6));
    assert!(result.has_toxic_endpoint);
}

#[test]
fn a_manual_lfl_makes_the_whole_flow_fuel() {
    let result = gas_dispersion(
        &names(&["methane"]),
        &[1.0],
        kilograms_per_second(1.0),
        Some(kelvins(288.15)),
        Some(pascals(101_325.0)),
        Some(meters_per_second(5.0)),
        Some("D"),
        Some(meters(0.1)),
        Some(false),
        Some(meters(0.0)),
        None,
        Some(0.05),
        None,
        None,
        None,
    )
    .expect("the captured state is valid");
    close(
        result.lower_flammable_limit.expect("the override"),
        0.05,
        "the override",
    );
    close(result.fuel_mole_fraction, 1.0, "the whole flow is fuel");
    close(
        result.distance_to_lfl.expect("an LFL distance").value,
        19.813_106_617_839_95,
        "the LFL distance",
    );
    close(
        result.flammable_cloud_volume.value,
        26.096_525_115_661_36,
        "the cloud volume",
    );
}

#[test]
fn a_forced_model_runs_against_the_auto_selection() {
    let heavy = gas_dispersion(
        &names(&["methane"]),
        &[1.0],
        kilograms_per_second(1.0),
        Some(kelvins(288.15)),
        Some(pascals(101_325.0)),
        Some(meters_per_second(5.0)),
        Some("D"),
        Some(meters(0.1)),
        Some(false),
        Some(meters(0.0)),
        Some("heavy_gas"),
        None,
        None,
        None,
        None,
    )
    .expect("the captured state is valid");
    assert_eq!(heavy.selected_model, "heavy_gas");
    close(
        heavy.distance_to_lfl.expect("an LFL distance").value,
        9.765_355_077_286_93,
        "the dense LFL distance",
    );

    let passive = gas_dispersion(
        &names(&["CO2"]),
        &[1.0],
        kilograms_per_second(5.0),
        Some(kelvins(288.15)),
        Some(pascals(101_325.0)),
        Some(meters_per_second(5.0)),
        Some("D"),
        Some(meters(0.1)),
        Some(false),
        Some(meters(0.0)),
        Some("gaussian_plume"),
        None,
        None,
        None,
        None,
    )
    .expect("the captured state is valid");
    assert_eq!(passive.selected_model, "gaussian_plume");
    assert!(passive.distance_to_lfl.is_none());
}

#[test]
fn the_benchmark_rows_hold() {
    let propane = gas_dispersion(
        &names(&["propane"]),
        &[1.0],
        kilograms_per_second(1.0),
        Some(kelvins(288.15)),
        Some(pascals(101_325.0)),
        Some(meters_per_second(5.0)),
        Some("D"),
        Some(meters(0.1)),
        Some(false),
        Some(meters(0.0)),
        None,
        None,
        None,
        None,
        None,
    )
    .expect("the captured state is valid");
    assert_eq!(propane.selected_model, "heavy_gas");
    close(
        propane.distance_to_lfl.expect("an LFL distance").value,
        11.290_953_671_181_775,
        "propane's dense LFL distance",
    );

    let co2 = gas_dispersion(
        &names(&["CO2"]),
        &[1.0],
        kilograms_per_second(10.0),
        Some(kelvins(283.15)),
        Some(pascals(101_325.0)),
        Some(meters_per_second(2.0)),
        Some("F"),
        Some(meters(0.1)),
        Some(false),
        Some(meters(0.0)),
        None,
        None,
        None,
        Some("CO2"),
        Some(50_000.0e-6),
    )
    .expect("the captured state is valid");
    close(
        co2.toxic_distance.expect("a toxic distance").value,
        19.186_450_541_681_06,
        "CO2's toxic distance",
    );
}

#[test]
fn the_source_term_rate_is_the_input() {
    let result = gas_dispersion(
        &names(&["methane"]),
        &[1.0],
        kilograms_per_second(0.437_386_632_671_170_23),
        Some(kelvins(288.15)),
        Some(pascals(101_325.0)),
        Some(meters_per_second(5.0)),
        Some("D"),
        Some(meters(0.1)),
        Some(false),
        Some(meters(0.0)),
        None,
        None,
        None,
        None,
        None,
    )
    .expect("the captured state is valid");
    close(
        result.mass_release_rate.value,
        0.437_386_632_671_170_23,
        "the peak rate",
    );
    close(
        result.distance_to_lfl.expect("an LFL distance").value,
        13.927_835_165_078_967,
        "the LFL distance",
    );
}

#[test]
fn the_guards_are_the_class_and_the_spec() {
    let zero = gas_dispersion(
        &names(&["methane"]),
        &[1.0],
        kilograms_per_second(0.0),
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    );
    assert!(zero.is_err(), "a non-positive rate is refused");

    let override_zero = gas_dispersion(
        &names(&["methane"]),
        &[1.0],
        kilograms_per_second(1.0),
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        Some(0.0),
        None,
        None,
        None,
    );
    assert!(override_zero.is_err(), "a stated LFL at zero is refused");

    let fraction = gas_dispersion(
        &names(&["methane"]),
        &[1.0],
        kilograms_per_second(1.0),
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        Some(1.5),
        None,
        None,
    );
    assert!(
        fraction.is_err(),
        "an endpoint fraction above one is refused"
    );

    let model = gas_dispersion(
        &names(&["methane"]),
        &[1.0],
        kilograms_per_second(1.0),
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        Some("plume"),
        None,
        None,
        None,
        None,
    );
    assert!(model.is_err(), "an unknown model word is refused");

    let cold = gas_dispersion(
        &names(&["methane"]),
        &[1.0],
        kilograms_per_second(1.0),
        Some(kelvins(0.0)),
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    );
    assert!(
        cold.is_err(),
        "a non-positive ambient temperature is refused"
    );
}

#[test]
fn a_stated_toxic_component_without_a_threshold_reports_no_endpoint() {
    let result = gas_dispersion(
        &names(&["methane", "H2S"]),
        &[0.98, 0.02],
        kilograms_per_second(2.0),
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        Some("H2S"),
        None,
    )
    .expect("the state is valid");
    assert_eq!(result.toxic_component, "");
    assert!(result.toxic_threshold.is_none());
    assert!(result.toxic_distance.is_none());
    assert!(!result.has_toxic_endpoint);
}
