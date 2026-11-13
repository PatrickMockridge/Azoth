//! `process.gaussian_plume` against `validation/neqsim/captures/process_gaussian_plume.tsv`.

use azoth_core::units::{
    kilograms_per_cubic_meter, kilograms_per_second, meters, meters_per_second,
};
use azoth_process::models::gaussian_plume::{GaussianPlumeResult, gaussian_plume};

fn close(actual: f64, expected: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= 1.0e-9,
        "{what}: {actual} against {expected}, relative error {error:.3e}"
    );
}

#[allow(clippy::too_many_arguments)] // one parameter per declared input
fn plume(
    emission_rate: f64,
    stack_height: f64,
    wind_speed: f64,
    stability: &str,
    terrain: &str,
    distance: f64,
    crosswind: f64,
    height: f64,
    threshold: f64,
) -> GaussianPlumeResult {
    gaussian_plume(
        kilograms_per_second(emission_rate),
        meters(stack_height),
        meters_per_second(wind_speed),
        stability,
        terrain,
        meters(distance),
        meters(crosswind),
        meters(height),
        kilograms_per_cubic_meter(threshold),
    )
    .expect("the captured state is valid")
}

#[test]
fn the_rural_d_profile_reproduces_the_capture() {
    let result = plume(1.0, 10.0, 5.0, "D", "rural", 100.0, 0.0, 0.0, 0.001);
    close(
        result.sigma_y.value,
        7.960_297_521_679_914,
        "sigma_y at 100 m",
    );
    close(
        result.sigma_z.value,
        5.595_028_849_441_882,
        "sigma_z at 100 m",
    );
    close(
        result.concentration.value,
        2.893_901_168_599_378_4e-4,
        "concentration at 100 m",
    );
    close(
        result.centerline_ground_concentration.value,
        2.893_901_168_599_378_4e-4,
        "centreline at 100 m",
    );
    // The 1e-3 threshold is never reached on this state.
    assert!(result.distance_to_threshold.is_none());

    let far = plume(1.0, 10.0, 5.0, "D", "rural", 1000.0, 0.0, 0.0, 0.001);
    close(
        far.sigma_y.value,
        76.277_007_139_647_38,
        "sigma_y at 1000 m",
    );
    close(
        far.sigma_z.value,
        37.947_331_922_020_55,
        "sigma_z at 1000 m",
    );
    close(
        far.concentration.value,
        2.124_347_514_885_570_2e-5,
        "concentration at 1000 m",
    );
    assert!(
        far.concentration.value < result.concentration.value,
        "the ground concentration falls with distance"
    );
}

#[test]
fn the_urban_b_coefficients_are_theirs() {
    let result = plume(1.0, 5.0, 3.0, "B", "urban", 500.0, 0.0, 0.0, 1.0e-4);
    close(
        result.sigma_y.value,
        146.059_348_668_044_3,
        "sigma_y, urban B",
    );
    close(
        result.sigma_z.value,
        97.979_589_711_327_14,
        "sigma_z, urban B",
    );
    close(
        result.distance_to_threshold.expect("a crossing").value,
        121.271_019_705_978_8,
        "the crossing distance",
    );
}

#[test]
fn the_two_image_terms_move_an_off_centre_receptor() {
    let result = plume(1.0, 10.0, 5.0, "D", "rural", 500.0, 25.0, 2.0, 1.0e-4);
    close(
        result.concentration.value,
        5.298_641_956_857_312e-5,
        "the off-centre concentration",
    );
    close(
        result.centerline_ground_concentration.value,
        6.525_134_622_135_576e-5,
        "the ground centreline",
    );
    assert!(result.concentration.value < result.centerline_ground_concentration.value);
}

#[test]
fn the_scan_reaches_the_threshold_and_stops_at_its_ends() {
    let crossing = plume(10.0, 0.0, 4.0, "D", "rural", 100.0, 0.0, 0.0, 0.001);
    close(
        crossing.distance_to_threshold.expect("a crossing").value,
        470.876_693_844_951_83,
        "the ground-level crossing",
    );

    let never_above = plume(1.0, 0.0, 5.0, "D", "rural", 100.0, 0.0, 0.0, 1.0e6);
    assert!(never_above.distance_to_threshold.is_none());

    let never_below = plume(10.0, 0.0, 4.0, "D", "rural", 100.0, 0.0, 0.0, 1.0e-30);
    close(
        never_below.distance_to_threshold.expect("the cap").value,
        100_000.0,
        "the scan's upper end",
    );
}

#[test]
fn a_zero_emission_answers_zero() {
    let result = plume(0.0, 0.0, 5.0, "D", "rural", 100.0, 0.0, 0.0, 1.0e-6);
    assert_eq!(result.concentration.value, 0.0);
    assert_eq!(result.centerline_ground_concentration.value, 0.0);
    assert!(result.distance_to_threshold.is_none());
}

#[test]
fn the_guards_are_the_class_constructors() {
    let wind = gaussian_plume(
        kilograms_per_second(1.0),
        meters(0.0),
        meters_per_second(0.0),
        "D",
        "rural",
        meters(100.0),
        meters(0.0),
        meters(0.0),
        kilograms_per_cubic_meter(1.0e-3),
    );
    assert!(wind.is_err(), "a non-positive wind speed is refused");

    let stability = gaussian_plume(
        kilograms_per_second(1.0),
        meters(0.0),
        meters_per_second(5.0),
        "G",
        "rural",
        meters(100.0),
        meters(0.0),
        meters(0.0),
        kilograms_per_cubic_meter(1.0e-3),
    );
    assert!(stability.is_err(), "an unknown stability class is refused");

    let threshold = gaussian_plume(
        kilograms_per_second(1.0),
        meters(0.0),
        meters_per_second(5.0),
        "D",
        "rural",
        meters(100.0),
        meters(0.0),
        meters(0.0),
        kilograms_per_cubic_meter(0.0),
    );
    assert!(threshold.is_err(), "a non-positive threshold is refused");
}
