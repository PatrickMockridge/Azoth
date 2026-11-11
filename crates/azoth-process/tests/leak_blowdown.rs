//! `process.leak_blowdown` against the states its cases state.

use azoth_core::units::{kelvins, meters, pascals, seconds};
use azoth_process::{LeakBlowdownResult, leak_blowdown};

#[allow(clippy::too_many_arguments)] // one argument per input the helper states
fn run(
    components: &[&str],
    z: &[f64],
    t: f64,
    p: f64,
    duration: f64,
    volume: f64,
    diameter: f64,
    coefficient: f64,
    back: f64,
    step: f64,
) -> LeakBlowdownResult {
    let components: Vec<String> = components.iter().map(|name| (*name).to_string()).collect();
    leak_blowdown(
        &components,
        z,
        kelvins(t),
        pascals(p),
        seconds(duration),
        Some(azoth_core::units::cubic_meters(volume)),
        Some(meters(diameter)),
        Some(coefficient),
        Some(pascals(back)),
        Some(seconds(step)),
        Some("Test"),
    )
    .expect("the state runs")
}

#[test]
fn the_first_point_is_the_class_s_scalar_oracle() {
    let result = run(
        &["methane"],
        &[1.0],
        300.0,
        5_000_000.0,
        10.0,
        1.0,
        0.01,
        0.62,
        101_325.0,
        1.0,
    );
    assert_eq!(result.status, "valid");
    assert_eq!(result.recorded_points, 11);
    let first = &result.points[0];
    let error =
        (first.mass_flow_rate.value - 0.430_505_623_641_902_6).abs() / 0.430_505_623_641_902_6;
    assert!(error <= 2.0e-5, "the first rate differs by {error:.3e}");
    assert_eq!(first.time.value, 0.0);
    assert_eq!(first.temperature.value, 300.0);
    assert_eq!(first.pressure.value, 5_000_000.0);
    assert_eq!(
        result.peak_mass_flow_rate.expect("a peak").value,
        first.mass_flow_rate.value
    );
}

#[test]
fn the_update_uses_the_initial_mass_denominator() {
    let result = run(
        &["methane"],
        &[1.0],
        300.0,
        5_000_000.0,
        2.0,
        1.0,
        0.01,
        0.62,
        101_325.0,
        1.0,
    );
    // The captured second state: `(mass/initialMass)^gamma`, not `(m2/m1)^gamma`.
    let second = &result.points[1];
    assert!((second.pressure.value - 4_910_169.899_738_68).abs() / 4_910_169.899_738_68 <= 2.0e-5);
    assert!(second.pressure.value < result.points[0].pressure.value);
}

#[test]
fn the_pressure_floor_is_a_zero_flow_plateau() {
    let result = run(
        &["methane"],
        &[1.0],
        300.0,
        5_000_000.0,
        60.0,
        1.0,
        0.01,
        0.62,
        101_325.0,
        1.0,
    );
    assert_eq!(result.recorded_points, 61);
    let last = &result.points[60];
    assert_eq!(last.pressure.value, 101_325.0);
    assert_eq!(last.mass_flow_rate.value, 0.0);
    assert!(
        result
            .diagnostic_codes
            .iter()
            .any(|code| code == "PRESSURE_FLOOR_REACHED")
    );
}

#[test]
fn a_single_step_depletes_the_vessel_and_reports_the_tail() {
    let result = run(
        &["methane"],
        &[1.0],
        300.0,
        5_000_000.0,
        100.0,
        1.0,
        0.01,
        0.62,
        101_325.0,
        100.0,
    );
    assert_eq!(result.number_of_points, 2);
    assert_eq!(result.recorded_points, 1);
    let total = result.total_mass_released.expect("a total").value;
    assert!(total > 35.0, "the unclamped total is {total} kg");
    assert!(
        result
            .diagnostic_codes
            .iter()
            .any(|code| code == "MASS_FLOOR_REACHED")
    );
}

#[test]
fn a_single_phase_liquid_reports_no_vapour_and_no_droplet() {
    let result = run(
        &["propane"],
        &[1.0],
        300.0,
        1_100_000.0,
        10.0,
        0.01,
        0.01,
        0.62,
        101_325.0,
        1.0,
    );
    assert_eq!(result.status, "valid");
    assert_eq!(result.points[0].vapour_mole_fraction, 0.0);
    assert_eq!(result.points[0].droplet_smd.value, 0.0);
}

#[test]
fn no_forward_flow_records_a_negative_momentum() {
    let result = run(
        &["methane"],
        &[1.0],
        300.0,
        500_000.0,
        3.0,
        1.0,
        0.01,
        0.62,
        600_000.0,
        1.0,
    );
    assert_eq!(result.total_mass_released.expect("a total").value, 0.0);
    assert!(result.points[0].jet_momentum.value < 0.0);
    assert!(
        result
            .diagnostic_codes
            .iter()
            .any(|code| code == "NO_FORWARD_FLOW")
    );
}

#[test]
fn a_two_phase_initial_state_carries_the_weighted_gamma_and_a_droplet() {
    let result = run(
        &["methane", "ethane"],
        &[0.5, 0.5],
        200.0,
        1_000_000.0,
        1.0,
        1.0,
        0.01,
        0.62,
        100_000.0,
        1.0,
    );
    assert_eq!(result.status, "valid");
    let first = &result.points[0];
    assert!(
        first.vapour_mole_fraction > 0.0 && first.vapour_mole_fraction < 1.0,
        "beta={}",
        first.vapour_mole_fraction
    );
    assert!(first.droplet_smd.value > 0.0);
}

#[test]
fn a_bad_composition_or_step_is_refused() {
    let components = vec!["methane".to_string()];
    let error = leak_blowdown(
        &components,
        &[0.5, 0.5],
        kelvins(300.0),
        pascals(5_000_000.0),
        seconds(1.0),
        Some(azoth_core::units::cubic_meters(1.0)),
        Some(meters(0.01)),
        Some(0.62),
        Some(pascals(101_325.0)),
        Some(seconds(1.0)),
        None,
    )
    .expect_err("two fractions for one component");
    assert_eq!(error.field(), Some("initial_z"));

    let error = leak_blowdown(
        &components,
        &[1.0],
        kelvins(300.0),
        pascals(5_000_000.0),
        seconds(1.0),
        Some(azoth_core::units::cubic_meters(1.0)),
        Some(meters(0.01)),
        Some(0.62),
        Some(pascals(101_325.0)),
        Some(seconds(0.0)),
        None,
    )
    .expect_err("a zero step");
    assert_eq!(error.field(), Some("time_step"));
}
