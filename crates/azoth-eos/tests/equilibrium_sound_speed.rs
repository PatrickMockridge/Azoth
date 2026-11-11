//! `eos.equilibrium_sound_speed` against `validation/neqsim/captures/process_equilibrium_sound_speed.tsv`.

use azoth_core::units::{joules_per_kilogram_kelvin, kelvins, pascals};
use azoth_eos::equilibrium_sound_speed::equilibrium_sound_speed;

fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= tolerance,
        "{what}: {actual} against {expected}, relative error {error:.3e} > {tolerance:.3e}"
    );
}

#[test]
fn the_dense_single_phase_agrees_with_the_class() {
    let result = equilibrium_sound_speed(
        &["CO2".to_string()],
        kelvins(298.15),
        pascals(15_000_000.0),
        &[1.0],
        None,
        None,
    )
    .expect("the state runs");

    assert_eq!(result.status, "converged");
    assert_eq!(result.stencil, "central");
    let speed = result.sound_speed.expect("a converged row has a speed");
    // The capture's `dense_co2` row; the class's own frozen phase value is the same number.
    close(speed, 456.22605621637035, 2.0e-3, "the dense speed");
    close(
        result.specific_entropy,
        -1537.9517199354361,
        1.0e-4,
        "the target entropy",
    );
    assert_eq!(
        result.samples.len(),
        3,
        "a central stencil carries three states"
    );
}

#[test]
fn the_pure_two_phase_state_is_a_seventh_of_the_phase_speed() {
    let result = equilibrium_sound_speed(
        &["CO2".to_string()],
        kelvins(278.5666675376772),
        pascals(4_000_000.0),
        &[1.0],
        None,
        Some(joules_per_kilogram_kelvin(-1537.9517199354361)),
    )
    .expect("the state runs");

    assert_eq!(result.status, "converged");
    let speed = result.sound_speed.expect("a converged row has a speed");
    close(speed, 57.08170902262859, 2.0e-3, "the equilibrium speed");
    // The class's frozen phase value on the same state is 408.758 m/s.
    assert!(speed < 408.758 / 5.0, "equilibrium is not the phase speed");
    for sample in &result.samples {
        assert_eq!(sample.phase.as_str(), "two_phase");
        assert!(sample.vapour_fraction.is_some());
    }
}

#[test]
fn a_step_outside_the_class_s_band_is_refused() {
    let error = equilibrium_sound_speed(
        &["CO2".to_string()],
        kelvins(298.15),
        pascals(15_000_000.0),
        &[1.0],
        Some(0.2),
        None,
    )
    .expect_err("the class refuses a relative step above 0.05");
    assert_eq!(error.field(), Some("relative_pressure_step"));
}
