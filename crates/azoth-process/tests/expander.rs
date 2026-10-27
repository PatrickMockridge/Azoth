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
