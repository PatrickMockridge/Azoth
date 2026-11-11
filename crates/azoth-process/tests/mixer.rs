//! `process.mixer`'s pressure-mismatch flag against
//! `validation/neqsim/captures/process_mixer.tsv`.
//!
//! The flag is a boolean and has no slot in a `TestCase`; the model case pins the spread,
//! and this holds the flag to the capture's three rows.

use azoth_core::units::{kelvins, pascals};
use azoth_process::models::mixer::mixer;

fn components() -> Vec<String> {
    ["n-butane", "n-pentane"]
        .iter()
        .map(|name| (*name).to_string())
        .collect()
}

fn run(first: f64, second: f64, specified: Option<f64>) -> bool {
    mixer(
        &components(),
        &[1.0, 2.0],
        &[vec![0.6, 0.4], vec![0.3, 0.7]],
        &[pascals(first), pascals(second)],
        &[kelvins(300.0), kelvins(300.0)],
        specified.map(pascals),
    )
    .expect("the captured feeds mix")
    .pressure_mismatch
}

#[test]
fn the_pressure_mismatch_flag_is_the_class_tolerance() {
    // 4 bar apart, no specified outlet: flagged.
    assert!(run(1.0e6, 0.6e6, None));
    // The same feeds with 8 bara specified: suppressed.
    assert!(!run(1.0e6, 0.6e6, Some(0.8e6)));
    // 0.4 bar apart, below the class's 0.5 bar: clear.
    assert!(!run(1.0e6, 0.96e6, None));
}
