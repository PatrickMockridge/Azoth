//! `process.probit` against `validation/neqsim/captures/process_probit.tsv`.

use azoth_process::models::probit::probit;

fn close(actual: f64, expected: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= 1.0e-9,
        "{what}: {actual} against {expected}, relative error {error:.3e}"
    );
}

#[test]
fn five_is_half_and_the_class_erf_says_so() {
    let result = probit(5.0, 0.0, 1.0, Some(1.0), None).expect("a dose is stated");
    close(result.probit.expect("a probit"), 5.0, "the probit");
    // Abramowitz-Stegun 7.1.26 answers 1e-9 at zero, not zero.
    close(result.probability, 0.500_000_000_5, "the probability");
}

#[test]
fn a_stated_probit_is_the_probability_argument() {
    let result = probit(5.0, 0.0, 1.0, Some(1.0), Some(5.0)).expect("both are stated");
    close(result.probability, 0.500_000_000_5, "probability(Y)");
}

#[test]
fn the_factory_rows_reproduce_the_capture() {
    let h2s = probit(-31.42, 3.008, 1.43, Some(100.0f64.powf(1.43) * 10.0), None)
        .expect("a dose is stated");
    close(
        h2s.probit.expect("a probit"),
        6.811_106_062_496_066,
        "H2S probit",
    );
    close(h2s.probability, 0.964_937_836_684_984_8, "H2S probability");

    let thermal = probit(
        -14.9,
        2.56,
        1.0,
        Some(60.0 * 5000.0f64.powf(4.0 / 3.0) / 1.0e4),
        None,
    )
    .expect("a dose is stated");
    close(
        thermal.probit.expect("a probit"),
        1.075_070_153_730_305_1,
        "thermal probit",
    );
    close(
        thermal.probability,
        4.339_551_543_125_486e-5,
        "thermal probability",
    );

    let blast = probit(-77.1, 6.91, 1.0, Some(20_000.0), None).expect("a dose is stated");
    close(
        blast.probit.expect("a probit"),
        -8.666_901_011_975_355,
        "blast probit",
    );
    assert_eq!(blast.probability, 0.0, "the CDF saturates at zero");
}

#[test]
fn the_non_positive_dose_is_refused() {
    assert!(probit(5.0, 0.0, 1.0, Some(0.0), None).is_err());
    assert!(probit(5.0, 0.0, 1.0, Some(-1.0), None).is_err());
    assert!(
        probit(5.0, 0.0, 1.0, None, None).is_err(),
        "neither a dose nor a probit value is a state"
    );
}
