//! `process.filter`'s clean-drop models and capture against
//! `validation/neqsim/captures/process_filter_models.tsv`.
//!
//! The interpolation itself is held in `kernels/filter.rs`'s own tests. What is here is the
//! model's wiring: the Ergun bed on a real state, and the beta-ratio capture's three answers.

use azoth_core::units::{
    kelvins, kilogram_per_hour, meters, pascal_seconds, pascals, square_meters,
};
use azoth_process::models::filter::{FilterBetaPoint, filter};

fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= tolerance,
        "{what}: {actual} against {expected}, relative error {error:.3e} > {tolerance:.3e}"
    );
}

#[test]
fn the_ergun_bed_matches_the_capture() {
    let result = filter(
        &["methane".to_string()],
        17.314578182246326,
        &[1.0],
        pascals(2_000_000.0),
        kelvins(298.15),
        pascals(10_000.0),
        Some("ergun"),
        None,
        None,
        None,
        Some(square_meters(1.0)),
        Some(meters(1.5)),
        Some(meters(0.002)),
        Some(0.40),
        Some(pascal_seconds(1.1529270713459012e-5)),
        None,
        None,
        None,
        None,
    )
    .expect("the captured state filters");

    close(
        result.clean_drop.value,
        145.7628909180396,
        5.0e-4,
        "clean_drop",
    );
    close(
        result.applied_drop.value,
        145.7628909180396,
        5.0e-4,
        "applied_drop",
    );
    close(
        result.outlet_p.value,
        1_999_854.237109082,
        5.0e-4,
        "outlet_p",
    );
    close(result.outlet_h.value, 518.0600145934741, 5.0e-4, "outlet_h");
}

#[test]
fn the_beta_curve_captures_particles() {
    let curve = [
        FilterBetaPoint {
            particle_size: meters(5.0e-6),
            beta_ratio: 2.0,
        },
        FilterBetaPoint {
            particle_size: meters(1.0e-5),
            beta_ratio: 100.0,
        },
        FilterBetaPoint {
            particle_size: meters(2.0e-5),
            beta_ratio: 1000.0,
        },
    ];
    let result = filter(
        &["methane".to_string()],
        17.314578182246326,
        &[1.0],
        pascals(2_000_000.0),
        kelvins(298.15),
        pascals(1000.0),
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        Some(&curve),
        Some(meters(1.0e-5)),
        None,
        Some(1.0e-4),
    )
    .expect("the captured state filters");

    close(
        result.nominal_removal_efficiency,
        0.99,
        1.0e-12,
        "efficiency",
    );
    close(
        result.outlet_particle_concentration,
        1.0e-6,
        1.0e-9,
        "outlet_particle_concentration",
    );
    close(
        result.captured_rate.get::<kilogram_per_hour>(),
        0.099,
        1.0e-9,
        "captured_rate",
    );
}
