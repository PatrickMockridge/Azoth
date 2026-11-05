//! `process.filter`'s clean-drop models and capture against
//! `validation/neqsim/captures/process_filter_models.tsv`.
//!
//! The interpolation itself is held in `kernels/filter.rs`'s own tests. What is here is the
//! model's wiring: the Ergun bed on a real state, and the beta-ratio capture's three answers.

use azoth_core::units::{
    cubic_meters, kelvins, kilogram, kilogram_per_hour, kilograms, kilograms_per_hour, meters,
    pascal_seconds, pascals, seconds, square_meters,
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
        None,
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

/// The PR methane feed the loading capture uses, with the steady loading inputs.
#[allow(clippy::too_many_arguments)]
fn loaded_filter(
    drop: f64,
    nominal: Option<f64>,
    concentration: Option<f64>,
    solids: Option<f64>,
    capacity: Option<f64>,
    increase: Option<f64>,
    start: Option<f64>,
    cracking: Option<f64>,
    terminal: Option<f64>,
) -> azoth_process::models::FilterResult {
    filter(
        &["methane".to_string()],
        17.314578182246326,
        &[1.0],
        pascals(2_000_000.0),
        kelvins(298.15),
        pascals(drop),
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
        nominal,
        concentration,
        solids.map(kilograms),
        capacity.map(kilograms),
        increase.map(pascals),
        start,
        cracking.map(pascals),
        terminal.map(pascals),
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    )
    .expect("the captured state filters")
}

/// One `runTransient` step on the PR methane feed; the active rates are the ones stated.
#[allow(clippy::too_many_arguments)]
fn dynamic_filter(
    drop: f64,
    solids: f64,
    capacity: f64,
    increase: f64,
    start: f64,
    increment: f64,
    holdup: Option<f64>,
    rate: Option<f64>,
    nominal: Option<f64>,
    concentration: Option<f64>,
    backwash: Option<f64>,
    regeneration: Option<f64>,
    cracking: Option<f64>,
) -> azoth_process::models::FilterResult {
    filter(
        &["methane".to_string()],
        17.314578182246326,
        &[1.0],
        pascals(2_000_000.0),
        kelvins(298.15),
        pascals(drop),
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
        nominal,
        concentration,
        Some(kilograms(solids)),
        Some(kilograms(capacity)),
        Some(pascals(increase)),
        Some(start),
        cracking.map(pascals),
        None,
        Some(seconds(increment)),
        holdup.map(cubic_meters),
        rate.map(kilograms_per_hour),
        backwash.map(kilograms_per_hour),
        regeneration.map(kilograms_per_hour),
        backwash.map(|_| true),
        regeneration.map(|_| true),
    )
    .expect("the captured state filters")
}

#[test]
fn the_loaded_element_matches_the_capture() {
    let result = loaded_filter(
        20_000.0,
        None,
        None,
        Some(9.0),
        Some(10.0),
        Some(100_000.0),
        Some(0.5),
        None,
        None,
    );
    close(
        result.unrestricted_drop.value,
        110_000.0,
        1.0e-12,
        "unrestricted_drop",
    );
    close(
        result.applied_drop.value,
        110_000.0,
        1.0e-12,
        "applied_drop",
    );
    close(result.loading_fraction, 0.9, 1.0e-12, "loading_fraction");
    close(
        result.breakthrough_fraction,
        0.8,
        1.0e-12,
        "breakthrough_fraction",
    );
    close(
        result.differential_pressure_utilization,
        0.22,
        1.0e-12,
        "differential_pressure_utilization",
    );
    assert!(!result.replacement_required);
    close(result.outlet_p.value, 1_890_000.0, 1.0e-12, "outlet_p");
}

#[test]
fn the_bypass_carries_half_the_flow() {
    let result = loaded_filter(
        20_000.0,
        Some(0.99),
        Some(1.0e-4),
        Some(5.0),
        Some(10.0),
        Some(360_000.0),
        Some(0.5),
        Some(50_000.0),
        None,
    );
    close(
        result.unrestricted_drop.value,
        200_000.0,
        1.0e-12,
        "unrestricted_drop",
    );
    close(result.applied_drop.value, 50_000.0, 1.0e-12, "applied_drop");
    close(result.bypass_fraction, 0.5, 1.0e-12, "bypass_fraction");
    close(
        result.current_removal_efficiency,
        0.495,
        1.0e-12,
        "current_removal_efficiency",
    );
    close(
        result.captured_rate.get::<kilogram_per_hour>(),
        0.0495,
        1.0e-9,
        "captured_rate",
    );
    close(
        result.outlet_particle_concentration,
        5.05e-5,
        1.0e-9,
        "outlet_particle_concentration",
    );
    close(result.outlet_p.value, 1_950_000.0, 1.0e-12, "outlet_p");
}

#[test]
fn no_loading_inputs_is_the_class_defaults() {
    let result = loaded_filter(20_000.0, None, None, None, None, None, None, None, None);
    assert_eq!(result.loading_fraction, 0.0);
    assert_eq!(result.breakthrough_fraction, 0.0);
    assert_eq!(result.bypass_fraction, 0.0);
    close(
        result.unrestricted_drop.value,
        20_000.0,
        1.0e-12,
        "unrestricted_drop",
    );
    close(
        result.differential_pressure_utilization,
        0.04,
        1.0e-12,
        "differential_pressure_utilization",
    );
    assert!(!result.replacement_required);
}

#[test]
fn a_non_positive_capacity_is_the_sentinel() {
    let result = loaded_filter(
        20_000.0,
        None,
        None,
        Some(9.0),
        Some(0.0),
        Some(100_000.0),
        Some(0.5),
        None,
        None,
    );
    assert_eq!(result.loading_fraction, 0.0);
    close(
        result.unrestricted_drop.value,
        20_000.0,
        1.0e-12,
        "unrestricted_drop",
    );
}

#[test]
fn replacement_required_reaches_both_branches() {
    let terminal = loaded_filter(
        20_000.0,
        None,
        None,
        Some(5.0),
        Some(10.0),
        Some(360_000.0),
        Some(0.5),
        Some(50_000.0),
        Some(50_000.0),
    );
    close(
        terminal.differential_pressure_utilization,
        4.0,
        1.0e-12,
        "terminal utilization",
    );
    assert!(terminal.replacement_required);

    let loading = loaded_filter(
        20_000.0,
        None,
        None,
        Some(11.0),
        Some(10.0),
        Some(100_000.0),
        Some(0.5),
        None,
        None,
    );
    close(loading.loading_fraction, 1.1, 1.0e-12, "loading fraction");
    assert!(loading.replacement_required);

    let boundary = loaded_filter(
        20_000.0,
        None,
        None,
        Some(5.0),
        Some(10.0),
        Some(100_000.0),
        Some(0.5),
        None,
        None,
    );
    assert!(!boundary.replacement_required);
}

#[test]
fn the_first_step_accumulates_the_generic_rate() {
    let result = dynamic_filter(
        10_000.0,
        0.0,
        10.0,
        100_000.0,
        0.5,
        3600.0,
        Some(1.0),
        Some(5.0),
        None,
        None,
        None,
        None,
        None,
    );
    close(
        result.solids_loading.get::<kilogram>(),
        5.0,
        1.0e-12,
        "loading",
    );
    close(result.loading_fraction, 0.5, 1.0e-12, "fraction");
    close(result.applied_drop.value, 60_000.0, 1.0e-12, "applied");
    close(
        result.holdup_residence_time.expect("a holdup").value,
        48.475_315_062_552_13,
        1.0e-6,
        "residence",
    );
}

#[test]
fn the_second_step_starts_from_the_carried_loading() {
    let result = dynamic_filter(
        10_000.0,
        5.0,
        10.0,
        100_000.0,
        0.5,
        1800.0,
        None,
        Some(5.0),
        None,
        None,
        None,
        None,
        None,
    );
    close(
        result.solids_loading.get::<kilogram>(),
        7.5,
        1.0e-12,
        "loading",
    );
    close(result.breakthrough_fraction, 0.5, 1.0e-12, "breakthrough");
    close(result.applied_drop.value, 85_000.0, 1.0e-12, "applied");
}

#[test]
fn the_backwash_and_regeneration_remove_loading() {
    let result = dynamic_filter(
        20_000.0,
        9.0,
        10.0,
        100_000.0,
        0.5,
        1800.0,
        None,
        Some(0.0),
        None,
        None,
        Some(4.0),
        Some(2.0),
        None,
    );
    close(
        result.solids_loading.get::<kilogram>(),
        6.0,
        1.0e-12,
        "loading",
    );
    close(result.loading_fraction, 0.6, 1.0e-12, "fraction");
    close(result.breakthrough_fraction, 0.2, 1.0e-12, "breakthrough");
    close(result.applied_drop.value, 80_000.0, 1.0e-12, "applied");
}

#[test]
fn the_concentration_accretion_uses_the_pre_step_efficiency() {
    let result = dynamic_filter(
        10_000.0,
        7.9,
        10.0,
        100_000.0,
        0.5,
        28_800.0,
        None,
        None,
        Some(0.99),
        Some(1.0e-4),
        None,
        None,
        None,
    );
    close(
        result.solids_loading.get::<kilogram>(),
        8.232_64,
        1.0e-12,
        "loading",
    );
    close(result.loading_fraction, 0.823_264, 1.0e-12, "fraction");
    close(
        result.breakthrough_fraction,
        0.646_528,
        1.0e-12,
        "breakthrough",
    );
    close(
        result.current_removal_efficiency,
        0.349_937_28,
        1.0e-12,
        "current efficiency",
    );
    close(
        result.captured_rate.get::<kilogram_per_hour>(),
        0.034_993_728,
        1.0e-9,
        "captured",
    );
}

#[test]
fn the_pre_step_bypass_is_in_the_accretion() {
    let result = dynamic_filter(
        20_000.0,
        5.0,
        10.0,
        360_000.0,
        0.8,
        3600.0,
        None,
        None,
        Some(0.99),
        Some(1.0e-4),
        None,
        None,
        Some(50_000.0),
    );
    close(
        result.solids_loading.get::<kilogram>(),
        5.049_5,
        1.0e-12,
        "loading",
    );
    close(
        result.bypass_fraction,
        0.502_212_724_399_558_3,
        1.0e-12,
        "bypass",
    );
    close(
        result.current_removal_efficiency,
        0.492_809_402_844_437_25,
        1.0e-12,
        "current efficiency",
    );
    close(
        result.captured_rate.get::<kilogram_per_hour>(),
        0.049_280_940_284_443_73,
        1.0e-9,
        "captured",
    );
}

#[test]
fn a_zero_flow_has_no_residence_time() {
    let result = filter(
        &["methane".to_string()],
        0.0,
        &[1.0],
        pascals(2_000_000.0),
        kelvins(298.15),
        pascals(10_000.0),
        // 1-11: the model selector, reference flow, exponent, curve, four media, viscosity,
        // the performance curve and the particle size.
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
        // 12-13: the nominal efficiency and the inlet concentration.
        None,
        None,
        // 14-19: loading, capacity, increase, breakthrough start, cracking, terminal.
        None,
        None,
        None,
        None,
        None,
        None,
        // 20-26: the step, holdup, generic rate, backwash, regeneration and the two flags.
        Some(seconds(3600.0)),
        Some(cubic_meters(1.0)),
        Some(kilograms_per_hour(1.0)),
        None,
        None,
        None,
        None,
    )
    .expect("the captured state filters");
    assert!(result.holdup_residence_time.is_none());
}

#[test]
fn both_capture_sources_are_refused() {
    let result = filter(
        &["methane".to_string()],
        17.314578182246326,
        &[1.0],
        pascals(2_000_000.0),
        kelvins(298.15),
        pascals(10_000.0),
        // 1-11: the model selector, reference flow, exponent, curve, four media, viscosity,
        // the performance curve and the particle size.
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
        // 12-13: the nominal efficiency and the inlet concentration.
        None,
        Some(1.0e-4),
        // 14-19: loading, capacity, increase, breakthrough start, cracking, terminal.
        None,
        None,
        None,
        None,
        None,
        None,
        // 20-26: the step, holdup, generic rate, backwash, regeneration and the two flags.
        None,
        None,
        Some(kilograms_per_hour(1.0)),
        None,
        None,
        None,
        None,
    );
    assert!(result.is_err());
}
