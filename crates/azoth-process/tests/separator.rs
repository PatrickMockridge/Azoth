//! `process.separator`'s capacity surface against
//! `validation/neqsim/captures/process_separator_capacity.tsv`.
//!
//! The flash itself is held by `python/tests/models/test_model_cases.py` and the cross-impl
//! comparison; what is here is the second fold on it - the five getters, the two density rules
//! and the states the class answers but this port refuses.

use azoth_core::units::{kelvins, meters, meters_per_second, pascals};
use azoth_process::models::separator::separator;

const METHANE_BUTANE: [&str; 2] = ["methane", "n-butane"];

fn names() -> Vec<String> {
    METHANE_BUTANE
        .iter()
        .map(|name| (*name).to_string())
        .collect()
}

#[allow(clippy::too_many_arguments)] // one parameter per state the rows sweep
fn run(
    components: &[String],
    feed_n: f64,
    z: &[f64],
    p: f64,
    t: f64,
    diameter: Option<f64>,
    load_factor: Option<f64>,
    orientation: Option<&str>,
    level_fraction: Option<f64>,
) -> azoth_process::models::separator::SeparatorResult {
    separator(
        components,
        feed_n,
        z,
        pascals(p),
        kelvins(t),
        pascals(0.0),
        0.0,
        None,
        diameter.map(meters),
        load_factor.map(meters_per_second),
        orientation,
        level_fraction,
    )
    .expect("the captured state separates")
}

fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= tolerance,
        "{what}: {actual} against {expected}, relative error {error:.3e} > {tolerance:.3e}"
    );
}

#[test]
fn the_vertical_capacity_matches_the_capture() {
    let result = run(
        &names(),
        1.0,
        &[0.7, 0.3],
        2_000_000.0,
        300.0,
        Some(1.0),
        Some(0.11),
        Some("vertical"),
        Some(0.8),
    );
    close(
        result.capacity_utilization.expect("utilisation"),
        0.002_126_127_241_504_300_3,
        2.0e-4,
        "capacity_utilization",
    );
    close(
        result
            .max_allowable_gas_velocity
            .expect("max velocity")
            .value,
        0.559_861_417_115_024,
        2.0e-4,
        "max_allowable_gas_velocity",
    );
    close(
        result
            .gas_superficial_velocity
            .expect("superficial velocity")
            .value,
        0.001_190_336_610_395_454_6,
        2.0e-4,
        "gas_superficial_velocity",
    );
    close(
        result.gas_load_factor.expect("load factor").value,
        2.338_739_965_654_731e-4,
        2.0e-4,
        "gas_load_factor",
    );
    assert_eq!(result.is_overloaded, Some(false));
}

/// Issue #3154: the vertical area is the full cross-section, and the horizontal one is a fifth
/// of it. The sizing velocity does not move, because it reads densities and the design `K`.
#[test]
fn the_orientation_derates_the_area_and_not_the_velocity() {
    let vertical = run(
        &names(),
        1.0,
        &[0.7, 0.3],
        2_000_000.0,
        300.0,
        Some(1.0),
        Some(0.11),
        Some("vertical"),
        Some(0.8),
    );
    let horizontal = run(
        &names(),
        1.0,
        &[0.7, 0.3],
        2_000_000.0,
        300.0,
        Some(1.0),
        Some(0.11),
        Some("horizontal"),
        Some(0.8),
    );
    close(
        vertical.max_allowable_gas_velocity.expect("vertical").value,
        horizontal
            .max_allowable_gas_velocity
            .expect("horizontal")
            .value,
        1.0e-12,
        "max velocity is the orientation's to move",
    );
    close(
        horizontal
            .gas_superficial_velocity
            .expect("horizontal")
            .value,
        5.0 * vertical.gas_superficial_velocity.expect("vertical").value,
        1.0e-12,
        "horizontal superficial velocity",
    );
    close(
        horizontal.capacity_utilization.expect("horizontal"),
        5.0 * vertical.capacity_utilization.expect("vertical"),
        1.0e-12,
        "horizontal utilisation",
    );
}

#[test]
fn a_vertical_vessel_ignores_the_liquid_level_fraction() {
    let low = run(
        &names(),
        1.0,
        &[0.7, 0.3],
        2_000_000.0,
        300.0,
        Some(1.0),
        Some(0.11),
        Some("vertical"),
        Some(0.2),
    );
    let high = run(
        &names(),
        1.0,
        &[0.7, 0.3],
        2_000_000.0,
        300.0,
        Some(1.0),
        Some(0.11),
        Some("vertical"),
        Some(0.8),
    );
    assert_eq!(low.capacity_utilization, high.capacity_utilization);
    assert_eq!(low.gas_load_factor, high.gas_load_factor);
}

#[test]
fn a_dry_gas_uses_the_sizing_fallback() {
    let result = run(
        &["methane".to_string()],
        1.0,
        &[1.0],
        5_000_000.0,
        300.0,
        Some(1.0),
        Some(0.11),
        Some("vertical"),
        Some(0.8),
    );
    close(
        result
            .max_allowable_gas_velocity
            .expect("max velocity")
            .value,
        0.575_421_321_890_263_7,
        2.0e-4,
        "max_allowable_gas_velocity",
    );
    close(
        result.gas_load_factor.expect("load factor").value,
        1.107_583_857_069_589_1e-4,
        2.0e-4,
        "gas_load_factor",
    );
}

#[test]
fn the_overload_flag_is_crossed() {
    let result = run(
        &names(),
        50.0,
        &[0.7, 0.3],
        2_000_000.0,
        300.0,
        Some(0.2),
        Some(0.11),
        Some("vertical"),
        Some(0.8),
    );
    close(
        result.capacity_utilization.expect("utilisation"),
        2.657_659_051_880_372_7,
        2.0e-4,
        "capacity_utilization",
    );
    assert_eq!(result.is_overloaded, Some(true));
}

#[test]
fn a_liquid_only_feed_uses_no_capacity() {
    let result = run(
        &["n-butane".to_string()],
        1.0,
        &[1.0],
        500_000.0,
        250.0,
        Some(1.0),
        Some(0.11),
        Some("vertical"),
        None,
    );
    assert_eq!(result.capacity_utilization, Some(0.0));
    assert_eq!(result.is_overloaded, Some(false));
    assert!(result.max_allowable_gas_velocity.is_none());
    assert!(result.gas_superficial_velocity.is_none());
    assert!(result.gas_load_factor.is_none());
}

#[test]
fn half_a_vessel_is_refused_both_ways() {
    let refused = |diameter: Option<f64>, load_factor: Option<f64>| {
        separator(
            &names(),
            1.0,
            &[0.7, 0.3],
            pascals(2_000_000.0),
            kelvins(300.0),
            pascals(0.0),
            0.0,
            None,
            diameter.map(meters),
            load_factor.map(meters_per_second),
            None,
            None,
        )
    };
    assert!(refused(Some(1.0), None).is_err());
    assert!(refused(None, Some(0.11)).is_err());
    assert!(refused(None, None).is_ok());
}

#[test]
fn the_gas_area_inputs_need_the_vessel() {
    let refused = |orientation: Option<&str>, level: Option<f64>| {
        separator(
            &names(),
            1.0,
            &[0.7, 0.3],
            pascals(2_000_000.0),
            kelvins(300.0),
            pascals(0.0),
            0.0,
            None,
            None,
            None,
            orientation,
            level,
        )
    };
    assert!(refused(Some("vertical"), None).is_err());
    assert!(refused(None, Some(0.5)).is_err());
}

#[test]
fn an_orientation_that_is_neither_is_refused() {
    let result = separator(
        &names(),
        1.0,
        &[0.7, 0.3],
        pascals(2_000_000.0),
        kelvins(300.0),
        pascals(0.0),
        0.0,
        None,
        Some(meters(1.0)),
        Some(meters_per_second(0.11)),
        Some("inclined"),
        None,
    );
    assert!(result.is_err());
}
