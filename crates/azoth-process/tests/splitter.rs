//! The splitter's capacity map: the class's lazy initialization and its two constraints.

use azoth_core::units::{meters, meters_per_second, pascals};
use azoth_process::Stream;
use azoth_process::kernels::splitter::{MAX_UTILIZATION, capacity};

fn feed(pressure_pa: f64) -> Stream {
    Stream::from_pt(
        vec!["n-butane".into(), "n-pentane".into()],
        vec![0.6, 0.4],
        1.0,
        pascals(pressure_pa),
        azoth_core::units::kelvins(300.0),
    )
    .expect("the butane / pentane feed resolves and flashes")
}

fn close(a: f64, b: f64) {
    let scale = 1.0 + a.abs() + b.abs();
    assert!((a - b).abs() < 1e-9 * scale, "{a} vs {b}");
}

#[test]
fn a_fresh_splitter_has_no_constraints() {
    let report = capacity(&feed(1.0e6), None, None, None, false).expect("capacity");
    assert!(report.pressure_drop_value.is_none());
    assert!(report.pressure_drop_design.is_none());
    assert!(report.pressure_drop_utilization.is_none());
    assert!(report.velocity_value.is_none());
    assert!(report.velocity_design.is_none());
    assert!(report.velocity_utilization.is_none());
    assert_eq!(report.max_utilization, 0.0);
    assert!(!report.capacity_exceeded);
    assert!(report.bottleneck_constraint.is_none());
}

#[test]
fn the_header_diameter_alone_does_not_initialize_the_map() {
    let report = capacity(&feed(1.0e6), None, None, Some(meters(0.02)), false).expect("capacity");
    assert!(report.velocity_value.is_none());
    assert_eq!(report.max_utilization, 0.0);
    assert!(!report.capacity_exceeded);
}

#[test]
fn enabling_analysis_initializes_the_class_defaults() {
    let report = capacity(&feed(1.0e6), None, None, None, true).expect("capacity");
    assert_eq!(report.pressure_drop_value.expect("present").value, 0.0);
    assert_eq!(report.pressure_drop_design.expect("present").value, 5000.0);
    assert_eq!(report.pressure_drop_utilization.expect("present"), 0.0);
    assert_eq!(report.velocity_value.expect("present").value, 0.0);
    assert_eq!(report.velocity_design.expect("present").value, 30.0);
    assert_eq!(report.velocity_utilization.expect("present"), 0.0);
    assert_eq!(report.max_utilization, 0.0);
    assert!(report.bottleneck_constraint.is_none());
}

#[test]
fn a_sized_header_reports_the_velocity_and_the_bottleneck() {
    let report = capacity(
        &feed(1.0e6),
        Some(pascals(2000.0)),
        None,
        Some(meters(0.02)),
        true,
    )
    .expect("capacity");
    close(
        report.velocity_value.expect("present").value,
        0.34597673984807986,
    );
    close(
        report.velocity_utilization.expect("present"),
        0.011532557994935995,
    );
    assert_eq!(report.bottleneck_constraint, Some("velocity"));
    assert!(!report.capacity_exceeded);
}

#[test]
fn a_non_positive_limit_suppresses_only_its_own_constraint() {
    let no_drop = capacity(
        &feed(1.0e6),
        Some(pascals(0.0)),
        None,
        Some(meters(0.02)),
        false,
    )
    .expect("capacity");
    assert!(no_drop.pressure_drop_value.is_none());
    assert!(no_drop.velocity_value.is_some());
    assert_eq!(no_drop.bottleneck_constraint, Some("velocity"));

    let no_velocity = capacity(
        &feed(1.0e6),
        None,
        Some(meters_per_second(0.0)),
        Some(meters(0.02)),
        false,
    )
    .expect("capacity");
    assert!(no_velocity.velocity_value.is_none());
    assert!(no_velocity.pressure_drop_value.is_some());
    assert_eq!(no_velocity.max_utilization, 0.0);
    assert!(no_velocity.bottleneck_constraint.is_none());
}

#[test]
fn a_tight_limit_is_exceeded_even_with_analysis_disabled() {
    let report = capacity(
        &feed(1.0e6),
        None,
        Some(meters_per_second(0.1)),
        Some(meters(0.02)),
        false,
    )
    .expect("capacity");
    close(
        report.velocity_utilization.expect("present"),
        3.4597673984807984,
    );
    assert!(report.capacity_exceeded);
    assert_eq!(report.bottleneck_constraint, Some("velocity"));
}

#[test]
fn the_utilization_is_capped_at_9_99() {
    let report = capacity(
        &feed(1.0e6),
        None,
        Some(meters_per_second(0.01)),
        Some(meters(0.02)),
        false,
    )
    .expect("capacity");
    assert_eq!(
        report.velocity_utilization.expect("present"),
        MAX_UTILIZATION
    );
    assert_eq!(report.max_utilization, MAX_UTILIZATION);
    assert!(report.capacity_exceeded);
}

#[test]
fn the_velocity_reads_the_two_phase_feeds_system_flow() {
    let two_phase = feed(1.5e5);
    assert!(
        two_phase
            .vapour_fraction
            .is_some_and(|beta| beta > 0.0 && beta < 1.0),
        "the 1.5 bara feed is two phases"
    );
    let report = capacity(&two_phase, None, None, Some(meters(0.05)), true).expect("capacity");
    close(
        report.velocity_value.expect("present").value,
        4.968249374085196,
    );
    close(
        report.velocity_utilization.expect("present"),
        0.16560831246950652,
    );
}
