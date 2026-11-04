//! The splitter's capacity map and its mechanical design.

use azoth_core::units::{kelvins, meters, meters_per_second, pascals};
use azoth_process::Stream;
use azoth_process::kernels::splitter::{MAX_UTILIZATION, capacity};
use azoth_process::models::splitter as run_splitter;

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

/// A flash-derived value: the feed's density is 6.2e-8 relative from NeqSim's, the measured
/// difference the capacity capture records, so the oracle values are held at 1e-6.
fn near(a: f64, b: f64) {
    let scale = 1.0 + a.abs() + b.abs();
    assert!((a - b).abs() < 1e-6 * scale, "{a} vs {b}");
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

fn model(
    mechanical_design: Option<bool>,
    capacity_analysis_enabled: Option<bool>,
    design_velocity: Option<azoth_core::units::Velocity>,
    header_diameter: Option<azoth_core::units::Length>,
) -> azoth_core::Result<azoth_process::models::SplitterResult> {
    run_splitter(
        &["n-butane".into(), "n-pentane".into()],
        1.0,
        &[0.6, 0.4],
        pascals(1.0e6),
        kelvins(300.0),
        &[0.3, 0.7],
        None,
        None,
        header_diameter,
        capacity_analysis_enabled,
        mechanical_design,
        design_velocity,
        None,
        None,
        None,
        None,
    )
}

#[test]
fn the_design_outputs_are_absent_without_their_flag() {
    let result = model(None, None, None, None).expect("the model runs");
    assert!(result.design_header_diameter.is_none());
    assert!(result.design_branch_diameter.is_none());
    assert!(result.design_header_length.is_none());
    assert!(result.design_wall_thickness.is_none());
    assert!(result.manifold_volume.is_none());
    assert!(result.design_total_pressure_drop.is_none());
    assert!(result.header_weight.is_none());
    assert!(result.branch_weight.is_none());
    assert!(result.total_weight.is_none());
}

#[test]
fn the_design_is_the_kernels_answer_on_the_record() {
    let result = model(Some(true), None, None, None).expect("the design runs");
    assert_eq!(
        result.design_header_diameter.expect("present").value,
        0.0127
    );
    assert_eq!(
        result.design_branch_diameter.expect("present").value,
        0.0127
    );
    assert_eq!(result.design_header_length.expect("present").value, 1.5);
    assert_eq!(result.design_wall_thickness.expect("present").value, 0.003);
    close(
        result.manifold_volume.expect("present").value,
        1.9001530466156163e-4,
    );
    near(
        result.design_total_pressure_drop.expect("present").value,
        107.92248156311979,
    );
    close(
        result.header_weight.expect("present").value,
        1.7423351396625337,
    );
    close(
        result.branch_weight.expect("present").value,
        1.1615567597750225,
    );
    close(
        result.total_weight.expect("present").value,
        2.903891899437556,
    );
}

#[test]
fn the_designed_header_is_what_the_capacity_reads() {
    let result = model(Some(true), Some(true), None, None).expect("the design and map run");
    near(
        result.velocity_value.expect("present").value,
        0.8580239773996624,
    );
    assert_eq!(result.velocity_design.expect("present").value, 30.0);
    near(
        result.velocity_utilization.expect("present"),
        0.028600799246655415,
    );
    assert_eq!(result.bottleneck_constraint.as_deref(), Some("velocity"));
}

#[test]
fn a_design_input_without_the_flag_is_refused() {
    let error = model(None, None, Some(meters_per_second(10.0)), None)
        .expect_err("design_velocity without mechanical_design");
    assert!(error.to_string().contains("design_velocity"));
}

#[test]
fn a_stated_header_with_the_flag_is_refused() {
    let error =
        model(Some(true), None, None, Some(meters(0.02))).expect_err("the design sizes the header");
    assert!(error.to_string().contains("header_diameter"));
}
