//! `process.pfp_demand` against `validation/neqsim/captures/process_pfp_demand.tsv`.

use azoth_core::units::{kelvins, meters, seconds, watts_per_meter_kelvin, watts_per_square_meter};
use azoth_process::models::pfp_demand::{PfpDemandResult, pfp_demand};

fn close(actual: f64, expected: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= 1.0e-9,
        "{what}: {actual} against {expected}, relative error {error:.3e}"
    );
}

#[allow(clippy::too_many_arguments)]
fn demand(
    flux: f64,
    thickness: f64,
    survival: f64,
    alpha: Option<f64>,
    initial: Option<f64>,
    critical: Option<f64>,
    fire: Option<f64>,
    conductivity: Option<f64>,
    fire_type: Option<&str>,
) -> PfpDemandResult {
    pfp_demand(
        watts_per_square_meter(flux),
        meters(thickness),
        seconds(survival),
        alpha,
        initial.map(kelvins),
        critical.map(kelvins),
        fire.map(kelvins),
        conductivity.map(watts_per_meter_kelvin),
        fire_type,
    )
    .expect("the captured state is valid")
}

fn pool(flux: f64, thickness: f64, survival: f64) -> PfpDemandResult {
    demand(
        flux, thickness, survival, None, None, None, None, None, None,
    )
}

#[test]
fn the_pool_fire_state_reproduces_the_capture() {
    let result = pool(100.0e3, 0.012, 3600.0);
    close(result.steel_mass_per_area, 94.2, "the steel mass per area");
    close(
        result.bare_steel_time_to_critical.value,
        197.453_666_666_666_66,
        "the bare-steel time",
    );
    close(
        result.required_survival_time.value,
        3600.0,
        "the echoed survival time",
    );
    assert!(result.pfp_required);
    close(
        result.required_pfp_thickness.value,
        0.028_057_215_110_380_333,
        "the thickness",
    );
    close(
        result.required_pfp_thickness_mm.value,
        0.028_057_215_110_380_334,
        "the thickness as a length",
    );
    assert_eq!(result.rating, "h60");
}

#[test]
fn the_jet_fire_state_is_j120() {
    let result = demand(
        250.0e3,
        0.012,
        7200.0,
        None,
        None,
        None,
        None,
        None,
        Some("jet"),
    );
    close(
        result.bare_steel_time_to_critical.value,
        78.981_466_666_666_66,
        "the bare-steel time",
    );
    assert!(result.pfp_required);
    close(
        result.required_pfp_thickness.value,
        0.056_114_430_220_760_67,
        "the thickness",
    );
    assert_eq!(result.rating, "j120");
}

#[test]
fn the_thick_wall_needs_no_protection() {
    let result = pool(100.0e3, 0.20, 600.0);
    close(
        result.steel_mass_per_area,
        1570.0,
        "the steel mass per area",
    );
    close(
        result.bare_steel_time_to_critical.value,
        3290.894_444_444_444_6,
        "the bare-steel time",
    );
    assert!(!result.pfp_required);
    assert_eq!(result.required_pfp_thickness.value, 0.0);
    assert_eq!(result.required_pfp_thickness_mm.value, 0.0);
    assert_eq!(result.rating, "none");
}

#[test]
fn an_equal_survival_time_is_not_required() {
    let result = pool(100.0e3, 0.012, 197.453_666_666_666_66);
    assert_eq!(
        result.bare_steel_time_to_critical.value, 197.453_666_666_666_66,
        "the bare time is the capture's double, bit for bit"
    );
    assert!(!result.pfp_required, "the class's test is strict");
    assert_eq!(result.rating, "none");
    assert_eq!(result.required_pfp_thickness.value, 0.0);
}

#[test]
fn the_rating_is_h60_at_exactly_the_threshold() {
    let result = pool(100.0e3, 0.012, 3600.0 + 1.0e-6);
    assert!(result.pfp_required);
    close(
        result.required_pfp_thickness.value,
        0.028_057_215_118_174,
        "the thickness",
    );
    assert_eq!(result.rating, "h60");
}

#[test]
fn the_rating_is_h120_just_past_the_threshold() {
    let result = pool(100.0e3, 0.012, (3600.0 + 1.0e-6_f64).next_up());
    assert!(result.pfp_required);
    close(
        result.required_pfp_thickness.value,
        0.028_057_215_118_174_005,
        "the thickness",
    );
    assert_eq!(result.rating, "h120");
}

#[test]
fn a_short_jet_fire_is_j60() {
    let result = demand(
        250.0e3,
        0.012,
        3600.0,
        None,
        None,
        None,
        None,
        None,
        Some("jet"),
    );
    assert_eq!(result.rating, "j60");
}

#[test]
fn the_driving_temperature_is_floored_at_one() {
    let result = demand(
        100.0e3,
        0.012,
        3600.0,
        None,
        None,
        None,
        Some(481.0),
        None,
        None,
    );
    assert!(result.pfp_required);
    close(
        result.required_pfp_thickness.value,
        4.051_583_409_441_203e-05,
        "the floored thickness",
    );
    assert_eq!(result.rating, "h60");
}

#[test]
fn the_overrides_move_the_bare_time_and_the_thickness() {
    let result = demand(
        100.0e3,
        0.012,
        3600.0,
        Some(0.5),
        Some(293.15),
        Some(573.15),
        Some(973.15),
        Some(0.4),
        None,
    );
    close(
        result.bare_steel_time_to_critical.value,
        258.4848,
        "the bare-steel time",
    );
    close(
        result.required_pfp_thickness.value,
        0.060_166_013_630_201_85,
        "the thickness",
    );
    assert_eq!(result.rating, "h60");
}

#[test]
fn the_port_refuses_a_critical_temperature_at_or_below_the_initial() {
    for critical in [288.15, 273.15] {
        let result = pfp_demand(
            watts_per_square_meter(100.0e3),
            meters(0.012),
            seconds(3600.0),
            None,
            Some(kelvins(288.15)),
            Some(kelvins(critical)),
            None,
            None,
            None,
        );
        assert!(
            result.is_err(),
            "a critical temperature of {critical} K is not a wall state"
        );
    }
}

#[test]
fn the_range_guards_are_the_class_constructor_and_the_captures() {
    assert!(
        pfp_demand(
            watts_per_square_meter(0.0),
            meters(0.012),
            seconds(3600.0),
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .is_err()
    );
    assert!(
        pfp_demand(
            watts_per_square_meter(100.0e3),
            meters(0.0),
            seconds(3600.0),
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .is_err()
    );
    assert!(
        pfp_demand(
            watts_per_square_meter(100.0e3),
            meters(0.012),
            seconds(0.0),
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .is_err()
    );
    assert!(
        pfp_demand(
            watts_per_square_meter(100.0e3),
            meters(0.012),
            seconds(3600.0),
            Some(0.0),
            None,
            None,
            None,
            None,
            None,
        )
        .is_err()
    );
    assert!(
        pfp_demand(
            watts_per_square_meter(100.0e3),
            meters(0.012),
            seconds(3600.0),
            Some(1.5),
            None,
            None,
            None,
            None,
            None,
        )
        .is_err()
    );
    assert!(
        pfp_demand(
            watts_per_square_meter(100.0e3),
            meters(0.012),
            seconds(3600.0),
            None,
            None,
            None,
            None,
            Some(watts_per_meter_kelvin(0.0)),
            None,
        )
        .is_err()
    );
    assert!(
        pfp_demand(
            watts_per_square_meter(f64::NAN),
            meters(0.012),
            seconds(3600.0),
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .is_err()
    );
    assert!(
        pfp_demand(
            watts_per_square_meter(100.0e3),
            meters(0.012),
            seconds(3600.0),
            None,
            None,
            None,
            Some(kelvins(f64::INFINITY)),
            None,
            None,
        )
        .is_err()
    );
}

#[test]
fn an_unknown_fire_type_is_refused() {
    let result = pfp_demand(
        watts_per_square_meter(100.0e3),
        meters(0.012),
        seconds(3600.0),
        None,
        None,
        None,
        None,
        None,
        Some("pool_fire"),
    );
    assert!(result.is_err());
}
