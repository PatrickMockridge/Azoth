//! `process.pool_fire` against `validation/neqsim/captures/process_pool_fire.tsv`.

use azoth_core::units::{kilograms_per_square_meter_second, meters, watts_per_square_meter};
use azoth_process::models::pool_fire::{PoolFireResult, pool_fire};

fn close(actual: f64, expected: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= 1.0e-9,
        "{what}: {actual} against {expected}, relative error {error:.3e}"
    );
}

fn fire(
    diameter: f64,
    burn: f64,
    air: Option<f64>,
    distance: Option<f64>,
    target: Option<f64>,
) -> PoolFireResult {
    pool_fire(
        meters(diameter),
        kilograms_per_square_meter_second(burn),
        air,
        distance.map(meters),
        target.map(watts_per_square_meter),
    )
    .expect("the captured state is valid")
}

fn check_state(distance: f64, view: f64, flux: f64) {
    let result = fire(10.0, 0.06, None, Some(distance), None);
    close(
        result.total_burn_rate.value,
        4.712_388_980_384_69,
        "the burn rate",
    );
    close(
        result.flame_height.value,
        16.678_690_349_069_907,
        "Thomas's flame height",
    );
    close(
        result.surface_emissive_power.value,
        56_143.305_429_464_26,
        "Mudan's SEP",
    );
    close(
        result.view_factor.expect("a distance was stated"),
        view,
        "the view factor",
    );
    close(
        result
            .incident_heat_flux
            .expect("a distance was stated")
            .value,
        flux,
        "the flux",
    );
}

#[test]
fn the_test_state_reproduces_the_capture() {
    check_state(5.1, 0.208_570_099_341_240_66, 11_709.814_790_768_976);
    check_state(10.0, 0.120_848_316_125_664_49, 6_784.823_922_879_632);
    check_state(20.0, 0.094_972_274_802_021_87, 5_332.057_431_540_926);
    check_state(50.0, 0.070_863_485_946_023_75, 3_978.510_335_264_159);
}

#[test]
fn the_inverse_is_the_class_bisection() {
    let result = fire(10.0, 0.06, None, None, Some(10_000.0));
    close(
        result.distance_to_flux.expect("a reachable target").value,
        5.781_778_829_395_185,
        "the distance to 10 kW/m2",
    );
}

#[test]
fn an_unreachable_target_is_absent() {
    // The class answers its near bound, 5.1 m, with a flux of 11.7 kW/m2.
    assert!(
        fire(10.0, 0.06, None, None, Some(12_500.0))
            .distance_to_flux
            .is_none()
    );
    assert!(
        fire(10.0, 0.06, None, None, Some(1.0e6))
            .distance_to_flux
            .is_none()
    );
}

#[test]
fn zero_distance_is_the_radius_clearance() {
    let result = fire(10.0, 0.06, None, Some(0.0), None);
    close(
        result.view_factor.expect("a distance was stated"),
        0.225_098_106_768_600_7,
        "the view factor at the centre",
    );
    close(
        result
            .incident_heat_flux
            .expect("a distance was stated")
            .value,
        12_637.751_759_903_706,
        "the flux at the centre",
    );
}

#[test]
fn the_guards_are_the_class_constructor_and_the_captures() {
    assert!(
        pool_fire(
            meters(0.0),
            kilograms_per_square_meter_second(0.06),
            None,
            None,
            None,
        )
        .is_err(),
        "the constructor refuses a zero diameter"
    );
    assert!(
        pool_fire(
            meters(10.0),
            kilograms_per_square_meter_second(0.0),
            None,
            None,
            None
        )
        .is_err(),
        "the constructor refuses a zero burn rate"
    );
    assert!(
        pool_fire(
            meters(10.0),
            kilograms_per_square_meter_second(0.06),
            Some(0.0),
            None,
            None,
        )
        .is_err(),
        "the flame height divides by the air density"
    );
    assert!(
        pool_fire(
            meters(10.0),
            kilograms_per_square_meter_second(0.06),
            None,
            None,
            Some(watts_per_square_meter(0.0)),
        )
        .is_err(),
        "a non-positive target is not a contour"
    );
    assert!(
        pool_fire(
            meters(10.0),
            kilograms_per_square_meter_second(f64::NAN),
            None,
            None,
            None,
        )
        .is_err(),
        "a non-finite input is refused"
    );
}
