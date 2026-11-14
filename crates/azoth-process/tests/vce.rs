//! `process.vce` against `validation/neqsim/captures/process_vce.tsv`.

use azoth_core::units::{kilograms, meters, pascals};
use azoth_process::models::vce::{VceResult, vce};

fn close(actual: f64, expected: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= 1.0e-9,
        "{what}: {actual} against {expected}, relative error {error:.3e}"
    );
}

fn explosion(
    mass: f64,
    heat_of_combustion: f64,
    strength_class: f64,
    ambient: Option<f64>,
    distance: Option<f64>,
    target: Option<f64>,
) -> VceResult {
    vce(
        kilograms(mass),
        heat_of_combustion,
        strength_class,
        ambient.map(pascals),
        distance.map(meters),
        target.map(pascals),
    )
    .expect("the captured state is valid")
}

#[test]
fn the_test_state_reproduces_the_capture() {
    let near = explosion(1000.0, 46.0e6, 4.0, None, Some(20.0), None);
    close(
        near.combustion_energy.value,
        46_000_000_000.0,
        "the combustion energy",
    );
    close(
        near.scaled_distance.expect("a distance was stated"),
        0.260_225_342_756_867_5,
        "the scaled distance at 20 m",
    );
    close(
        near.dimensionless_overpressure
            .expect("a distance was stated"),
        0.1,
        "the class-4 plateau",
    );
    close(
        near.overpressure.expect("a distance was stated").value,
        10_132.5,
        "the overpressure at 20 m",
    );

    let far = explosion(1000.0, 46.0e6, 4.0, None, Some(200.0), None);
    close(
        far.scaled_distance.expect("a distance was stated"),
        2.602_253_427_568_675,
        "the scaled distance at 200 m",
    );
    close(
        far.dimensionless_overpressure
            .expect("a distance was stated"),
        0.038_428_232_600_477_94,
        "the power-law value",
    );
    close(
        far.overpressure.expect("a distance was stated").value,
        3_893.740_668_243_427,
        "the overpressure at 200 m",
    );
}

#[test]
fn the_inverse_is_the_class_bisection() {
    let result = explosion(1000.0, 46.0e6, 4.0, None, None, Some(5000.0));
    close(
        result
            .distance_to_overpressure
            .expect("a reachable target")
            .value,
        155.749_626_729_737_03,
        "the distance to 5 kPa",
    );
}

#[test]
fn the_ten_plateaus_are_the_capture() {
    let plateaus = [
        (1.0, 0.01, 1013.25),
        (2.0, 0.02, 2026.5),
        (3.0, 0.05, 5066.25),
        (4.0, 0.10, 10_132.5),
        (5.0, 0.20, 20_265.0),
        (6.0, 0.50, 50_662.5),
        (7.0, 1.00, 101_325.0),
        (8.0, 2.00, 202_650.0),
        (9.0, 5.00, 506_625.0),
        (10.0, 10.0, 1_013_250.0),
    ];
    for (class, dimensionless, overpressure) in plateaus {
        let result = explosion(1000.0, 46.0e6, class, None, Some(20.0), None);
        close(
            result
                .dimensionless_overpressure
                .expect("a distance was stated"),
            dimensionless,
            "the plateau",
        );
        close(
            result.overpressure.expect("a distance was stated").value,
            overpressure,
            "the overpressure",
        );
    }
}

#[test]
fn an_unreachable_target_is_absent() {
    // The class answers its near bound, 0.1 m, with 10.1 kPa for class 4 and 1.0 kPa for
    // class 1.
    assert!(
        explosion(1000.0, 46.0e6, 4.0, None, None, Some(20_700.0))
            .distance_to_overpressure
            .is_none()
    );
    assert!(
        explosion(1000.0, 46.0e6, 1.0, None, None, Some(1.0e6))
            .distance_to_overpressure
            .is_none()
    );
}

#[test]
fn the_guards_are_the_class_constructor_and_the_captures() {
    assert!(
        vce(
            kilograms(1000.0),
            46.0e6,
            4.5,
            None,
            Some(meters(20.0)),
            None,
        )
        .is_err(),
        "a non-integral strength class is refused rather than truncated"
    );
    assert!(
        vce(
            kilograms(1000.0),
            46.0e6,
            11.0,
            None,
            Some(meters(20.0)),
            None,
        )
        .is_err(),
        "the constructor refuses a class above ten"
    );
    assert!(
        vce(
            kilograms(1000.0),
            46.0e6,
            4.0,
            None,
            Some(meters(0.0)),
            None,
        )
        .is_err(),
        "the class answers an infinite overpressure at zero distance"
    );
    assert!(
        vce(
            kilograms(1000.0),
            46.0e6,
            4.0,
            None,
            None,
            Some(pascals(0.0)),
        )
        .is_err(),
        "the class throws on a non-positive target"
    );
    assert!(
        vce(kilograms(1000.0), f64::NAN, 4.0, None, None, None,).is_err(),
        "a non-finite input is refused"
    );
}
