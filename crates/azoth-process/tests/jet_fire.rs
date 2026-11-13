//! `process.jet_fire` against `validation/neqsim/captures/process_jet_fire.tsv`.

use azoth_core::units::{kilograms_per_second, meters, watts_per_square_meter};
use azoth_process::models::jet_fire::{JetFireResult, jet_fire};

fn close(actual: f64, expected: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= 1.0e-9,
        "{what}: {actual} against {expected}, relative error {error:.3e}"
    );
}

fn fire(
    mass: f64,
    heat_of_combustion: f64,
    radiative_fraction: f64,
    humidity: Option<f64>,
    distance: Option<f64>,
    target: Option<f64>,
) -> JetFireResult {
    jet_fire(
        kilograms_per_second(mass),
        heat_of_combustion,
        radiative_fraction,
        humidity,
        distance.map(meters),
        target.map(watts_per_square_meter),
    )
    .expect("the captured state is valid")
}

#[test]
fn the_test_state_reproduces_the_capture() {
    let near = fire(2.0, 50.0e6, 0.20, None, Some(20.0), None);
    close(
        near.total_radiative_power.value,
        2.0e7,
        "the radiated power",
    );
    close(
        near.transmissivity.expect("a distance was stated"),
        0.532_667_632_681_730_9,
        "transmissivity at 20 m",
    );
    close(
        near.incident_heat_flux
            .expect("a distance was stated")
            .value,
        2119.417_169_158_887,
        "the flux at 20 m",
    );
    assert!(near.distance_to_flux.is_none(), "no target was stated");

    let far = fire(2.0, 50.0e6, 0.20, None, Some(80.0), None);
    close(
        far.transmissivity.expect("a distance was stated"),
        0.396_408_908_804_208_2,
        "transmissivity at 80 m",
    );
    close(
        far.incident_heat_flux.expect("a distance was stated").value,
        98.578_808_315_397_01,
        "the flux at 80 m",
    );
}

#[test]
fn the_inverse_is_the_class_fixed_point() {
    let twelve = fire(2.0, 50.0e6, 0.20, None, None, Some(12500.0));
    close(
        twelve.distance_to_flux.expect("a target was stated").value,
        8.780_123_032_005_108,
        "the distance to 12.5 kW/m2",
    );
    assert!(twelve.transmissivity.is_none(), "no distance was stated");

    let small = fire(2.0, 50.0e6, 0.20, None, None, Some(10.0));
    close(
        small.distance_to_flux.expect("a target was stated").value,
        214.477_057_488_732,
        "the distance to 10 W/m2, inside the 25-iteration cap",
    );
}

#[test]
fn humidity_moves_the_transmissivity() {
    let half = fire(2.0, 50.0e6, 0.20, Some(0.5), Some(50.0), None);
    close(
        half.transmissivity.expect("a distance was stated"),
        0.444_528_982_429_716_4,
        "tau at half humidity",
    );
    close(
        half.incident_heat_flux
            .expect("a distance was stated")
            .value,
        282.995_939_605_198_6,
        "the flux at half humidity",
    );

    let saturated = fire(2.0, 50.0e6, 0.20, Some(1.0), Some(50.0), None);
    close(
        saturated.transmissivity.expect("a distance was stated"),
        0.372_872_136_397_780_5,
        "tau at saturation",
    );
    close(
        saturated
            .incident_heat_flux
            .expect("a distance was stated")
            .value,
        237.377_774_595_768_76,
        "the flux at saturation",
    );
}

#[test]
fn inside_one_metre_the_transmissivity_is_evaluated_at_one() {
    let result = fire(2.0, 50.0e6, 0.20, Some(0.5), Some(0.5), None);
    close(
        result.transmissivity.expect("a distance was stated"),
        0.768_486_697_041_438_2,
        "tau at the one-metre floor",
    );
    close(
        result
            .incident_heat_flux
            .expect("a distance was stated")
            .value,
        4_892_338.261_380_348,
        "the flux at half a metre",
    );
}

#[test]
fn the_guards_are_the_class_constructor_and_the_captures() {
    assert!(
        jet_fire(
            kilograms_per_second(2.0),
            50.0e6,
            0.20,
            None,
            Some(meters(0.0)),
            None,
        )
        .is_err(),
        "the class answers an infinite flux at zero distance"
    );
    assert!(
        jet_fire(
            kilograms_per_second(2.0),
            50.0e6,
            0.20,
            None,
            None,
            Some(watts_per_square_meter(0.0)),
        )
        .is_err(),
        "the class throws on a non-positive target"
    );
    assert!(
        jet_fire(
            kilograms_per_second(2.0),
            50.0e6,
            0.20,
            Some(0.0),
            Some(meters(50.0)),
            None,
        )
        .is_err(),
        "zero humidity makes the class's transmissivity NaN"
    );
    assert!(
        jet_fire(
            kilograms_per_second(2.0),
            50.0e6,
            1.5,
            None,
            Some(meters(20.0)),
            None,
        )
        .is_err(),
        "the constructor refuses a radiative fraction above one"
    );
    assert!(
        jet_fire(kilograms_per_second(2.0), 0.0, 0.20, None, None, None).is_err(),
        "the constructor refuses a non-positive heat of combustion"
    );
    assert!(
        jet_fire(kilograms_per_second(2.0), f64::NAN, 0.20, None, None, None).is_err(),
        "a non-finite input is refused"
    );
}
