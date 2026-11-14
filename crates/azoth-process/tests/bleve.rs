//! `process.bleve` against `validation/neqsim/captures/process_bleve.tsv`.

use azoth_core::units::{kilograms, meters, watts_per_square_meter};
use azoth_process::models::bleve::{BleveResult, bleve};

fn close(actual: f64, expected: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= 1.0e-9,
        "{what}: {actual} against {expected}, relative error {error:.3e}"
    );
}

fn fireball(
    mass: f64,
    heat_of_combustion: f64,
    radiative_fraction: f64,
    humidity: Option<f64>,
    distance: Option<f64>,
    target: Option<f64>,
) -> BleveResult {
    bleve(
        kilograms(mass),
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
    let result = fireball(1000.0, 46.0e6, 0.30, None, None, None);
    close(
        result.fireball_diameter.value,
        58.0,
        "the fireball diameter",
    );
    close(result.fireball_duration.value, 4.5, "the fireball duration");
    close(
        result.surface_emissive_power.value,
        290_175.480_865_128_3,
        "the surface emissive power",
    );
    assert!(result.transmissivity.is_none(), "no distance was stated");
    assert!(
        result.incident_heat_flux.is_none(),
        "no distance was stated"
    );
}

#[test]
fn the_flux_at_one_hundred_metres() {
    let result = fireball(1000.0, 46.0e6, 0.30, None, Some(100.0), None);
    close(
        result.transmissivity.expect("a distance was stated"),
        0.372_872_136_397_780_5,
        "the transmissivity",
    );
    close(
        result
            .incident_heat_flux
            .expect("a distance was stated")
            .value,
        9099.481_359_504_47,
        "the flux",
    );
}

#[test]
fn the_transmissivity_is_the_jet_fires() {
    // `captures/process_jet_fire.tsv`'s `jet_fire_humidity_05` value at the same distance.
    let result = fireball(1000.0, 46.0e6, 0.30, Some(0.5), Some(50.0), None);
    close(
        result.transmissivity.expect("a distance was stated"),
        0.444_528_982_429_716_4,
        "the shared Wayne transmissivity",
    );
    close(
        result
            .incident_heat_flux
            .expect("a distance was stated")
            .value,
        43_392.710_739_463_79,
        "the flux at 50 m",
    );
}

#[test]
fn the_contours_are_the_class_bisections() {
    let twelve = fireball(1000.0, 46.0e6, 0.30, None, None, Some(12_500.0));
    close(
        twelve.distance_to_flux.expect("a reachable target").value,
        86.991_797_885_879_6,
        "the distance to 12.5 kW/m2",
    );
    let small = fireball(1000.0, 46.0e6, 0.30, None, None, Some(10.0));
    close(
        small.distance_to_flux.expect("a reachable target").value,
        1305.503_413_614_169_4,
        "the distance to 10 W/m2",
    );
}

#[test]
fn an_unreachable_target_is_absent() {
    // The class answers its near bound, 29.1 m, with a 143 kW/m2 flux.
    assert!(
        fireball(1000.0, 46.0e6, 0.30, None, None, Some(1.0e6))
            .distance_to_flux
            .is_none()
    );
}

#[test]
fn the_guards_include_the_fraction_the_class_omits() {
    assert!(
        bleve(
            kilograms(1000.0),
            46.0e6,
            -0.30,
            None,
            Some(meters(100.0)),
            None,
        )
        .is_err(),
        "a negative radiative fraction answers a negative SEP in the class"
    );
    assert!(
        bleve(
            kilograms(1000.0),
            46.0e6,
            1.50,
            None,
            Some(meters(100.0)),
            None,
        )
        .is_err(),
        "a fraction above one is refused"
    );
    assert!(
        bleve(
            kilograms(1000.0),
            46.0e6,
            0.30,
            Some(0.0),
            Some(meters(50.0)),
            None,
        )
        .is_err(),
        "zero humidity makes the shared transmissivity NaN"
    );
    assert!(
        bleve(
            kilograms(1000.0),
            46.0e6,
            0.30,
            None,
            Some(meters(0.0)),
            None,
        )
        .is_err(),
        "the class answers an infinite flux at zero distance"
    );
    assert!(
        bleve(
            kilograms(1000.0),
            46.0e6,
            0.30,
            None,
            None,
            Some(watts_per_square_meter(0.0)),
        )
        .is_err(),
        "a non-positive target is not a contour"
    );
    assert!(
        bleve(kilograms(1000.0), f64::NAN, 0.30, None, None, None).is_err(),
        "a non-finite input is refused"
    );
}
