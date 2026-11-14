//! `process.fire_water_demand` against
//! `validation/neqsim/captures/process_fire_water_demand.tsv`.

use azoth_core::units::{seconds, square_meters};
use azoth_process::models::fire_water_demand::{
    FireWaterDemandResult, ProtectedObject, fire_water_demand,
};

fn close(actual: f64, expected: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= 1.0e-9,
        "{what}: {actual} against {expected}, relative error {error:.3e}"
    );
}

fn lpm(value: f64) -> f64 {
    value * 60_000.0
}

fn m3h(value: f64) -> f64 {
    value * 3600.0
}

fn demand(
    area: f64,
    rate: Option<f64>,
    duration_s: f64,
    foam: Option<f64>,
    factor: Option<f64>,
    rows: &[(f64, f64)],
) -> FireWaterDemandResult {
    let objects: Vec<ProtectedObject> = rows
        .iter()
        .map(|(surface, row_rate)| ProtectedObject {
            surface_area: square_meters(*surface),
            application_rate: *row_rate,
        })
        .collect();
    fire_water_demand(
        square_meters(area),
        rate,
        Some(seconds(duration_s)),
        foam,
        factor,
        Some(&objects),
    )
    .expect("the captured state is valid")
}

#[test]
fn the_area_demand_follows_the_application_rate() {
    let result = demand(510.0, Some(10.0), 1800.0, None, None, &[]);
    close(lpm(result.area_demand.value), 5100.0, "the area demand");
    assert_eq!(result.object_demand.value, 0.0);
    close(lpm(result.total_demand.value), 5100.0, "the total demand");
    close(
        m3h(result.total_demand_m3_per_hour.value),
        306.0,
        "the hourly rate",
    );
    close(result.water_volume.value, 153.0, "the water volume");
    assert_eq!(result.foam_concentrate_volume.value, 0.0);
    assert_eq!(result.object_to_area_demand_ratio, Some(0.0));
}

#[test]
fn the_wellhead_rate_doubles_the_process_area_rate() {
    let result = demand(510.0, Some(20.0), 1800.0, None, None, &[]);
    close(lpm(result.area_demand.value), 10200.0, "the area demand");
    close(
        m3h(result.total_demand_m3_per_hour.value),
        612.0,
        "the hourly rate",
    );
    close(result.water_volume.value, 306.0, "the water volume");
}

#[test]
fn the_horizontal_vessel_surface_is_shell_and_ends() {
    let result = demand(
        0.0,
        Some(10.0),
        1800.0,
        None,
        None,
        &[(15.756_720_785_932_68, 10.0)],
    );
    close(
        lpm(result.object_demand.value),
        157.567_207_859_326_8,
        "the object demand",
    );
    close(
        lpm(result.total_demand.value),
        157.567_207_859_326_8,
        "the total demand",
    );
    close(
        m3h(result.total_demand_m3_per_hour.value),
        9.454_032_471_559_607,
        "the hourly rate",
    );
    close(
        result.water_volume.value,
        4.727_016_235_779_804,
        "the water volume",
    );
    assert_eq!(result.object_to_area_demand_ratio, None);
}

#[test]
fn the_selective_object_protection_ratio_is_small() {
    let result = demand(
        510.0,
        Some(10.0),
        1800.0,
        None,
        None,
        &[(15.756_720_785_932_68, 10.0)],
    );
    close(
        lpm(result.total_demand.value),
        5257.567_207_859_327,
        "the total demand",
    );
    close(
        result
            .object_to_area_demand_ratio
            .expect("an area demand exists"),
        0.030_895_530_952_809_176,
        "the object-to-area ratio",
    );
}

#[test]
fn the_foam_concentrate_scales_with_the_water_volume() {
    let result = demand(100.0, Some(10.0), 1800.0, Some(3.0), None, &[]);
    close(result.water_volume.value, 30.0, "the water volume");
    close(result.foam_concentrate_volume.value, 0.9, "the foam volume");
    close(
        m3h(result.total_demand_m3_per_hour.value),
        60.0,
        "the hourly rate",
    );
}

#[test]
fn the_simultaneous_factor_multiplies_the_demand() {
    let result = demand(200.0, Some(10.0), 1800.0, None, Some(2.0), &[]);
    close(lpm(result.area_demand.value), 2000.0, "the area demand");
    close(lpm(result.total_demand.value), 4000.0, "the total demand");
    close(
        m3h(result.total_demand_m3_per_hour.value),
        240.0,
        "the hourly rate",
    );
    close(result.water_volume.value, 120.0, "the water volume");
}

#[test]
fn a_zero_area_demand_has_no_ratio() {
    let result = demand(0.0, Some(10.0), 1800.0, None, None, &[]);
    assert_eq!(result.area_demand.value, 0.0);
    assert_eq!(result.total_demand.value, 0.0);
    assert_eq!(result.object_to_area_demand_ratio, None);
}

#[test]
fn the_defaults_are_the_class_ones() {
    let result = demand(510.0, None, 1800.0, None, None, &[]);
    close(lpm(result.area_demand.value), 5100.0, "the NORSOK rate");
    close(result.water_volume.value, 153.0, "the 30-minute duration");
}

#[test]
fn the_port_refuses_a_negative_object_row() {
    let result = fire_water_demand(
        square_meters(0.0),
        None,
        None,
        None,
        None,
        Some(&[ProtectedObject {
            surface_area: square_meters(-1.0),
            application_rate: 10.0,
        }]),
    );
    assert!(result.is_err(), "a negative surface is not an object");
}

#[test]
fn the_range_guards_are_the_class_constructor_and_setters() {
    assert!(fire_water_demand(square_meters(-1.0), None, None, None, None, None).is_err());
    assert!(fire_water_demand(square_meters(510.0), Some(-1.0), None, None, None, None).is_err());
    assert!(
        fire_water_demand(
            square_meters(510.0),
            None,
            Some(seconds(0.0)),
            None,
            None,
            None
        )
        .is_err()
    );
    assert!(fire_water_demand(square_meters(510.0), None, None, Some(11.0), None, None).is_err());
    assert!(fire_water_demand(square_meters(510.0), None, None, None, Some(0.5), None).is_err());
    assert!(fire_water_demand(square_meters(f64::NAN), None, None, None, None, None).is_err());
    assert!(
        fire_water_demand(
            square_meters(510.0),
            Some(f64::INFINITY),
            None,
            None,
            None,
            None
        )
        .is_err()
    );
}
