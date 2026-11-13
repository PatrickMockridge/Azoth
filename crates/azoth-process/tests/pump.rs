//! `process.pump`'s performance curve against `validation/neqsim/captures/process_pump_curve.tsv`.
//!
//! The reduced fit itself is held row by row in `pump_chart.rs`'s own tests, which read the
//! capture's fluid-free rows. **What is here is the branch `Pump.run` takes**: the chart's head
//! setting the pressure, the efficiency setting the shaft power, and the NPSH pair on a running
//! pump. The eight-speed run row is pinned by the model case and the cross-impl test; this is
//! the two-speed NPSH table, which is the state that also carries a `cavitating` answer.

use azoth_core::units::{cubic_meters_per_hour, kelvins, meters, pascals};
use azoth_process::models::pump::{PumpCurvePoint, pump};

fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= tolerance,
        "{what}: {actual} against {expected}, relative error {error:.3e} > {tolerance:.3e}"
    );
}

/// `PumpNPSHCurveTest`'s two-speed table and its NPSH column, as the capture states it.
fn npsh_curve() -> Vec<PumpCurvePoint> {
    let speeds = [1000.0, 1500.0];
    let flows = [
        [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0],
        [15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0],
    ];
    let heads = [
        [120.0, 118.0, 115.0, 110.0, 103.0, 94.0, 83.0, 70.0],
        [270.0, 265.5, 258.8, 247.5, 231.8, 211.5, 186.8, 157.5],
    ];
    let efficiencies = [
        [0.60, 0.70, 0.78, 0.82, 0.81, 0.76, 0.68, 0.55],
        [0.62, 0.71, 0.79, 0.83, 0.82, 0.77, 0.69, 0.56],
    ];
    let npsh = [
        [2.0, 2.2, 2.5, 3.0, 3.8, 4.8, 6.2, 8.0],
        [4.5, 4.95, 5.625, 6.75, 8.55, 10.8, 13.95, 18.0],
    ];
    let mut points = Vec::new();
    for row in 0..2 {
        for column in 0..8 {
            points.push(PumpCurvePoint {
                speed: speeds[row],
                flow: cubic_meters_per_hour(flows[row][column]),
                head: meters(heads[row][column]),
                efficiency: efficiencies[row][column],
                npsh_required: Some(meters(npsh[row][column])),
            });
        }
    }
    points
}

#[test]
fn the_npsh_curve_run_matches_the_capture() {
    let curve = npsh_curve();
    let result = pump(
        &["water".to_string()],
        615.228050698492,
        &[1.0],
        pascals(200_000.0),
        kelvins(298.15),
        None,
        0.75,
        Some(&curve),
        Some(1000.0),
        None,
        None,
        None,
        Some(true),
        Some(1.3),
        None,
        None,
        None,
        None,
        None,
    )
    .expect("the captured state pumps");

    close(result.outlet_n, 615.228050698492, 1.0e-12, "outlet_n");
    close(
        result.outlet_p.value,
        1_259_620.8540335423,
        1.0e-4,
        "outlet_p",
    );
    close(
        result.outlet_t.value,
        298.21808354632856,
        1.0e-5,
        "outlet_t",
    );
    close(
        result.outlet_h.value,
        -44_859.52897068571,
        1.0e-5,
        "outlet_h",
    );
    close(result.head.value, 109.65976180825011, 1.0e-6, "head");
    close(result.efficiency, 0.820389061190195, 1.0e-6, "efficiency");
    close(
        result.npsh_required.value,
        3.012675168303445,
        1.0e-6,
        "npsh_required",
    );
    close(
        result.npsh_available.value,
        20.42087595770533,
        1.0e-4,
        "npsh_available",
    );
    close(result.power.value, 14_528.414709493518, 1.0e-4, "power");
    assert!(!result.cavitating, "the capture's margin is not exceeded");
}

#[test]
fn a_pump_with_neither_a_curve_nor_a_pressure_is_refused() {
    let error = pump(
        &["water".to_string()],
        615.228050698492,
        &[1.0],
        pascals(200_000.0),
        kelvins(298.15),
        None,
        0.75,
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
    .expect_err("there is no pressure to reach");
    assert_eq!(error.field(), Some("outlet_pressure"));
}

#[test]
fn the_captured_plain_design_is_reproduced() {
    let result = pump(
        &["n-butane".to_string()],
        1.0,
        &[1.0],
        pascals(500_000.0),
        kelvins(250.0),
        Some(pascals(2_000_000.0)),
        0.75,
        None,
        Some(1000.0),
        None,
        None,
        None,
        None,
        None,
        Some(true),
        None,
        None,
        None,
        None,
    )
    .expect("the design computes");

    assert_eq!(result.design_pump_type.as_deref(), Some("overhung"));
    assert_eq!(result.design_api610_pump_type.as_deref(), Some("oh2"));
    assert_eq!(
        result.design_seal_type.as_deref(),
        Some("single_mechanical")
    );
    close(
        result
            .design_rated_flow
            .expect("a flow")
            .get::<azoth_core::units::cubic_meter_per_hour>(),
        0.339_772_084_194_310_9,
        1e-6,
        "rated flow",
    );
    close(
        result.design_rated_head.expect("a head").value,
        248.290_111_165_232_93,
        1e-6,
        "rated head",
    );
    close(
        result.design_hydraulic_power.expect("a power").value,
        141.571_701_747_629_54,
        1e-6,
        "hydraulic power",
    );
    close(result.design_rated_speed.unwrap(), 1000.0, 1e-12, "speed");
    close(
        result.design_specific_speed.unwrap(),
        0.155_318_443_441_676_7,
        1e-6,
        "specific speed",
    );
    close(
        result.design_suction_specific_speed.unwrap(),
        298.332_347_349_396_6,
        1e-6,
        "suction specific speed",
    );
    assert_eq!(result.design_number_of_stages, Some(1.0));
    assert!(result.design_bep_flow.is_none());
    assert_eq!(result.design_bep_source.as_deref(), Some("not_available"));
    assert_eq!(
        result.design_npshr_source.as_deref(),
        Some("screening_estimate")
    );
    assert!(result.design_shutoff_head.is_none());
    close(
        result.design_pressure.expect("a pressure").value,
        2_150_000.0,
        1e-6,
        "design pressure",
    );
    close(
        result.design_temperature.expect("a temperature").value,
        280.762_151_467_612_08,
        1e-9,
        "design temperature",
    );
    close(
        result.design_impeller_diameter.expect("a diameter").value,
        1.0,
        1e-12,
        "impeller diameter",
    );
    close(
        result.design_impeller_width.expect("a width").value,
        0.010,
        1e-12,
        "impeller width",
    );
    close(
        result.design_shaft_diameter.expect("a shaft").value,
        0.020,
        1e-12,
        "shaft diameter",
    );
    close(
        result.design_casing_wall_thickness.expect("a wall").value,
        0.012_006_957_852_172_963,
        1e-6,
        "wall",
    );
    close(
        result
            .design_max_allowable_working_pressure
            .expect("a MAWP")
            .value,
        2_150_000.0,
        1e-6,
        "MAWP",
    );
    close(
        result.design_driver_power.expect("a driver").value,
        750.0,
        1e-9,
        "driver power",
    );
    close(result.design_driver_margin.unwrap(), 1.25, 1e-12, "margin");
    close(
        result.design_suction_nozzle_size.expect("a nozzle").value,
        0.0254,
        1e-9,
        "suction nozzle",
    );
    close(
        result
            .design_minimum_flow
            .expect("a flow")
            .get::<azoth_core::units::cubic_meter_per_hour>(),
        0.220_851_854_726_302_03,
        1e-6,
        "minimum flow",
    );
    close(
        result.design_casing_weight.expect("a weight").value,
        75.697_916_278_539_41,
        1e-6,
        "casing weight",
    );
    close(
        result.design_motor_weight.expect("a weight").value,
        22.25,
        1e-9,
        "motor weight",
    );
    close(
        result.design_baseplate_weight.expect("a weight").value,
        58.768_749_767_123_644,
        1e-6,
        "baseplate weight",
    );
    close(
        result.design_inner_diameter.expect("a diameter").value,
        1.15,
        1e-9,
        "inner diameter",
    );
    close(
        result.design_outer_diameter.expect("a diameter").value,
        1.174_013_915_704_345_9,
        1e-6,
        "outer diameter",
    );
    close(
        result.design_tantan_length.expect("a length").value,
        0.01,
        1e-9,
        "tantan",
    );
    close(
        result.design_weight_total.expect("a weight").value,
        261.977_128_283_269,
        1e-6,
        "total weight",
    );
    close(
        result.design_module_length.expect("a length").value,
        3.343_301_270_189_222,
        1e-6,
        "module length",
    );
    close(
        result.design_module_width.expect("a length").value,
        2.5,
        1e-9,
        "module width",
    );
    close(
        result.design_module_height.expect("a length").value,
        1.3,
        1e-9,
        "module height",
    );
}

#[test]
fn the_captured_multistage_design_is_reproduced() {
    let result = pump(
        &["water".to_string()],
        2775.9,
        &[1.0],
        pascals(100_000.0),
        kelvins(293.15),
        Some(pascals(6_000_000.0)),
        0.75,
        None,
        Some(3000.0),
        None,
        None,
        None,
        None,
        None,
        Some(true),
        None,
        None,
        None,
        None,
    )
    .expect("the design computes");

    assert_eq!(result.design_pump_type.as_deref(), Some("between_bearings"));
    assert_eq!(result.design_api610_pump_type.as_deref(), Some("bb3"));
    assert_eq!(result.design_number_of_stages, Some(4.0));
    close(
        result
            .design_rated_flow
            .expect("a flow")
            .get::<azoth_core::units::cubic_meter_per_hour>(),
        182.025_774_256_735_05,
        1e-6,
        "rated flow",
    );
    close(
        result.design_rated_head.expect("a head").value,
        608.100_425_325_190_7,
        1e-6,
        "rated head",
    );
    close(
        result.design_hydraulic_power.expect("a power").value,
        298_320.018_920_760_2,
        1e-6,
        "hydraulic power",
    );
    close(
        result.design_driver_power.expect("a driver").value,
        630_000.0,
        1e-9,
        "driver power",
    );
    close(
        result.design_impeller_diameter.expect("a diameter").value,
        0.347_686_063_065_166_8,
        1e-6,
        "impeller",
    );
    close(
        result.design_shaft_diameter.expect("a shaft").value,
        0.070,
        1e-9,
        "shaft",
    );
    close(
        result.design_casing_wall_thickness.expect("a wall").value,
        0.012_692_762_920_695_419,
        1e-6,
        "wall",
    );
    close(
        result.design_suction_nozzle_size.expect("a nozzle").value,
        0.2032,
        1e-9,
        "suction nozzle",
    );
    close(
        result.design_discharge_nozzle_size.expect("a nozzle").value,
        0.1524,
        1e-9,
        "discharge nozzle",
    );
    close(
        result.design_weight_total.expect("a weight").value,
        2617.792_198_958_328,
        1e-6,
        "total weight",
    );
}

#[test]
fn the_captured_purchaser_design_is_reproduced() {
    let result = pump(
        &["n-butane".to_string()],
        1.0,
        &[1.0],
        pascals(500_000.0),
        kelvins(250.0),
        Some(pascals(2_000_000.0)),
        0.75,
        None,
        Some(1000.0),
        None,
        None,
        None,
        None,
        None,
        Some(true),
        Some("oh2".to_string()),
        Some("single_mechanical".to_string()),
        Some(pascals(600_000.0)),
        Some(meters(400.0)),
    )
    .expect("the design computes");

    assert_eq!(result.design_api610_pump_type.as_deref(), Some("oh2"));
    assert_eq!(
        result.design_shutoff_source.as_deref(),
        Some("purchaser_input")
    );
    close(
        result.design_shutoff_head.expect("a head").value,
        400.0,
        1e-12,
        "shutoff",
    );
    close(
        result.design_pressure.expect("a pressure").value,
        3_016_527.976_825_907,
        1e-6,
        "design pressure",
    );
    close(
        result.design_casing_wall_thickness.expect("a wall").value,
        0.015_661_068_359_330_123,
        1e-6,
        "wall",
    );
}

#[test]
fn a_sealless_pump_is_the_magnetic_drive() {
    let result = pump(
        &["n-butane".to_string()],
        1.0,
        &[1.0],
        pascals(500_000.0),
        kelvins(250.0),
        Some(pascals(2_000_000.0)),
        0.75,
        None,
        Some(1000.0),
        None,
        None,
        None,
        None,
        None,
        Some(true),
        None,
        Some("magnetic_drive".to_string()),
        None,
        None,
    )
    .expect("the design computes");
    assert_eq!(result.design_seal_type.as_deref(), Some("magnetic_drive"));
}

#[test]
fn the_captured_vendor_npsh_design_is_reproduced() {
    let curve = npsh_curve();
    let result = pump(
        &["water".to_string()],
        615.228050698492,
        &[1.0],
        pascals(200_000.0),
        kelvins(298.15),
        None,
        0.75,
        Some(&curve),
        Some(1000.0),
        None,
        None,
        None,
        Some(true),
        Some(1.3),
        Some(true),
        None,
        None,
        None,
        None,
    )
    .expect("the design computes");

    assert_eq!(result.design_bep_source.as_deref(), Some("vendor_curve"));
    assert_eq!(result.design_npshr_source.as_deref(), Some("vendor_curve"));
    assert_eq!(
        result.design_shutoff_source.as_deref(),
        Some("vendor_curve")
    );
    close(
        result
            .design_bep_flow
            .expect("a BEP flow")
            .get::<azoth_core::units::cubic_meter_per_hour>(),
        43.6,
        1e-6,
        "BEP flow",
    );
    close(
        result.design_bep_head.expect("a BEP head").value,
        107.635_665_925_925_89,
        1e-6,
        "BEP head",
    );
    close(
        result.design_shutoff_head.expect("a head").value,
        119.707_870_370_370_34,
        1e-6,
        "shutoff",
    );
    close(
        result.design_rated_head.expect("a head").value,
        109.622_314_285_104_57,
        1e-6,
        "rated head",
    );
    close(
        result.design_driver_power.expect("a driver").value,
        18_500.0,
        1e-9,
        "driver power",
    );
    close(
        result.design_specific_speed.unwrap(),
        3.293_251_078_281_470_3,
        1e-6,
        "specific speed",
    );
    close(
        result
            .design_minimum_flow
            .expect("a flow")
            .get::<azoth_core::units::cubic_meter_per_hour>(),
        28.34,
        1e-9,
        "minimum flow",
    );
    close(
        result
            .design_maximum_flow
            .expect("a flow")
            .get::<azoth_core::units::cubic_meter_per_hour>(),
        50.14,
        1e-9,
        "maximum flow",
    );
}
