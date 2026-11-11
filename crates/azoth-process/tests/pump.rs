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
    )
    .expect_err("there is no pressure to reach");
    assert_eq!(error.field(), Some("outlet_pressure"));
}
