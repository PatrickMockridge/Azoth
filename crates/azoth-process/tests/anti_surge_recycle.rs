//! `process.anti_surge_recycle`: the branch the loop must not enter, and its cap.

use azoth_core::units::{kelvins, pascals};
use azoth_process::compressor_chart::{CompressorChart, CurvePoint};
use azoth_process::kernels::anti_surge_recycle::{RecycleSettings, anti_surge_recycle};
use azoth_process::stream::Stream;
use azoth_process::surge_curve::SurgeCurve;

/// A machine whose inlet is a few m3/hr, on a compact three-point map. The measured flows are
/// in the inlet's own range, so the reduced curve is interpolated rather than extrapolated.
fn machine(surge_flow: f64) -> (Stream, CompressorChart) {
    let feed = Stream::from_pt(
        vec!["methane".to_string(), "n-butane".to_string()],
        vec![0.9, 0.1],
        1.0,
        pascals(3_000_000.0),
        kelvins(320.0),
    )
    .expect("the captured state is a stream");

    let points = [
        CurvePoint {
            speed: 3000.0,
            flow: 2.0,
            head: 5000.0,
            efficiency: 0.8,
        },
        CurvePoint {
            speed: 3000.0,
            flow: 3.0,
            head: 4500.0,
            efficiency: 0.8,
        },
        CurvePoint {
            speed: 3000.0,
            flow: 4.0,
            head: 4000.0,
            efficiency: 0.8,
        },
    ];
    let mut chart = CompressorChart::new(&points).expect("the compact map fits");
    chart.set_surge_curve(
        SurgeCurve::new(
            &[surge_flow, surge_flow * 1.2, surge_flow * 1.4],
            &[0.0, 10_000.0, 20_000.0],
        )
        .expect("the boundary fits"),
    );
    (feed, chart)
}

#[test]
fn a_natural_flow_above_the_control_line_opens_nothing() {
    // The surge flow is below the inlet's, so the control line is too: the branch the loop
    // must not enter is the one that answers without a pass.
    let (feed, chart) = machine(0.5);
    let outcome = anti_surge_recycle(&feed, &chart, 3000.0, &RecycleSettings::default())
        .expect("the natural point is solvable");

    assert!(!outcome.recycle_active);
    assert!(outcome.converged);
    assert_eq!(outcome.iterations, 0);
    assert_eq!(outcome.recycle_mass_flow, 0.0);
    assert!(
        outcome.control_flow < outcome.inlet_volume_flow,
        "control {} against inlet {}",
        outcome.control_flow,
        outcome.inlet_volume_flow
    );
    assert!(outcome.distance_to_surge > 0.0);
}

#[test]
fn an_exhausted_cap_reports_the_pass_it_ran() {
    // The surge flow is above the inlet's, so a recycle opens; a cap of one pass and a gate no
    // pass can meet stop the loop with `converged` false and the one pass it ran reported.
    let (feed, chart) = machine(3.0);
    let settings = RecycleSettings {
        tolerance: 1.0e-12,
        max_iterations: 1,
        ..RecycleSettings::default()
    };
    let outcome =
        anti_surge_recycle(&feed, &chart, 3000.0, &settings).expect("one pass is solvable");

    assert!(outcome.recycle_active);
    assert!(!outcome.converged);
    assert_eq!(outcome.iterations, 1);
    assert!(outcome.recycle_mass_flow > 0.0);
    // The pass ran and did not meet the gate, so the reported recycle mass is the update that
    // followed it - the state the class's own `break`-before-update rule reports here.
    assert!(
        (outcome.inlet_volume_flow - outcome.control_flow).abs() > 0.0,
        "a one-pass solve reported a converged inlet"
    );
}
