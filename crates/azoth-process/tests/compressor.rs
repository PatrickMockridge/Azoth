//! `process.compressor`'s power and its chart-type refusals.
//!
//! The chart's own rows and the chart branch are held by `compressor_chart.rs`'s tests and
//! the model case; this is the isentropic power the record now publishes and the named
//! refusals the spec declares.

use azoth_core::units::{kelvins, pascals};
use azoth_process::models::compressor::compressor;

#[test]
fn the_isentropic_power_is_the_enthalpy_step() {
    let result = compressor(
        &["methane".to_string(), "n-butane".to_string()],
        1.0,
        &[0.9, 0.1],
        pascals(3_000_000.0),
        kelvins(320.0),
        Some(pascals(6_000_000.0)),
        0.75,
        None,
        None,
        None,
    )
    .expect("the captured state compresses");

    // `process_compressor.tsv`'s `power_kW=2.4188668614585818`.
    let expected = 2418.8668614585818;
    let error = (result.power.value - expected).abs() / expected;
    assert!(
        error < 1.0e-5,
        "power: {} against {expected}, relative error {error:.3e}",
        result.power.value
    );
    assert_eq!(result.polytropic_head.value, 0.0);
    assert_eq!(result.polytropic_efficiency, 0.0);
}

#[test]
fn the_other_chart_types_are_refused_by_name() {
    for chart_type in ["interpolate", "interpolate_and_extrapolate", "khader_2015"] {
        let error = compressor(
            &["methane".to_string()],
            1.0,
            &[1.0],
            pascals(3_000_000.0),
            kelvins(320.0),
            Some(pascals(6_000_000.0)),
            0.75,
            None,
            None,
            Some(chart_type),
        )
        .expect_err("the chart type is a declared refusal");
        assert_eq!(error.field(), Some("chart_type"), "{chart_type}");
    }
}
