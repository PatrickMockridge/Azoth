//! Spec-driven tests for `eos.valve_scale_drift`.

use azoth_core::units::{meters, meters_per_second, seconds};
use azoth_eos::model_gen;
use azoth_eos::results::ValveScaleDriftResult;
use azoth_eos::valve_scale_drift::valve_scale_drift;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.valve_scale_drift";

fn call(case: &azoth_core::spec::TestCase) -> ValveScaleDriftResult {
    valve_scale_drift(
        meters(common::input(case, "port_diameter")),
        meters_per_second(common::input(case, "growth_rate")),
        seconds(common::input(case, "elapsed_time")),
        common::input(case, "clean_opening_percent"),
        common::input(case, "valve_kv"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

fn run(growth: f64, elapsed: f64, clean: f64) -> ValveScaleDriftResult {
    valve_scale_drift(
        meters(0.05),
        meters_per_second(growth),
        seconds(elapsed),
        clean,
        100.0,
    )
    .expect("the state runs")
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        for (field, got) in [
            ("deposit_thickness", result.deposit_thickness.value),
            ("elapsed_time", result.elapsed_time.value),
            ("open_area_fraction", result.open_area_fraction),
            ("fouling_fraction", result.fouling_fraction),
            ("effective_kv", result.effective_kv),
            ("time_to_plug", result.time_to_plug.value),
            ("opening_percent", result.opening_percent),
            ("time_to_pin", result.time_to_pin.value),
        ] {
            if let Some(want) = case.expected_value(field) {
                common::assert_close(got, want, case.tolerance, &format!("{context} ({field})"));
            }
        }
        common::assert_consistent(&result, context);
    }
}

/// **The cap is the class's `0.999999`**, reached when the deposit closes the port, and the
/// effective Kv is then zero rather than negative.
#[test]
fn a_closed_port_caps_the_fouling() {
    let result = run(3.168808781402895e-11, 864_000_000.0, 50.0);
    assert_eq!(result.open_area_fraction, 0.0);
    assert_eq!(result.fouling_fraction, 0.999_999);
    assert_eq!(result.effective_kv, 0.0);
    assert_eq!(result.opening_percent, 100.0);
    assert_eq!(result.time_to_plug.value, 0.0);
    assert_eq!(result.time_to_pin.value, 0.0);
}

/// **Zero growth or no port gives an infinite time to plug and to pin**, and a deposit still
/// accumulates against an unset port while the area stays one.
#[test]
fn infinite_times_need_no_deposit() {
    let still = run(0.0, 31_557_600.0, 50.0);
    assert!(still.time_to_plug.value.is_infinite());
    assert!(still.time_to_pin.value.is_infinite());
    assert_eq!(still.deposit_thickness.value, 0.0);

    let no_port = valve_scale_drift(
        meters(0.0),
        meters_per_second(3.168808781402895e-11),
        seconds(31_557_600.0),
        50.0,
        100.0,
    )
    .expect("the state runs");
    assert!(no_port.time_to_plug.value.is_infinite());
    assert!(no_port.time_to_pin.value.is_infinite());
    assert_eq!(no_port.open_area_fraction, 1.0);
    assert_eq!(no_port.deposit_thickness.value, 0.001);
}

/// **The opening pins earlier than the port plugs**, because a valve holding half its
/// opening runs out of controllability at a third of the thickness a full closure needs.
#[test]
fn the_opening_pins_before_the_port_plugs() {
    let result = run(3.168808781402895e-11, 0.0, 50.0);
    assert!(result.time_to_pin.value < result.time_to_plug.value);
}
