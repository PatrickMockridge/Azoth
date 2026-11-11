//! `process.heat_exchanger`'s `runTransient`, held to
//! `validation/neqsim/captures/process_heat_exchanger_transient.tsv`.
//!
//! The wall trajectory is flow-independent while no holdup is stated, so the upstream test's
//! ten steps are reproduced at the capture's cased flow. The CSTR states and the outlets
//! inherit the flash, so they carry a tolerance; the wall does not.

use azoth_core::AzothError;
use azoth_core::units::{kelvins, pascals};
use azoth_process::Stream;
use azoth_process::kernels::heat_exchanger::{TransientSettings, heat_exchanger_transient};

/// 100 000 kg/hr of methane, the capture's cased flow.
const N: f64 = 1731.4578182246325;

fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= tolerance,
        "{what}: {actual} against {expected}, relative error {error:.3e} > {tolerance:.3e}"
    );
}

fn methane(temperature: f64, pressure: f64, flow: f64) -> Stream {
    Stream::from_pt(
        vec!["methane".to_string()],
        vec![1.0],
        flow,
        pascals(pressure),
        kelvins(temperature),
    )
    .expect("methane flashes")
}

fn settings() -> TransientSettings {
    TransientSettings {
        time_step_s: 1.0,
        wall_mass_kg: 200.0,
        wall_cp_j_per_kg_k: 500.0,
        wall_temperature_k: None,
        shell_holdup_volume_m3: 0.0,
        tube_holdup_volume_m3: 0.0,
        shell_side_htc_w_per_m2_k: 500.0,
        tube_side_htc_w_per_m2_k: 800.0,
        heat_transfer_area_m2: 20.0,
        shell_fluid_temperature_k: None,
        tube_fluid_temperature_k: None,
        tube_is_hot: false,
    }
}

#[test]
fn the_first_step_seeds_the_wall_to_the_mean_of_the_inlets() {
    let hot = methane(353.15, 1.0e6, N);
    let cold = methane(293.15, 1.0e6, N);
    let outcome = heat_exchanger_transient(&hot, &cold, &settings()).expect("the step runs");
    let state = outcome
        .transient
        .expect("the dynamic branch carries a state");
    close(state.wall_temperature.value, 321.35, 1e-12, "wall seed");
    assert!(state.shell_fluid_temperature.is_none());
    assert!(state.tube_fluid_temperature.is_none());
    close(
        outcome.hot_out.t.value,
        348.6848889155141,
        1e-5,
        "hot outlet",
    );
    close(
        outcome.cold_out.t.value,
        300.6934680836265,
        1e-5,
        "cold outlet",
    );
    close(outcome.duty.value, 300_000.0, 1e-5, "duty");
}

#[test]
fn ten_steps_of_the_wall_ode_follow_the_capture() {
    let hot = methane(353.15, 1.0e6, N);
    let cold = methane(293.15, 1.0e6, N);
    let expected = [
        321.35,
        320.018,
        319.03232,
        318.3029168,
        317.763158432,
        317.36373723968,
        317.0681655573632,
        316.8494425124488,
        316.68758745921207,
        316.5678147198169,
    ];
    let mut settings = settings();
    for (step, want) in expected.iter().enumerate() {
        let outcome = heat_exchanger_transient(&hot, &cold, &settings).expect("the step runs");
        let state = outcome
            .transient
            .expect("the dynamic branch carries a state");
        close(
            state.wall_temperature.value,
            *want,
            1e-9,
            &format!("step {}", step + 1),
        );
        settings.wall_temperature_k = Some(state.wall_temperature.value);
    }
}

#[test]
fn the_holdup_cstrs_carry_the_fluid_temperatures() {
    let hot = methane(423.15, 3.0e6, N);
    let cold = methane(303.15, 3.0e6, N);
    let mut settings = settings();
    settings.wall_mass_kg = 400.0;
    settings.heat_transfer_area_m2 = 2.0;
    settings.shell_holdup_volume_m3 = 40.0;
    settings.tube_holdup_volume_m3 = 30.0;

    let first = heat_exchanger_transient(&hot, &cold, &settings).expect("the step runs");
    let state = first.transient.expect("the dynamic branch carries a state");
    close(
        state.wall_temperature.value,
        362.9703196616055,
        1e-9,
        "wall 1",
    );
    close(
        state.shell_fluid_temperature.expect("shell holdup").value,
        423.10975376724326,
        1e-6,
        "shell CSTR 1",
    );
    close(
        state.tube_fluid_temperature.expect("tube holdup").value,
        303.2151115961644,
        1e-6,
        "tube CSTR 1",
    );

    settings.wall_temperature_k = Some(state.wall_temperature.value);
    settings.shell_fluid_temperature_k = state.shell_fluid_temperature.map(|t| t.value);
    settings.tube_fluid_temperature_k = state.tube_fluid_temperature.map(|t| t.value);
    let second = heat_exchanger_transient(&hot, &cold, &settings).expect("the step runs");
    let state = second
        .transient
        .expect("the dynamic branch carries a state");
    close(
        state.wall_temperature.value,
        362.79327843795187,
        1e-9,
        "wall 2",
    );
    close(
        state.shell_fluid_temperature.expect("shell holdup").value,
        423.0714363325721,
        1e-6,
        "shell CSTR 2",
    );
    close(
        state.tube_fluid_temperature.expect("tube holdup").value,
        303.27696878554656,
        1e-6,
        "tube CSTR 2",
    );
}

#[test]
fn the_outlet_shift_does_not_scale_with_dt() {
    let hot = methane(353.15, 1.0e6, N);
    let cold = methane(293.15, 1.0e6, N);
    let mut first: Option<(f64, f64, f64)> = None;
    for (dt, wall) in [(0.0, 323.15), (0.5, 322.25), (1.0, 321.35)] {
        let mut settings = settings();
        settings.time_step_s = dt;
        settings.wall_temperature_k = Some(323.15);
        let outcome = heat_exchanger_transient(&hot, &cold, &settings).expect("the step runs");
        let state = outcome
            .transient
            .expect("the dynamic branch carries a state");
        close(state.wall_temperature.value, wall, 1e-12, "wall");
        let here = (
            outcome.hot_out.t.value,
            outcome.cold_out.t.value,
            outcome.duty.value,
        );
        if let Some(before) = first {
            close(here.0, before.0, 1e-12, "hot outlet");
            close(here.1, before.1, 1e-12, "cold outlet");
            close(here.2, before.2, 1e-12, "duty");
        } else {
            first = Some(here);
        }
    }
}

#[test]
fn a_hot_tube_maps_the_duty_to_the_hot_side() {
    let cold = methane(293.15, 1.0e6, N);
    let hot = methane(353.15, 1.0e6, N);
    let mut settings = settings();
    settings.tube_is_hot = true;
    let outcome = heat_exchanger_transient(&hot, &cold, &settings).expect("the step runs");
    let state = outcome
        .transient
        .expect("the dynamic branch carries a state");
    close(state.wall_temperature.value, 324.95, 1e-12, "wall");
    // The class's own duty is `qShellToWall = -300 000`, on the cold shell; this id's duty is
    // the hot side's release, `-tube_duty = 480 000`.
    close(
        outcome.duty.value,
        479_999.9989547981,
        1e-6,
        "hot-side duty",
    );
}

#[test]
fn a_zero_flow_side_is_refused() {
    let hot = methane(353.15, 1.0e6, 0.0);
    let cold = methane(293.15, 1.0e6, N);
    let error = heat_exchanger_transient(&hot, &cold, &settings())
        .err()
        .expect("zero flow");
    assert!(
        matches!(&error, AzothError::InvalidInput { field, .. } if field.as_str() == "hot_in_n"),
        "{error:?}"
    );
}

#[test]
fn a_two_phase_side_is_refused() {
    let two_phase = Stream::from_pt(
        vec!["methane".to_string(), "n-butane".to_string()],
        vec![0.5, 0.5],
        N,
        pascals(500_000.0),
        kelvins(250.0),
    )
    .expect("the mixture flashes");
    let cold = methane(293.15, 500_000.0, N);
    let error = heat_exchanger_transient(&two_phase, &cold, &settings())
        .err()
        .expect("two phases");
    assert!(
        matches!(&error, AzothError::InvalidInput { field, .. } if field.as_str() == "stream"),
        "{error:?}"
    );
}

#[test]
fn a_non_positive_wall_cp_is_refused() {
    let hot = methane(353.15, 1.0e6, N);
    let cold = methane(293.15, 1.0e6, N);
    let mut settings = settings();
    settings.wall_cp_j_per_kg_k = 0.0;
    let error = heat_exchanger_transient(&hot, &cold, &settings)
        .err()
        .expect("singular wall");
    assert!(
        matches!(&error, AzothError::InvalidInput { field, .. } if field.as_str() == "wall_cp"),
        "{error:?}"
    );
}

#[test]
fn a_negative_step_is_refused() {
    let hot = methane(353.15, 1.0e6, N);
    let cold = methane(293.15, 1.0e6, N);
    let mut settings = settings();
    settings.time_step_s = -1.0;
    let error = heat_exchanger_transient(&hot, &cold, &settings)
        .err()
        .expect("backwards step");
    assert!(
        matches!(&error, AzothError::InvalidInput { field, .. } if field.as_str() == "time_step"),
        "{error:?}"
    );
}

#[test]
fn a_non_finite_carried_temperature_is_refused() {
    let hot = methane(353.15, 1.0e6, N);
    let cold = methane(293.15, 1.0e6, N);
    let mut settings = settings();
    settings.wall_temperature_k = Some(f64::INFINITY);
    let error = heat_exchanger_transient(&hot, &cold, &settings)
        .err()
        .expect("no such state");
    assert!(
        matches!(&error, AzothError::InvalidInput { field, .. } if field.as_str() == "wall_temperature"),
        "{error:?}"
    );
}

#[test]
fn a_nan_wall_temperature_is_the_first_step_seed() {
    let hot = methane(353.15, 1.0e6, N);
    let cold = methane(293.15, 1.0e6, N);
    let mut settings = settings();
    settings.wall_temperature_k = Some(f64::NAN);
    let outcome = heat_exchanger_transient(&hot, &cold, &settings).expect("the step runs");
    let state = outcome
        .transient
        .expect("the dynamic branch carries a state");
    close(state.wall_temperature.value, 321.35, 1e-12, "wall seed");
}
