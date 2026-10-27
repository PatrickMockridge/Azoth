//! `process.pipe`'s fittings, flow solve and Beggs-and-Brills option against
//! `validation/neqsim/captures/process_pipe_modes.tsv`.
//!
//! The adiabatic equations themselves are held by the model's cases; this is the new
//! wiring: the effective length, the bisection, and the march.

use azoth_core::units::{kelvins, meters, pascals};
use azoth_process::models::pipe::{PipeFitting, pipe};

fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= tolerance,
        "{what}: {actual} against {expected}, relative error {error:.3e} > {tolerance:.3e}"
    );
}

fn run(
    fittings: Option<&[PipeFitting]>,
    outlet_pressure: Option<f64>,
    correlation: Option<&str>,
) -> azoth_process::models::PipeResult {
    pipe(
        &["methane".to_string(), "ethane".to_string()],
        4894.970603565765,
        &[0.95, 0.05],
        pascals(10_000_000.0),
        kelvins(298.15),
        meters(50_000.0),
        meters(0.5),
        meters(1.0e-5),
        fittings,
        outlet_pressure.map(pascals),
        correlation,
        Some(5),
        None,
    )
    .expect("the captured line runs")
}

#[test]
fn the_fittings_add_their_equivalent_length() {
    let fittings = [PipeFitting {
        ld_ratio: 30.0,
        count: 2.0,
    }];
    let result = run(Some(&fittings), None, None);
    close(
        result.effective_length.value,
        50_030.0,
        1.0e-12,
        "effective_length",
    );
    // The plain line's drop is 10.539 bar; two elbows add 0.0068.
    close(
        result.pressure_drop.value,
        1_054_592.850707958,
        5.0e-4,
        "pressure_drop",
    );
    close(
        result.outlet_p.value,
        8_945_407.149292042,
        5.0e-4,
        "outlet_p",
    );
}

#[test]
fn the_stated_outlet_pressure_solves_the_flow() {
    let result = run(None, Some(8_000_000.0), None);
    close(result.outlet_n, 6656.090851893372, 2.0e-3, "outlet_n");
    close(result.outlet_p.value, 8_000_000.0, 1.0e-12, "outlet_p");
    close(
        result.pressure_drop.value,
        2_000_000.0,
        1.0e-12,
        "pressure_drop",
    );
}

#[test]
fn the_beggs_and_brills_march_runs_over_the_same_geometry() {
    let result = run(None, None, Some("beggs_and_brills"));
    close(
        result.outlet_p.value,
        8_942_233.52757637,
        5.0e-4,
        "outlet_p",
    );
    close(
        result.pressure_drop.value,
        1_057_766.47242363,
        5.0e-4,
        "pressure_drop",
    );
}
