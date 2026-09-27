//! `ColumnPumparound`: a liquid draw that **comes back to another tray**.
//!
//! The draw itself this port already carries - `SimpleTray.setLiquidPumparoundDrawFraction` is
//! one of the three fraction vectors, and a tray hands its pumparound stream out beside its
//! liquid. What the draw alone does **not** carry is the return: `addLiquidPumparound(name,
//! drawTray, returnTray, fraction, drop)` cools the withdrawn liquid by `drop` kelvin and feeds
//! it in at another tray, so the liquid leaves a stage and arrives at a different one.
//!
//! **That is a recycle, and the class converges it the ordinary way.** Each outer pass runs the
//! whole column, rebuilds every return from the draw it has just produced -
//! `updateReturnStream`'s own clone, restate at `T - drop`, re-flash - and stops when the
//! return's *flow* stops moving. There are no candidates here and no acceptance rule: unlike the
//! side-draw flow tear, a pumparound return is a stream the column already knows how to take,
//! and the iteration is a fixed point on its flow.
//!
//! **The oracle is `side_draw_pumparound_return_binary`**, the binary column this port
//! reproduces exactly with `addLiquidPumparound("PA", 1, 3, 0.10, 5.0)` on it. Measured there:
//! the return is the draw at exactly `T - 5 K` (`334.5849285970497` -> `329.5849285970497`),
//! the flow settles to a `2.8e-5` relative change in 19 inner iterations against the class's
//! own `1e-4` tolerance, the cooler takes `-438.78378433523903` W, and the profile moves - tray
//! 1 lands at `334.58` K against the ideal column's `336.15`.

use azoth_core::units::kelvins;
use azoth_core::{AzothError, Result};

use crate::kernels::distillation_column::{ColumnOutcome, ColumnSetup, PumparoundDiagnostics};
use crate::stream::Stream;

/// `DEFAULT_PUMPAROUND_TOLERANCE`: the class's own field initialiser.
const DEFAULT_TOLERANCE: f64 = 1.0e-4;

/// `maxPumparoundIterations`' own initialiser.
const DEFAULT_ITERATIONS: usize = 12;

/// The floor the relative change's scale takes, as `updateReturnStream` takes it.
const SCALE_FLOOR: f64 = 1.0e-12;

/// **One liquid pumparound**: `ColumnPumparound` minus its cached streams and its name.
///
/// The draw's *fraction* is stated here and written into the tray's pumparound vector, which is
/// the same fraction vector the id already declares - so a pumparound with a return and a
/// pumparound without one are one mechanism with an extra inlet, and not two.
#[derive(Debug, Clone, PartialEq)]
pub struct PumparoundReturn {
    /// The tray the liquid leaves, `addLiquidPumparound`'s `drawTrayNumber`.
    pub draw_tray: usize,
    /// The tray it arrives at, `returnTrayNumber`.
    pub return_tray: usize,
    /// The fraction of the draw tray's liquid withdrawn.
    pub fraction: f64,
    /// The temperature drop from the draw to the return, K.
    pub temperature_drop: f64,
}

/// **The pumparound returns, iterated to a fixed point on their own flow.**
///
/// `solve` is the inner solve - `distillation_column`'s own body - passed in because the returns
/// are an outer loop over it. Each pass carries the returns built from the previous pass's
/// draws, which is what makes it a recycle rather than a one-shot join.
///
/// # Errors
/// [`AzothError::invalid_input`] for a draw tray or a return tray the column does not have, for
/// a fraction outside `[0, 1]`, for a non-finite drop, and for two pumparounds drawing from one
/// tray - each of which `addLiquidPumparound` refuses; and whatever the last pass's solve raised.
pub fn solve_with_returns(
    setup: &ColumnSetup,
    solve: &dyn Fn(&ColumnSetup) -> Result<ColumnOutcome>,
) -> Result<ColumnOutcome> {
    let tray_count =
        setup.number_of_stages + usize::from(setup.has_reboiler) + usize::from(setup.has_condenser);
    for pumparound in &setup.pumparound_returns {
        validate(setup, tray_count, pumparound)?;
    }

    let tolerance = setup.pumparound_tolerance.unwrap_or(DEFAULT_TOLERANCE);
    let iterations_cap = setup
        .pumparound_max_iterations
        .unwrap_or(DEFAULT_ITERATIONS);

    // **The draw fractions are stated on the trays**, which is the half the id already carries:
    // `addLiquidPumparound` writes `drawFraction` onto the draw tray's own setter.
    let mut active = setup.clone();
    active.pumparound_returns = Vec::new();
    {
        let mut fractions = active
            .pumparound_fractions
            .clone()
            .unwrap_or_else(|| vec![0.0; tray_count]);
        if fractions.len() < tray_count {
            fractions.resize(tray_count, 0.0);
        }
        for pumparound in &setup.pumparound_returns {
            fractions[pumparound.draw_tray] = pumparound.fraction;
        }
        active.pumparound_fractions = Some(fractions);
    }

    let mut returns: Vec<Option<Stream>> = vec![None; tray_count];
    let mut iterations = 0_usize;
    let mut relative_change = f64::INFINITY;
    let mut outcome = solve(&active)?;

    for iteration in 0..iterations_cap {
        iterations = iteration + 1;
        let mut changed = 0.0_f64;
        for pumparound in &setup.pumparound_returns {
            let draw = outcome.pumparounds[pumparound.draw_tray].clone().ok_or_else(|| {
                AzothError::invalid_input(
                    "pumparound_returns",
                    format!(
                        "the pumparound drawing from tray {} withdrew nothing, so its return has \
                         no stream to carry to tray {}",
                        pumparound.draw_tray, pumparound.return_tray
                    ),
                )
            })?;
            let returned = cooled(&draw, pumparound.temperature_drop)?;
            changed = changed.max(relative_flow_change(
                returns[pumparound.return_tray].as_ref(),
                &returned,
            ));
            returns[pumparound.return_tray] = Some(returned);
        }
        relative_change = changed;
        active.pumparound_returns = Vec::new();
        active.pumparound_inlets = returns.clone();
        let next = solve(&active)?;
        outcome = next;
        if relative_change <= tolerance {
            break;
        }
    }

    let (flows, duties) = returns_of(&outcome, setup)?;
    outcome.pumparound = Some(PumparoundDiagnostics {
        iterations,
        relative_change,
        converged: relative_change <= tolerance,
        return_n: flows,
        duty: duties,
    });
    Ok(outcome)
}

fn validate(setup: &ColumnSetup, tray_count: usize, pumparound: &PumparoundReturn) -> Result<()> {
    for (name, tray) in [
        ("pumparound_returns_draw_tray", pumparound.draw_tray),
        ("pumparound_returns_return_tray", pumparound.return_tray),
    ] {
        if tray >= tray_count {
            return Err(AzothError::invalid_input(
                name,
                format!(
                    "tray {tray} of a column with {tray_count} tray(s), which is not a tray it has"
                ),
            ));
        }
    }
    if !pumparound.fraction.is_finite() || !(0.0..=1.0).contains(&pumparound.fraction) {
        return Err(AzothError::invalid_input(
            "pumparound_returns_fraction",
            format!(
                "a draw fraction of {} is outside [0, 1], which `addLiquidPumparound` refuses",
                pumparound.fraction
            ),
        ));
    }
    if !pumparound.temperature_drop.is_finite() {
        return Err(AzothError::invalid_input(
            "pumparound_returns_temperature_drop",
            format!(
                "a temperature drop of {} is not finite, which `addLiquidPumparound` refuses",
                pumparound.temperature_drop
            ),
        ));
    }
    let duplicate = setup
        .pumparound_returns
        .iter()
        .filter(|other| other.draw_tray == pumparound.draw_tray)
        .count();
    if duplicate > 1 {
        return Err(AzothError::invalid_input(
            "pumparound_returns_draw_tray",
            format!(
                "tray {} owns more than one pumparound: `addLiquidPumparound` refuses a second \
                 draw from the same tray, because one tray has one liquid pumparound stream",
                pumparound.draw_tray
            ),
        ));
    }
    Ok(())
}

/// **`updateReturnStream`**: the draw, cooled by the stated drop and re-flashed.
fn cooled(draw: &Stream, temperature_drop: f64) -> Result<Stream> {
    let temperature = draw.t.value - temperature_drop;
    if !temperature.is_finite() || temperature <= 0.0 {
        return Err(AzothError::invalid_input(
            "pumparound_returns_temperature_drop",
            format!(
                "cooling the draw by {temperature_drop} K puts the return at {temperature} K, \
                 which is not a temperature: `updateReturnStream` raises on the same state"
            ),
        ));
    }
    Stream::from_pt(
        draw.components.clone(),
        draw.z.clone(),
        draw.n,
        draw.p,
        kelvins(temperature),
    )
}

/// `updateReturnStream`'s own relative change: the two flows against the larger of them.
fn relative_flow_change(previous: Option<&Stream>, current: &Stream) -> f64 {
    let previous = previous.map_or(0.0, |stream| stream.n);
    (previous - current.n).abs() / previous.max(current.n).max(SCALE_FLOOR)
}

/// The published returns: each one's flow and its cooler's duty.
///
/// **The duty is the return's enthalpy less the draw's**, which is the class's own subtraction
/// and is negative for a cooler.
fn returns_of(outcome: &ColumnOutcome, setup: &ColumnSetup) -> Result<(Vec<f64>, Vec<f64>)> {
    let mut flows = Vec::new();
    let mut duties = Vec::new();
    for pumparound in &setup.pumparound_returns {
        let draw = outcome.pumparounds[pumparound.draw_tray]
            .clone()
            .ok_or_else(|| {
                AzothError::invalid_input(
                    "pumparound_returns",
                    format!(
                        "the pumparound drawing from tray {} published no draw",
                        pumparound.draw_tray
                    ),
                )
            })?;
        let returned = cooled(&draw, pumparound.temperature_drop)?;
        flows.push(returned.n);
        duties.push(returned.n * (returned.h.value - draw.h.value));
    }
    Ok((flows, duties))
}
