//! `DistillationColumn.solveSingleSideDrawFlowSpecification`: a draw's **flow** as a tear.
//!
//! Every other draw this port carries is a *fraction* of a stage's own phase, decided once.
//! A flow specification turns that fraction into a **controlled variable**: the column states
//! the mass flow a draw must deliver and moves its fraction until it does.
//!
//! **The class's search is over cold candidates, and each is a whole column solve.** NeqSim
//! deep-copies the column per candidate (`getSingleSideDrawCandidateTemplate`), sets the
//! fraction, and solves it from scratch, and only a candidate whose inner solve reports
//! `RIGOROUS_CONVERGED` or `RECONCILED_PRODUCTS` may update the controller or become the
//! published state - so a `FALLBACK_PRODUCTS` state cannot leak into the products. **The copy
//! is the one part this port does not need**: a kernel is a pure function of its setup, so a
//! cold candidate is another call. Acceptance is the remaining half, and it reduces to the
//! *same rule* - `Ok` here is an accepted state, because this kernel refuses a solve it cannot
//! reach rather than publishing a fallback for it.
//!
//! **Three things are refused rather than reduced**, each by name: more than one specification,
//! a specification beside a pumparound, and a specification on an end. The first two are
//! `solveWithColumnTearVariables`' **coordinated** problem - several tear variables updated
//! together, with its own convergence test - and this module implements the *independent*
//! single-variable search the class uses when there is exactly one.
//!
//! **The one-shot continuation retry is not ported, and that is a measured gap rather than an
//! omission.** `createSingleSideDrawContinuationCandidate` re-solves a rejected fraction *warm*,
//! from a copy of the nearest accepted state with `setDoInitializion(false)`, and this kernel
//! has no warm start: `sweep` always seeds from `seed_network`. A candidate the class rescues
//! that way is simply rejected here, so a state that needs the retry takes more candidates and
//! may exhaust the cap. On the captured oracle the retry never fires (zero rejections in two
//! candidates); on the multistage rows in the same capture the class rejects 18 of 30 and does
//! not converge either, so what the retry buys there is unmeasured.

use azoth_core::{AzothError, Result};

use crate::kernels::distillation_column::{
    ColumnOutcome, ColumnSetup, SideDrawFlow, SideDrawPhase, TearDiagnostics,
};

/// `SIDE_DRAW_CANDIDATE_SCAN_STEP`: the grid the bounded scan walks.
const CANDIDATE_SCAN_STEP: f64 = 5.0e-3;

/// `addSideDrawFlowSpecification`'s seed: a flow target on an uncontrolled draw starts the
/// fraction at five per cent rather than at zero, which would be a column that draws nothing.
const SEED_FRACTION: f64 = 0.05;

/// The two fractions are the same number within this, which is `wasSideDrawFractionAttempted`'s
/// own identity tolerance.
const FRACTION_IDENTITY: f64 = 1.0e-12;

/// The floor the residual's scale takes, so a target of zero is a defined question.
const SCALE_FLOOR: f64 = 1.0e-12;

/// **One side-draw flow target, solved by the class's own candidate search.**
///
/// `solve` is the inner solve - `distillation_column`'s own body - passed in rather than called,
/// because the tear is an outer loop over it and the two modules would otherwise be circular.
/// Each candidate is a full cold solve at one fraction.
///
/// # Errors
/// [`AzothError::invalid_input`] for a second specification, for a specification beside a
/// pumparound, and for a specification on an end; and whatever the last accepted candidate's
/// solve raised, where no candidate was ever accepted.
pub fn solve_single_flow(
    setup: &ColumnSetup,
    specification: &SideDrawFlow,
    solve: &dyn Fn(&ColumnSetup) -> Result<ColumnOutcome>,
) -> Result<ColumnOutcome> {
    let tray_count =
        setup.number_of_stages + usize::from(setup.has_reboiler) + usize::from(setup.has_condenser);
    refuse_coordinated(setup, tray_count)?;
    let tray = specification.tray;
    if tray >= tray_count {
        return Err(AzothError::invalid_input(
            "side_draw_flow_tray",
            format!(
                "a side-draw flow specification names tray {tray} of a column with {tray_count} \
                 tray(s), which is not a tray it has"
            ),
        ));
    }
    if is_end(setup, tray_count, tray) {
        return Err(AzothError::invalid_input(
            "side_draw_flow_tray",
            format!(
                "a side-draw flow specification names {tray}, which is this port's reboiler or \
                 condenser: the ends are `column::reboiler` and `column::condenser` rather than \
                 stages, so there is no draw for the specification to move"
            ),
        ));
    }

    let maximum = maximum_fraction(setup, tray, specification.phase);
    let mut candidate_fraction = fraction_of(setup, tray, specification.phase);
    if candidate_fraction <= 0.0 {
        candidate_fraction = SEED_FRACTION;
    }
    let mut attempted: Vec<f64> = Vec::new();
    let mut accepted_fractions: Vec<f64> = Vec::new();
    let mut accepted_flows: Vec<f64> = Vec::new();
    let mut history = String::new();
    let mut best_residual = f64::INFINITY;
    let mut accepted_any = false;
    let mut accepted: Option<ColumnOutcome> = None;
    let mut accepted_flow = f64::NAN;
    let mut rejected = 0_usize;
    let mut rollbacks = 0_usize;
    let mut inner_iterations = 0_u32;
    let mut iterations = 0_usize;
    let mut converged = false;

    for iteration in 0..specification.max_iterations {
        if !candidate_fraction.is_finite() {
            break;
        }
        iterations = iteration + 1;
        // numerics-ok: `solveSingleSideDrawFlowSpecification` writes this clamp itself -
        // `candidateFraction = Math.max(0.0, Math.min(maximumFraction, candidateFraction))` -
        // and it is a *search* bound, not a physical quantity: the fraction was proposed by the
        // selection above and a value outside the draw's available range is a candidate the
        // class steps back inside rather than a state to refuse.
        candidate_fraction = candidate_fraction.clamp(0.0, maximum);
        if was_attempted(&attempted, candidate_fraction) {
            candidate_fraction = select_next(
                specification,
                &accepted_fractions,
                &accepted_flows,
                &attempted,
                maximum,
            );
            if !candidate_fraction.is_finite() {
                break;
            }
        }
        attempted.push(candidate_fraction);

        let candidate_setup = with_fraction(setup, tray, specification.phase, candidate_fraction);
        let candidate = solve(&candidate_setup);
        let (outcome, flow) = match candidate {
            Ok(outcome) => {
                inner_iterations += outcome.iterations;
                let flow = draw_mass_flow(&outcome, tray, specification.phase);
                match flow {
                    Some(flow) if flow.is_finite() => (Some(outcome), flow),
                    _ => (None, f64::NAN),
                }
            }
            Err(_) => (None, f64::NAN),
        };
        append_history(
            &mut history,
            iterations,
            candidate_fraction,
            flow,
            outcome.is_some(),
        );

        if outcome.is_none() {
            rejected += 1;
            if accepted_any {
                rollbacks += 1;
            }
            candidate_fraction = select_next(
                specification,
                &accepted_fractions,
                &accepted_flows,
                &attempted,
                maximum,
            );
            continue;
        }

        accepted_any = true;
        accepted_fractions.push(candidate_fraction);
        accepted_flows.push(flow);
        let residual =
            (flow - specification.target).abs() / specification.target.abs().max(SCALE_FLOOR);
        if residual < best_residual {
            best_residual = residual;
            accepted_flow = flow;
            accepted = outcome;
        }
        if residual <= specification.tolerance {
            converged = true;
            break;
        }
        candidate_fraction = select_next(
            specification,
            &accepted_fractions,
            &accepted_flows,
            &attempted,
            maximum,
        );
    }

    let mut published = accepted.ok_or_else(|| {
        AzothError::invalid_input(
            "side_draw_flow_target",
            format!(
                "every candidate for a {} draw of {} kg/s on tray {tray} was rejected by the \
                 column solve, so the class's own `finalizeColumnTearConvergenceStatus` has \
                 nothing to publish - NeqSim reports the same state as \
                 `All single-side-draw flow candidates were rejected by the inner column solver`",
                specification.phase.name(),
                specification.target
            ),
        )
    })?;
    published.tear = Some(TearDiagnostics {
        iterations,
        residual: best_residual,
        converged,
        rejected_candidates: rejected,
        rollbacks,
        inner_iterations,
        history,
        actual_flow: accepted_flow,
        fraction: last_accepted_fraction(&accepted_fractions),
    });
    Ok(published)
}

/// A second specification, or one beside a pumparound, is the class's **coordinated** tear.
fn refuse_coordinated(setup: &ColumnSetup, _tray_count: usize) -> Result<()> {
    if setup.side_draw_flows.len() > 1 {
        return Err(AzothError::invalid_input(
            "side_draw_flow_target",
            format!(
                "{} side-draw flow specifications are stated: `solveWithColumnTearVariables` \
                 updates several tear variables together under its own convergence test, and \
                 this port carries the *independent* single-variable search \
                 `solveSingleSideDrawFlowSpecification` when exactly one is stated",
                setup.side_draw_flows.len()
            ),
        ));
    }
    if setup.pumparound_fractions.is_some() {
        return Err(AzothError::invalid_input(
            "side_draw_flow_target",
            "a side-draw flow specification is stated beside a pumparound: the class solves the \
             two together as coordinated tear variables, and this port carries the independent \
             single-variable search - which `process_column_tear.tsv` measures as the right \
             call, because the coordinated loop does not converge",
        ));
    }
    Ok(())
}

fn is_end(setup: &ColumnSetup, tray_count: usize, tray: usize) -> bool {
    (tray == 0 && setup.has_reboiler) || (tray + 1 == tray_count && setup.has_condenser)
}

/// `getMaximumSideDrawFraction`: all of a stage's vapour, or the liquid a pumparound left.
fn maximum_fraction(setup: &ColumnSetup, tray: usize, phase: SideDrawPhase) -> f64 {
    match phase {
        SideDrawPhase::Gas => 1.0,
        SideDrawPhase::Liquid => {
            let pumparound = setup
                .pumparound_fractions
                .as_ref()
                .and_then(|fractions| fractions.get(tray).copied())
                .unwrap_or(0.0);
            (1.0 - pumparound).max(0.0)
        }
    }
}

fn fraction_of(setup: &ColumnSetup, tray: usize, phase: SideDrawPhase) -> f64 {
    let vector = match phase {
        SideDrawPhase::Gas => setup.gas_side_draw_fractions.as_ref(),
        SideDrawPhase::Liquid => setup.liquid_side_draw_fractions.as_ref(),
    };
    vector
        .and_then(|fractions| fractions.get(tray).copied())
        .unwrap_or(0.0)
}

/// The same setup at one tray's stated fraction, the vector grown to the tray count where the
/// caller stated none - which is what `setSideDrawFraction` does to the tray.
fn with_fraction(
    setup: &ColumnSetup,
    tray: usize,
    phase: SideDrawPhase,
    fraction: f64,
) -> ColumnSetup {
    let tray_count =
        setup.number_of_stages + usize::from(setup.has_reboiler) + usize::from(setup.has_condenser);
    let mut next = setup.clone();
    let mut vector = match phase {
        SideDrawPhase::Gas => setup.gas_side_draw_fractions.clone(),
        SideDrawPhase::Liquid => setup.liquid_side_draw_fractions.clone(),
    }
    .unwrap_or_else(|| vec![0.0; tray_count]);
    if vector.len() < tray_count {
        vector.resize(tray_count, 0.0);
    }
    vector[tray] = fraction;
    match phase {
        SideDrawPhase::Gas => next.gas_side_draw_fractions = Some(vector),
        SideDrawPhase::Liquid => next.liquid_side_draw_fractions = Some(vector),
    }
    // **The tear is not itself a tear**: the inner solve must not recurse into this search.
    next.side_draw_flows = Vec::new();
    next
}

/// The mass flow the named draw carries, kg/s - `getSideDrawStream(...).getFlowRate(unit)`.
fn draw_mass_flow(outcome: &ColumnOutcome, tray: usize, phase: SideDrawPhase) -> Option<f64> {
    let draw = match phase {
        SideDrawPhase::Gas => outcome.gas_side_draws.get(tray),
        SideDrawPhase::Liquid => outcome.liquid_side_draws.get(tray),
    }
    .and_then(|draw| draw.as_ref())?;
    draw.mass_flow().ok()
}

fn last_accepted_fraction(accepted: &[f64]) -> f64 {
    accepted.last().copied().unwrap_or(f64::NAN)
}

/// `wasSideDrawFractionAttempted`, with the class's own identity tolerance.
fn was_attempted(attempted: &[f64], fraction: f64) -> bool {
    if !fraction.is_finite() {
        return true;
    }
    attempted
        .iter()
        .any(|prior| (prior - fraction).abs() <= FRACTION_IDENTITY)
}

/// **`selectNextSingleSideDrawCandidate`**: a secant through the accepted probes where they
/// bracket the target, else the multiplicative proposal, else a bounded deterministic grid scan.
///
/// **Only accepted probes may propose.** A rejected candidate's flow is not a measurement, so
/// the class interpolates on the accepted pair alone and steps *past* a rejected fraction
/// rather than through it - which is what keeps a `FALLBACK_PRODUCTS` state out of the
/// controller.
fn select_next(
    specification: &SideDrawFlow,
    accepted_fractions: &[f64],
    accepted_flows: &[f64],
    attempted: &[f64],
    maximum: f64,
) -> f64 {
    if accepted_fractions.is_empty() {
        let origin = attempted.last().copied().unwrap_or(0.0);
        return next_grid_fraction(origin, attempted, maximum);
    }

    let target = specification.target;
    let mut best_index = 0;
    let mut best_residual = f64::INFINITY;
    let mut lower: Option<usize> = None;
    let mut lower_residual = f64::INFINITY;
    let mut upper: Option<usize> = None;
    let mut upper_residual = f64::INFINITY;
    for (index, flow) in accepted_flows.iter().enumerate() {
        let residual = (flow - target).abs();
        if residual < best_residual {
            best_residual = residual;
            best_index = index;
        }
        if *flow <= target && residual < lower_residual {
            lower_residual = residual;
            lower = Some(index);
        }
        if *flow >= target && residual < upper_residual {
            upper_residual = residual;
            upper = Some(index);
        }
    }

    if let (Some(lower), Some(upper)) = (lower, upper) {
        if lower != upper {
            let lower_fraction = accepted_fractions[lower];
            let upper_fraction = accepted_fractions[upper];
            let denominator = accepted_flows[upper] - accepted_flows[lower];
            let interpolated = if denominator.abs() > SCALE_FLOOR {
                lower_fraction
                    + (target - accepted_flows[lower]) * (upper_fraction - lower_fraction)
                        / denominator
            } else {
                0.5 * (lower_fraction + upper_fraction)
            };
            let interpolated = interpolated.clamp(0.0, maximum);
            if !was_attempted(attempted, interpolated) {
                return interpolated;
            }
        }
    }

    let best_fraction = accepted_fractions[best_index];
    let multiplicative =
        multiplicative_candidate(best_fraction, target, accepted_flows[best_index], maximum);
    if !was_attempted(attempted, multiplicative) {
        return multiplicative;
    }
    let nearest_rejected = nearest_rejected_fraction(best_fraction, attempted, accepted_fractions);
    if nearest_rejected.is_finite()
        && (nearest_rejected - best_fraction).abs() < CANDIDATE_SCAN_STEP
    {
        let opposite =
            best_fraction + (5.0 * CANDIDATE_SCAN_STEP).copysign(best_fraction - nearest_rejected);
        let opposite = opposite.clamp(0.0, maximum);
        if !was_attempted(attempted, opposite) {
            return opposite;
        }
    }
    next_grid_fraction(best_fraction, attempted, maximum)
}

/// `calculateSideDrawMultiplicativeCandidate`, bounded.
fn multiplicative_candidate(fraction: f64, target: f64, actual: f64, maximum: f64) -> f64 {
    if target <= SCALE_FLOOR {
        return 0.0;
    }
    let candidate = if actual.abs() <= SCALE_FLOOR {
        fraction + CANDIDATE_SCAN_STEP
    } else {
        fraction * target / actual
    };
    candidate.clamp(0.0, maximum)
}

/// `findNearestRejectedSideDrawFraction`.
fn nearest_rejected_fraction(origin: f64, attempted: &[f64], accepted: &[f64]) -> f64 {
    let mut nearest = f64::NAN;
    let mut nearest_distance = f64::INFINITY;
    for fraction in attempted {
        if was_attempted(accepted, *fraction) {
            continue;
        }
        let distance = (fraction - origin).abs();
        if distance < nearest_distance {
            nearest_distance = distance;
            nearest = *fraction;
        }
    }
    nearest
}

/// `nextUntriedSideDrawGridFractionAround`: the deterministic scan, above then below.
fn next_grid_fraction(origin: f64, attempted: &[f64], maximum: f64) -> f64 {
    let first_upper =
        ((origin + FRACTION_IDENTITY) / CANDIDATE_SCAN_STEP).ceil() * CANDIDATE_SCAN_STEP;
    let first_lower =
        ((origin - FRACTION_IDENTITY) / CANDIDATE_SCAN_STEP).floor() * CANDIDATE_SCAN_STEP;
    let maximum_steps = (maximum / CANDIDATE_SCAN_STEP).ceil() as usize + 1;
    for step in 0..=maximum_steps {
        let offset = step as f64 * CANDIDATE_SCAN_STEP;
        let upper = first_upper + offset;
        if upper <= maximum + FRACTION_IDENTITY {
            let upper = upper.clamp(0.0, maximum);
            if !was_attempted(attempted, upper) {
                return upper;
            }
        }
        let lower = first_lower - offset;
        if lower >= -FRACTION_IDENTITY {
            let lower = lower.clamp(0.0, maximum);
            if !was_attempted(attempted, lower) {
                return lower;
            }
        }
    }
    f64::NAN
}

/// `appendSideDrawCandidateHistory`: the public trace, one entry per candidate.
fn append_history(
    history: &mut String,
    iteration: usize,
    fraction: f64,
    flow: f64,
    accepted: bool,
) {
    if !history.is_empty() {
        history.push_str("; ");
    }
    history.push_str(&format!(
        "#{iteration} fraction={fraction}, flow={flow}, accepted={accepted}"
    ));
}
