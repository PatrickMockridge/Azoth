//! The counter-current profile: the class's `solveFixedPointProfile`.

use std::f64::consts::PI;

use azoth_core::Result;

use super::fallbacks::Fallbacks;
use super::step::{SegmentComputation, SegmentResult, calculate_segment, component_moles};
use super::transport::SnapshotSettings;
use crate::stream::Stream;

/// One pass's answer, as the class's `CounterCurrentSolution` holds it.
#[derive(Debug, Clone)]
struct CounterCurrentSolution {
    gas_outlet: Stream,
    liquid_outlet: Stream,
    liquid_leaving: Vec<Stream>,
    segments: Vec<SegmentResult>,
}

/// What a column's profile solve answered.
#[derive(Debug, Clone)]
pub struct ColumnOutcome {
    /// The gas leaving the top segment.
    pub gas_outlet: Stream,
    /// The liquid leaving the **bottom** segment, which is the column's liquid outlet.
    pub liquid_outlet: Stream,
    /// One record per segment, bottom first.
    pub segments: Vec<SegmentResult>,
    /// The iterations taken.
    pub iterations: usize,
    /// The outlet residual at the last iteration, in mol/s.
    pub convergence_residual: f64,
    /// The sum of the per-component transfers' magnitudes, mol/s.
    pub total_absolute_molar_transfer: f64,
    /// The per-component transfer totals, in the order the segments produced them.
    pub component_transfer_totals: Vec<(String, f64)>,
    /// Whether the gate was met - and a bed of no height is a *success*, not a miss.
    pub converged: bool,
    /// Which of the class's constants stood in for a missing property.
    pub fallbacks: Fallbacks,
}

/// The loop's parameters.
#[derive(Debug, Clone)]
pub struct ProfileSettings {
    /// The packed height, in metres.
    pub packed_height: f64,
    /// The axial slices.
    pub number_of_segments: usize,
    /// The convergence gate, mol/s.
    pub tolerance: f64,
    /// The iteration cap.
    pub max_iterations: usize,
}

/// `solveFixedPointProfile`.
///
/// **On exhaustion the class accepts its last iterate silently**, and so does this: the loop
/// ends, `acceptSolution` runs, and the answer is whatever the last pass produced. What makes
/// that visible rather than quiet is the outcome's `converged`, `iterations` and
/// `convergence_residual`, which the result carries.
pub fn solve_fixed_point_profile(
    gas_in: &Stream,
    liquid_in: &Stream,
    profile: &ProfileSettings,
    settings: &SnapshotSettings,
    transfer_components: &[String],
    matrix_model: bool,
) -> Result<ColumnOutcome> {
    let segment_height = profile.packed_height / profile.number_of_segments as f64;
    let segment_volume =
        PI * settings.column_diameter * settings.column_diameter / 4.0 * segment_height;

    // `initializeLiquidProfile`: every segment's entering liquid is the feed.
    let mut liquid_entering = vec![liquid_in.clone(); profile.number_of_segments];
    let mut previous_gas: Option<Stream> = None;
    let mut previous_liquid: Option<Stream> = None;
    let mut solution = None;
    let mut iterations = 0;
    let mut residual = f64::INFINITY;
    let mut converged = false;

    for iteration in 1..=profile.max_iterations {
        let pass = run_one_profile_iteration(
            gas_in,
            liquid_in,
            &liquid_entering,
            segment_height,
            segment_volume,
            transfer_components,
            settings,
            matrix_model,
        )?;
        residual = outlet_residual(
            previous_gas.as_ref(),
            &pass.gas_outlet,
            previous_liquid.as_ref(),
            &pass.liquid_outlet,
        );
        iterations = iteration;
        // **A bed of no height is a converged solve, not a missed gate**: the class accepts it
        // on the first pass, and warning there would fail its own zero-height test.
        if residual <= profile.tolerance || profile.packed_height == 0.0 {
            converged = true;
            solution = Some(pass);
            break;
        }
        previous_gas = Some(pass.gas_outlet.clone());
        previous_liquid = Some(pass.liquid_outlet.clone());
        liquid_entering = update_liquid_profile(liquid_in, &pass.liquid_leaving);
        solution = Some(pass);
    }

    let solution = solution.expect("the cap is at least one iteration, so a pass ran");
    Ok(accept_solution(solution, iterations, residual, converged))
}

/// `runOneProfileIteration`: one gas pass, bottom to top, against the current liquid profile.
#[allow(clippy::too_many_arguments)]
fn run_one_profile_iteration(
    gas_in: &Stream,
    liquid_in: &Stream,
    liquid_entering: &[Stream],
    segment_height: f64,
    segment_volume: f64,
    transfer_components: &[String],
    settings: &SnapshotSettings,
    matrix_model: bool,
) -> Result<CounterCurrentSolution> {
    let _ = liquid_in;
    let mut gas_current = gas_in.clone();
    let mut liquid_leaving = Vec::with_capacity(liquid_entering.len());
    let mut segments = Vec::with_capacity(liquid_entering.len());
    for (index, entering) in liquid_entering.iter().enumerate() {
        let computation: SegmentComputation = calculate_segment(
            index,
            &gas_current,
            entering,
            segment_height,
            segment_volume,
            transfer_components,
            settings,
            matrix_model,
        )?;
        gas_current = computation.gas;
        liquid_leaving.push(computation.liquid);
        segments.push(computation.result);
    }
    // **The column's liquid outlet is the *bottom* segment's liquid** - `liquidLeaving.get(0)`
    // - and not the last one. It is the easiest line in the class to invert, and inverting it
    // converges to a different column.
    let liquid_outlet = liquid_leaving
        .first()
        .ok_or_else(|| {
            azoth_core::AzothError::invalid_input(
                "number_of_segments",
                "a column of no segments has no profile".to_string(),
            )
        })?
        .clone();
    Ok(CounterCurrentSolution {
        gas_outlet: gas_current,
        liquid_outlet,
        liquid_leaving,
        segments,
    })
}

/// `updateLiquidProfile`: the counter-current shift, with the feed entering at the top.
fn update_liquid_profile(liquid_in: &Stream, liquid_leaving: &[Stream]) -> Vec<Stream> {
    let count = liquid_leaving.len();
    (0..count)
        .map(|segment| {
            if segment + 1 == count {
                liquid_in.clone()
            } else {
                liquid_leaving[segment + 1].clone()
            }
        })
        .collect()
}

/// `calculateOutletResidual`: the largest per-component change between successive outlets.
///
/// **The two sides are compared to their own kind, never across.** The class takes the max of
/// `|gas_before − gas_after|` and `|liquid_before − liquid_after|` per component, each read by
/// name and zero where the system does not carry it. Comparing one pass's *gas* against the
/// other's *liquid* instead is a real defect and not a detail: CO2 in the gas is `1.47` mol/s
/// where the same component in the liquid is `0.0035`, so the residual settles at their
/// difference and the loop runs to its cap on a state that has already converged.
fn outlet_residual(
    previous_gas: Option<&Stream>,
    current_gas: &Stream,
    previous_liquid: Option<&Stream>,
    current_liquid: &Stream,
) -> f64 {
    let (Some(previous_gas), Some(previous_liquid)) = (previous_gas, previous_liquid) else {
        return f64::INFINITY;
    };
    let mut residual: f64 = 0.0;
    for component in components_of(current_gas, current_liquid) {
        residual = residual
            .max((moles_of(previous_gas, &component) - moles_of(current_gas, &component)).abs());
        residual = residual.max(
            (moles_of(previous_liquid, &component) - moles_of(current_liquid, &component)).abs(),
        );
    }
    residual
}

/// The union of two systems' component names, in the class's own order.
fn components_of(gas: &Stream, liquid: &Stream) -> Vec<String> {
    let mut names = gas.components.clone();
    for name in &liquid.components {
        if !names.contains(name) {
            names.push(name.clone());
        }
    }
    names
}

/// A component's molar flow in a system, by name, zero where it is absent.
fn moles_of(stream: &Stream, component: &str) -> f64 {
    stream
        .components
        .iter()
        .position(|name| name == component)
        .map_or(0.0, |index| component_moles(stream)[index])
}

/// `acceptSolution`: fold the pass's per-segment transfers into the column's totals.
fn accept_solution(
    solution: CounterCurrentSolution,
    iterations: usize,
    convergence_residual: f64,
    converged: bool,
) -> ColumnOutcome {
    let mut totals: Vec<(String, f64)> = Vec::new();
    let mut total_absolute = 0.0;
    let mut fallbacks = Fallbacks::default();
    for segment in &solution.segments {
        fallbacks.merge(segment.fallbacks);
        for (component, transfer) in &segment.component_transfer {
            total_absolute += transfer.abs();
            match totals.iter_mut().find(|(name, _)| name == component) {
                Some((_, total)) => *total += transfer,
                None => totals.push((component.clone(), *transfer)),
            }
        }
    }
    ColumnOutcome {
        gas_outlet: solution.gas_outlet,
        liquid_outlet: solution.liquid_outlet,
        segments: solution.segments,
        iterations,
        convergence_residual,
        total_absolute_molar_transfer: total_absolute,
        component_transfer_totals: totals,
        converged,
        fallbacks,
    }
}
