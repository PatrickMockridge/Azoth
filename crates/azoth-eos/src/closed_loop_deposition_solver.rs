//! `eos.closed_loop_deposition_solver` - NeqSim's `ClosedLoopDepositionSolver`.
//!
//! The pipe is re-run at a shrinking effective diameter and the accumulator evaluated at each
//! run; the worst segment's thickness sets the next diameter until the diameter stops moving,
//! the iteration cap is reached, or the line blocks.

use azoth_core::units::{
    Length, MassDensity, MassFlux, Pressure, ThermodynamicTemperature, Time, meters,
};
use azoth_core::{Result, apply_checks};

use crate::model_gen;
use crate::pipe_beggs_and_brills::pipe_beggs_and_brills;
use crate::results::ClosedLoopDepositionSolverResult;
use crate::scale_deposition_accumulator::scale_deposition_accumulator;

/// The lowest tolerance the class's setter admits, in metres.
const MINIMUM_TOLERANCE: f64 = 1.0e-7;

/// The closed-loop deposited-diameter solve over a Beggs-and-Brills line.
///
/// # Errors
/// * Whatever the pipe or the accumulator refuses.
///
/// # Example
/// ```
/// use azoth_core::units::{
///     kilograms_per_cubic_meter, kelvins, kilograms_per_square_meter_second, meters, pascals,
///     seconds,
/// };
/// use azoth_eos::closed_loop_deposition_solver::closed_loop_deposition_solver;
///
/// let r = closed_loop_deposition_solver(
///     &["methane".to_owned(), "CO2".to_owned(), "water".to_owned()],
///     80.0,
///     &[0.85, 0.05, 0.10],
///     pascals(1.0e7),
///     kelvins(333.15),
///     meters(1000.0),
///     meters(0.1),
///     meters(0.0),
///     10,
///     meters(1.0e-5),
///     kilograms_per_cubic_meter(1.0),
///     kilograms_per_cubic_meter(0.5),
///     kilograms_per_cubic_meter(0.1),
///     kilograms_per_cubic_meter(0.01),
///     kilograms_per_cubic_meter(20.0),
///     kilograms_per_cubic_meter(50.0),
///     6.5,
///     seconds(31_557_600.0),
///     kilograms_per_cubic_meter(2700.0),
///     kilograms_per_square_meter_second(1.5844043908606944e-8),
///     0.0,
///     meters(1.0e-4),
///     10,
///     Some("srk"),
/// )?;
/// assert_eq!(r.iterations_taken, 1);
/// assert!(r.converged);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn closed_loop_deposition_solver(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
    length: Length,
    diameter: Length,
    elevation: Length,
    number_of_increments: usize,
    roughness: Length,
    calcium: MassDensity,
    bicarbonate: MassDensity,
    sulphate: MassDensity,
    barium: MassDensity,
    sodium: MassDensity,
    total_dissolved_solids: MassDensity,
    ph: f64,
    service_time: Time,
    scale_density: MassDensity,
    base_rate: MassFlux,
    inhibitor_efficiency: f64,
    tolerance: Length,
    max_iterations: usize,
    eos: Option<&str>,
) -> Result<ClosedLoopDepositionSolverResult> {
    let spec = &model_gen::CLOSED_LOOP_DEPOSITION_SOLVER_SPEC;
    let mut warnings = Vec::new();

    apply_checks(spec.input_checks(), |_| None, &mut warnings)?;

    let tolerance = tolerance.value.max(MINIMUM_TOLERANCE);
    let max_iterations = max_iterations.max(1);
    let original = diameter.value;

    let mut current = original;
    let mut previous = original;
    let mut converged = false;
    let mut iterations_taken = 0;
    let mut final_effective_diameter = 0.0;
    let mut diameter_history = Vec::with_capacity(max_iterations);
    let mut max_thickness_history = Vec::with_capacity(max_iterations);
    let mut velocity_history = Vec::with_capacity(max_iterations);

    for index in 0..max_iterations {
        iterations_taken = u32::try_from(index + 1).unwrap_or(u32::MAX);
        let run = pipe_beggs_and_brills(
            components,
            inlet_n,
            inlet_z,
            inlet_p,
            inlet_t,
            length,
            meters(current),
            elevation,
            number_of_increments,
            roughness,
            eos,
        )?;
        let deposition = scale_deposition_accumulator(
            calcium,
            bicarbonate,
            sulphate,
            barium,
            sodium,
            total_dissolved_solids,
            ph,
            meters(current),
            length,
            &run.temperature_profile,
            &run.pressure_profile,
            service_time,
            scale_density,
            base_rate,
            inhibitor_efficiency,
        )?;

        let thickness = deposition.max_thickness.value;
        diameter_history.push(current);
        max_thickness_history.push(thickness);
        let velocity = &run.mixture_superficial_velocity_profile;
        velocity_history.push(if velocity.is_empty() {
            0.0
        } else {
            velocity.iter().sum::<f64>() / velocity.len() as f64
        });

        let new_diameter = original - 2.0 * thickness;
        if new_diameter <= 0.0 {
            final_effective_diameter = 0.0;
            converged = true;
            break;
        }
        if (new_diameter - previous).abs() < tolerance {
            final_effective_diameter = new_diameter;
            converged = true;
            break;
        }
        previous = current;
        current = new_diameter;
    }
    if !converged {
        final_effective_diameter = current;
    }

    Ok(ClosedLoopDepositionSolverResult {
        iterations_taken,
        converged,
        final_effective_diameter: meters(final_effective_diameter),
        diameter_history,
        max_thickness_history,
        velocity_history,
        warnings,
    })
}
