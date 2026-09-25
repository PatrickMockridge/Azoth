//! The tear's fixed point: `Recycle`'s four residuals, its `solved()` and its acceleration.
//!
//! **A recycle is a declared tear and this is what converging one means.** The class's `Recycle`
//! mixes its inlets, flashes, and compares the result against the *previous* iteration's - so what
//! it converges is the difference between two successive states of the same stream, not a
//! residual in the equations. Four quantities are compared, each with its own tolerance, and all
//! four default to `1e-2`.
//!
//! **The four are not in the same units and the class's own javadoc says so.** `flowBalanceCheck`
//! returns an absolute `kg/s` difference for a stream below `1 kg/s` and a *percentage* above it,
//! so one tolerance means a different thing on a small loop and a large one. That is why
//! `absoluteFlowTolerance` exists as an alternative criterion - and why it is off by default.
//!
//! **The class's own numbers are carried rather than invented**: the four tolerances, the
//! iteration cap, the flow floor and the Wegstein bounds are all `Recycle`'s defaults, and
//! `specs/flowsheets/`'s `[[recycles]]` entries carry them so a flowsheet's tear and a case's tear
//! are the same machine.
//!
//! **The low-flow cutoff is `deactivateOnLowFlow`, and it is a switch rather than a tolerance.** A
//! tear whose inlet carries less than `minimumFlow` kg/hr is marked **inactive**, its four
//! residuals are reported as exactly zero, and `solved()` returns true - so a loop that carries
//! nothing ends after one pass rather than running until the zero-flow floor closes it. That is
//! what the shipped `demo.toml` does, and it is why the NeqSim capture beside it records
//! `recycle_iterations=1`.
//!
//! **The two accelerations are carried to different depths, and the difference is the class's own
//! write-back.** `apply_composition` is the port of `applyStreamValues`, and what it writes is the
//! **composition and nothing else** - so an acceleration can only move a composition, and a loop
//! whose slow variable is a flow is not accelerated at all. Beyond that, the class writes through
//! `Component.setx` on **both phases**, and a two-phase system does not read back what it was
//! written (measured in `AccelerationProbe`: `[1.0, 0.0]` written, `[1.0, 0.4]` read). Wegstein's
//! step is bounded - `q` is clamped to `[-5, 0]` - so that gap moves it by `1.9e-11` relative, and
//! it is wired and held to the flowsheet capture. Broyden's step is unbounded, so the same gap
//! moves it by two orders of magnitude, and the declaration is **refused** by
//! `Recycle::unsupported_acceleration` rather than run. The arithmetic of both is transcribed and
//! held to `captures/process_acceleration.tsv`, step by step.
//!
//! **What is not ported**: `runTransient`, the adaptive-acceleration ladder, and the
//! `RecycleController`'s priority levels. The first is a transient role a steady-state flowsheet
//! never takes; the second's auto-upgrade path is reached by `applyAutoAdaptiveAcceleration`, which
//! nothing in the sequential path calls; the third matters only for chained loops.

use azoth_core::{AzothError, Result};

use crate::stream::Stream;

/// `Recycle`'s own defaults, as its fields declare them.
pub const FLOW_TOLERANCE: f64 = 1e-2;
/// `Recycle.compositionTolerance`.
pub const COMPOSITION_TOLERANCE: f64 = 1e-2;
/// `Recycle.temperatureTolerance`.
pub const TEMPERATURE_TOLERANCE: f64 = 1e-2;
/// `Recycle.pressureTolerance`.
pub const PRESSURE_TOLERANCE: f64 = 1e-2;
/// `Recycle.maxIterations`.
pub const MAX_ITERATIONS: u32 = 10;
/// `Recycle.minimumFlow`, kg/hr. A loop below this carries nothing to converge.
pub const MINIMUM_FLOW_KG_PER_HR: f64 = 1e-20;
/// `Recycle.wegsteinQMin` - the floor on the q-factor, and the value a slope of one takes.
pub const WEGSTEIN_Q_MIN: f64 = -5.0;
/// `Recycle.wegsteinQMax`.
pub const WEGSTEIN_Q_MAX: f64 = 0.0;
/// `Recycle.wegsteinDelayIterations`: acceleration starts after this many passes.
pub const WEGSTEIN_DELAY_ITERATIONS: u32 = 2;

/// Which of `AccelerationMethod`'s three runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Acceleration {
    /// `x_{n+1} = g(x_n)`. **The class's default**, and the slowest.
    #[default]
    DirectSubstitution,
    /// The secant extrapolation, applied after `wegstein_delay_iterations`.
    Wegstein,
    /// Broyden's quasi-Newton update, through `BroydenAccelerator`.
    Broyden,
}

impl Acceleration {
    /// The acceleration a declaration's name selects.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] on a name that is none of the three, which the class's own
    /// enum makes impossible - and a declaration can carry anything, so it is refused here.
    pub fn named(name: &str) -> Result<Self> {
        match name {
            "direct_substitution" => Ok(Self::DirectSubstitution),
            "wegstein" => Ok(Self::Wegstein),
            "broyden" => Ok(Self::Broyden),
            other => Err(AzothError::invalid_input(
                "acceleration_method",
                format!("`{other}` is not direct_substitution, wegstein or broyden"),
            )),
        }
    }
}

/// A tear's convergence settings, as the class defaults them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RecycleSettings {
    /// `flowTolerance`, in the mixed unit `flow_residual` returns.
    pub flow_tolerance: f64,
    /// `compositionTolerance`, a sum of absolute mole-fraction differences.
    pub composition_tolerance: f64,
    /// `temperatureTolerance`, a sum of percentage changes.
    pub temperature_tolerance: f64,
    /// `pressureTolerance`, a sum of percentage changes.
    pub pressure_tolerance: f64,
    /// `maxIterations`.
    pub max_iterations: u32,
    /// `minimumFlow`, kg/hr.
    pub minimum_flow_kg_per_hr: f64,
    /// Which acceleration runs.
    pub acceleration: Acceleration,
}

impl Default for RecycleSettings {
    fn default() -> Self {
        Self {
            flow_tolerance: FLOW_TOLERANCE,
            composition_tolerance: COMPOSITION_TOLERANCE,
            temperature_tolerance: TEMPERATURE_TOLERANCE,
            pressure_tolerance: PRESSURE_TOLERANCE,
            max_iterations: MAX_ITERATIONS,
            minimum_flow_kg_per_hr: MINIMUM_FLOW_KG_PER_HR,
            acceleration: Acceleration::DirectSubstitution,
        }
    }
}

/// What one pass left behind: the four residuals and the absolute flow change beside them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Residuals {
    /// `flowBalanceCheck`: an absolute kg/s difference below `1 kg/s`, a percentage at or above.
    pub flow: f64,
    /// `absoluteFlowChange`, kg/hr - the quantity `absoluteFlowTolerance` would compare.
    pub absolute_flow_change_kg_per_hr: f64,
    /// `compositionBalanceCheck`: the summed absolute mole-fraction difference.
    pub composition: f64,
    /// `temperatureBalanceCheck`: the summed percentage change over the phases.
    pub temperature: f64,
    /// `pressureBalanceCheck`: the same over the phases.
    pub pressure: f64,
}

impl Residuals {
    /// The residuals of a tear the low-flow cutoff has switched off - all five zero.
    ///
    /// **`deactivateOnLowFlow` sets them, and it sets them to zero rather than measuring
    /// them.** A loop below `minimumFlow` carries no physically meaningful inventory, so the class
    /// reports nothing to converge rather than reporting a residual on a negligible stream. That
    /// is why the numbers here are exactly zero and not merely small: an empty tear's true
    /// residuals are zero too, and the two are different statements - one is measured, this one is
    /// declared.
    #[must_use]
    pub fn deactivated() -> Self {
        Self {
            flow: 0.0,
            absolute_flow_change_kg_per_hr: 0.0,
            composition: 0.0,
            temperature: 0.0,
            pressure: 0.0,
        }
    }
}

/// The four residuals between this pass's mixed stream and the previous pass's.
///
/// **The port carries one phase where the class sums over `getNumberOfPhases()`.** A `Stream` here
/// is a single state; the class's is a flashed system whose temperature and pressure checks add a
/// term per phase. For a one-phase stream the sums are the terms, which is the case every
/// flowsheet in `specs/flowsheets/` is in - and a two-phase tear would need the flash the session
/// does not yet carry, which is named as owed rather than guessed at.
///
/// # Errors
/// [`AzothError::InvalidInput`] if the two streams carry different component sets, which the class
/// answers with the constant `10.0` rather than an error.
pub fn residuals(mixed: &Stream, previous: &Stream) -> Result<Residuals> {
    if mixed.components.len() != previous.components.len() {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "a recycle compares two states of one stream, and these carry {} and {} components",
                mixed.components.len(),
                previous.components.len()
            ),
        ));
    }

    let mixed_kg_per_hr = mass_flow_kg_per_hr(mixed)?;
    let previous_kg_per_hr = mass_flow_kg_per_hr(previous)?;
    let absolute_flow_change_kg_per_hr = (mixed_kg_per_hr - previous_kg_per_hr).abs();

    // **The mixed unit, which is the class's own.** Below one kg/s the residual is the absolute
    // difference; at or above it, the percentage of the current flow.
    let mixed_kg_per_s = mixed_kg_per_hr / 3600.0;
    let flow = if mixed_kg_per_s < 1.0 {
        (mixed_kg_per_s - previous_kg_per_hr / 3600.0).abs()
    } else {
        (mixed_kg_per_s - previous_kg_per_hr / 3600.0).abs() / mixed_kg_per_s * 100.0
    };

    let composition: f64 = mixed
        .z
        .iter()
        .zip(&previous.z)
        .map(|(a, b)| (a - b).abs())
        .sum();
    let temperature = relative_percent(mixed.t.value, previous.t.value);
    let pressure = relative_percent(mixed.p.value, previous.p.value);

    Ok(Residuals {
        flow,
        absolute_flow_change_kg_per_hr,
        composition,
        temperature,
        pressure,
    })
}

/// `|a - b| / b * 100`, the class's own percentage.
fn relative_percent(value: f64, previous: f64) -> f64 {
    if previous == 0.0 {
        return if value == 0.0 { 0.0 } else { f64::INFINITY };
    }
    ((value - previous) / previous).abs() * 100.0
}

/// A stream's mass flow in kg/hr, which is the unit both the floor and the change are in.
///
/// # Errors
/// [`AzothError::PropertyUnavailable`] if a component carries no molar mass - the class would
/// divide by a zero it never reads, and this refuses instead.
pub fn mass_flow_kg_per_hr(stream: &Stream) -> Result<f64> {
    Ok(stream.n * 3600.0 * stream.molar_mass()?.value)
}

/// `Recycle.solved()`.
///
/// **All four residuals, and `iterations > 1`.** The last is the class's "legacy second
/// observation": one pass compares a state against itself, since the previous-iteration stream
/// starts as the current one, so a first pass would always read as converged. It is not a
/// formality - it is what stops a recycle from declaring victory before it has run.
///
/// **The zero-flow floor is checked first and needs both streams.** A loop carrying less than
/// `max(minimum_flow, 1e-20)` kg/hr has nothing to converge; the class requires *both* the outlet
/// and the previous state below the floor, so a loop collapsing from a real flow is not mistaken
/// for an empty one.
///
/// **An inactive tear is solved whatever its residuals say**, and the test is
/// `!active && iterations > 0` rather than `!active` - a recycle that has been switched off
/// *before it ever ran* has nothing to report and is not solved. `RecycleController.solvedAll()`
/// skips a deactivated recycle for the same reason, so the loop's exit condition and the tear's
/// own answer agree about a loop that carries nothing.
#[must_use]
pub fn solved(
    residuals: &Residuals,
    settings: &RecycleSettings,
    outlet_kg_per_hr: f64,
    previous_kg_per_hr: f64,
    iterations: u32,
    active: bool,
) -> bool {
    if !active && iterations > 0 {
        return true;
    }

    let zero_flow_floor = settings.minimum_flow_kg_per_hr.max(1e-20);
    if outlet_kg_per_hr < zero_flow_floor && previous_kg_per_hr < zero_flow_floor && iterations > 1
    {
        return true;
    }

    let flow_converged = residuals.flow.abs() < settings.flow_tolerance;
    let composition_converged = residuals.composition.abs() < settings.composition_tolerance;
    let temperature_converged = residuals.temperature.abs() < settings.temperature_tolerance;
    let pressure_converged = residuals.pressure.abs() < settings.pressure_tolerance;

    flow_converged
        && composition_converged
        && temperature_converged
        && pressure_converged
        && iterations > 1
}

/// `Recycle.extractStreamValues`: the tear's state as the acceleration vector sees it.
///
/// `[T, P, n_mol_per_s, x_0 … x_{n-1}]`. **The order is the class's and it is load-bearing**:
/// `apply_composition` reads the fractions from index three, so a vector in another order would
/// accelerate the temperature into the composition slot.
#[must_use]
pub fn extract(stream: &Stream) -> Vec<f64> {
    let mut values = Vec::with_capacity(3 + stream.z.len());
    values.push(stream.t.value);
    values.push(stream.p.value);
    values.push(stream.n);
    values.extend_from_slice(&stream.z);
    values
}

/// `Recycle.applyStreamValues`: write an accelerated vector back onto a stream.
///
/// **The composition is the only field that moves, and that is the class's own choice** - its
/// comment beside the method says "T, P, and flow are handled elsewhere". So the acceleration can
/// only accelerate a *composition*: a loop whose slow variable is a flow, a temperature or a
/// pressure is not affected by it at all. That is a real limit rather than a tidying, and it is
/// the reason this port can apply the class's acceleration faithfully and still show no
/// difference on a flowsheet whose tear accumulates in its flow.
///
/// Fractions below zero are **clamped rather than refused**, and a vector whose fractions sum to
/// at most `1e-15` is **skipped rather than normalised** - both the class's own guards, and both
/// reachable from an accelerated step that overshoots past a physical composition.
#[must_use]
pub fn apply_composition(stream: &Stream, values: &[f64]) -> Stream {
    let mut stream = stream.clone();
    let components = stream.z.len();
    if values.len() < 3 + components {
        return stream;
    }

    let fractions: Vec<f64> = values[3..3 + components]
        .iter()
        .map(|fraction| fraction.max(0.0))
        .collect();
    let sum: f64 = fractions.iter().sum();
    if sum > 1e-15 {
        stream.z = fractions.iter().map(|fraction| fraction / sum).collect();
    }
    stream
}

/// **Wegstein's acceleration**, per variable.
///
/// `x_{n+1} = q g(x_n) + (1 - q) x_n` with `q = s / (s - 1)` from the secant slope
/// `s = Δoutput / Δinput`, clamped to `[wegstein_q_min, wegstein_q_max]` - `[-5, 0]` by the
/// class's defaults, so `q` is never positive and the step never overshoots the fixed point.
///
/// **Three of the class's guards are load-bearing and none is tidied.** A slope of exactly one is
/// the *diverging* case and takes `q = wegstein_q_min` rather than the division's infinity; a
/// `Δinput` below `1e-15` is a slope of zero rather than a division by it; and the first pass has
/// no previous values, so `q` is all zeros and the acceleration is the identity - which is what
/// makes the delay unnecessary on the first pass and is why the class's own guard is a *check*
/// (`iterations > wegsteinDelayIterations`) rather than the only thing preventing it.
///
/// **A flat secant gives `q = 0`, and at `q = 0` the formula returns the input** - the step does
/// not move. The class's comment beside `wegsteinQFactors` says an array of zeros is "direct
/// substitution", which is true of the *array* and not of the formula: `q = 0` is a stalled step,
/// and the path that actually substitutes directly is the early return above.
#[must_use]
pub fn wegstein(
    input: &[f64],
    output: &[f64],
    previous_input: Option<&[f64]>,
    previous_output: Option<&[f64]>,
    q_min: f64,
    q_max: f64,
) -> WegsteinStep {
    let identity = |output: &[f64]| WegsteinStep {
        values: output.to_vec(),
        q_factors: vec![0.0; output.len()],
    };

    let (Some(previous_input), Some(previous_output)) = (previous_input, previous_output) else {
        return identity(output);
    };
    if previous_input.len() != input.len() || previous_output.len() != output.len() {
        return identity(output);
    }

    let mut values = Vec::with_capacity(input.len());
    let mut q_factors = Vec::with_capacity(input.len());
    for i in 0..input.len() {
        let delta_input = input[i] - previous_input[i];
        let delta_output = output[i] - previous_output[i];
        let slope = if delta_input.abs() > 1e-15 {
            delta_output / delta_input
        } else {
            0.0
        };
        let q = if (slope - 1.0).abs() > 1e-10 {
            slope / (slope - 1.0)
        } else {
            q_min
        };
        let q = q.clamp(q_min, q_max);
        q_factors.push(q);
        values.push(q * output[i] + (1.0 - q) * input[i]);
    }
    WegsteinStep { values, q_factors }
}

/// What one Wegstein call answers: the step, and the `q` it took.
///
/// **The `q` factors are published because the class publishes them** (`getWegsteinQFactors`),
/// and because they are the only thing that makes the method checkable step by step: the step
/// itself is one arithmetic expression, while the `q` it used records *which* slope produced it -
/// so a port that took the branch one call early, or clamped at the wrong place, shows up in the
/// factors rather than in a value that merely looks reasonable.
#[derive(Debug, Clone, PartialEq)]
pub struct WegsteinStep {
    /// `q * g + (1 - q) * x`, per variable.
    pub values: Vec<f64>,
    /// `wegsteinQFactors` after the call: all zeros on the identity path.
    pub q_factors: Vec<f64>,
}

/// `BroydenAccelerator`: Broyden's "good" method on the fixed-point residual.
///
/// **The inverse Jacobian is the state, and it is carried across passes.** Broyden's method
/// approximates the Jacobian of the fixed-point map and applies rank-one corrections from the
/// secant condition, which avoids recomputing derivatives while still giving Newton-like steps -
/// and because the correction is applied to the *inverse*, each step is a matrix-vector product
/// rather than a solve.
///
/// **The delay is inside the accelerator and not at the call site**, unlike Wegstein's. The first
/// two calls return the map's own output and store the pair, which is what seeds the secant; the
/// first real step is the third call. The class's comment says the initial `-I` "corresponds to
/// assuming J ≈ -I initially (typical for convergent iterations)", which is exactly what makes
/// the first step a direct substitution.
///
/// **The three guards are the class's and none is tidied**: a `delta_x` whose norm is at or below
/// `1e-15` skips the update rather than dividing by it; a Sherman-Morrison denominator below it
/// skips it too; and a step longer than `max_step_size` is scaled down rather than refused.
#[derive(Debug, Clone, PartialEq)]
pub struct BroydenAccelerator {
    /// The inverse-Jacobian approximation, `-I` until the first update.
    inverse_jacobian: Vec<Vec<f64>>,
    /// The input vector of the previous call, and its residual `g(x) - x`.
    previous: Option<(Vec<f64>, Vec<f64>)>,
    /// Calls made, which is what the delay counts.
    iteration_count: u32,
    /// `delayIterations`: calls that substitute directly before the first step.
    pub delay_iterations: u32,
    /// `relaxationFactor`, applied to the Newton step.
    pub relaxation_factor: f64,
    /// `maxStepSize`, the ceiling on the step's length.
    pub max_step_size: f64,
}

impl Default for BroydenAccelerator {
    fn default() -> Self {
        Self {
            inverse_jacobian: Vec::new(),
            previous: None,
            iteration_count: 0,
            delay_iterations: 2,
            relaxation_factor: 1.0,
            max_step_size: f64::MAX,
        }
    }
}

/// The smallest `delta_x` norm and Sherman-Morrison denominator the update survives.
const BROYDEN_EPSILON: f64 = 1e-15;

impl BroydenAccelerator {
    /// A fresh accelerator, its dimension unknown until the first call.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// How many calls have been made, which the delay is measured in.
    #[must_use]
    pub fn iteration_count(&self) -> u32 {
        self.iteration_count
    }

    /// The fixed-point function's output, or an accelerated next iterate.
    ///
    /// `current_x` is the previous pass's published vector and `function_output` is `g(x)` for
    /// it. The returned vector is what the class then writes onto the tear's composition.
    pub fn accelerate(&mut self, current_x: &[f64], function_output: &[f64]) -> Vec<f64> {
        let n = current_x.len();

        // `initialize(n)` when the dimension changes, which the class reaches on the first call
        // because it creates the accelerator with dimension zero.
        if self.inverse_jacobian.len() != n {
            self.inverse_jacobian = vec![vec![0.0; n]; n];
            for (i, row) in self.inverse_jacobian.iter_mut().enumerate() {
                row[i] = -1.0;
            }
            self.previous = None;
            self.iteration_count = 0;
        }

        self.iteration_count += 1;

        // `f(x) = g(x) - x`, the residual the secant condition is written on.
        let current_f: Vec<f64> = function_output
            .iter()
            .zip(current_x)
            .map(|(output, x)| output - x)
            .collect();

        // The delay: substitute directly and keep the pair for the next call's slope.
        if self.iteration_count <= self.delay_iterations || self.previous.is_none() {
            self.previous = Some((current_x.to_vec(), current_f));
            return function_output.to_vec();
        }

        let (previous_x, previous_f) = self.previous.as_ref().expect("the delay period set it");
        let delta_x: Vec<f64> = current_x
            .iter()
            .zip(previous_x)
            .map(|(now, before)| now - before)
            .collect();
        let delta_f: Vec<f64> = current_f
            .iter()
            .zip(previous_f)
            .map(|(now, before)| now - before)
            .collect();

        if norm(&delta_x) > BROYDEN_EPSILON {
            self.update_inverse_jacobian(&delta_x, &delta_f);
        }

        // The Newton step, `-B^-1 f`. **The sign is in the matrix**: the matrix-vector product is
        // taken as-is and added, because the inverse starts at `-I` and the update preserves the
        // negation.
        let mut step: Vec<f64> = multiply(&self.inverse_jacobian, &current_f)
            .iter()
            .map(|entry| entry * self.relaxation_factor)
            .collect();

        let step_norm = norm(&step);
        if step_norm > self.max_step_size {
            let scale = self.max_step_size / step_norm;
            for entry in &mut step {
                *entry *= scale;
            }
        }

        self.previous = Some((current_x.to_vec(), current_f));
        current_x
            .iter()
            .zip(&step)
            .map(|(x, step)| x + step)
            .collect()
    }

    /// Broyden's "good" method, through Sherman-Morrison.
    fn update_inverse_jacobian(&mut self, delta_x: &[f64], delta_f: &[f64]) {
        let n = self.inverse_jacobian.len();
        let binv_delta_f = multiply(&self.inverse_jacobian, delta_f);
        let denominator = dot(delta_x, &binv_delta_f);
        if denominator.abs() < BROYDEN_EPSILON {
            return;
        }

        let numerator: Vec<f64> = delta_x
            .iter()
            .zip(&binv_delta_f)
            .map(|(x, binv)| x - binv)
            .collect();

        let mut delta_x_transpose_binv = vec![0.0; n];
        for (j, entry) in delta_x_transpose_binv.iter_mut().enumerate() {
            for (i, x) in delta_x.iter().enumerate() {
                *entry += x * self.inverse_jacobian[i][j];
            }
        }

        for (i, row) in self.inverse_jacobian.iter_mut().enumerate() {
            for (j, entry) in row.iter_mut().enumerate() {
                *entry += numerator[i] * delta_x_transpose_binv[j] / denominator;
            }
        }
    }
}

/// A matrix times a vector.
fn multiply(matrix: &[Vec<f64>], vector: &[f64]) -> Vec<f64> {
    matrix
        .iter()
        .map(|row| dot(row, vector))
        .collect::<Vec<f64>>()
}

/// A dot product.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

/// The Euclidean norm.
fn norm(v: &[f64]) -> f64 {
    dot(v, v).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn settings() -> RecycleSettings {
        RecycleSettings::default()
    }

    fn residual(flow: f64, composition: f64, temperature: f64, pressure: f64) -> Residuals {
        Residuals {
            flow,
            absolute_flow_change_kg_per_hr: 0.0,
            composition,
            temperature,
            pressure,
        }
    }

    /// **All four, and the second observation.** One residual inside its tolerance is not
    /// convergence, and neither is all four on the first pass - the class compares a state against
    /// itself there.
    #[test]
    fn convergence_needs_four_residuals_and_a_second_pass() {
        let inside = residual(0.0, 0.0, 0.0, 0.0);
        assert!(solved(&inside, &settings(), 1000.0, 1000.0, 2, true));
        assert!(
            !solved(&inside, &settings(), 1000.0, 1000.0, 1, true),
            "the first pass has nothing to compare against"
        );
        for outside in [
            residual(1.0, 0.0, 0.0, 0.0),
            residual(0.0, 1.0, 0.0, 0.0),
            residual(0.0, 0.0, 1.0, 0.0),
            residual(0.0, 0.0, 0.0, 1.0),
        ] {
            assert!(
                !solved(&outside, &settings(), 1000.0, 1000.0, 5, true),
                "one residual of 1 exceeds the class's 1e-2 tolerance"
            );
        }
    }

    /// **The zero-flow floor needs both streams below it**, and it is checked before the
    /// residuals - so an empty loop is solved and a collapsing one is not mistaken for it.
    #[test]
    fn an_empty_loop_is_solved_and_a_collapsing_one_is_not() {
        let wild = residual(50.0, 50.0, 50.0, 50.0);
        assert!(solved(&wild, &settings(), 1e-30, 1e-30, 5, true));
        assert!(!solved(&wild, &settings(), 1e-30, 10.0, 5, true));
        assert!(!solved(&wild, &settings(), 10.0, 1e-30, 5, true));
    }

    /// Wegstein's three guards, each a case the division would not survive.
    #[test]
    fn wegstein_clamps_the_q_factor_and_survives_its_own_edges() {
        // A slope of one is the diverging case: q takes the floor, not infinity.
        let one = wegstein(&[0.0], &[1.0], Some(&[-1.0]), Some(&[0.0]), -5.0, 0.0);
        assert_eq!(one.q_factors, vec![-5.0]);
        assert_eq!(one.values, vec![-5.0 * 1.0 + 6.0 * 0.0]);

        // No previous values is the identity, which is direct substitution - and the q factors
        // are all zeros, which is the array the class's own comment calls "direct substitution".
        let first = wegstein(&[1.0, 2.0], &[3.0, 4.0], None, None, -5.0, 0.0);
        assert_eq!(first.values, vec![3.0, 4.0]);
        assert_eq!(first.q_factors, vec![0.0, 0.0]);

        // **A slope of zero gives q = 0, and at q = 0 the formula returns the *input*.**
        // `x = q*g + (1-q)*x` is `x` there, so a flat secant stalls rather than substituting - it
        // is not the "all zeros = direct substitution" the class's own comment beside
        // `wegsteinQFactors` claims, and that comment is about the arrays and not the formula.
        let flat = wegstein(&[1.0], &[2.0], Some(&[0.0]), Some(&[2.0]), -5.0, 0.0);
        assert_eq!(flat.values, vec![1.0], "q = 0 does not move");
        assert_eq!(flat.q_factors, vec![0.0]);

        // And the clamp: a slope of 0.5 gives q = 0.5/(0.5 - 1) = -1, which is inside the
        // bounds, so the step is `q*g + (1-q)*x` = `-1.5 + 2.0` = `0.5` - a point the original
        // never visited, which is what acceleration is for.
        let half = wegstein(&[1.0], &[1.5], Some(&[0.0]), Some(&[1.0]), -5.0, 0.0);
        assert!(
            (half.values[0] - 0.5).abs() < 1e-12,
            "got {}",
            half.values[0]
        );
        assert!(
            (half.q_factors[0] + 1.0).abs() < 1e-12,
            "got {}",
            half.q_factors[0]
        );
    }

    /// The percentage form of the temperature and pressure residuals, and the mixed unit of the
    /// flow one - which is the class's own and the reason one tolerance means two things.
    #[test]
    fn the_percentage_residuals_are_percentages() {
        assert!((relative_percent(101.0, 100.0) - 1.0).abs() < 1e-12);
        assert!((relative_percent(100.0, 100.0)).abs() < 1e-12);
        assert!((relative_percent(99.0, 100.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn the_acceleration_names_are_the_class_s_own() {
        assert_eq!(
            Acceleration::named("direct_substitution").expect("known"),
            Acceleration::DirectSubstitution
        );
        assert_eq!(
            Acceleration::named("wegstein").expect("known"),
            Acceleration::Wegstein
        );
        assert_eq!(
            Acceleration::named("broyden").expect("known"),
            Acceleration::Broyden
        );
        assert!(Acceleration::named("newton").is_err());
        assert_eq!(Acceleration::default(), Acceleration::DirectSubstitution);
    }
}
