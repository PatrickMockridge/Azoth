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
) -> Vec<f64> {
    let (Some(previous_input), Some(previous_output)) = (previous_input, previous_output) else {
        return output.to_vec();
    };
    if previous_input.len() != input.len() || previous_output.len() != output.len() {
        return output.to_vec();
    }

    (0..input.len())
        .map(|i| {
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
            q * output[i] + (1.0 - q) * input[i]
        })
        .collect()
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
        assert_eq!(one, vec![-5.0 * 1.0 + 6.0 * 0.0]);

        // No previous values is the identity, which is direct substitution.
        assert_eq!(
            wegstein(&[1.0, 2.0], &[3.0, 4.0], None, None, -5.0, 0.0),
            vec![3.0, 4.0]
        );

        // **A slope of zero gives q = 0, and at q = 0 the formula returns the *input*.**
        // `x = q*g + (1-q)*x` is `x` there, so a flat secant stalls rather than substituting - it
        // is not the "all zeros = direct substitution" the class's own comment beside
        // `wegsteinQFactors` claims, and that comment is about the arrays and not the formula.
        let flat = wegstein(&[1.0], &[2.0], Some(&[0.0]), Some(&[2.0]), -5.0, 0.0);
        assert_eq!(flat, vec![1.0], "q = 0 does not move");

        // And the clamp: a slope of 0.5 gives q = 0.5/(0.5 - 1) = -1, which is inside the
        // bounds, so the step is `q*g + (1-q)*x` = `-1.5 + 2.0` = `0.5` - a point the original
        // never visited, which is what acceleration is for.
        let half = wegstein(&[1.0], &[1.5], Some(&[0.0]), Some(&[1.0]), -5.0, 0.0);
        assert!((half[0] - 0.5).abs() < 1e-12, "got {}", half[0]);
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
