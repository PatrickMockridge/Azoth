//! The reactive PH flash, from `flashops/reactiveflash/ReactiveMultiphasePHflash.java`.
//!
//! Given a pressure and a total enthalpy, it looks for the temperature: an outer loop on `T`
//! wraps the reactive TP flash, and each pass asks the inner flash what the fluid's enthalpy
//! is at that temperature and steps towards the one that was specified.
//!
//! # The loop is not the one its docstring describes
//!
//! **The class's own docstring says the outer loop is a Newton on `1/T` "following Michelsen
//! 1987". The code is not that.** `solveEnthalpySpec` is a secant on `T` with a bisection
//! fallback, bracketed on `[50, 5000] K`: the first pass - and any pass whose predecessor
//! moved less than `1e-12` - takes a heat-capacity Newton step, and every other pass
//! interpolates between the last two `(T, error)` pairs. The `1/T` variable appears nowhere in
//! the file. The docstring is a comment about an intention; this port follows the code, and
//! the capture records what the code does: three outer passes for the water-gas shift at 600 K.
//!
//! # The enthalpy is thermochemical, and the formation inventory is what makes it so
//!
//! NeqSim's process-stream enthalpy is a *sensible* one and excludes the ideal-gas formation
//! enthalpies. A reactive calculation cannot: the composition moves, and the heat that moving
//! it releases is part of the balance. So both sides of the comparison carry the inventory
//! `sum_i n_i dHf_i` - the specification adds it once at construction, and every trial state
//! adds its own - and the residual is the difference of the two thermochemical enthalpies over
//! the magnitude of the specification. The inventories do not cancel, because the composition
//! they are built from is not the same at the trial temperature as it was at the specification.
//!
//! # The inner flash is the caller's
//!
//! A pass needs three things from the state at a trial temperature: the flash's own pass count,
//! its **thermochemical** enthalpy and its heat capacity. All three come back through
//! [`InnerFlash`], so this module carries neither an equation of state nor the reactive flash
//! itself - [`crate::reactive_tp_flash`] and `azoth-eos`' enthalpy are the caller's business.

use azoth_core::Result;

/// The outer loop's pass cap, from `MAX_OUTER_ITER`.
pub const MAX_OUTER_ITERATIONS: usize = 200;

/// The tolerance on the normalised enthalpy residual, from `TOL`.
pub const TOL: f64 = 1.0e-8;

/// The largest temperature step a pass may take, from `MAX_T_STEP`.
pub const MAX_TEMPERATURE_STEP: f64 = 50.0;

/// The lowest temperature the loop will try, from `T_MIN`.
pub const T_MIN: f64 = 50.0;

/// The highest, from `T_MAX`.
pub const T_MAX: f64 = 5000.0;

/// The bracket width below which the loop declares convergence, from `(tHigh - tLow) < 1e-6`.
pub const BRACKET_TOLERANCE: f64 = 1.0e-6;

/// What one pass of the inner flash answers with.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhState {
    /// The flash's own passes at this temperature, which the loop reports summed.
    pub iterations: u32,
    /// The state's **thermochemical** enthalpy: its sensible enthalpy plus the formation
    /// inventory at this composition.
    pub thermochemical_enthalpy: f64,
    /// The state's heat capacity at constant pressure, which the first pass steps by.
    pub cp: f64,
}

/// The inner flash: a temperature in, the solved state's own numbers out.
pub type InnerFlash<'a> = &'a mut dyn FnMut(f64) -> Result<PhState>;

/// What the PH flash answers with.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhFlashOutcome {
    /// The temperature the loop stopped at - `getEquilibriumTemperature`.
    pub temperature: f64,
    /// `isConverged`. **Also true where the *bracket* closed rather than the residual**: the
    /// loop accepts a bracket narrower than [`BRACKET_TOLERANCE`] as an answer.
    pub converged: bool,
    /// `getOuterIterations`.
    pub outer_iterations: u32,
    /// `getTotalInnerIterations`: every inner flash's passes, summed over the outer loop.
    pub total_inner_iterations: u32,
}

/// `ReactiveMultiphasePHflash.solveEnthalpySpec`: the secant-with-bisection loop on `T`.
///
/// `initial_temperature` is where the loop starts - the system's temperature as the class
/// finds it, which its own test perturbs to make the search a search - and
/// `thermochemical_enthalpy_spec` is the specification **with the formation inventory already
/// added**, which is what the class's constructor builds.
///
/// # Errors
/// Whatever `inner` raises.
pub fn reactive_ph_flash(
    initial_temperature: f64,
    thermochemical_enthalpy_spec: f64,
    inner: InnerFlash<'_>,
) -> Result<PhFlashOutcome> {
    let mut absolute_spec = thermochemical_enthalpy_spec.abs();
    if absolute_spec < 1.0e-10 {
        // The class's own guard: a near-zero enthalpy would divide by nothing.
        absolute_spec = 1.0;
    }

    // Step 1: the inner flash at the starting temperature.
    let first = inner(initial_temperature)?;
    let mut total_inner_iterations = first.iterations;
    let mut temperature = initial_temperature;
    let mut error = (first.thermochemical_enthalpy - thermochemical_enthalpy_spec) / absolute_spec;
    let mut state = first;

    if error.abs() < TOL {
        return Ok(PhFlashOutcome {
            temperature,
            converged: true,
            outer_iterations: 0,
            total_inner_iterations,
        });
    }

    // The bracket, tracked on the sign of the residual: `H(T)` rises with `T`, so a positive
    // residual means the temperature is too high.
    let mut low = T_MIN;
    let mut high = T_MAX;
    if error > 0.0 {
        high = temperature;
    } else {
        low = temperature;
    }

    let mut previous_temperature = temperature;
    let mut previous_error = f64::NAN;
    let mut converged = false;
    let mut outer_iterations = 0_u32;

    for iteration in 0..MAX_OUTER_ITERATIONS {
        outer_iterations = iteration as u32 + 1;

        let has_bracket = low > T_MIN && high < T_MAX;
        let newton = |state: &PhState, temperature: f64| {
            let cp = if state.cp.abs() < 1.0e-20 {
                100.0
            } else {
                state.cp
            };
            temperature - (state.thermochemical_enthalpy - thermochemical_enthalpy_spec) / cp
        };

        let mut next = if iteration == 0
            || previous_error.is_nan()
            || (temperature - previous_temperature).abs() < 1.0e-12
        {
            // The first pass, and any pass that did not move: a heat-capacity Newton step.
            newton(&state, temperature)
        } else {
            // The secant step, which is what captures `dH/dT` including the reaction's own
            // contribution - the heat the shifting equilibrium absorbs or releases.
            let slope = (error - previous_error) / (temperature - previous_temperature);
            if slope.abs() > 1.0e-30 {
                temperature - error / slope
            } else {
                newton(&state, temperature)
            }
        };

        // A step larger than the cap is trimmed to it, and a step that leaves the bracket is
        // replaced by the bracket's midpoint - the guaranteed-convergence fallback.
        if (next - temperature).abs() > MAX_TEMPERATURE_STEP {
            next = temperature + (next - temperature).signum() * MAX_TEMPERATURE_STEP;
        }
        if has_bracket && (next <= low || next >= high) {
            next = 0.5 * (low + high);
        }
        next = next.clamp(T_MIN, T_MAX);

        // A step that went nowhere: bisect if there is a bracket, else move a kelvin the way
        // the residual points.
        if (next - temperature).abs() < 1.0e-12 {
            next = if has_bracket {
                0.5 * (low + high)
            } else {
                temperature + if error < 0.0 { 1.0 } else { -1.0 }
            };
        }

        previous_temperature = temperature;
        previous_error = error;

        state = inner(next)?;
        total_inner_iterations += state.iterations;
        temperature = next;
        error = (state.thermochemical_enthalpy - thermochemical_enthalpy_spec) / absolute_spec;

        if error > 0.0 && temperature < high {
            high = temperature;
        } else if error < 0.0 && temperature > low {
            low = temperature;
        }

        if error.abs() < TOL {
            converged = true;
            break;
        }
        let has_bracket = low > T_MIN && high < T_MAX;
        if has_bracket && (high - low) < BRACKET_TOLERANCE {
            converged = true;
            break;
        }
    }

    Ok(PhFlashOutcome {
        temperature,
        converged,
        outer_iterations,
        total_inner_iterations,
    })
}
