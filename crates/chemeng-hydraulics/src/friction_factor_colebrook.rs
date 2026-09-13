//! `hydraulics.friction_factor_colebrook` - the implicit Colebrook-White
//! friction factor.
//!
//! ```text
//! 1 / sqrt(f) = -2 * log10(relative_roughness / 3.7 + 2.51 / (Re * sqrt(f)))
//! ```
//!
//! Colebrook, C. F. (1939). "Turbulent flow in pipes, with particular reference
//! to the transition region between the smooth and rough pipe laws." J. Inst.
//! Civil Engineers 11(4), 133-156. DOI 10.1680/ijoti.1939.13150
//!
//! Spec: `specs/calcs/hydraulics/friction_factor_colebrook.yaml`
//!
//! The equation is implicit in `f`, so it is solved rather than evaluated.
//! Substituting `x = 1/sqrt(f)` makes it a fixed point in `x`, which converges
//! quickly and - crucially - converges the *same way* in both languages, since
//! the scheme is fixed by the spec rather than chosen per implementation.

use crate::results::ColebrookResult;
use crate::solver::{Convergence, fixed_point, require_converged};
use crate::spec_gen;
use chemeng_core::{Result, apply_checks};

/// Solve the Colebrook-White equation for the Darcy friction factor.
///
/// `relative_roughness` is `epsilon / D`.
///
/// The returned [`ColebrookResult`] carries the solver's own report. That is not
/// decoration: when `converged` is false, `f` is the last iterate of an
/// iteration that did not meet tolerance, and it does not solve the equation.
/// This function returns it only when it converged; a non-converged run is an
/// error, because a number that solves nothing should not be returned as if it
/// were an answer.
///
/// # Errors
/// * [`chemeng_core::ChemEngError::OutOfRange`] if `Re <= 0` (the `2.51/(Re*sqrt(f))`
///   term is singular) or if `relative_roughness < 0` (unphysical).
/// * [`chemeng_core::ChemEngError::SolverNotConverged`] if the iteration hits its
///   cap. With the spec's tolerance and cap this does not happen for physical
///   inputs, and a test asserts so.
///
/// Below `Re = 4000` the result is out of the range the correlation was fitted
/// to; the value is still returned, carrying an `OutOfValidRange` warning. Use
/// `64/Re` for laminar flow - a different equation, deliberately not part of
/// this calc.
///
/// # Example
/// ```
/// use chemeng_hydraulics::friction_factor_colebrook;
///
/// // Commercial steel: epsilon = 0.046 mm in a 100 mm pipe.
/// let r = friction_factor_colebrook(100_000.0, 4.6e-4)?;
/// assert!((r.f - 0.02016203).abs() < 1e-6);
/// assert!(r.converged);
/// # Ok::<(), chemeng_core::ChemEngError>(())
/// ```
pub fn friction_factor_colebrook(re: f64, relative_roughness: f64) -> Result<ColebrookResult> {
    let spec = &spec_gen::FRICTION_FACTOR_COLEBROOK_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "re" => Some(re),
            "relative_roughness" => Some(relative_roughness),
            _ => None,
        },
        &mut warnings,
    )?;

    let solver = spec
        .solver
        .expect("friction_factor_colebrook is implicit and must declare a solver in its spec");
    let convergence = Convergence::parse(solver.convergence)?;

    // x = 1/sqrt(f). The initial guess is a friction factor, so its reciprocal
    // square root is the starting x.
    let x0 = 1.0 / solver.initial_guess.sqrt();
    let outcome = fixed_point(
        x0,
        solver.tolerance,
        solver.max_iterations,
        convergence,
        |x| -2.0 * (relative_roughness / 3.7 + 2.51 * x / re).log10(),
    );
    let outcome = require_converged(outcome, solver.tolerance)?;

    let f = 1.0 / (outcome.x * outcome.x);

    // Range checks on the output. Nothing here depends on an optional input for
    // this calc, so no check can be skipped.
    apply_checks(
        spec.derived_checks(),
        |q| (q == "f").then_some(f),
        &mut warnings,
    )?;

    Ok(ColebrookResult {
        f,
        iterations: outcome.iterations,
        converged: outcome.converged,
        residual: outcome.residual,
        warnings,
    })
}

/// The fully-rough asymptote of the Colebrook equation, which is explicit.
///
/// As `Re -> infinity` the roughness term dominates and
/// `1/sqrt(f) -> -2*log10(relative_roughness/3.7)`. Exposed because it is a
/// useful sanity bound, and because a test uses it to confirm the solver
/// converges to the correct asymptote rather than to a nearby fixed point.
#[must_use]
pub fn fully_rough_limit(relative_roughness: f64) -> f64 {
    let x = -2.0 * (relative_roughness / 3.7).log10();
    1.0 / (x * x)
}
