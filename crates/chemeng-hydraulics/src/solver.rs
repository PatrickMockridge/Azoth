//! The iterative solver used by implicit equations.
//!
//! This module is small and deliberately not general. It exists because
//! `hydraulics.friction_factor_colebrook` is implicit in `f`, and the Python and
//! Rust implementations have to agree numerically. Two implementations running
//! the *same* algorithm converge to the same floating-point value; two
//! implementations running different-but-equally-valid algorithms (Newton vs
//! fixed point, say) agree only to within their difference, which is far larger
//! than a tight tolerance would allow.
//!
//! So the scheme is fixed by the spec - kind, tolerance, iteration cap, initial
//! guess and convergence rule - and both languages implement exactly it. Nothing
//! here is a general-purpose numerical library and it should not grow into one.

use chemeng_core::{ChemEngError, Result};

/// How successive iterates are compared to decide convergence.
///
/// Named in the spec, because the same tolerance means different things under
/// the two rules and the two implementations must agree on which is meant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Convergence {
    /// `|x_{k+1} - x_k| <= tolerance`.
    Absolute,
    /// `|x_{k+1} - x_k| <= tolerance * |x_{k+1}|`.
    Relative,
}

impl Convergence {
    /// Parse the spec's spelling.
    ///
    /// # Errors
    /// Unknown names are an error rather than a silent default: guessing which
    /// convergence rule was meant would change the answer.
    pub fn parse(name: &str) -> Result<Self> {
        match name {
            "absolute" => Ok(Self::Absolute),
            "relative" => Ok(Self::Relative),
            other => Err(ChemEngError::invalid_input(
                "convergence",
                format!("unknown convergence rule `{other}`; expected `absolute` or `relative`"),
            )),
        }
    }
}

/// What the iteration produced.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolverOutcome {
    /// The final iterate. Note this is the last value computed, which when
    /// `converged` is false is not a solution to anything.
    pub x: f64,
    /// Iterations performed.
    pub iterations: u32,
    /// Whether the stopping rule was met.
    pub converged: bool,
    /// `|x_k - x_{k-1}|` at the final step.
    pub residual: f64,
}

/// Iterate `x -> f(x)` from `x0` until the stopping rule is met or the cap is
/// reached.
///
/// Iteration counts are deterministic: the same inputs always perform the same
/// number of steps, which is what makes the Python and Rust results comparable
/// rather than merely close.
///
/// A `NaN` iterate is not converged - and cannot be, since every comparison with
/// `NaN` is false - so it falls through to the iteration cap rather than
/// masquerading as a solution.
#[must_use]
pub fn fixed_point(
    x0: f64,
    tolerance: f64,
    max_iterations: u32,
    convergence: Convergence,
    f: impl Fn(f64) -> f64,
) -> SolverOutcome {
    let mut x = x0;
    let mut residual = f64::INFINITY;

    for iteration in 1..=max_iterations {
        let next = f(x);
        residual = (next - x).abs();
        let close_enough = match convergence {
            Convergence::Absolute => residual <= tolerance,
            Convergence::Relative => residual <= tolerance * next.abs(),
        };
        x = next;
        if close_enough {
            return SolverOutcome {
                x,
                iterations: iteration,
                converged: true,
                residual,
            };
        }
    }

    SolverOutcome {
        x,
        iterations: max_iterations,
        converged: false,
        residual,
    }
}

/// Turn a non-converged outcome into an error.
///
/// A caller that gets `f` from a run that never converged has a number that
/// solves nothing, so this is an error rather than a warning - unlike an
/// out-of-range input, where the number is real and merely untrustworthy.
///
/// # Errors
/// Returns [`ChemEngError::SolverNotConverged`] when the outcome did not
/// converge.
pub fn require_converged(outcome: SolverOutcome, tolerance: f64) -> Result<SolverOutcome> {
    if outcome.converged {
        return Ok(outcome);
    }
    Err(ChemEngError::SolverNotConverged {
        iterations: outcome.iterations,
        residual: outcome.residual,
        tolerance,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn converges_on_a_contraction() {
        // x = (x + 2/x) / 2 has fixed point sqrt(2).
        let out = fixed_point(1.0, 1e-15, 100, Convergence::Relative, |x| {
            (x + 2.0 / x) / 2.0
        });
        assert!(out.converged);
        assert!(
            (out.x - std::f64::consts::SQRT_2).abs() < 1e-14,
            "{}",
            out.x
        );
    }

    #[test]
    fn iteration_count_is_deterministic() {
        let run = || {
            fixed_point(1.0, 1e-12, 100, Convergence::Relative, |x| {
                (x + 2.0 / x) / 2.0
            })
        };
        assert_eq!(run().iterations, run().iterations);
        assert_eq!(run().x.to_bits(), run().x.to_bits());
    }

    #[test]
    fn a_divergent_iteration_hits_the_cap_and_says_so() {
        // x -> 2x diverges; it must stop at the cap rather than run forever,
        // and must not claim to have converged.
        let out = fixed_point(1.0, 1e-12, 10, Convergence::Absolute, |x| 2.0 * x);
        assert!(!out.converged);
        assert_eq!(out.iterations, 10);
        assert!(require_converged(out, 1e-12).is_err());
    }

    #[test]
    fn nan_never_converges() {
        // NaN would otherwise slip through every comparison as "false", i.e. as
        // not converged, which is correct - but assert it explicitly so a future
        // refactor to a different stopping rule cannot silently accept NaN.
        let out = fixed_point(1.0, 1e-12, 5, Convergence::Relative, |_| f64::NAN);
        assert!(!out.converged);
        assert!(out.x.is_nan());
        assert!(require_converged(out, 1e-12).is_err());
    }

    #[test]
    fn relative_and_absolute_differ_when_the_magnitude_is_large() {
        // |dx| = 1e-6. Under an absolute tolerance of 1e-9 that is not converged;
        // under a relative tolerance of 1e-9 applied to x = 1e4 (threshold 1e-5)
        // it is. This is why the spec has to name the rule.
        let f = |x: f64| x + 1e-6;
        let abs = fixed_point(1e4, 1e-9, 5, Convergence::Absolute, f);
        let rel = fixed_point(1e4, 1e-9, 5, Convergence::Relative, f);
        assert!(!abs.converged);
        assert!(rel.converged);
    }

    #[test]
    fn convergence_names_are_parsed_strictly() {
        assert_eq!(
            Convergence::parse("relative").unwrap(),
            Convergence::Relative
        );
        assert_eq!(
            Convergence::parse("absolute").unwrap(),
            Convergence::Absolute
        );
        assert!(Convergence::parse("Relative").is_err());
        assert!(Convergence::parse("rel").is_err());
    }
}
