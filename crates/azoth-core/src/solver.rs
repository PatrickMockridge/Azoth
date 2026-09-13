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
//! guess and convergence rule - and both languages implement exactly it.
//!
//! This is not a general-purpose numerical library and must not become one. The
//! rule that governs what may be added here: **it holds the named schemes that
//! `specs/schema/calc.schema.json`'s `solver.kind` permits, and it grows only
//! when that enum grows** - each addition bringing the scheme in both languages
//! and the contract test that holds the enum to them, in one commit.
//!
//! It lives in `azoth-core` rather than in a namespace crate because every
//! namespace may need it and no namespace may depend on another. It was in
//! `azoth-hydraulics` while hydraulics was the only domain with an implicit
//! equation; the first calc outside it made that untenable.

use crate::{AzothError, Result};

/// The named solution schemes the spec schema permits.
///
/// A cross-language contract in the same way [`Convergence`] and
/// `WarningCode` are: the schema's `solver.kind` enum, the Python [`SolverKind`]
/// in `azoth/core/solver.py`, and [`SolverKind::ALL`] here must name the same
/// set, and `python/tests/test_solver_contract.py` asserts that rather than
/// trusting three hand-edited lists to stay in step. [`SolverKind::ALL`] is
/// exposed to Python through `azoth._core.solver_kinds`.
///
/// [`SolverKind`]: https://docs.rs/azoth/latest/azoth/core/solver
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolverKind {
    /// `x -> f(x)`, iterated from the spec's declared initial guess.
    FixedPoint,
    /// The real roots of a cubic, formed analytically and then polished.
    CubicRoots,
}

impl SolverKind {
    /// Every kind this crate implements, in the schema's spelling.
    ///
    /// Adding a name here is a claim that *both* implementations run it. The
    /// schema's own description of `solver.kind` says why that matters: a spec
    /// naming a kind neither can run "would describe a calculation neither
    /// implementation can run, which is the 'looks like validation, does
    /// nothing' failure this schema exists to catch".
    pub const ALL: &'static [SolverKind] = &[SolverKind::FixedPoint, SolverKind::CubicRoots];

    /// The schema's spelling of this kind.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::FixedPoint => "fixed_point",
            Self::CubicRoots => "cubic_roots",
        }
    }

    /// Parse the spec's spelling.
    ///
    /// # Errors
    /// Unknown names are an error rather than a silent default: guessing which
    /// scheme was meant would change the answer, and silently falling back to
    /// the only implemented kind would run a scheme the spec did not ask for.
    pub fn parse(name: &str) -> Result<Self> {
        Self::ALL
            .iter()
            .copied()
            .find(|kind| kind.as_str() == name)
            .ok_or_else(|| {
                let known: Vec<&str> = Self::ALL.iter().map(|k| k.as_str()).collect();
                AzothError::invalid_input(
                    "solver.kind",
                    format!("unknown solver kind `{name}`; expected one of {known:?}"),
                )
            })
    }
}

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
            other => Err(AzothError::invalid_input(
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

/// What [`cubic_roots`] produced.
///
/// Not a [`SolverOutcome`], because a cubic has up to three answers rather than
/// one and the caller - not this module - decides which of them is the physical
/// one. `roots` is returned whole so that decision stays with the calc that knows
/// what its roots mean.
#[derive(Debug, Clone, PartialEq)]
pub struct CubicRootsOutcome {
    /// The real roots, ascending, after polishing.
    ///
    /// Length 1 when the discriminant is positive and 3 when it is not, which
    /// includes the degenerate cases: a double or triple root is returned as
    /// repeated values rather than collapsed, so a caller that counts roots gets
    /// the algebraic count and not a de-duplicated one.
    pub roots: Vec<f64>,
    /// Newton steps taken, summed over the roots.
    pub iterations: u32,
    /// Whether every root met the stopping rule.
    pub converged: bool,
    /// The largest `|x_k - x_{k-1}|` seen on the final step of any root.
    pub residual: f64,
}

/// The real roots of `x**3 + c2*x**2 + c1*x + c0`, formed analytically then polished.
///
/// # Why analytic first
///
/// An iterative root-finder started from a fixed guess does not reliably find the
/// root a caller wants. Measured on the Peng-Robinson cubic for propane at
/// `Tr = 0.8, Pr = 0.25` - roots 0.0368, 0.1481, 0.7908 - Newton started at 0.30
/// converges to the *middle* root, at 0.50 to the liquid root and at 0.70 to the
/// vapour one. Which root comes back depends on where the search began, so an
/// implementation that guesses cannot be told what it found. Forming the roots
/// analytically and ordering them removes the guess.
///
/// # Why the polish is not optional
///
/// Three reasons, and all of them are measured rather than anticipated:
///
/// * **Cardano cancels catastrophically when the roots are close.** The
///   subtraction `cbrt(-q/2 + sqrt(d)) - cbrt(-q/2 - sqrt(d))` loses precision as
///   the two terms approach each other, which is exactly the near-triple-root case.
/// * **The discriminant's sign is not reliably computable near zero.** With the
///   rounded Peng-Robinson constants at the critical point it is `1.8e-12` -
///   positive, but small enough that a different rounding puts it on the other
///   side, selecting the three-real-root branch for a cubic whose roots are mostly
///   complex.
/// * **`cbrt`, `acos` and `cos` are not correctly rounded** in any libm, so the
///   analytic value is only close to begin with. The polish is what makes the two
///   languages agree to a tolerance worth asserting.
///
/// A Newton step on the polynomial is cheap, is exactly representable arithmetic,
/// and converges quadratically once it is near - and the analytic formation has
/// already put it near.
///
/// # Degenerate cases
///
/// A double or triple root is returned as repeated values, not collapsed. When the
/// depressed cubic's `p` and `q` are both zero the trigonometric branch would
/// divide by zero, so that case is answered directly as a triple root at `-c2/3`.
#[must_use]
pub fn cubic_roots(
    c2: f64,
    c1: f64,
    c0: f64,
    tolerance: f64,
    max_iterations: u32,
    convergence: Convergence,
) -> CubicRootsOutcome {
    // Depress x**3 + c2*x**2 + c1*x + c0 to w**3 + p*w + q by x = w - c2/3.
    let p = c1 - c2 * c2 / 3.0;
    let q = 2.0 * c2 * c2 * c2 / 27.0 - c2 * c1 / 3.0 + c0;
    let half_q = q / 2.0;
    let third_p = p / 3.0;
    let discriminant = half_q * half_q + third_p * third_p * third_p;

    let mut roots: Vec<f64> = if discriminant > 0.0 {
        // One real root and two complex conjugates. Cardano, with `cbrt` carrying
        // the sign so a negative argument does not become NaN.
        let root_d = discriminant.sqrt();
        let u = (-half_q + root_d).cbrt();
        let v = (-half_q - root_d).cbrt();
        vec![u + v - c2 / 3.0]
    } else {
        let radius = (-third_p * third_p * third_p).sqrt();
        if radius == 0.0 {
            // p = q = 0: the cubic is (w)**3, a triple root. The general form would
            // divide by `radius` here.
            vec![-c2 / 3.0; 3]
        } else {
            // Clamped because rounding can push the ratio a hair outside [-1, 1],
            // and `acos` of that is NaN rather than a slightly wrong angle.
            let cosine = (-half_q / radius).clamp(-1.0, 1.0);
            let angle = cosine.acos();
            let scale = 2.0 * (-third_p).sqrt();
            let mut three: Vec<f64> = (0..3)
                .map(|k| {
                    let two_pi_k = 2.0 * std::f64::consts::PI * f64::from(k);
                    scale * ((angle + two_pi_k) / 3.0).cos() - c2 / 3.0
                })
                .collect();
            three.sort_by(|a, b| a.partial_cmp(b).expect("no NaN before polishing"));
            three
        }
    };

    // Newton polish, on the polynomial rather than the depressed form: that is the
    // equation the caller wrote, and a step on it needs no back-substitution.
    let mut iterations = 0_u32;
    let mut residual = 0.0_f64;
    let mut converged = true;
    for root in &mut roots {
        let mut x = *root;
        let mut step_residual = f64::INFINITY;
        let mut root_converged = false;
        for _ in 0..max_iterations {
            let f = ((x + c2) * x + c1) * x + c0;
            let df = (3.0 * x + 2.0 * c2) * x + c1;
            if df == 0.0 {
                // A stationary point: Newton cannot step, and the analytic value is
                // the best available. Leave it and let the stopping rule decide.
                break;
            }
            let next = x - f / df;
            let delta = (next - x).abs();
            step_residual = delta;
            x = next;
            iterations += 1;
            let close_enough = match convergence {
                Convergence::Absolute => delta <= tolerance,
                Convergence::Relative => delta <= tolerance * x.abs(),
            };
            if close_enough {
                root_converged = true;
                break;
            }
        }
        if !root_converged {
            converged = false;
        }
        residual = residual.max(step_residual);
        *root = x;
    }

    CubicRootsOutcome {
        roots,
        iterations,
        converged,
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
/// Returns [`AzothError::SolverNotConverged`] when the outcome did not
/// converge.
pub fn require_converged(outcome: SolverOutcome, tolerance: f64) -> Result<SolverOutcome> {
    if outcome.converged {
        return Ok(outcome);
    }
    Err(AzothError::SolverNotConverged {
        iterations: outcome.iterations,
        residual: outcome.residual,
        tolerance,
    })
}

/// Turn a cubic outcome that did not converge into an error.
///
/// The same rule [`require_converged`] applies and for the same reason: a root
/// that did not meet tolerance is the last iterate of an iteration that did not
/// finish, and is not a root of anything. Separate from [`require_converged`]
/// rather than generic over it because the two outcomes are different shapes -
/// one answer against up to three - and a trait to unify them would buy one call
/// site and cost a reader a definition to chase.
///
/// # Errors
/// Returns [`AzothError::SolverNotConverged`] when the polish did not converge.
pub fn require_cubic_converged(
    outcome: CubicRootsOutcome,
    tolerance: f64,
) -> Result<CubicRootsOutcome> {
    if outcome.converged {
        return Ok(outcome);
    }
    Err(AzothError::SolverNotConverged {
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

    #[test]
    fn every_solver_kind_round_trips_through_its_name() {
        for kind in SolverKind::ALL {
            assert_eq!(SolverKind::parse(kind.as_str()).unwrap(), *kind);
        }
    }

    #[test]
    fn solver_kind_names_are_unique_and_lowercase() {
        let names: Vec<&str> = SolverKind::ALL.iter().map(|k| k.as_str()).collect();
        let unique: std::collections::BTreeSet<&str> = names.iter().copied().collect();
        assert_eq!(unique.len(), names.len(), "ALL contains a duplicate");
        for name in names {
            assert_eq!(name, name.to_lowercase(), "`{name}` is not lowercase");
        }
    }

    #[test]
    fn cubic_roots_finds_three_well_separated_roots_in_order() {
        // (x - 1)(x - 2)(x - 3) = x**3 - 6x**2 + 11x - 6
        let out = cubic_roots(-6.0, 11.0, -6.0, 1e-14, 50, Convergence::Relative);
        assert!(out.converged, "residual {:e}", out.residual);
        assert_eq!(out.roots.len(), 3);
        for (got, want) in out.roots.iter().zip([1.0, 2.0, 3.0]) {
            assert!((got - want).abs() < 1e-12, "got {got}, want {want}");
        }
    }

    #[test]
    fn cubic_roots_returns_one_root_when_the_other_two_are_complex() {
        // x**3 + x + 1 has one real root, near -0.682327803828.
        let out = cubic_roots(0.0, 1.0, 1.0, 1e-14, 50, Convergence::Relative);
        assert_eq!(out.roots.len(), 1);
        assert!(out.converged);
        let root = out.roots[0];
        let residual = ((root) * root + 1.0) * root + 1.0;
        assert!(residual.abs() < 1e-12, "not a root: f = {residual:e}");
    }

    #[test]
    fn cubic_roots_handles_a_triple_root_without_dividing_by_zero() {
        // (x - 2)**3 = x**3 - 6x**2 + 12x - 8. This is the degenerate case the
        // trigonometric branch cannot take, because its radius is zero.
        let out = cubic_roots(-6.0, 12.0, -8.0, 1e-14, 50, Convergence::Relative);
        assert_eq!(out.roots.len(), 3, "a triple root is three roots, not one");
        for got in &out.roots {
            assert!((got - 2.0).abs() < 1e-9, "got {got}, want 2");
        }
    }

    #[test]
    fn cubic_roots_is_deterministic() {
        // The same property `fixed_point` is held to, and for the same reason: the
        // two languages are compared against each other, so a scheme that lands on
        // different roots from identical input would make that comparison
        // meaningless rather than merely noisy.
        let run = || cubic_roots(-6.0, 11.0, -6.0, 1e-14, 50, Convergence::Relative);
        let first = run();
        let second = run();
        assert_eq!(first.iterations, second.iterations);
        for (a, b) in first.roots.iter().zip(&second.roots) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }

    #[test]
    fn cubic_roots_answers_a_polynomial_it_was_not_shaped_for() {
        // x**3 - 2x**2 - 5x + 6 = (x - 1)(x + 2)(x - 3), roots -2, 1, 3. Included
        // because its roots straddle zero, so an implementation that assumed
        // positive roots - which a Z factor always is - fails here.
        let out = cubic_roots(-2.0, -5.0, 6.0, 1e-14, 50, Convergence::Relative);
        for (got, want) in out.roots.iter().zip([-2.0, 1.0, 3.0]) {
            assert!((got - want).abs() < 1e-12, "got {got}, want {want}");
        }
    }

    #[test]
    fn an_unknown_solver_kind_is_an_error_naming_what_is_known() {
        let err = SolverKind::parse("bisection").unwrap_err();
        let message = err.to_string();
        assert!(message.contains("bisection"), "{message}");
        assert!(
            message.contains("fixed_point"),
            "the error should say what IS implemented: {message}"
        );
    }
}
