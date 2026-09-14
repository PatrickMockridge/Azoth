//! `eos.pr_z_factor` - the Peng-Robinson compressibility factor.
//!
//! ```text
//! z**3 - (1 - B)*z**2 + (A - 3*B**2 - 2*B)*z - (A*B - B**2 - B**3) = 0
//! ```
//!
//! Spec: `specs/calcs/eos/pr_z_factor.toml`, which carries the provenance, the proof
//! that the admissible root count is one or three, and what the answer does near the
//! critical point.
//!
//! This calc returns the outermost two roots that are admissible - `z > B` - and
//! discards the middle one, which lies on the unstable branch between the spinodals.

use azoth_core::solver::{Convergence, cubic_roots, require_cubic_converged};
use azoth_core::{Result, apply_checks};

use crate::results::{PrZFactorResult, RootStructure};
use crate::spec_gen;

/// The Peng-Robinson compressibility factor, for one state.
///
/// `a_reduced` and `b_reduced` are the cubic's dimensionless parameters, from
/// [`crate::pr_alpha_ab`]. They already carry the composition, so this calc cannot
/// tell whether it was given a pure component or a mixture's pseudo-parameters, and
/// does not apply a mixing rule.
///
/// Returns the smallest and largest *admissible* roots, and how many there were.
/// The middle root is deliberately not returned - see the module documentation.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `b_reduced <= 0` (the cubic
///   degenerates to a trivial double root at zero pressure) or `a_reduced < 0`
///   (which no Peng-Robinson state produces).
/// * [`azoth_core::AzothError::SolverNotConverged`] if the polish hits its cap.
///
/// # Example
/// ```
/// use azoth_eos::{pr_alpha_ab, pr_kappa, pr_z_factor};
///
/// let kappa = pr_kappa(0.152)?.kappa;
/// let ab = pr_alpha_ab(kappa, 0.8, 0.25)?;
/// let z = pr_z_factor(ab.a_reduced, ab.b_reduced)?;
/// assert!((z.z_min - 0.036771130006003336).abs() < 1e-12);
/// assert!((z.z_max - 0.7907792374136378).abs() < 1e-12);
/// assert_eq!(z.root_structure, azoth_eos::RootStructure::Three);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn pr_z_factor(a_reduced: f64, b_reduced: f64) -> Result<PrZFactorResult> {
    let spec = &spec_gen::PR_Z_FACTOR_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "a_reduced" => Some(a_reduced),
            "b_reduced" => Some(b_reduced),
            _ => None,
        },
        &mut warnings,
    )?;

    let solver = spec
        .solver
        .expect("pr_z_factor is implicit and must declare a solver in its spec");
    let convergence = Convergence::parse(solver.convergence)?;

    // The monic cubic z**3 + c2*z**2 + c1*z + c0, coefficients straight from the
    // published form so a reader can check them against the equation above.
    let c2 = -(1.0 - b_reduced);
    let c1 = a_reduced - 3.0 * b_reduced * b_reduced - 2.0 * b_reduced;
    let c0 = -(a_reduced * b_reduced - b_reduced * b_reduced - b_reduced * b_reduced * b_reduced);

    let outcome = cubic_roots(
        c2,
        c1,
        c0,
        solver.tolerance,
        solver.max_iterations,
        convergence,
    );
    let outcome = require_cubic_converged(outcome, solver.tolerance)?;

    // Admissible means `z > B`: `z = B` is the zero-volume limit, and a root below it
    // makes `ln(z - B)` the logarithm of a negative number.
    //
    // For `B > 0` there is always at least one: the polynomial equals `-2*B**2` at
    // `z = B` and tends to positive infinity, so it crosses zero above `B`. An empty
    // set would mean `B <= 0`, which the input checks have already refused.
    let admissible: Vec<f64> = outcome
        .roots
        .iter()
        .copied()
        .filter(|z| *z > b_reduced)
        .collect();
    debug_assert!(
        !admissible.is_empty(),
        "no admissible root for B = {b_reduced}, which the bounds should have refused"
    );

    // `f64::min`/`max` rather than `first()`/`last()`: the roots come back sorted,
    // but relying on that here would make the admissibility filter's correctness
    // depend on an ordering the solver happens to guarantee.
    let (z_min, z_max) = admissible
        .iter()
        .fold((f64::MAX, f64::MIN), |(lo, hi), z| (lo.min(*z), hi.max(*z)));

    // One or three, and nothing else - see `RootStructure`. Not a panic for the
    // impossible case: this is library code, and the project keeps panics for a
    // malformed *spec*, which is a programmer's error, rather than for a state
    // that arithmetic is merely expected never to reach. The assertion is in the
    // test build, where the property sweep that establishes the theorem runs.
    debug_assert!(
        admissible.len() == 1 || admissible.len() == 3,
        "the admissible root count is 1 or 3, never {}",
        admissible.len()
    );
    let root_structure = if admissible.len() == 1 {
        RootStructure::One
    } else {
        RootStructure::Three
    };

    apply_checks(
        spec.derived_checks(),
        |name| match name {
            "z_min" => Some(z_min),
            "z_max" => Some(z_max),
            _ => None,
        },
        &mut warnings,
    )?;

    Ok(PrZFactorResult {
        z_min,
        z_max,
        root_structure,
        iterations: outcome.iterations,
        converged: outcome.converged,
        residual: outcome.residual,
        warnings,
    })
}
