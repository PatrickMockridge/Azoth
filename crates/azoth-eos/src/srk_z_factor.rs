//! `eos.srk_z_factor` - the Soave-Redlich-Kwong compressibility factor.
//!
//! ```text
//! z**3 - z**2 + (a_reduced - b_reduced - b_reduced**2)*z - a_reduced*b_reduced = 0
//! ```
//!
//! Spec: `specs/calcs/eos/srk_z_factor.toml`, which carries the provenance and the same
//! admissible-root theorem as the Peng-Robinson form.
//!
//! This calc returns the outermost two roots that are admissible - `z > B` - and
//! discards the middle one, exactly as [`crate::pr_z_factor`] does; only the cubic's
//! coefficients differ.

use azoth_core::solver::{Convergence, cubic_roots, require_cubic_converged};
use azoth_core::{Result, apply_checks};

use crate::cubic::Cubic;
use crate::results::{RootStructure, SrkZFactorResult};
use crate::spec_gen;

/// The Soave-Redlich-Kwong compressibility factor, for one state.
///
/// `a_reduced` and `b_reduced` are the cubic's dimensionless parameters, from
/// [`crate::srk_alpha_ab`], already carrying any composition.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `b_reduced <= 0` or `a_reduced < 0`.
/// * [`azoth_core::AzothError::SolverNotConverged`] if the polish hits its cap.
///
/// # Example
/// ```
/// use azoth_eos::{srk_alpha_ab, srk_kappa, srk_z_factor};
///
/// let kappa = srk_kappa(0.152)?.kappa;
/// let ab = srk_alpha_ab(kappa, 0.8, 0.25)?;
/// let z = srk_z_factor(ab.a_reduced, ab.b_reduced)?;
/// assert!((z.z_min - 0.04171316701447374).abs() < 1e-12);
/// assert!((z.z_max - 0.8019551972557889).abs() < 1e-12);
/// assert_eq!(z.root_structure, azoth_eos::RootStructure::Three);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn srk_z_factor(a_reduced: f64, b_reduced: f64) -> Result<SrkZFactorResult> {
    let spec = &spec_gen::SRK_Z_FACTOR_SPEC;
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
        .expect("srk_z_factor is implicit and must declare a solver in its spec");
    let convergence = Convergence::parse(solver.convergence)?;

    // The monic cubic z**3 + c2*z**2 + c1*z + c0, coefficients from the cubic's
    // geometry so a reader can check them against the equation above.
    let (c2, c1, c0) = Cubic::Srk.z_coefficients(a_reduced, b_reduced);

    let outcome = cubic_roots(
        c2,
        c1,
        c0,
        solver.tolerance,
        solver.max_iterations,
        convergence,
    );
    let outcome = require_cubic_converged(outcome, solver.tolerance)?;

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

    let (z_min, z_max) = admissible
        .iter()
        .fold((f64::MAX, f64::MIN), |(lo, hi), z| (lo.min(*z), hi.max(*z)));

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

    Ok(SrkZFactorResult {
        z_min,
        z_max,
        root_structure,
        iterations: outcome.iterations,
        converged: outcome.converged,
        residual: outcome.residual,
        warnings,
    })
}
