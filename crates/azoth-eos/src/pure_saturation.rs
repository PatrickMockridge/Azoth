//! `eos.pure_saturation` - pure-component saturation pressure.
//!
//! The pressure at which a pure component's vapour and liquid roots have equal
//! fugacity, found by bisection. A *model* rather than a calculation: what the spec
//! pins down is the procedure, not an equation, and this module reads the procedure
//! from the generated table rather than choosing it.
//!
//! Spec: `specs/models/eos/pure_saturation.yaml`
//!
//! # This composes the kernels and adds a search, and nothing else
//!
//! Every number it computes comes from a registered calculation:
//! [`crate::pr_kappa`] for the attraction coefficient, [`crate::pr_alpha_ab`] for the
//! reduced parameters at a trial pressure, [`crate::pr_z_factor`] for the roots and
//! their admissibility, and [`crate::pr_departure`] for the two fugacity coefficients.
//! The only thing here that is not a kernel is the loop that searches for the
//! pressure where the two agree.
//!
//! That is deliberate, and it is the model layer's whole contract: a second
//! implementation of the Peng-Robinson equation living in the model layer would be a
//! *third* implementation of it - untested, uncited, and cross-checked by nothing -
//! while the two that exist still claimed to be the independent pair.
//!
//! # The algorithm
//!
//! 1. **Bracket.** Scan the reduced pressure upward from the spec's `lower` to its
//!    `upper` in `steps` points, and take the *last* one at which the cubic still has
//!    three admissible roots. That is the spinodal; above it there is one root, no
//!    liquid branch, and nothing to equate.
//!
//!    The window is **one-sided**, which is worth stating because it is easy to
//!    assume otherwise: a cubic has three real roots at every pressure below the
//!    spinodal, including pressures so low that the "liquid" root describes a molar
//!    volume no liquid could have. The bracket's lower end is therefore the scan's
//!    first point and is arbitrary on purpose - the residual is positive and
//!    monotonically decreasing below the saturation pressure, so the bracket only has
//!    to straddle the root, not be tight.
//!
//! 2. **Bisect** on the reduced pressure until the bracket's *width* meets the
//!    tolerance, relatively. Not until the residual is small: the residual is the
//!    thing being solved for, so stopping on it would be circular whenever the two
//!    fugacities disagree for a reason other than the pressure.
//!
//! 3. **Return** the bracket's midpoint times `Pc`.

use azoth_core::units::{Pressure, ThermodynamicTemperature, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::{PureSaturationResult, RootStructure};
use crate::{model_gen, pr_alpha_ab, pr_departure, pr_kappa, pr_z_factor};

/// The saturation pressure of a pure component at a temperature.
///
/// `Tc`, `Pc` and `omega` are the caller's - this library ships no component
/// databank - and `T` must be below `Tc`.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `T`, `Tc` or `Pc` is not positive, or if
///   `T >= Tc` - a state with no saturation pressure. The bound is on the ratio
///   `T/Tc`, so the error names what is actually wrong rather than picking one of the
///   two inputs arbitrarily.
/// * [`AzothError::OutOfRange`] with a different detail if the scan finds no pressure
///   at which the cubic has a liquid branch, which a below-critical temperature
///   should never produce and which means the scan's bounds are wrong.
/// * [`AzothError::SolverNotConverged`] if the bisection hits its cap.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, pascals};
/// use azoth_eos::pure_saturation;
///
/// let r = pure_saturation(kelvins(369.83), pascals(4_248_000.0), 0.1523, kelvins(300.0))?;
/// assert!((r.p_sat.value / 1e5 - 9.9767).abs() < 1e-3);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tc`, `Pc` and `T` are the symbols in the chemistry
pub fn pure_saturation(
    Tc: ThermodynamicTemperature,
    Pc: Pressure,
    omega: f64,
    T: ThermodynamicTemperature,
) -> Result<PureSaturationResult> {
    let spec = &model_gen::PURE_SATURATION_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "Tc" => Some(Tc.value),
            "Pc" => Some(Pc.value),
            "omega" => Some(omega),
            "T" => Some(T.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let reduced_temperature = T.value / Tc.value;
    apply_checks(
        spec.derived_checks(),
        |quantity| (quantity == "t_over_tc").then_some(reduced_temperature),
        &mut warnings,
    )?;

    let algorithm = spec.algorithm;
    let Some(bracket) = algorithm.bracket else {
        // The spec's scheme is `saturation_pressure_bisection`, which searches a
        // bracket. `bracket` is optional in the schema because a *nested* scheme
        // need not have one - Rachford-Rice's interval comes from its K-values -
        // so a top-level scheme that needs one and has none is a generator bug.
        // Returned rather than panicked, per this library's no-panic rule.
        return Err(AzothError::InvalidInput {
            field: "algorithm.bracket".to_string(),
            reason: format!(
                "scheme `{}` searches for a bracket but the spec declares none",
                algorithm.scheme
            ),
        });
    };
    let kappa = pr_kappa(omega)?.kappa;

    // The two extreme fugacities at a trial reduced pressure, or `None` when the
    // cubic has no liquid branch there - which is exactly what `pr_z_factor`'s
    // `root_structure` reports, so this does not second-guess its admissibility rule.
    let pair_at = |pr: f64| -> Option<(f64, f64)> {
        let ab = pr_alpha_ab(kappa, reduced_temperature, pr).ok()?;
        let roots = pr_z_factor(ab.a_reduced, ab.b_reduced).ok()?;
        if roots.root_structure != RootStructure::Three {
            return None;
        }
        let liquid = pr_departure(
            ab.a_reduced,
            ab.b_reduced,
            roots.z_min,
            kappa,
            reduced_temperature,
        )
        .ok()?;
        let vapour = pr_departure(
            ab.a_reduced,
            ab.b_reduced,
            roots.z_max,
            kappa,
            reduced_temperature,
        )
        .ok()?;
        Some((liquid.ln_phi, vapour.ln_phi))
    };

    // 1. Bracket. Linear and with the spec's step count, so both implementations
    //    land on the same point - the scan's resolution is the one thing that moves
    //    the answer, and it is in the spec for that reason.
    let mut upper = None;
    for step in 0..bracket.steps {
        let fraction = f64::from(step) / f64::from(bracket.steps - 1);
        let pr = bracket.lower + (bracket.upper - bracket.lower) * fraction;
        if pair_at(pr).is_some() {
            upper = Some(pr);
        }
    }
    let Some(mut hi) = upper else {
        return Err(AzothError::OutOfRange {
            field: "t_over_tc".to_string(),
            value: reduced_temperature,
            detail: format!(
                "no pressure between Pr = {} and Pr = {} has a liquid branch, so there \
                 is nothing to equate. A below-critical temperature should always have \
                 one, so this means the scan's bounds are wrong rather than the state",
                bracket.lower, bracket.upper
            ),
        });
    };
    let mut lo = bracket.lower;
    let convergence = azoth_core::solver::Convergence::parse(algorithm.convergence)?;

    // 2. Bisect on the bracket's width.
    let mut iterations = 0;
    for step in 1..=algorithm.max_iterations {
        iterations = step;
        let mid = 0.5 * (lo + hi);
        let width = hi - lo;
        let close_enough = match convergence {
            azoth_core::solver::Convergence::Absolute => width <= algorithm.tolerance,
            azoth_core::solver::Convergence::Relative => width <= algorithm.tolerance * mid.abs(),
        };
        if close_enough {
            break;
        }
        let Some((liquid, vapour)) = pair_at(mid) else {
            // The bracket's upper end is the highest three-root pressure, so a
            // midpoint always has one; breaking rather than panicking keeps a
            // surprise here from taking the interpreter down.
            break;
        };
        if liquid - vapour > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    let reduced_pressure = 0.5 * (lo + hi);
    let residual = (hi - lo) / (2.0 * reduced_pressure.abs());
    let Some((liquid, _)) = pair_at(reduced_pressure) else {
        return Err(AzothError::SolverNotConverged {
            iterations,
            residual,
            tolerance: algorithm.tolerance,
        });
    };

    Ok(PureSaturationResult {
        p_sat: pascals(reduced_pressure * Pc.value),
        ln_phi: liquid,
        iterations,
        residual,
        warnings,
    })
}
