//! `eos.pure_saturation` - pure-component saturation pressure.
//!
//! The pressure at which a pure component's vapour and liquid roots have equal
//! fugacity, found by bisection on reduced pressure.
//!
//! Spec: `specs/models/eos/pure_saturation.yaml`, which carries the procedure - the
//! bracketing rule, the tolerance and the cap - since a procedure that differs between
//! two implementations reaches a slightly different answer.
//!
//! Every number is a registered calc's: [`crate::pr_kappa`], [`crate::pr_alpha_ab`],
//! [`crate::pr_z_factor`] and [`crate::pr_departure`]. Only the search is here.

use azoth_core::units::{Pressure, ThermodynamicTemperature, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::algorithm_of;
use crate::model_gen;
use crate::results::{PureSaturationResult, RootStructure};
use crate::{pr_alpha_ab, pr_departure, pr_kappa, pr_z_factor};

/// The saturation pressure of a pure component at a temperature.
///
/// `Tc`, `Pc` and `omega` are the caller's. A caller who wants them from a name
/// takes them off `crate::databank::entry`, the same lookup `databank::mixture_of`
/// performs on the mixture path; this function is the scalar kernel underneath and
/// has no opinion about where its numbers came from. `T` must be below `Tc`.
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
/// assert!((r.p_sat.value / 1e5 - 9.9791).abs() < 1e-3);
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

    // `Tc` and `Pc` are refused here rather than by a bound in the spec. They used to
    // be declared `valid_range` entries: the model took the constants as arguments, so
    // a bound on them was a bound on an input. The model names a *component* now and
    // the constants come from the databank, which leaves no input for a bound to
    // guard - and a bound the generator emits against nothing is a check that never
    // fires, which is worse than no check at all. `Component::new` refuses a
    // non-positive pair on the mixture path; this is the same refusal on the scalar
    // one, and it is what keeps `Tc` from being a divisor by zero.
    // A NaN is refused alongside zero and the negatives: `NaN > 0.0` is false, so the
    // comparison is written as a refusal rather than as an acceptance, and a NaN
    // critical temperature would otherwise propagate silently through every root.
    if Tc.value.is_nan() || Tc.value <= 0.0 {
        return Err(AzothError::out_of_range(
            "Tc",
            Tc.value,
            "a critical temperature is absolute and positive by definition; it is also \
             a divisor in the reduced temperature",
        ));
    }
    if Pc.value.is_nan() || Pc.value <= 0.0 {
        return Err(AzothError::out_of_range(
            "Pc",
            Pc.value,
            "a critical pressure is positive by definition, and it converts the reduced \
             answer back to pascals",
        ));
    }

    apply_checks(
        spec.input_checks(),
        |quantity| (quantity == "T").then_some(T.value),
        &mut warnings,
    )?;

    let reduced_temperature = T.value / Tc.value;
    apply_checks(
        spec.derived_checks(),
        |quantity| (quantity == "t_over_tc").then_some(reduced_temperature),
        &mut warnings,
    )?;

    let algorithm = algorithm_of(spec)?;
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

    // 1. Bracket. Linear and with the spec's step count, so both implementations land
    //    on the same point.
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
