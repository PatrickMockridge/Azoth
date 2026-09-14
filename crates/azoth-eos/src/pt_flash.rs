//! `eos.pt_flash` - the isothermal two-phase flash.
//!
//! Spec: `specs/models/eos/pt_flash.yaml`
//!
//! The two-phase split of a mixture at a fixed temperature and pressure: Wilson
//! K-value estimates, then successive substitution, with Rachford-Rice bisected each
//! iteration. A *model* rather than a calculation - what the spec pins down is the
//! procedure, and this module reads the procedure from the generated table rather than
//! choosing it.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result, apply_checks};

use crate::algorithm_of;
use crate::mixture::{Mixture, RootSide, wilson_k};
use crate::model_gen;

use crate::results::{Phase, PtFlashResult};

/// The interval on which Rachford-Rice has its physical root, or `None` if it has
/// none.
///
/// `g(beta) = sum_i z_i (K_i - 1) / (1 + beta (K_i - 1))` has poles at `1/(1 - K_i)`,
/// so the root the flash wants is the interval on which `1 + beta (K_i - 1) > 0` for
/// every `i` - the one that keeps `x_i = z_i / (1 + beta (K_i - 1))` and `y_i = K_i x_i`
/// non-negative:
///
/// ```text
/// max over {i : K_i > 1} of 1/(1 - K_i)  <  beta  <  min over {i : K_i < 1} of 1/(1 - K_i)
/// ```
///
/// It exists only when the K-values straddle one. `None` says they do not, which is a
/// proof that the feed has no two-phase solution at these K-values rather than a
/// numerical failure: `sum_i y_i = sum_i K_i x_i = 1` alongside `sum_i x_i = 1` needs
/// every `K_i` below one and above one at once.
fn rachford_rice_bounds(k: &[f64]) -> Option<(f64, f64)> {
    let mut lo = f64::NEG_INFINITY;
    let mut hi = f64::INFINITY;
    for &value in k {
        if value > 1.0 {
            lo = lo.max(1.0 / (1.0 - value));
        } else if value < 1.0 {
            hi = hi.min(1.0 / (1.0 - value));
        } else {
            // `K_i = 1` exactly puts a pole at infinity and makes `g` degenerate.
            // The caller detects this as the trivial solution before asking.
            return None;
        }
    }
    (lo.is_finite() && hi.is_finite()).then_some((lo, hi))
}

/// The vapour fraction that solves Rachford-Rice, by bisection.
///
/// `bounds` must be the interval [`rachford_rice_bounds`] returned, on which `g` is
/// continuous, strictly decreasing, and positive at the lower end.
///
/// Returns `(beta, iterations)`.
fn rachford_rice(
    z: &[f64],
    k: &[f64],
    bounds: (f64, f64),
    tolerance: f64,
    max_iterations: u32,
) -> (f64, u32) {
    let (mut lo, mut hi) = bounds;
    let g = |beta: f64| -> f64 {
        z.iter()
            .zip(k)
            .map(|(&zi, &ki)| zi * (ki - 1.0) / (1.0 + beta * (ki - 1.0)))
            .sum()
    };

    let mut iterations = 0;
    for step in 1..=max_iterations {
        iterations = step;
        let mid = 0.5 * (lo + hi);
        if hi - lo <= tolerance {
            break;
        }
        if g(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (0.5 * (lo + hi), iterations)
}

/// The compositions of the two phases at a vapour fraction.
fn compositions(z: &[f64], k: &[f64], beta: f64) -> (Vec<f64>, Vec<f64>) {
    let x: Vec<f64> = z
        .iter()
        .zip(k)
        .map(|(&zi, &ki)| zi / (1.0 + beta * (ki - 1.0)))
        .collect();
    let y = k.iter().zip(&x).map(|(&ki, &xi)| ki * xi).collect();
    (x, y)
}

/// The rms change in `ln K` across one iteration.
fn rms_delta(ln_k_new: &[f64], k: &[f64]) -> f64 {
    let n = ln_k_new.len();
    let total: f64 = ln_k_new
        .iter()
        .zip(k)
        .map(|(&new, &old)| {
            let delta = new - old.ln();
            delta * delta
        })
        .sum();
    (total / n as f64).sqrt()
}

/// The `ln K` below which the iteration has found the trivial solution.
///
/// `|ln K_i| < 1e-8` for every `i` means the two phases have converged onto the
/// feed. It is compared against the *K-values* and never against `beta`, which is
/// indeterminate there - see the spec's correction 2 for the feed a `beta` test would
/// have mislabelled.
const TRIVIAL_TOLERANCE: f64 = 1.0e-08;

/// The isothermal two-phase flash of a mixture at a temperature and pressure.
///
/// `z` is the overall composition, in mole fractions, and is **checked rather than
/// renormalised** - silently rescaling a caller's composition would make their error
/// invisible in a way that changes every number downstream.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
/// * [`AzothError::InvalidInput`] if `z` is the wrong length, has a negative entry,
///   or does not sum to one.
/// * [`AzothError::SolverNotConverged`] if the iteration hits its cap.
/// * Propagates the kernels' range checks.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, pascals};
/// use azoth_eos::{databank, pt_flash};
///
/// let mixture = databank::mixture_of(&["methane", "n-butane"])
///     .expect("the pair resolves")
///     .0;
/// let r = pt_flash(&mixture, kelvins(330.0), pascals(2_500_000.0), &[0.6, 0.4])?;
/// assert_eq!(r.phase, azoth_eos::Phase::TwoPhase);
/// assert!((r.beta.expect("a split has a vapour fraction") - 0.8422055475803881).abs() < 1e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn pt_flash(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<PtFlashResult> {
    let spec = &model_gen::PT_FLASH_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "P" => Some(p.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let n = mixture.len();
    if z.len() != n {
        return Err(AzothError::invalid_input(
            "z",
            format!("a feed for {n} components has {} entries", z.len()),
        ));
    }
    if let Some(bad) = z.iter().position(|&value| value < 0.0) {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "z[{bad}] is {} but a mole fraction cannot be negative",
                z[bad]
            ),
        ));
    }
    let sum: f64 = z.iter().sum();
    if (sum - 1.0).abs() > 1.0e-09 {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "the feed's mole fractions sum to {sum}, not to one. Renormalising it \
                 here would make a composition error invisible in every number \
                 downstream, so it is refused instead"
            ),
        ));
    }

    let min_t_over_tc = mixture
        .components()
        .iter()
        .map(|c| t.value / c.tc.value)
        .fold(f64::INFINITY, f64::min);
    apply_checks(
        spec.derived_checks(),
        |quantity| (quantity == "min_t_over_tc").then_some(min_t_over_tc),
        &mut warnings,
    )?;

    let reduced = mixture.reduced_parameters(t, p)?;
    warnings.extend(reduced.warnings.iter().cloned());

    let algorithm = algorithm_of(spec)?;
    let inner = algorithm.inner.ok_or_else(|| AzothError::InvalidInput {
        field: "algorithm.inner".to_string(),
        reason: format!(
            "scheme `{}` solves Rachford-Rice every iteration but the spec declares \
             no inner scheme",
            algorithm.scheme
        ),
    })?;

    let mut k = wilson_k(mixture, t, p);
    let mut iterations = 0;
    let mut residual = f64::NAN;

    // How the loop finished. `None` means it is still iterating, so a value that
    // is still `None` after the loop hit its cap is a genuine failure to converge.
    let mut settled: Option<Outcome> = None;

    for step in 1..=algorithm.max_iterations {
        iterations = step;

        // Both guards run *before* the Rachford-Rice solve, because both are cases
        // where the bracket is a division by zero: `K_i = 1` puts a pole at
        // infinity, and K-values all on one side of one put a bracket end there.
        if is_trivial(&k) {
            settled = Some(Outcome::Trivial);
            break;
        }
        // A K-value that is not finite and positive is the iteration diverging
        // rather than a state to diagnose. Guarded because every path below treats
        // the K-values as a composition ratio, and `NaN > 1.0` being false would
        // otherwise report a diverged iteration as a single-phase feed.
        if k.iter().any(|value| !value.is_finite() || *value <= 0.0) {
            return Err(AzothError::SolverNotConverged {
                iterations,
                residual,
                tolerance: algorithm.tolerance,
            });
        }
        let Some(bounds) = rachford_rice_bounds(&k) else {
            settled = Some(Outcome::SinglePhase(if k.iter().all(|&v| v > 1.0) {
                Phase::AllVapour
            } else {
                Phase::AllLiquid
            }));
            break;
        };

        let (beta, _) = rachford_rice(z, &k, bounds, inner.tolerance, inner.max_iterations);
        let (x, y) = compositions(z, &k, beta);
        let liquid = mixture.phase_state(&reduced, &x, RootSide::Liquid)?;
        let vapour = mixture.phase_state(&reduced, &y, RootSide::Vapour)?;

        let ln_k_new: Vec<f64> = liquid
            .ln_phi
            .iter()
            .zip(&vapour.ln_phi)
            .map(|(&l, &v)| l - v)
            .collect();
        residual = rms_delta(&ln_k_new, &k);
        k = ln_k_new.iter().map(|value| value.exp()).collect();

        if residual <= algorithm.tolerance {
            break;
        }
    }

    // Classify what the loop settled on. Re-evaluated once more at the converged
    // `K`, so the returned `beta`, `x`, `y` and `K` are one consistent state rather
    // than one state and its predecessor's vapour fraction.
    let (beta, x, y, phase) = match settled {
        // The trivial solution is stated *exactly* - `x = y = z` and `K = 1` -
        // rather than as the last iterate that approached it. That iterate is a
        // function of where the bisection stopped on an identically-zero function
        // and is not reproducible; this is. See the spec's correction 2.
        Some(Outcome::Trivial) => {
            k = vec![1.0; n];
            warnings.push(trivial_warning());
            (None, z.to_vec(), z.to_vec(), Phase::Trivial)
        }
        // No root at all: the feed is single phase by proof, and `beta` is absent
        // because there is no vapour fraction to report.
        Some(Outcome::SinglePhase(phase)) => {
            warnings.push(single_phase_warning(phase));
            (None, z.to_vec(), z.to_vec(), phase)
        }
        None if residual > algorithm.tolerance => {
            return Err(AzothError::SolverNotConverged {
                iterations,
                residual,
                tolerance: algorithm.tolerance,
            });
        }
        None => {
            if is_trivial(&k) {
                k = vec![1.0; n];
                warnings.push(trivial_warning());
                (None, z.to_vec(), z.to_vec(), Phase::Trivial)
            } else {
                let bounds = rachford_rice_bounds(&k).ok_or(AzothError::SolverNotConverged {
                    iterations,
                    residual,
                    tolerance: algorithm.tolerance,
                })?;
                let (beta, _) = rachford_rice(z, &k, bounds, inner.tolerance, inner.max_iterations);
                let (x, y) = compositions(z, &k, beta);
                let phase = match beta {
                    value if value < 0.0 => Phase::AllLiquid,
                    value if value > 1.0 => Phase::AllVapour,
                    _ => Phase::TwoPhase,
                };
                if phase != Phase::TwoPhase {
                    warnings.push(negative_flash_warning(phase, beta));
                }
                (Some(beta), x, y, phase)
            }
        }
    };

    let liquid = mixture.phase_state(&reduced, &x, RootSide::Liquid)?;
    let vapour = mixture.phase_state(&reduced, &y, RootSide::Vapour)?;

    Ok(PtFlashResult {
        beta,
        x,
        y,
        k,
        ln_phi_liquid: liquid.ln_phi,
        ln_phi_vapour: vapour.ln_phi,
        z_liquid: liquid.z,
        z_vapour: vapour.z,
        min_t_over_tc,
        phase,
        iterations,
        residual,
        warnings,
    })
}

/// How the iteration finished, when it finished somewhere other than convergence.
#[derive(Debug, Clone, Copy)]
enum Outcome {
    /// Converged to `x = y = z`, where the two phases have met.
    Trivial,
    /// No Rachford-Rice root exists for these K-values, so the feed is single phase.
    SinglePhase(Phase),
}

/// Whether the K-values have converged onto the feed.
///
/// Compared against the *K-values* and never against `beta`, which is indeterminate
/// there.
fn is_trivial(k: &[f64]) -> bool {
    k.iter().all(|value| value.ln().abs() < TRIVIAL_TOLERANCE)
}

/// The caveat every trivial solution carries.
fn trivial_warning() -> azoth_core::Warning {
    azoth_core::Warning::new(
        azoth_core::WarningCode::TrivialSolution,
        "the iteration converged to x = y = z, so the feed is single phase and there \
         is no vapour fraction. Which single phase it is, this model does not say - \
         that needs a stability analysis it does not perform. `beta` is absent rather \
         than zero, and `phase` is `trivial`.",
    )
}

/// The caveat a feed carries when no Rachford-Rice root exists at all.
fn single_phase_warning(phase: Phase) -> azoth_core::Warning {
    let which = match phase {
        Phase::AllVapour => "vapour",
        _ => "liquid",
    };
    azoth_core::Warning::new(
        azoth_core::WarningCode::TrivialSolution,
        format!(
            "every K-value is on the same side of one, so the Rachford-Rice equation \
             has no root and the feed has no two-phase solution at this temperature \
             and pressure. The feed is single-phase {which}, and `beta` is absent \
             rather than zero because there is no vapour fraction to report."
        ),
    )
}

/// The caveat a converged-but-out-of-range vapour fraction carries.
fn negative_flash_warning(phase: Phase, beta: f64) -> azoth_core::Warning {
    let which = match phase {
        Phase::AllVapour => "superheated vapour",
        _ => "subcooled liquid",
    };
    azoth_core::Warning::new(
        azoth_core::WarningCode::OutOfValidRange,
        format!(
            "the vapour fraction is {beta}, outside [0, 1], so the feed is \
             single-phase {which}. It is reported because the negative flash is a \
             real reading - it is the amount of the absent phase that would have to \
             be added to bring the feed to saturation - but it is not a phase split."
        ),
    )
}
