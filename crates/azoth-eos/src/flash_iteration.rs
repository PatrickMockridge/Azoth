//! The Rachford-Rice and successive-substitution scaffold the isothermal flashes
//! share.
//!
//! [`crate::pt_flash`] and [`crate::ge_nrtl_flash`] differ in exactly one thing: where
//! the K-value update's two fugacity coefficients come from. Peng-Robinson gives both
//! from the cubic; the gamma-phi flash takes the liquid's from an activity-coefficient
//! phase and the vapour's from SRK. Everything else - the physical bracket of the
//! Rachford-Rice root, the bisection on it, the convergence measure, what a trivial
//! solution is and why its vapour fraction is absent rather than zero - is the same
//! for both, and subtle enough that writing it twice would invite the two copies to
//! disagree.

use crate::results::Phase;

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
pub(crate) fn rachford_rice_bounds(k: &[f64]) -> Option<(f64, f64)> {
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
pub(crate) fn rachford_rice(
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
pub(crate) fn compositions(z: &[f64], k: &[f64], beta: f64) -> (Vec<f64>, Vec<f64>) {
    let x: Vec<f64> = z
        .iter()
        .zip(k)
        .map(|(&zi, &ki)| zi / (1.0 + beta * (ki - 1.0)))
        .collect();
    let y = k.iter().zip(&x).map(|(&ki, &xi)| ki * xi).collect();
    (x, y)
}

/// The rms change in `ln K` across one iteration.
pub(crate) fn rms_delta(ln_k_new: &[f64], k: &[f64]) -> f64 {
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
/// `|ln K_i| < 1e-6` for every `i` means the two phases have converged onto the
/// feed. It is compared against the *K-values* and never against `beta`, which is
/// indeterminate there.
///
/// **The threshold must not sit on the scale of the answer.** The second-order
/// scheme lands on `|ln K|` within `1e-8` of one at a state where the feed is single
/// phase (methane/n-butane 355 K, 120 bar), and `1e-8` made the trivial declaration a
/// coin flip on the last bit of the fugacity coefficient - a one-ulp change in
/// `Cubic::eos_z_minus_one` moved it from `trivial` at 128 iterations to no
/// convergence at all. A split has `|ln K|` of order one, so three further decades
/// cost nothing.
pub(crate) const TRIVIAL_TOLERANCE: f64 = 1.0e-06;

/// Whether the K-values have converged onto the feed.
///
/// Compared against the *K-values* and never against `beta`, which is indeterminate
/// there.
pub(crate) fn is_trivial(k: &[f64]) -> bool {
    k.iter().all(|value| value.ln().abs() < TRIVIAL_TOLERANCE)
}

/// How the iteration finished, when it finished somewhere other than convergence.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Outcome {
    /// Converged to `x = y = z`, where the two phases have met.
    Trivial,
    /// No Rachford-Rice root exists for these K-values, so the feed is single phase.
    SinglePhase(Phase),
}

/// The caveat every trivial solution carries.
pub(crate) fn trivial_warning() -> azoth_core::Warning {
    azoth_core::Warning::new(
        azoth_core::WarningCode::TrivialSolution,
        "the iteration converged to x = y = z, so the feed is single phase and there \
         is no vapour fraction. Which single phase it is, this model does not say - \
         that needs a stability analysis it does not perform. `beta` is absent rather \
         than zero, and `phase` is `trivial`.",
    )
}

/// The caveat a feed carries when no Rachford-Rice root exists at all.
pub(crate) fn single_phase_warning(phase: Phase) -> azoth_core::Warning {
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
pub(crate) fn negative_flash_warning(phase: Phase, beta: f64) -> azoth_core::Warning {
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
