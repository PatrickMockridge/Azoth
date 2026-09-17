//! `eos.rachford_rice` - the vapour fraction that solves the Rachford-Rice equation.
//!
//! ```text
//! g(beta) = sum_i z_i (K_i - 1) / (1 + beta (K_i - 1))  =  0
//! ```
//!
//! A *model* rather than a calculation, and `kind = "procedure"`: it is a loop, so its
//! spec fixes the scheme, the convergence rule, the tolerance and the cap, and both
//! kernels run those rather than choosing them.
//!
//! Spec: `specs/models/eos/rachford_rice.toml`. NeqSim's `RachfordRice` is the port
//! source. The class carries two solvers behind a static `method` field whose default is
//! `Nielsen2023`, and **that default is the one ported**: nothing in NeqSim 3.20.0 calls
//! `setMethod`, so every flash in the library runs Nielsen's, and `calcBetaMichelsen2001`
//! is reachable only by calling it directly. The alternative is named in the spec's
//! assumptions with that measurement rather than ported as a branch nothing takes.
//!
//! Nielsen & Lia's formulation is worth a sentence because the arithmetic is not the
//! textbook one. It solves the same equation in a rescaled pair of unknowns -
//! `a = (alpha - alpha_min)/(alpha_max - alpha)` and `b = 1/(alpha - alpha_min)` for
//! `alpha = 1/(1 - K)` - so that a mixture spanning many orders of magnitude in `K` does
//! not lose precision to cancellation in `g` itself. When the root is above one half the
//! equation is solved in `1 / K` instead, which turns the search back around so the
//! arithmetic always runs on the interval with the root below one half.
//!
//! Two of NeqSim's returns are **not** ported, and the spec names both. `calcBeta`
//! declines to iterate when the root lies outside `[0, 1]` - it answers `g(0) < 0` with
//! the lower clamp and `g(1) > 0` with the upper one - and clamps what it does iterate to
//! `[1e-12, 1 - 1e-12]`. Both are the same move: they report *which single phase the feed
//! is* where azoth reports *what the vapour fraction is*. A negative vapour fraction is a
//! reading, not a failure - it is the amount of the absent phase that would have to be
//! added to reach saturation, which is what `eos.pt_flash` reports it as - so the root is
//! returned as the equation gives it and the clamp answers only the feeds where no root
//! exists at all.
//!
//! **`eos.pt_flash` does not call this**, and the reason is measured rather than
//! aesthetic: NeqSim stops at `1e-10` on its two rescaled residuals, which is about
//! `1e-11` in `beta`, and `eos.pt_flash` asserts the Rachford-Rice residual of every
//! answer at `1e-12`. A flash that took this solver would fail its own material-balance
//! invariant. The two procedures solve one equation and stand side by side.

use azoth_core::{AzothError, Result, apply_checks};

use crate::algorithm_of;
use crate::model_gen;
use crate::results::RachfordRiceResult;

/// NeqSim's `ThermodynamicModelSettings.phaseFractionMinimumLimit`.
///
/// The answer NeqSim returns for a feed that cannot split, at whichever end the feed
/// lies: a phase fraction of `1e-12` rather than an exact zero, so that a caller taking
/// a logarithm of it gets a number.
const PHASE_FRACTION_MINIMUM_LIMIT: f64 = 1.0e-12;

/// A K-value below this is an ion, which stays in the liquid and does not take part in
/// the split. NeqSim's own threshold, and its own reason.
const ION_THRESHOLD: f64 = 1.0e-30;

/// Which side of one the K-values lie on, and therefore whether a split exists.
enum Sides {
    /// Some `K` above one and some below: the equation has a root, and this returns it.
    Splits,
    /// No `K` above one. `g` has no root at all, and the feed is a liquid.
    AllLiquid,
    /// No `K` below one.
    AllVapour,
}

/// Whether the K-values straddle one.
///
/// `g(beta) = sum_i z_i (K_i - 1) / (1 + beta (K_i - 1))` has poles at `1/(1 - K_i)` and
/// is strictly decreasing between them, so it crosses zero exactly once on the interval
/// that keeps every `1 + beta (K_i - 1)` positive - and only if the K-values straddle
/// one, since the sum is then a mix of increasing and decreasing terms in `beta`. K-values
/// exactly one contribute `z_i * 0` to every term and are ignored here as they are in the
/// solver, which drops them by the same guard NeqSim uses.
///
/// It is not the same test as NeqSim's. NeqSim asks the sign of `g` at the two ends, which
/// answers "is the root inside `[0, 1]`"; this asks whether a root exists, which is the
/// weaker question and the one that leaves the root outside `[0, 1]` visible.
fn sides(k: &[f64]) -> Sides {
    let mut above = false;
    let mut below = false;
    for &value in k {
        if value < ION_THRESHOLD {
            continue;
        }
        if value > 1.0 {
            above = true;
        } else if value < 1.0 {
            below = true;
        }
    }
    match (above, below) {
        (true, true) => Sides::Splits,
        (true, false) => Sides::AllVapour,
        // Every K exactly one is the trivial solution, where `x = y = z`; the lower
        // clamp is an arbitrary answer to a question with none, and `eos.pt_flash`
        // detects that state before it asks.
        _ => Sides::AllLiquid,
    }
}

/// `g(beta) = sum_i z_i (K_i - 1) / (1 + beta (K_i - 1))`, ions excluded.
fn residual(z: &[f64], k: &[f64], beta: f64) -> f64 {
    z.iter()
        .zip(k)
        .filter(|&(_, &ki)| ki >= ION_THRESHOLD)
        .map(|(&zi, &ki)| zi * (ki - 1.0) / (1.0 + beta * (ki - 1.0)))
        .sum()
}

/// The root of the Rachford-Rice equation for a feed and a set of K-values.
///
/// `z` is the overall composition and `k` the K-values, `K_i = y_i / x_i`, at the same
/// temperature and pressure. A component with `K_i < 1e-30` is an ion and is skipped, as
/// NeqSim skips it.
///
/// The root is returned as the equation gives it and is **not clamped to `[0, 1]`**: a
/// negative vapour fraction is the negative flash, and a caller that wants NeqSim's
/// single-phase dispatch has to make it. A feed whose K-values do not straddle one has no
/// root at all, and that is not an error - the answer is NeqSim's `1e-12` at the end the
/// feed lies towards.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `z` and `k` differ in length, or `z` is empty.
/// * [`AzothError::OutOfRange`] if any K-value is not positive.
/// * [`AzothError::SolverNotConverged`] if the iteration hits its cap.
///
/// # Example
/// ```
/// use azoth_eos::rachford_rice::rachford_rice;
///
/// let r = rachford_rice(&[0.6, 0.4], &[7.304244305324782, 0.33749596785762953])?;
/// assert!((r.beta - 0.8422055475803871).abs() < 1e-12);
///
/// // Every K below one: no root, so NeqSim's clamp comes back rather than an error.
/// let single = rachford_rice(&[0.5, 0.5], &[0.2, 0.3])?;
/// assert!((single.beta - 1e-12).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn rachford_rice(z: &[f64], k: &[f64]) -> Result<RachfordRiceResult> {
    let spec = &model_gen::RACHFORD_RICE_SPEC;
    let algorithm = algorithm_of(spec)?;
    let mut warnings = Vec::new();

    if z.len() != k.len() {
        return Err(AzothError::invalid_input(
            "K",
            format!(
                "the feed has {} component(s) and there are {} K-value(s)",
                z.len(),
                k.len()
            ),
        ));
    }
    if z.is_empty() {
        return Err(AzothError::invalid_input(
            "z",
            "a feed of zero components has no vapour fraction",
        ));
    }
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "K" => k.iter().copied().reduce(f64::min),
            _ => None,
        },
        &mut warnings,
    )?;

    let beta = match sides(k) {
        Sides::AllLiquid => PHASE_FRACTION_MINIMUM_LIMIT,
        Sides::AllVapour => 1.0 - PHASE_FRACTION_MINIMUM_LIMIT,
        Sides::Splits => nielsen_2023(z, k, algorithm.tolerance, algorithm.max_iterations)?,
    };

    if !beta.is_finite() {
        return Err(AzothError::SolverNotConverged {
            iterations: 0,
            residual: f64::NAN,
            tolerance: algorithm.tolerance,
        });
    }

    Ok(RachfordRiceResult { beta, warnings })
}

/// NeqSim's `calcBetaNielsen2023`: the same equation solved in Nielsen & Lia's rescaled
/// pair, so a mixture spanning many orders of magnitude in `K` keeps its precision.
///
/// Transcribed from NeqSim's method with the two single-phase returns removed, so it is
/// called only where a root exists and it always finds one.
fn nielsen_2023(z: &[f64], k: &[f64], tolerance: f64, max_iterations: u32) -> Result<f64> {
    // `h` is `g` at the starting point, and its sign selects which unknown pair is used.
    // Solving in `1 / K` when the root is above one half turns the search back around, so
    // the arithmetic always runs on the interval where the root is below one half.
    let h = residual(z, k, 0.5);

    let work_k: Vec<f64> = if h > 0.0 {
        k.iter()
            .map(|&ki| if ki < ION_THRESHOLD { ki } else { 1.0 / ki })
            .collect()
    } else {
        k.to_vec()
    };

    let mut k_max = 0.0;
    let mut k_min = f64::MAX;
    let mut found = false;
    for &ki in &work_k {
        if ki < ION_THRESHOLD {
            continue;
        }
        if !found {
            k_max = ki;
            k_min = ki;
            found = true;
        } else if ki < k_min {
            k_min = ki;
        } else if ki > k_max {
            k_max = ki;
        }
    }
    if !found {
        // Every component is an ion, so there is no split to find.
        return Ok(PHASE_FRACTION_MINIMUM_LIMIT);
    }

    let alpha_min = 1.0 / (1.0 - k_max);
    let alpha_max = 1.0 / (1.0 - k_min);

    let alpha = 0.5;
    let mut a = (alpha - alpha_min) / (alpha_max - alpha);
    let mut b = 1.0 / (alpha - alpha_min);

    // The per-component constants of the rescaled form. `c` is the reciprocal of the
    // distance from `K = 1` and `d` rescales it into the `alpha` interval. A component
    // with `K` exactly one has no distance to rescale, and NeqSim's denominator floor
    // sends its term to zero in `funk` and in `hb` alike - which is what the equation
    // says it contributes.
    let mut c = vec![0.0; work_k.len()];
    let mut d = vec![0.0; work_k.len()];
    for (i, &ki) in work_k.iter().enumerate() {
        if ki < ION_THRESHOLD {
            continue;
        }
        let ki = ki.clamp(1.0e-25, 1.0e25);
        let mut denom = 1.0 - ki;
        if denom.abs() < 1.0e-25 {
            denom = if denom < 0.0 { -1.0e-25 } else { 1.0e-25 };
        }
        c[i] = 1.0 / denom;
        d[i] = (alpha_min - c[i]) / (alpha_max - alpha_min);
    }

    let mut a_max = alpha;
    let mut b_max = 1.0e20;
    let mut a_min = 0.0;
    let mut b_min = 1.0 / (alpha_max - alpha_min);

    let mut iterations: u32 = 0;
    let (funk, hb) = loop {
        iterations += 1;
        let mut f = 0.0;
        let mut f_deriv = 0.0;
        let mut h = 0.0;
        let mut h_deriv = 0.0;

        for (i, &ki) in work_k.iter().enumerate() {
            if ki < ION_THRESHOLD {
                continue;
            }
            f -= z[i] * a * (1.0 + a) / (d[i] + a * (1.0 + d[i]));
            f_deriv -=
                z[i] * (a * a + (1.0 + a) * (1.0 + a) * d[i]) / (d[i] + a * (1.0 + d[i])).powi(2);
            h += z[i] * b / (1.0 + b * (alpha_min - c[i]));
            h_deriv += z[i] / (1.0 + b * (alpha_min - c[i])).powi(2);
        }

        if (f.abs() < tolerance && h.abs() < tolerance) || iterations >= max_iterations {
            break (f, h);
        }

        if f > 0.0 {
            a_max = a;
        } else {
            a_min = a;
        }
        if h > 0.0 {
            b_max = b;
        } else {
            b_min = b;
        }

        a -= f / f_deriv;
        if a > a_max || a < a_min {
            a = 0.5 * (a_max + a_min);
        }
        b -= h / h_deriv;
        if b > b_max || b < b_min {
            b = 0.5 * (b_max + b_min);
        }
    };

    if funk.abs() >= tolerance || hb.abs() >= tolerance {
        return Err(AzothError::SolverNotConverged {
            iterations,
            residual: funk.abs().max(hb.abs()),
            tolerance,
        });
    }

    let mut v = -((1.0 / b) / a - alpha_max);
    if h > 0.0 {
        v = 1.0 - v;
    }
    Ok(v)
}
