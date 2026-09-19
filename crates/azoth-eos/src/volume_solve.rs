//! The volume solve every equation of state in this crate reaches its state through.
//!
//! NeqSim's `PhasePCSAFTRahmat.calcVolume` is a damped Newton on the volume:
//! `h = P - P_calc`, `dh = -dP_calc/dv`, and `v += 0.9 h/dh`, stopping when the step is a
//! relative `1e-10` or after a hundred of them. Two things about it are the same whichever
//! equation is behind `P_calc`, and one is not:
//!
//! - **The vapour branch Newtons from the ideal gas and the liquid branch cannot.** At a
//!   pressure where the isotherm has one root that root is the vapour one, so a dilute seed
//!   converges there and stops; the liquid root is a separate zero the seed never points at,
//!   and it is found by walking up from the floor instead.
//! - **The floor is where the model's own terms diverge**, which each equation states for
//!   itself: for a packing-fraction model it is `(pi/6) N_A md3`, and below it there is no
//!   state rather than a small volume.
//!
//! So the arithmetic lives here once and each model passes its own residual, its own floor
//! and its own way of saying that a branch has no root.

use azoth_core::{AzothError, Result};

use crate::mixture::RootSide;

/// Solve `P_calc(v) = p` for the molar volume, returning `(v, Newton steps)`.
///
/// `residual(v)` is `(P_calc(v) - p, dP_calc/dv)` in SI - the value and the slope the
/// Newton step divides by. `ideal` is the ideal-gas volume `R T/p`, which is the vapour
/// seed and the top of the liquid walk's bracket; `floor` is the smallest volume the model
/// admits, below which the terms diverge.
///
/// The returned step count is `0` on the liquid branch, which is bracketed and bisected
/// rather than Newtoned, so a caller can tell the two paths apart - and, more to the point,
/// **see a vapour solve that did not converge**, which is otherwise invisible because the
/// walk answers and the walk finds the same root.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if the wanted branch has no zero at this state, which is a
///   real answer rather than a failure: the isotherm has no root on that side.
pub fn molar_volume(
    residual: impl Fn(f64) -> Result<(f64, f64)>,
    ideal: f64,
    floor: f64,
    side: RootSide,
    no_root: impl Fn(&str) -> AzothError,
) -> Result<(f64, usize)> {
    let value = |v: f64| -> Result<f64> { Ok(residual(v)?.0) };
    let slope = |v: f64| -> Result<f64> { Ok(residual(v)?.1) };

    // A bisection on a fixed count rather than on the residual: near the floor the residual
    // is a difference of large terms, and one that stopped at an absolute tolerance there
    // would stop early.
    fn bisect(value: &impl Fn(f64) -> Result<f64>, lo: f64, hi: f64) -> Result<f64> {
        let (mut lo, mut hi) = (lo, hi);
        let sign = value(lo)?.signum();
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if value(mid)?.signum() == sign {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Ok(0.5 * (lo + hi))
    }

    match side {
        RootSide::Vapour => {
            let mut v = ideal;
            let mut converged = None;
            for step in 0..100 {
                // `v += 0.9 (P - P_calc)/(dP_calc/dv)`, NeqSim's step: the residual is
                // `P_calc - P`, so the sign is Newton's with the residual negated.
                let delta = -0.9 * value(v)? / slope(v)?;
                if !delta.is_finite() || v + delta <= floor {
                    break;
                }
                let relative = delta.abs() / v;
                v += delta;
                if relative < 1.0e-10 {
                    converged = Some(step + 1);
                    break;
                }
            }
            match converged {
                Some(steps) => Ok((v, steps)),
                // Newton left the domain, a step was not a number, or it stopped moving.
                // The vapour root is the **last** zero in volume, so the walk is taken for
                // the last bracket rather than the first.
                None => {
                    let (lo, hi) = walk(&value, floor, ideal)?
                        .last()
                        .copied()
                        .ok_or_else(|| no_root("vapour"))?;
                    Ok((bisect(&value, lo, hi)?, 0))
                }
            }
        }
        RootSide::Liquid => {
            // **The lowest zero above the floor.** The walk is geometric because the root's
            // *ratio* to the floor is what is bounded, not its distance.
            let (lo, hi) = walk(&value, floor, ideal)?
                .first()
                .copied()
                .ok_or_else(|| no_root("liquid"))?;
            Ok((bisect(&value, lo, hi)?, 0))
        }
    }
}

/// The sign-change brackets of the residual between the floor and a volume far enough into
/// ideality that it can only be negative, in increasing volume.
///
/// The residual is `+inf` at the floor, so the walk starts a hair above it, and at the top
/// `P_calc v/(RT)` is within `1e-4` of one - below any pressure a caller asks for - so at
/// least one change is always present.
///
/// # Errors
/// As the residual.
pub fn walk(
    value: &impl Fn(f64) -> Result<f64>,
    floor: f64,
    ideal: f64,
) -> Result<Vec<(f64, f64)>> {
    let start = floor * (1.0 + 1.0e-9);
    let top = 1.0e4 * ideal;
    let steps = 128;
    let mut previous = start;
    let mut previous_value = value(previous)?;
    let mut out = Vec::new();
    for step in 1..=steps {
        let v = start * (top / start).powf(step as f64 / steps as f64);
        let current = value(v)?;
        if previous_value.signum() != current.signum() {
            out.push((previous, v));
        }
        previous = v;
        previous_value = current;
    }
    Ok(out)
}
