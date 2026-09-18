//! `eos.pt_phase_envelope` - the Michelsen natural-parameter phase envelope.
//!
//! Spec: `specs/models/eos/pt_phase_envelope.toml`
//!
//! The PT phase envelope is the locus where a mixture of composition `z` is at a phase
//! boundary. It is traced as two independent natural-parameter continuations - the
//! bubble branch at a tiny vapour fraction and the dew branch at one minus it - each
//! bootstrapped at a low pressure where its K-values are far from one, then continued
//! upward to the critical point where every K-value collapses to one. The `N+2`
//! unknowns `(ln K_i, ln T, ln P)` satisfy `N` isofugacity equations, one material
//! balance, and one specification equation that pins the most sensitive variable so the
//! trace can turn the corner at the cricondenbar and cricondentherm.
//!
//! The next point is predicted by a cubic through the last four converged points and
//! then polished by a damped Newton step. The Jacobian is taken by central difference,
//! the same choice the flash-property solver makes rather than porting NeqSim's analytic
//! fugacity-derivative surface.

use azoth_core::units::{Pressure, kelvins, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::mixture::{Mixture, RootSide};
use crate::model_gen;
use crate::phase_boundary::{self, WILSON_CONSTANT};
use crate::results::PtPhaseEnvelopeResult;

/// The vapour fraction on the bubble branch: `x_i ~= z_i`, so the feed is liquid.
///
/// Small but not below the central-difference step, so the numerical Jacobian can
/// resolve how the compositions depend on the K-values.
const BUBBLE_BETA: f64 = 1.0e-6;
/// The vapour fraction on the dew branch: `y_i ~= z_i`, so the feed is vapour.
const DEW_BETA: f64 = 1.0 - 1.0e-6;
/// Maximum temperature step per continuation point, in kelvin.
const D_TMAX: f64 = 10.0;
/// Maximum pressure step per continuation point, in pascal (10 bar).
const D_PMAX: f64 = 10.0e5;
/// The pressure ceiling, in pascal (1000 bar).
const MAX_PRESSURE: f64 = 1000.0e5;
/// The Newton tolerance.
const NEWTON_TOL: f64 = 1.0e-5;
/// The Newton iteration cap per point.
const NEWTON_MAX: u32 = 50;
/// The number of continuation points that pin pressure as the specification.
const PRESSURE_SPEC_POINTS: usize = 5;
/// How many converged points the cubic predictor fits through.
const HISTORY_LEN: usize = 4;

/// Solve `A x = b` by Gaussian elimination with partial pivoting, in place.
#[allow(clippy::needless_range_loop)] // the elimination is row-wise; indices read clearly
fn solve_linear(a: &mut [Vec<f64>], b: &mut [f64]) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        let mut pivot = col;
        for row in (col + 1)..n {
            if a[row][col].abs() > a[pivot][col].abs() {
                pivot = row;
            }
        }
        if a[pivot][col].abs() < 1.0e-15 {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in (col + 1)..n {
            let factor = a[row][col] / a[col][col];
            for j in col..n {
                a[row][j] -= factor * a[col][j];
            }
            b[row] -= factor * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let mut s = b[row];
        for j in (row + 1)..n {
            s -= a[row][j] * x[j];
        }
        x[row] = s / a[row][row];
    }
    Some(x)
}

/// The two phase compositions at a vapour fraction, from `K_i = y_i / x_i`.
fn compositions(z: &[f64], k: &[f64], beta: f64) -> (Vec<f64>, Vec<f64>) {
    let n = z.len();
    let mut x = vec![0.0; n];
    let mut y = vec![0.0; n];
    for i in 0..n {
        let d = 1.0 - beta + beta * k[i];
        x[i] = z[i] / d;
        y[i] = k[i] * x[i];
    }
    (x, y)
}

/// What the system's last row pins: an ordinary continuation point, or the critical point.
///
/// The critical row is `sum_i (ln K_i)^2 = 0`, which is zero exactly when every K-value is one.
/// Its Jacobian row is `2 ln K_i` and zero in the two state columns, which the central
/// difference below produces without being told.
#[derive(Debug, Clone, Copy)]
enum LastRow {
    /// A state variable pinned to a value: row `u[index] - value`.
    Pinned { index: usize, value: f64 },
    /// `sum_i (ln K_i)^2`, the criticality condition.
    Critical,
}

/// The `N+2` residuals: `[isofugacity; material balance; last row]`.
///
/// `u` is `[ln K_i; ln T; ln P]`.
fn residual(mixture: &Mixture, u: &[f64], beta: f64, z: &[f64], last: LastRow) -> Result<Vec<f64>> {
    let n = z.len();
    let k: Vec<f64> = u[..n].iter().map(|&lnk| lnk.exp()).collect();
    let t = u[n].exp();
    let p = u[n + 1].exp();
    let (x, y) = compositions(z, &k, beta);
    let reduced = mixture.reduced_parameters(kelvins(t), pascals(p))?;
    let liquid = mixture.phase_state(&reduced, &x, RootSide::Liquid)?;
    let vapour = mixture.phase_state(&reduced, &y, RootSide::Vapour)?;

    let mut f = vec![0.0; n + 2];
    for i in 0..n {
        // ln K_i = ln phi_i^L - ln phi_i^V.
        f[i] = u[i] - liquid.ln_phi[i] + vapour.ln_phi[i];
    }
    f[n] = y.iter().sum::<f64>() - x.iter().sum::<f64>();
    f[n + 1] = match last {
        LastRow::Pinned { index, value } => u[index] - value,
        LastRow::Critical => u[..n].iter().map(|value| value * value).sum(),
    };
    Ok(f)
}

/// The `N+2` by `N+2` Jacobian, by central difference over the state variables.
fn jacobian(
    mixture: &Mixture,
    u: &[f64],
    beta: f64,
    z: &[f64],
    last: LastRow,
) -> Result<Vec<Vec<f64>>> {
    let m = u.len();
    let f0 = residual(mixture, u, beta, z, last)?;
    let mut jac = vec![vec![0.0; m]; m];
    for j in 0..m {
        let delta = (1.0e-6 * u[j].abs()).max(1.0e-8);
        let mut up = u.to_vec();
        up[j] += delta;
        let fp = residual(mixture, &up, beta, z, last)?;
        for (i, row) in jac.iter_mut().enumerate() {
            row[j] = (fp[i] - f0[i]) / delta;
        }
    }
    Ok(jac)
}

fn norm2(v: &[f64]) -> f64 {
    v.iter().map(|&x| x * x).sum::<f64>().sqrt()
}

/// Clamp the predicted state back into the physical region. The cubic predictor can
/// overshoot a steep, nearly-critical segment into a state the root finder cannot
/// evaluate.
fn clamp_state(u: &mut [f64], n: usize) {
    for value in u.iter_mut().take(n) {
        if !value.is_finite() {
            *value = 0.0;
        }
        *value = value.clamp(-20.0, 20.0);
    }
    if !u[n].is_finite() {
        u[n] = 300.0_f64.ln();
    }
    if !u[n + 1].is_finite() {
        u[n + 1] = 1.0e5_f64.ln();
    }
    let t = u[n].exp().clamp(10.0, 2000.0);
    let p = u[n + 1].exp().clamp(1.0, MAX_PRESSURE * 1.5);
    u[n] = t.ln();
    u[n + 1] = p.ln();
}

/// One damped Newton solve at a fixed specification. Returns the solution and the
/// number of Newton iterations taken, which the step-size control reads.
fn newton(
    mixture: &Mixture,
    u: &[f64],
    beta: f64,
    z: &[f64],
    last: LastRow,
) -> Result<(Vec<f64>, u32)> {
    let mut u = u.to_vec();
    let mut iters = 0;
    for _ in 0..NEWTON_MAX {
        iters += 1;
        let f = residual(mixture, &u, beta, z, last)?;
        let norm = norm2(&f);
        if norm < NEWTON_TOL {
            return Ok((u, iters));
        }
        let mut jac = jacobian(mixture, &u, beta, z, last)?;
        let mut b = f.clone();
        let dx = solve_linear(&mut jac, &mut b).ok_or(AzothError::SolverNotConverged {
            iterations: iters,
            residual: norm,
            tolerance: NEWTON_TOL,
        })?;
        if dx.iter().any(|&d| !d.is_finite()) {
            return Err(AzothError::SolverNotConverged {
                iterations: iters,
                residual: norm,
                tolerance: NEWTON_TOL,
            });
        }
        let mut step = 1.0;
        let mut next = u.clone();
        for _ in 0..24 {
            for (i, &ui) in u.iter().enumerate() {
                next[i] = ui - step * dx[i];
            }
            clamp_state(&mut next, z.len());
            let fnext = residual(mixture, &next, beta, z, last)?;
            if norm2(&fnext) < norm {
                break;
            }
            step *= 0.5;
        }
        u = next;
    }
    Err(AzothError::SolverNotConverged {
        iterations: iters,
        residual: norm2(&residual(mixture, &u, beta, z, last)?),
        tolerance: NEWTON_TOL,
    })
}

/// The continuation sensitivity `du/d(specval)`: `J^-1` applied to a unit forcing on
/// the specification equation, the last row.
fn sensitivity(
    mixture: &Mixture,
    u: &[f64],
    beta: f64,
    z: &[f64],
    speceq: usize,
) -> Result<Vec<f64>> {
    let n = z.len();
    let mut jac = jacobian(
        mixture,
        u,
        beta,
        z,
        LastRow::Pinned {
            index: speceq,
            value: u[speceq],
        },
    )?;
    let mut e = vec![0.0; n + 2];
    e[n + 1] = 1.0;
    solve_linear(&mut jac, &mut e).ok_or(AzothError::SolverNotConverged {
        iterations: 0,
        residual: f64::NAN,
        tolerance: NEWTON_TOL,
    })
}

/// The index of the variable the continuation pins: `ln T` or `ln P`, whichever is
/// more sensitive. The K-values are excluded because they are not monotonic along the
/// envelope, and pinning one makes the trace oscillate.
fn find_spec(dxds: &[f64], n: usize) -> usize {
    if dxds[n].abs() >= dxds[n + 1].abs() {
        n
    } else {
        n + 1
    }
}

/// Fit a cubic through `history` (four states) as a function of `speceq`, and evaluate
/// every variable at `sny`.
fn cubic_predict(history: &[Vec<f64>], speceq: usize, sny: f64) -> Vec<f64> {
    let m = history[0].len();
    let s: Vec<f64> = history.iter().map(|h| h[speceq]).collect();
    let mut vandermonde = vec![vec![0.0; HISTORY_LEN]; HISTORY_LEN];
    for (i, &si) in s.iter().enumerate() {
        vandermonde[i][0] = 1.0;
        vandermonde[i][1] = si;
        vandermonde[i][2] = si * si;
        vandermonde[i][3] = si * si * si;
    }
    let mut result = vec![0.0; m];
    for j in 0..m {
        let v: Vec<f64> = history.iter().map(|h| h[j]).collect();
        let mut a = vandermonde.clone();
        let mut b = v.clone();
        let coeffs = solve_linear(&mut a, &mut b).unwrap_or(vec![v[3], 0.0, 0.0, 0.0]);
        result[j] = coeffs[0] + sny * (coeffs[1] + sny * (coeffs[2] + sny * coeffs[3]));
    }
    result
}

/// The most Newton steps the critical refinement takes, upstream's ten.
/// Upstream's ten, and it is a guard rather than a cap: measured, fifty steps let the Newton
/// wander to `355.78 K` with a branch gap of `21.0 K`.
const CRIT_MAX: u32 = 10;
/// Its convergence tolerance, upstream's.
///
/// Unreachable on this system in double precision, and that is a fact about the system rather
/// than about the constant: the refinement stops at its best iterate instead of converging, so
/// the tolerance is what it is aimed at rather than what it attains. Kept at upstream's value
/// because lowering it changed nothing - measured, `1e-7` returned the same points.
const CRIT_TOL: f64 = 1.0e-10;
/// The largest single correction it accepts, upstream's `|dx_j| > 0.5` abort.
const CRIT_STEP: f64 = 0.5;
/// How far a candidate critical point may sit from where the refinement started, relative.
///
/// Upstream compares each candidate against its polynomial's extrapolation and keeps it only
/// when the deviations are under `0.10` in temperature and `0.20` in pressure. The guard is
/// what stops a diverging Newton's wild iterate from winning the best-tracking on
/// `sum (ln K)^2` alone, which it otherwise does: measured, a step that jumped 373 K to 437 K
/// had the smallest sum of any iterate. The starting state is used here instead of the
/// polynomial, because the refinement starts where the trace crossed the K-window and the two
/// are within a few per cent of each other there.
const CRIT_T_BAND: f64 = 0.10;
const CRIT_P_BAND: f64 = 0.20;

/// The critical point, refined from a state the trace has brought close to it.
///
/// **This is what replaces a test with a computation.** The trace stops its branches where the
/// lightest component's K-value falls below `1.05` and the heaviest's rises above `0.95`, which
/// is a place a heuristic happened to fire: it put azoth's critical point 6.86 K and 4.54 bar
/// from NeqSim's. Upstream refines the same crossing with a Newton on `sum_i (ln K_i)^2 = 0`
/// alongside the `N` isofugacity rows and the material balance, which is the definition of the
/// critical point rather than a proxy for it. Measured on methane/n-butane 50/50, the refinement
/// moves the answer from `367.4573 K / 103.5078 bar` to `375.5008 K / 97.1374 bar` against
/// NeqSim's `374.3214 / 98.9683`, and closes the two branches' gap from 9.40 K to 2.64 K.
///
/// Returns `None` if it does not converge, which the caller reports as the unrefined point rather
/// than as an error: the trace's own answer is still a boundary. The **best** iterate is kept as
/// well as the last, upstream's rule and for its reason - the Newton can step past the point
/// where `sum (ln K)^2` is smallest and come back with a larger one.
fn calc_crit(mixture: &Mixture, u0: &[f64], beta: f64, z: &[f64]) -> Result<Option<(f64, f64)>> {
    let n = z.len();
    let mut u = u0.to_vec();
    let start = (u0[n].exp(), u0[n + 1].exp());
    let mut best: Option<(f64, f64, f64)> = None;

    for _ in 0..CRIT_MAX {
        let f = residual(mixture, &u, beta, z, LastRow::Critical)?;
        let sum_ln_k2: f64 = u[..n].iter().map(|value| value * value).sum();
        let (t, p) = (u[n].exp(), u[n + 1].exp());
        let inside = (t - start.0).abs() <= CRIT_T_BAND * start.0
            && (p - start.1).abs() <= CRIT_P_BAND * start.1;
        if inside
            && sum_ln_k2.is_finite()
            && best.is_none_or(|(smallest, _, _)| sum_ln_k2 < smallest)
        {
            best = Some((sum_ln_k2, t, p));
        }
        if norm2(&f) < CRIT_TOL {
            return Ok(Some((t, p)));
        }

        let mut jac = jacobian(mixture, &u, beta, z, LastRow::Critical)?;
        // **Levenberg-Marquardt, and the reason the refinement is deterministic.** The system
        // is *singular* at the critical point - that is what makes it critical - so the Newton
        // is ill-conditioned exactly where it is aimed, and it stops at an iterate chosen by a
        // knife-edge comparison. Two implementations of the same arithmetic then return points
        // 0.035 K apart, which is not a defect in either: measured, methane/n-butane gave
        // `375.9999` in one kernel and `376.0348` in the other before this. The regulariser is
        // `flash_newton.rs`'s, at the same magnitude, so the solve stays away from the
        // singularity and the iteration converges rather than stopping beside it.
        let trace: f64 = (0..n + 2).map(|i| jac[i][i].abs()).sum();
        let lambda = 1.0e-8 * trace / ((n + 2) as f64);
        for (i, row) in jac.iter_mut().enumerate() {
            row[i] += lambda;
        }
        let mut b = f.clone();
        let Some(dx) = solve_linear(&mut jac, &mut b) else {
            break;
        };
        if dx.iter().any(|d| !d.is_finite() || d.abs() > CRIT_STEP) {
            break;
        }
        if norm2(&dx) < CRIT_TOL {
            return Ok(Some((t, p)));
        }
        // **The step-halving line search the continuation Newton uses, and it is not optional
        // here.** Upstream takes the full step because it restarts the refinement from a
        // polynomial extrapolated to the `K = 1` point, which is already inside the critical
        // basin; this starts from wherever the trace crossed the `K`-window. Measured without
        // the search: `376.696 -> 375.501 -> 403.180`, the last outside the deviation band.
        let norm = norm2(&f);
        let mut step = 1.0;
        let mut accepted = false;
        for _ in 0..24 {
            let mut next = u.clone();
            for (i, &d) in dx.iter().enumerate() {
                next[i] = u[i] - step * d;
            }
            clamp_state(&mut next, n);
            let fnext = residual(mixture, &next, beta, z, LastRow::Critical)?;
            if norm2(&fnext) < norm {
                u = next;
                accepted = true;
                break;
            }
            step *= 0.5;
        }
        if !accepted {
            break;
        }
    }

    Ok(best.map(|(_, t, p)| (t, p)))
}

/// The indices of the lightest and heaviest component by critical temperature.
fn light_heavy(mixture: &Mixture) -> (usize, usize) {
    let mut lc = 0;
    let mut hc = 0;
    for (i, c) in mixture.components().iter().enumerate() {
        if c.tc.value < mixture.components()[lc].tc.value {
            lc = i;
        }
        if c.tc.value > mixture.components()[hc].tc.value {
            hc = i;
        }
    }
    (lc, hc)
}

/// Wilson's estimate of the bubble (`beta < 0.5`) or dew (`beta > 0.5`) temperature at
/// a pressure, iterated on the Wilson K-value sum until `sum z_i K_i = 1` (or
/// `sum z_i / K_i = 1`).
fn wilson_temperature(mixture: &Mixture, p: Pressure, z: &[f64], beta: f64) -> f64 {
    let n = z.len();
    let (lc, hc) = light_heavy(mixture);
    let idx = if beta < 0.5 { lc } else { hc };
    let c = &mixture.components()[idx];
    let ln_p_ratio = (p.value / c.pc.value).ln();
    let mut t = c.tc.value * WILSON_CONSTANT * (1.0 + c.omega)
        / (WILSON_CONSTANT * (1.0 + c.omega) - ln_p_ratio);
    let mut told = 0.0;
    for _ in 0..1000 {
        let psat = phase_boundary::wilson_psat(mixture, kelvins(t));
        let kwil: Vec<f64> = psat.iter().map(|&ps| ps / p.value).collect();
        let (mut s, mut dsdt) = (0.0, 0.0);
        for i in 0..n {
            let ci = &mixture.components()[i];
            let dlnkdt = WILSON_CONSTANT * (1.0 + ci.omega) * ci.tc.value / (t * t);
            if beta < 0.5 {
                s += z[i] * kwil[i];
                dsdt += z[i] * kwil[i] * dlnkdt;
            } else {
                s += z[i] / kwil[i];
                dsdt -= z[i] / kwil[i] * dlnkdt;
            }
        }
        s -= 1.0;
        if (s / dsdt).abs() > 0.1 * t {
            t -= 0.001 * s / dsdt;
        } else {
            t -= s / dsdt;
        }
        if (t - told).abs() < 1e-5 {
            break;
        }
        told = t;
    }
    t
}

/// One branch of the envelope.
struct BranchTrace {
    /// The traced temperatures.
    t: Vec<f64>,
    /// The traced pressures, parallel to `t`.
    p: Vec<f64>,
    /// The cricondenbar, `(T, P)`.
    cricondenbar: (f64, f64),
    /// The cricondentherm, `(T, P)`.
    cricondentherm: (f64, f64),
    /// The critical point this branch refined to, when the refinement converged.
    ///
    /// The two branches refine independently and should land on the same point, which is the
    /// evidence that they meet and what `residual` reports.
    critical: Option<(f64, f64)>,
}

/// One branch of the envelope, traced upward from a low pressure to the critical point.
fn trace_branch(mixture: &Mixture, p: Pressure, z: &[f64], beta: f64) -> Result<BranchTrace> {
    let spec = &model_gen::PT_PHASE_ENVELOPE_SPEC;
    let max_iterations = spec
        .algorithm
        .as_ref()
        .map(|a| a.max_iterations)
        .unwrap_or(9980) as usize;
    let n = z.len();
    let (lc, hc) = light_heavy(mixture);

    let start_t = wilson_temperature(mixture, p, z, beta);
    let psat = phase_boundary::wilson_psat(mixture, kelvins(start_t));
    let mut u = vec![0.0; n + 2];
    for (i, &ps) in psat.iter().enumerate() {
        u[i] = (ps / p.value).ln();
    }
    u[n] = start_t.ln();
    u[n + 1] = p.value.ln();

    let mut speceq = n + 1;
    let mut specval = u[speceq];
    let mut t_vec: Vec<f64> = Vec::new();
    let mut p_vec: Vec<f64> = Vec::new();
    let mut cricondenbar = (start_t, p.value);
    let mut cricondentherm = (start_t, p.value);
    let mut critical: Option<(f64, f64)> = None;
    let mut history: Vec<Vec<f64>> = Vec::new();
    let mut dxds: Vec<f64> = vec![0.0; n + 2];
    let mut ds = 0.1;
    let mut last_iters = 2;

    let (converged, _) = newton(
        mixture,
        &u,
        beta,
        z,
        LastRow::Pinned {
            index: speceq,
            value: specval,
        },
    )?;
    u = converged;

    for _ in 0..max_iterations {
        let t = u[n].exp();
        let pv = u[n + 1].exp();
        let k: Vec<f64> = u[..n].iter().map(|&lnk| lnk.exp()).collect();

        if t > cricondentherm.0 {
            cricondentherm = (t, pv);
        }
        if pv > cricondenbar.1 {
            cricondenbar = (t, pv);
        }

        // The critical point: every K-value approaches one. The crossing is what *detects* it
        // and the refinement above is what computes it.
        if k[lc] < 1.05 && k[hc] > 0.95 {
            critical = calc_crit(mixture, &u, beta, z)?;
            break;
        }

        t_vec.push(t);
        p_vec.push(pv);
        history.push(u.clone());
        if history.len() > HISTORY_LEN {
            history.remove(0);
        }

        // Exit when the pressure ceiling is reached.
        if pv > MAX_PRESSURE {
            break;
        }

        // Predict the next point.
        if t_vec.len() <= PRESSURE_SPEC_POINTS {
            // Pin pressure for the first few points: step `ln P` up by a fixed amount.
            speceq = n + 1;
            dxds = sensitivity(mixture, &u, beta, z, speceq)?;
            ds = 0.1 / dxds[n + 1];
            for (i, &d) in dxds.iter().enumerate() {
                u[i] += d * ds;
            }
            specval = u[n + 1];
        } else {
            // Pick the most sensitive variable, re-derive the sensitivity for it, and
            // predict by cubic extrapolation through the last four points.
            speceq = find_spec(&dxds, n);
            let sign = dxds[speceq].signum();
            ds = sign * ds.abs();
            dxds = sensitivity(mixture, &u, beta, z, speceq)?;
            if last_iters > 6 {
                ds *= 0.5;
            } else if last_iters < 3 {
                ds *= 1.1;
            } else if last_iters == 4 {
                ds *= 0.9;
            } else if last_iters > 4 {
                ds *= 0.7;
            }
            ds = ds.signum() * ds.abs();
            // Clamp to the temperature and pressure steps.
            let t_cur = u[n].exp();
            let p_cur = u[n + 1].exp();
            if dxds[n].abs() * ds.abs() > (1.0 + D_TMAX / t_cur).ln() {
                ds = ds.signum() * (1.0 + D_TMAX / t_cur).ln() / dxds[n].abs();
            }
            if dxds[n + 1].abs() * ds.abs() > (1.0 + D_PMAX / p_cur).ln() {
                ds = ds.signum() * (1.0 + D_PMAX / p_cur).ln() / dxds[n + 1].abs();
            }
            if history.len() == HISTORY_LEN {
                let sny = ds * dxds[speceq] + history[HISTORY_LEN - 1][speceq];
                u = cubic_predict(&history, speceq, sny);
                specval = sny;
            } else {
                for (i, &d) in dxds.iter().enumerate() {
                    u[i] += d * ds;
                }
                specval = u[speceq];
            }
        }
        clamp_state(&mut u, n);

        let predicted = u.clone();
        let mut solved_point = false;
        for _ in 0..8 {
            match newton(
                mixture,
                &u,
                beta,
                z,
                LastRow::Pinned {
                    index: speceq,
                    value: specval,
                },
            ) {
                Ok((solved, iters)) => {
                    u = solved;
                    last_iters = iters;
                    solved_point = true;
                    break;
                }
                Err(_) => {
                    ds *= 0.5;
                    let prev = history.last().cloned().unwrap_or(predicted.clone());
                    for (i, &d) in dxds.iter().enumerate() {
                        u[i] = prev[i] + d * ds;
                    }
                    specval = u[speceq];
                }
            }
        }
        if !solved_point {
            break;
        }
    }

    Ok(BranchTrace {
        t: t_vec,
        p: p_vec,
        cricondenbar,
        cricondentherm,
        critical,
    })
}

/// The PT phase envelope of a mixture of composition `z`, traced upward from a low
/// pressure to the critical point on both the bubble and dew branches.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the mixture has one component, or `z` is not a
///   composition.
/// * [`AzothError::OutOfRange`] if `P` is not positive.
/// * [`AzothError::SolverNotConverged`] if a branch's continuation hits its cap.
pub fn pt_phase_envelope(
    mixture: &Mixture,
    p: Pressure,
    z: &[f64],
) -> Result<PtPhaseEnvelopeResult> {
    let spec = &model_gen::PT_PHASE_ENVELOPE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| (quantity == "P").then_some(p.value),
        &mut warnings,
    )?;

    let n = mixture.len();
    if n < 2 {
        return Err(AzothError::invalid_input(
            "components",
            "a phase envelope needs two phases with different compositions, and one \
             component cannot have them. A pure component's envelope is its \
             vapour-pressure curve, which `eos.pure_saturation` computes.",
        ));
    }
    phase_boundary::check_composition(z, n, "z")?;

    let bubble = trace_branch(mixture, p, z, BUBBLE_BETA)?;
    let dew = trace_branch(mixture, p, z, DEW_BETA)?;

    let cricondenbar = if bubble.cricondenbar.1 >= dew.cricondenbar.1 {
        bubble.cricondenbar
    } else {
        dew.cricondenbar
    };
    let cricondentherm = if bubble.cricondentherm.0 >= dew.cricondentherm.0 {
        bubble.cricondentherm
    } else {
        dew.cricondentherm
    };

    // The critical point is where the two branches meet: take the last point of the
    // bubble branch, which stops at the K-value crossing.
    // The critical point is the refined one when either branch reached it, and either will do:
    // they are the same point and each branch computed it alone, so a disagreement is a finding
    // rather than a choice. Falling back to the branch endpoint keeps a trace that never got
    // near criticality reporting where it stopped rather than a `NaN`.
    let critical = bubble
        .critical
        .or(dew.critical)
        .or_else(|| {
            bubble
                .t
                .last()
                .zip(bubble.p.last())
                .map(|(&t, &pr)| (t, pr))
        })
        .unwrap_or((f64::NAN, f64::NAN));
    // The residual is the two branches' disagreement about the critical temperature: zero when
    // both refined to the same point, which is the statement that they meet. A branch that did
    // not refine reports the gap between where the branches stopped - the older, weaker claim.
    let residual = match (bubble.critical, dew.critical) {
        (Some((t_b, _)), Some((t_d, _))) => (t_b - t_d).abs(),
        _ => bubble
            .t
            .last()
            .zip(dew.t.last())
            .map(|(&t_b, &t_d)| (t_b - t_d).abs())
            .unwrap_or(f64::NAN),
    };

    let iterations = (bubble.t.len() + dew.t.len()) as u32;
    Ok(PtPhaseEnvelopeResult {
        dew_temperature: dew.t,
        dew_pressure: dew.p,
        bubble_temperature: bubble.t,
        bubble_pressure: bubble.p,
        cricondenbar_temperature: kelvins(cricondenbar.0),
        cricondenbar_pressure: pascals(cricondenbar.1),
        cricondentherm_temperature: kelvins(cricondentherm.0),
        cricondentherm_pressure: pascals(cricondentherm.1),
        critical_temperature: kelvins(critical.0),
        critical_pressure: pascals(critical.1),
        iterations,
        residual,
        warnings,
    })
}
