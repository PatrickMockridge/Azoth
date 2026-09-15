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

/// The `N+2` residuals: `[isofugacity; material balance; specification]`.
///
/// `u` is `[ln K_i; ln T; ln P]`; the specification is `u[speceq] - specval`.
fn residual(
    mixture: &Mixture,
    u: &[f64],
    beta: f64,
    z: &[f64],
    speceq: usize,
    specval: f64,
) -> Result<Vec<f64>> {
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
    f[n + 1] = u[speceq] - specval;
    Ok(f)
}

/// The `N+2` by `N+2` Jacobian, by central difference over the state variables.
fn jacobian(
    mixture: &Mixture,
    u: &[f64],
    beta: f64,
    z: &[f64],
    speceq: usize,
    specval: f64,
) -> Result<Vec<Vec<f64>>> {
    let m = u.len();
    let f0 = residual(mixture, u, beta, z, speceq, specval)?;
    let mut jac = vec![vec![0.0; m]; m];
    for j in 0..m {
        let delta = (1.0e-6 * u[j].abs()).max(1.0e-8);
        let mut up = u.to_vec();
        up[j] += delta;
        let fp = residual(mixture, &up, beta, z, speceq, specval)?;
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
    speceq: usize,
    specval: f64,
) -> Result<(Vec<f64>, u32)> {
    let mut u = u.to_vec();
    let mut iters = 0;
    for _ in 0..NEWTON_MAX {
        iters += 1;
        let f = residual(mixture, &u, beta, z, speceq, specval)?;
        let norm = norm2(&f);
        if norm < NEWTON_TOL {
            return Ok((u, iters));
        }
        let mut jac = jacobian(mixture, &u, beta, z, speceq, specval)?;
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
            let fnext = residual(mixture, &next, beta, z, speceq, specval)?;
            if norm2(&fnext) < norm {
                break;
            }
            step *= 0.5;
        }
        u = next;
    }
    Err(AzothError::SolverNotConverged {
        iterations: iters,
        residual: norm2(&residual(mixture, &u, beta, z, speceq, specval)?),
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
    let mut jac = jacobian(mixture, u, beta, z, speceq, u[speceq])?;
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
    let mut history: Vec<Vec<f64>> = Vec::new();
    let mut dxds: Vec<f64> = vec![0.0; n + 2];
    let mut ds = 0.1;
    let mut last_iters = 2;

    let (converged, _) = newton(mixture, &u, beta, z, speceq, specval)?;
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

        // The critical point: every K-value approaches one. Stop the branch here.
        if k[lc] < 1.05 && k[hc] > 0.95 {
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
            match newton(mixture, &u, beta, z, speceq, specval) {
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
    let critical = bubble
        .t
        .last()
        .zip(bubble.p.last())
        .map(|(&t, &pr)| (t, pr))
        .unwrap_or((f64::NAN, f64::NAN));
    // The two branches stop at the K-value crossing but not exactly at the same point;
    // the residual reports the temperature gap between their endpoints.
    let residual = bubble
        .t
        .last()
        .zip(dew.t.last())
        .map(|(&t_b, &t_d)| (t_b - t_d).abs())
        .unwrap_or(f64::NAN);

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
