//! `eos.pt_phase_envelope` - the Michelsen natural-parameter phase envelope.
//!
//! Spec: `specs/models/eos/pt_phase_envelope.toml`
//!
//! The PT phase envelope is the locus where a mixture of composition `z` is at a phase
//! boundary. It is traced as two independent natural-parameter continuations - the
//! bubble branch at a tiny vapour fraction and the dew branch at one minus it - each
//! bootstrapped at a low pressure where its K-values are far from one, then continued
//! upward until the two phases have nearly collapsed onto one another. The `N+2`
//! unknowns `(ln K_i, ln T, ln P)` satisfy `N` isofugacity equations, one material
//! balance, and one specification equation that pins the most sensitive variable so the
//! trace can turn the corner at the cricondenbar and cricondentherm.
//!
//! The next point is predicted by a cubic through the last four converged points and
//! then polished by a damped Newton step. **The Jacobian is analytic**, from
//! [`Mixture::phase_derivatives`] - the surface P6 item 2 ported, and the reason that item
//! exists. It was a central difference until this pass, which cost `N+3` state evaluations per
//! Newton step, each of them two cubic solves; NeqSim's own `setJac` is analytic too, so the
//! difference was a divergence as well as an expense.
//!
//! **The critical point is computed, not refined out of the trace.** At the critical point of
//! the isopleth the two phases have the feed's composition, which is exactly the problem
//! [`crate::critical_point`] solves by the Heidemann-Khalil construction; the other candidate,
//! Newtoning `sum_i (ln K_i)^2 = 0` beside the isofugacity rows from a state the trace has
//! brought close, does not work here and was measured not to: with every `K_i = 1` the two
//! phases have the same composition, so the isofugacity rows are `ln phi^V_i - ln phi^L_i`,
//! which is zero wherever the cubic has one root, and the system is therefore satisfied on a
//! two-parameter *family* rather than at a point. Measured on methane/n-butane 50/50 at three
//! separate states, all three residuals were exactly zero and the Newton never moved. Where it
//! does move it is not converging: across five compositions it returned `42.74 K` on one branch
//! against Heidemann-Khalil's `275.62 K`.
//!
//! **`RobustPhaseEnvelope` is named and not ported.** NeqSim carries a second envelope
//! driver under that name, 517 lines in `phaseenvelopeops/`, and **nothing constructs
//! it** - a repo-wide grep for `new RobustPhaseEnvelope(` returns zero sites. The
//! package's reachable drivers are the ones already dispositioned: this trace's source
//! `PTPhaseEnvelopeMichelsen` for the envelope itself, and `SysNewtonRhapsonPhaseEnvelope`
//! for its points.

use azoth_core::units::{Pressure, kelvins, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::mixture::{Mixture, RootSide, normalise};
use crate::model_gen;
use crate::phase_boundary::{self, WILSON_CONSTANT};
use crate::results::PtPhaseEnvelopeResult;

/// The vapour fraction on the bubble branch: `x_i ~= z_i`, so the feed is liquid.
///
/// Upstream's own value - `ThermodynamicOperations.calcPTphaseEnvelope` passes `phasefraction =
/// 1e-10` - and it was `1e-6` here until this pass. Measured, the difference is nil: the two
/// give the same envelope to every reported digit, because the trace's state is fixed by the
/// isofugacity equations rather than by how far from the boundary it is placed.
const BUBBLE_BETA: f64 = 1.0e-10;
/// The vapour fraction on the dew branch: `y_i ~= z_i`, so the feed is vapour.
const DEW_BETA: f64 = 1.0 - 1.0e-10;
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

/// A copy of a vector scaled to sum to one.
fn normalised(values: &[f64]) -> Vec<f64> {
    let mut out = values.to_vec();
    normalise(&mut out);
    out
}

/// The two phase compositions at a vapour fraction, from `K_i = y_i / x_i`.
///
/// **The vectors returned are mole *numbers*, not compositions**, and their sums are one only
/// when the material balance already holds - which is what the residual's last-but-one row
/// tests. That is deliberate: the row *is* their difference, the Rachford-Rice residual, and it
/// is identically zero if they are normalised. Callers that want a composition normalise.
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

/// What the system's last row pins: a state variable held at a value.
///
/// The continuation's specification equation. There was a second variant here, the criticality
/// row `sum_i (ln K_i)^2 = 0`, until the critical point stopped being refined out of the trace;
/// see the module note for why that refinement is not a computation.
#[derive(Debug, Clone, Copy)]
enum LastRow {
    /// A state variable pinned to a value: row `u[index] - value`.
    Pinned { index: usize, value: f64 },
}

/// The residuals, and on request their analytic Jacobian.
///
/// `u` is `[ln K_i; ln T; ln P]`; `last` is the system's final row.
struct Point {
    /// The `N+2` residuals.
    f: Vec<f64>,
    /// The `N+2` by `N+2` Jacobian, when it was asked for.
    jac: Option<Vec<Vec<f64>>>,
}

/// The residuals at a state, and their analytic Jacobian.
///
/// **The Jacobian is analytic**, from [`Mixture::phase_derivatives`] - the surface P6 item 2
/// ported. NeqSim's `setJac` assembles the same rows from `getdfugdx`, `getdfugdt` and
/// `getdfugdp`, and there is no finite difference anywhere in `SysNewtonRhapsonPhaseEnvelope`.
/// This replaces a central difference that cost `N+3` state evaluations per Newton step, each
/// of them two cubic solves, with one.
///
/// With `D_i = 1 - beta + beta K_i`, `x_i = z_i/D_i` and `y_i = K_i x_i`, and the compositions
/// **normalised** because that is the space `d ln phi_i/dn_j` is a derivative on:
///
/// ```text
/// d x_i / d ln K_m = x_i (beta K_m / D_m) (x_m - delta_im)
/// d y_i / d ln K_m = y_i [ (1 - beta K_i/D_i) delta_im - y_m (1 - beta K_m/D_m) ]
/// ```
///
/// Both sum to zero over `i`, which is the normalisation surviving the derivative: a row of the
/// Jacobian cannot move mass off the simplex. **Upstream's `setJac` omits this term**, contracting
/// `dfugdx` with the derivative of the *raw* vector. Ours has to be the exact one, because the
/// central difference it replaced differentiates the residual as written and the two must agree.
///
/// The material-balance row is the exception and uses the **raw** derivatives, because that row
/// differentiates the raw sums - using the simplex form there makes it identically zero and the
/// Rachford-Rice content of the system disappears.
fn point(
    mixture: &Mixture,
    u: &[f64],
    beta: f64,
    z: &[f64],
    last: LastRow,
    want_jacobian: bool,
) -> Result<Point> {
    let n = z.len();
    let k: Vec<f64> = u[..n].iter().map(|&lnk| lnk.exp()).collect();
    let t = u[n].exp();
    let p = u[n + 1].exp();
    let (x_raw, y_raw) = compositions(z, &k, beta);
    // **The states are evaluated at normalised compositions.** `K_i = y_i/x_i` fixes the
    // *ratio*, so `z_i/(1 - beta + beta K_i)` is a mole-number vector whose sum is one only
    // when the material balance happens to hold; at any other iterate it is not a composition
    // at all, and `phase_state` builds `a_mix` and `b_mix` from the vector it is given - which
    // the cubic is not invariant under. The raw sums are kept below because the material
    // balance row *is* their difference.
    let x = normalised(&x_raw);
    let y = normalised(&y_raw);
    let reduced = mixture.reduced_parameters(kelvins(t), pascals(p))?;
    let liquid = mixture.phase_state(&reduced, &x, RootSide::Liquid)?;
    let vapour = mixture.phase_state(&reduced, &y, RootSide::Vapour)?;

    let mut f = vec![0.0; n + 2];
    for i in 0..n {
        // ln K_i = ln phi_i^L - ln phi_i^V.
        f[i] = u[i] - liquid.ln_phi[i] + vapour.ln_phi[i];
    }
    f[n] = y_raw.iter().sum::<f64>() - x_raw.iter().sum::<f64>();
    let LastRow::Pinned { index, value } = last;
    f[n + 1] = u[index] - value;

    if !want_jacobian {
        return Ok(Point { f, jac: None });
    }

    let factor = |i: usize| beta * k[i] / (1.0 - beta + beta * k[i]);
    let m = n + 2;
    let mut jac = vec![vec![0.0; m]; m];
    let mut dx = vec![vec![0.0; n]; n];
    let mut dy = vec![vec![0.0; n]; n];
    for i in 0..n {
        for mm in 0..n {
            let diagonal = if i == mm { 1.0 } else { 0.0 };
            dx[i][mm] = x[i] * factor(mm) * (x[mm] - diagonal);
            dy[i][mm] = y[i] * ((1.0 - factor(i)) * diagonal - y[mm] * (1.0 - factor(mm)));
        }
    }

    let liquid_derivatives = mixture.phase_derivatives(&reduced, &x, liquid.z)?;
    let vapour_derivatives = mixture.phase_derivatives(&reduced, &y, vapour.z)?;

    for (i, row) in jac.iter_mut().enumerate().take(n) {
        for (j, value) in row.iter_mut().enumerate().take(n) {
            let diagonal = if i == j { 1.0 } else { 0.0 };
            let d_vapour: f64 = (0..n)
                .map(|mm| vapour_derivatives.d_ln_phi_dn[i][mm] * dy[mm][j])
                .sum();
            let d_liquid: f64 = (0..n)
                .map(|mm| liquid_derivatives.d_ln_phi_dn[i][mm] * dx[mm][j])
                .sum();
            *value = diagonal + d_vapour - d_liquid;
        }
        // The chain rule for the `ln T` and `ln P` unknowns: a derivative in `T` is `T` times a
        // derivative in `ln T`.
        row[n] = t * (vapour_derivatives.d_ln_phi_dt[i] - liquid_derivatives.d_ln_phi_dt[i]);
        row[n + 1] = p * (vapour_derivatives.d_ln_phi_dp[i] - liquid_derivatives.d_ln_phi_dp[i]);
    }
    // The material-balance row, in the raw derivatives: `X_i = z_i/D_i` and `Y_i = K_i X_i`, both
    // functions of `K_i` alone.
    for j in 0..n {
        jac[n][j] = y_raw[j] * (1.0 - factor(j)) + x_raw[j] * factor(j);
    }
    // The specification row: the pinned variable, and nothing else.
    let LastRow::Pinned { index, .. } = last;
    jac[n + 1][index] = 1.0;

    Ok(Point { f, jac: Some(jac) })
}

/// The residuals alone, for the callers that do not need the Jacobian.
fn residual(mixture: &Mixture, u: &[f64], beta: f64, z: &[f64], last: LastRow) -> Result<Vec<f64>> {
    Ok(point(mixture, u, beta, z, last, false)?.f)
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
        let Point { f, jac } = point(mixture, &u, beta, z, last, true)?;
        let norm = norm2(&f);
        if norm < NEWTON_TOL {
            return Ok((u, iters));
        }
        let mut jac = jac.expect("the Jacobian was asked for");
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
    let mut jac = point(
        mixture,
        u,
        beta,
        z,
        LastRow::Pinned {
            index: speceq,
            value: u[speceq],
        },
        true,
    )?
    .jac
    .expect("the Jacobian was asked for");
    let mut e = vec![0.0; n + 2];
    e[n + 1] = 1.0;
    solve_linear(&mut jac, &mut e).ok_or(AzothError::SolverNotConverged {
        iterations: 0,
        residual: f64::NAN,
        tolerance: NEWTON_TOL,
    })
}

/// The index of the variable the continuation pins: whichever of the `N+2` is most sensitive.
///
/// **The K-values are candidates**, and excluding them was wrong. They were excluded as "not
/// monotonic along the envelope, and pinning one makes the trace oscillate" - which is true of
/// the trace on its own, and is why this cannot land without the acceptance-and-retry in
/// `trace_branch`. Measured, widening the choice alone sends the dew branch to its 9,980-point
/// cap and pulls the cricondentherm from 384.0182 K to 376.4847 K.
///
/// Upstream's `findSpecEq` maximises `|dxds|` over all `N+2` and takes the winner outright. It
/// also computes a second sensitivity - the relative change in `exp(u_i)` since the last point -
/// and uses that one only when the two agree, so the maximum of `|dxds|` is the decider. That
/// second measure is not computed here: on the validated case the two agreed at every point, and
/// a quantity that never decides is a quantity nothing can check.
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

/// **The critical point of the envelope, computed rather than refined out of the trace.**
///
/// The critical point of the isopleth at `z` is the state where the two phases become one, and
/// there both have the feed's composition - which is exactly the mixture critical point that
/// [`crate::critical_point`] solves. The alternative, which this replaced, was to Newton
/// `sum_i (ln K_i)^2 = 0` beside the isofugacity rows from wherever the trace came closest, and
/// it does not compute anything: at every `K_i = 1` the two phases have the same composition, so
/// the isofugacity rows reduce to `ln phi^V_i - ln phi^L_i`, which vanishes wherever the cubic
/// has a single root. The system is satisfied on a two-parameter *family* near criticality, not
/// at a point; measured at three separate states, all three residuals were exactly zero and the
/// Newton returned its seed. Where the seed is off that family the iteration does move, and then
/// it is unreliable: measured across five compositions of methane/n-butane it returned `42.74 K`
/// on one branch, against `275.62 K` from the construction below.
///
/// Returns `None` if the construction does not converge, which the caller reports as the point
/// where the trace stopped rather than as an error.
fn critical(mixture: &Mixture, z: &[f64]) -> Option<(f64, f64)> {
    crate::critical_point::critical_point(mixture, z)
        .ok()
        .map(|point| (point.tc.value, point.pc.value))
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
    /// Where this branch stopped, when it ran into the critical region rather than into the
    /// pressure ceiling or a failure to converge.
    ///
    /// The two branches approach the critical point from opposite sides, so the two stopping
    /// points bracket it and their gap is what `residual` reports.
    stopped: Option<(f64, f64)>,
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
    let mut stopped: Option<(f64, f64)> = None;
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

        // The branch has run into the critical region: every K-value is approaching one. The
        // state is recorded as where this branch *stopped* - it is a converged point of the
        // continuation and the last thing this branch has to say about the critical point. What
        // the critical point *is* comes from `critical`, and not from here.
        if k[lc] < 1.05 && k[hc] > 0.95 {
            stopped = Some((t, pv));
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
        stopped,
    })
}

/// The PT phase envelope of a mixture of composition `z`, traced upward from a low
/// pressure towards the critical point on both the bubble and dew branches.
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

    // **The critical point is computed, not read off the trace**, and it is one point: the two
    // phases become one there, so a branch cannot have its own. It is computed only when a branch
    // actually ran into the critical region; a trace that never got near it reports where it
    // stopped rather than a critical point it never found.
    let critical = if bubble.stopped.is_some() || dew.stopped.is_some() {
        critical(mixture, z)
    } else {
        None
    }
    .or(bubble.stopped)
    .or(dew.stopped)
    .or_else(|| {
        bubble
            .t
            .last()
            .zip(bubble.p.last())
            .map(|(&t, &pr)| (t, pr))
    })
    .unwrap_or((f64::NAN, f64::NAN));

    // The residual is how far apart the two branches stopped. They approach the critical point
    // from opposite sides, so the gap is the width of the bracket they put around it: the trace's
    // own statement about whether it reached the critical point, and one the computed critical
    // point above cannot make for it. A branch that ran into the pressure ceiling instead of the
    // critical region has not met the other, and there is no gap to report.
    let residual = match (bubble.stopped, dew.stopped) {
        (Some((t_b, _)), Some((t_d, _))) => (t_b - t_d).abs(),
        _ => f64::NAN,
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
