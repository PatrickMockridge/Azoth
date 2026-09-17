//! The second-order step the slowly-converging flashes hand over to.
//!
//! NeqSim's `SysNewtonRhapsonTPflash`: a Newton step on the isofugacity residuals in
//! the *vapour mole numbers* `u_i = beta y_i`, which is the variable set that makes
//! the residuals' Jacobian the Hessian of Michelsen's reduced Gibbs energy
//!
//! ```text
//! Q(u) = sum_i u_i ln(y_i phi_i^V) + (z_i - u_i) ln(x_i phi_i^L)
//! ```
//!
//! and the step a descent direction for it. That is why the line search below is on
//! `Q` rather than on the residual: the residual's gradient *is* `Q`'s, so `Q` is the
//! function whose minimisation the equations describe, and a step that does not
//! decrease it is a step away from the answer.
//!
//! Unlike [`crate::flash_iteration`], which is a scheme both implementations run for
//! their whole solve, this runs one step per outer iteration: the caller alternates,
//! doing its own scheme until it stops making progress and this from then on. See
//! `eos.pt_flash`'s `algorithm.fallback`.
//!
//! Every derivative here is analytic, from
//! [`crate::mixture::Mixture::phase_derivatives`]. Nothing is differenced.

use azoth_core::{AzothError, ModelAlgorithm, Result};

/// The step's own linear solve failing is the one thing this cannot proceed past.
use crate::mixture::{Mixture, ReducedParameters, RootSide};
use crate::results::Phase;

/// One Newton step's outcome.
pub(crate) struct Step {
    /// The vapour mole numbers the step lands on.
    pub u: Vec<f64>,
    /// `max_i |g_i|` at the iterate the step landed on.
    ///
    /// The **stopping measure**, and it is the residual rather than NeqSim's
    /// `deviation`. NeqSim stops on the relative step, which is a different quantity
    /// from the one the scheme it replaced stops on - an rms change in `ln K` - so the
    /// two answers stop at different distances from the same root. Measured: that put
    /// the answer 5.2e-10 away from NeqSim's own on the SRK flash case, against a case
    /// tolerance of 1e-10, where the outer scheme alone reproduced it to 1e-12. The
    /// residual is what the equations actually say, so stopping on it makes the two
    /// schemes stop in the same place.
    pub residual: f64,
}

/// The two phases at a set of vapour mole numbers.
fn split(feed: &[f64], u: &[f64]) -> (f64, Vec<f64>, Vec<f64>) {
    let beta: f64 = u.iter().sum();
    let y: Vec<f64> = u.iter().map(|value| value / beta).collect();
    let x: Vec<f64> = feed
        .iter()
        .zip(u)
        .map(|(&zi, &ui)| (zi - ui) / (1.0 - beta))
        .collect();
    (beta, x, y)
}

/// Michelsen's reduced Gibbs energy at a set of vapour mole numbers.
///
/// Its gradient in `u` is the residual vector, which is what makes it the right thing
/// to take a line search on: a step that decreases `Q` decreases every residual.
fn gibbs_energy(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    feed: &[f64],
    u: &[f64],
) -> Result<f64> {
    let (_, x, y) = split(feed, u);
    let liquid = mixture.phase_state(reduced, &x, RootSide::Liquid)?;
    let vapour = mixture.phase_state(reduced, &y, RootSide::Vapour)?;
    let mut total = 0.0;
    for i in 0..u.len() {
        total += u[i] * (vapour.ln_phi[i] + y[i].ln());
        total += (feed[i] - u[i]) * (liquid.ln_phi[i] + x[i].ln());
    }
    Ok(total)
}

/// Whether a trial `u` describes two phases at all.
///
/// NeqSim's `isFeasible`: every `u_i` strictly between zero and `z_i`, and `beta`
/// strictly inside `(0, 1)`. A step that leaves this region has asked for a negative
/// mole number, and the compositions below would take its logarithm.
fn is_feasible(feed: &[f64], u: &[f64]) -> bool {
    const EDGE: f64 = 1.0e-15;
    if u.iter()
        .zip(feed)
        .any(|(&ui, &zi)| ui < EDGE || ui > zi - EDGE)
    {
        return false;
    }
    let beta: f64 = u.iter().sum();
    beta > 1.0e-12 && beta < 1.0 - 1.0e-12
}

/// One Newton step, from an iterate the caller's own scheme produced.
///
/// # Errors
/// * [`AzothError::SolverNotConverged`] if the linear system is singular at this
///   iterate, which is a state the step cannot improve rather than one to retry.
/// * Propagates the phase models' range checks, including a cubic that has no root
///   at a trial point.
pub(crate) fn step(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    feed: &[f64],
    u: &[f64],
    algorithm: &ModelAlgorithm,
) -> Result<Step> {
    let n = feed.len();
    let (_, x, y) = split(feed, u);
    let liquid = mixture.phase_state(reduced, &x, RootSide::Liquid)?;
    let vapour = mixture.phase_state(reduced, &y, RootSide::Vapour)?;
    let d_vapour = mixture.phase_derivatives(reduced, &y, vapour.z)?;
    let d_liquid = mixture.phase_derivatives(reduced, &x, liquid.z)?;
    let beta: f64 = u.iter().sum();

    // `g_i = ln(y_i phi_i^V) - ln(x_i phi_i^L)`, the isofugacity residual.
    let mut residual = vec![0.0; n];
    for i in 0..n {
        residual[i] = vapour.ln_phi[i] + y[i].ln() - liquid.ln_phi[i] - x[i].ln();
    }
    let current = gibbs_energy(mixture, reduced, feed, u)?;

    // NeqSim's Jacobian, and it is the same matrix the derivation gives. The two
    // `- 1` terms are the normalisation of each phase's fractions and the two
    // composition derivatives are `d ln phi_i / d n_j` at unit total, which is why
    // each is divided by its own phase's mole number.
    let mut jacobian = vec![vec![0.0; n]; n];
    for (i, row) in jacobian.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            let dij = if i == j { 1.0 } else { 0.0 };
            *entry = (d_vapour.d_ln_phi_dn[i][j] + dij / y[i] - 1.0) / beta
                + (d_liquid.d_ln_phi_dn[i][j] + dij / x[i] - 1.0) / (1.0 - beta);
        }
    }
    // Levenberg-Marquardt regularisation, NeqSim's and at NeqSim's magnitude: the
    // diagonal is made strictly dominant enough to keep the solve away from a
    // singular matrix at a point where the two phases have nearly met.
    let trace: f64 = (0..n).map(|i| jacobian[i][i].abs()).sum();
    let lambda = 1.0e-8 * trace / n as f64;
    for (i, row) in jacobian.iter_mut().enumerate() {
        row[i] += lambda;
    }

    let delta = solve(&jacobian, &residual, algorithm)?;
    // `slope = g . dx`, positive because `dx` is `J^-1 g` and `J` is `Q`'s Hessian.
    let slope: f64 = (0..n).map(|i| residual[i] * delta[i]).sum();

    // Armijo backtracking on `Q`, NeqSim's `c1 = 1e-4` and eight halvings. A trial the
    // cubic cannot even evaluate - a root that does not exist at that composition - is
    // treated as a failure to decrease `Q` and halved through.
    //
    // **The halved step is taken even when none of the eight satisfies Armijo**, which
    // is NeqSim's behaviour and not an oversight in it: `alpha` has been divided down
    // to `1/256` by then, so the step is small enough that its direction matters more
    // than its descent, and refusing it outright stops the solve over a line search
    // that has already done its job. Measured: refusing it costs the states that
    // converge to the *trivial* solution - where the Jacobian is degenerate and no
    // step descends - their convergence, because they spend the outer cap on refused
    // attempts that NeqSim would have spent inching forward.
    let mut alpha = 1.0;
    for _ in 0..8 {
        let trial: Vec<f64> = (0..n).map(|i| u[i] - alpha * delta[i]).collect();
        if is_feasible(feed, &trial) {
            if let Ok(trial_energy) = gibbs_energy(mixture, reduced, feed, &trial) {
                if trial_energy <= current - 1.0e-4 * alpha * slope {
                    break;
                }
            }
        }
        alpha *= 0.5;
    }

    let stepped: Vec<f64> = (0..n).map(|i| u[i] - alpha * delta[i]).collect();

    // The residual at the iterate the step landed on, which is what the solve stops
    // on. One more phase pair per step, and it buys the two schemes stopping in the
    // same place rather than at the same relative step.
    let (_, landed_x, landed_y) = split(feed, &stepped);
    let landed_liquid = mixture.phase_state(reduced, &landed_x, RootSide::Liquid)?;
    let landed_vapour = mixture.phase_state(reduced, &landed_y, RootSide::Vapour)?;
    let landed = (0..n)
        .map(|i| {
            (landed_vapour.ln_phi[i] + landed_y[i].ln()
                - landed_liquid.ln_phi[i]
                - landed_x[i].ln())
            .abs()
        })
        .fold(0.0, f64::max);

    Ok(Step {
        u: stepped,
        residual: landed,
    })
}

/// The `u` a converged flash hands over at, from the phase pair.
pub(crate) fn vapour_moles(beta: f64, y: &[f64]) -> Vec<f64> {
    y.iter().map(|value| beta * value).collect()
}

/// The two phases at a set of vapour mole numbers, for a caller reading the answer back.
pub(crate) fn split_at(feed: &[f64], u: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let (_, x, y) = split(feed, u);
    (x, y)
}

/// The phase a converged `u` describes.
pub(crate) fn phase_of(beta: f64) -> Phase {
    if beta < 0.0 {
        Phase::AllLiquid
    } else if beta > 1.0 {
        Phase::AllVapour
    } else {
        Phase::TwoPhase
    }
}

/// `J dx = g` by Gaussian elimination with partial pivoting.
///
/// Written out rather than pulled in: the system is `N x N` with `N` the component
/// count, the crate carries no linear algebra and adding one for a two-by-two would
/// be the larger dependency. NeqSim uses an LU decomposition; this is the same
/// arithmetic without the factorisation.
fn solve(j: &[Vec<f64>], g: &[f64], algorithm: &ModelAlgorithm) -> Result<Vec<f64>> {
    let n = g.len();
    let mut a: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut row = j[i].clone();
            row.push(g[i]);
            row
        })
        .collect();

    for column in 0..n {
        let pivot = (column..n)
            .max_by(|&p, &q| {
                a[p][column]
                    .abs()
                    .partial_cmp(&a[q][column].abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(column);
        a.swap(column, pivot);
        if a[column][column] == 0.0 {
            return Err(AzothError::SolverNotConverged {
                iterations: 0,
                residual: f64::NAN,
                tolerance: algorithm.tolerance,
            });
        }
        let pivot_value = a[column][column];
        let (upper, lower) = a.split_at_mut(column + 1);
        let pivot_row = &upper[column];
        for row in lower {
            let factor = row[column] / pivot_value;
            for (entry, &pivot_entry) in row.iter_mut().zip(pivot_row.iter()).skip(column) {
                *entry -= factor * pivot_entry;
            }
        }
    }

    // The condition number is not estimated: the step is a *trial*, and the line
    // search is what decides whether to take it, so a poorly-conditioned solve shows
    // up as a step that fails to decrease `Q` rather than as a wrong answer.
    let mut dx = vec![0.0; n];
    for row in (0..n).rev() {
        let sum: f64 = ((row + 1)..n).map(|k| a[row][k] * dx[k]).sum();
        dx[row] = (a[row][n] - sum) / a[row][row];
    }
    Ok(dx)
}
