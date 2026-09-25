//! Naphtali-Sandholm: the column's simultaneous MESH correction.
//!
//! The second solve. Where [`crate::kernels::distillation_column`]'s substitution core moves one
//! tray at a time and iterates the profile, this linearises the **whole column at once**: one
//! Newton step on the `N * (C + 2)` variables, with the block-tridiagonal structure of the MESH
//! Jacobian exploited by [`super::block_tridiagonal`].
//!
//! **The method is the source's own citation**: `NaphtaliSandholmSolver`'s javadoc names
//! Naphtali & Sandholm (1971), "Multicomponent Separation Calculations by Linearization",
//! AIChE Journal 17(1) 148-153, and the class is a faithful implementation of it -
//! `N` tray blocks of `C + 2` variables (`C` liquid component flows, the temperature, the
//! vapour flow), `C` component material balances, the energy balance, and the summation
//! equation `sum(K x) - 1`, in that order.
//!
//! # What is the class's, and what is this port's
//!
//! Everything numerical is the class's, quoted rather than improved: the residual and its
//! scalings, the finite-difference Jacobian (`h = max(|x| * 1e-4, 1e-8)`, one perturbed
//! variable at a time, only the perturbed tray's thermodynamics re-evaluated), the
//! block-tridiagonal solve, the trust region (`|dT| <= 10 K`, `|dV| <= 0.5 V + 0.05 flowScale`,
//! `|dliq| <= 0.5 liq + 1e-3 flowScale`), the Armijo backtracking line search (`c = 1e-4`,
//! `rho = 0.5`, 15 backtracks), and the acceptance gates.
//!
//! **The flows are mol/hr inside the solve, and mol/s outside it.** The class works in mol/hr,
//! and three of its constants are absolute rather than relative: `flowScale` itself (the
//! residual's divisor, `totalFeed / N`), the trust region's floors (`0.05 * flowScale` on the
//! vapour step and `1e-3 * flowScale` on a component step) and the perturbations' own floor
//! (`1e-8`). A solve scaled by one library's flow unit and stopped by the other's tolerance is
//! a different solve, so the port converts once at each boundary rather than carrying a
//! conversion through every arithmetic line - which is what `flow_scale`'s own comment says.
//!
//! **The seed is the class's warm start**, `initializeTrayStateFromColumn`: the MESH variables
//! are mapped from the column's own converged trays. NeqSim takes that path when the column has
//! been solved before; a port that owns both solves can always take it.
//!
//! **The class's cold seed is not ported, and `initialize_from_column` records the
//! measurement**: its per-component shoot leaves the top tray's light components on the
//! overhead and `runBostonSullivanRefinement` is the 464-line step that repairs it, without
//! which the first Newton step is `1e11` times the variables it moves.
//!
//! # What the port cannot be reached with
//!
//! * **A reboiler with no temperature pin** runs `solveBubblePointMethod` in the class, whose
//!   `V`-cascade is the basin-finding this port does not carry.
//! * **A product specification** re-runs the whole column inside
//!   `solveWithSpecificationTargets`, so the class's NS-with-a-specification is a second
//!   integration; this port's outer specification loop drives the substitution solve.

// The equations and the elimination are the class's own index arithmetic: the MESH residual is
// written over tray and component indices and the linear algebra over rows and columns, so the
// loops index rather than iterate. Renaming them to iterators would rewrite the port rather
// than clean it up - the precedent `azoth-eos`'s `gerg2008` sets.
#![allow(clippy::needless_range_loop)]

use azoth_core::units::{Pressure, kelvins, pascals};
use azoth_core::{AzothError, Result};
use azoth_eos::{RootSide, molar_enthalpy_entropy};

use super::block_tridiagonal::BlockTridiagonal;
use crate::kernels::distillation_column::{ColumnOutcome, ColumnSetup, SolverType, TrayProfile};
use crate::stream::Stream;

/// The class's own flow unit, which every absolute constant of this solver is written in.
const MOL_PER_HOUR: f64 = 3600.0;

/// Maximum non-descent line-search steps before the best state is restored.
const MAX_NON_DESCENT_LINE_SEARCH_STEPS: usize = 3;
/// Forced-root fugacity fixed-point sweeps for one tray evaluation.
const THERMO_K_VALUE_ITERATIONS: usize = 2;
/// Maximum forced-root fugacity sweeps when refining a retained column state.
#[allow(dead_code)] // the class's warm-start seed is not ported; the constant names its budget
const THERMO_WARM_START_K_VALUE_ITERATIONS: usize = 3;
/// Convergence tolerance for the largest absolute logarithmic K-value update.
const THERMO_K_VALUE_TOLERANCE: f64 = 1.0e-8;
/// Maximum Newton steps.
const MAX_ITERATIONS: usize = 80;
/// The class's own default convergence tolerance, and the one `solveNaphtaliSandholm` sets.
const TOLERANCE: f64 = 1.0e-8;
/// Relative perturbation for the numerical Jacobian, and its absolute floor.
const PERTURBATION: f64 = 1.0e-4;
const MIN_PERTURBATION: f64 = 1.0e-8;
/// The largest per-tray, per-component imbalance an early exit accepts, relative to the feed.
const COMPONENT_IMBALANCE_TOLERANCE: f64 = 5.0e-3;

/// The derived state a Jacobian build holds still: the K-values, the vapour component flows,
/// the liquid totals and the two enthalpies.
type DerivedState = (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<f64>, Vec<f64>, Vec<f64>);

/// The whole solve: the MESH state, its thermodynamics, and the Newton loop over both.
struct Mesh {
    /// Tray count, including the two ends.
    n: usize,
    /// Component count.
    c: usize,
    /// Variables per tray: `C` liquid component flows, the temperature, the vapour flow.
    m: usize,
    /// Tray temperatures, K.
    t: Vec<f64>,
    /// Tray pressures, Pa.
    p: Vec<f64>,
    /// Vapour leaving each tray, mol/s.
    v: Vec<f64>,
    /// Liquid leaving each tray, mol/s.
    l: Vec<f64>,
    /// Liquid component flows, mol/s.
    liq: Vec<Vec<f64>>,
    /// Vapour component flows, mol/s - **derived**, `K x V`, never a variable.
    vap: Vec<Vec<f64>>,
    /// K-values.
    k: Vec<Vec<f64>>,
    /// Liquid and vapour molar enthalpies, J/mol.
    hl: Vec<f64>,
    hv: Vec<f64>,
    /// The feed's split, per tray and component, mol/s.
    feed_liq: Vec<Vec<f64>>,
    feed_vap: Vec<Vec<f64>>,
    /// The feed's liquid and vapour enthalpies, J/mol.
    feed_hl: f64,
    feed_hv: f64,
    /// The feed's liquid and vapour totals **per tray**, mol/s: zero everywhere but the tray
    /// the feed enters, which is what the class's `feedLTotal`/`feedVTotal` arrays are.
    feed_l_total: Vec<f64>,
    feed_v_total: Vec<f64>,
    /// The fraction of a tray's traffic that continues inward - one everywhere, because side
    /// draws are not ported and the tray's own fractions are declared out of this pass.
    internal_vapour_fraction: Vec<f64>,
    internal_liquid_fraction: Vec<f64>,
    /// The pinned temperature where a tray has one, `NaN` where its energy equation stands.
    fixed_temperature: Vec<f64>,
    /// The residual's scaling, `totalFeed / N` floored at one.
    flow_scale: f64,
    /// The temperature residual's scaling, 100 K.
    temp_scale: f64,
    /// The mixture and its ideal-gas model, resolved once from the feed's components.
    mixture: azoth_eos::Mixture,
    ideal_gas: azoth_eos::IdealGasModel,
    /// The components' critical constants in the Wilson correlation's units: bara and kelvin.
    wilson: Vec<(f64, f64, f64)>,
    /// The feed's substances, by name, as the outcome's streams carry them.
    components: Vec<String>,
}

impl Mesh {
    /// Every variable's value, in the order the Jacobian is assembled in.
    fn variable(&self, tray: usize, k: usize) -> f64 {
        if k < self.c {
            self.liq[tray][k]
        } else if k == self.c {
            self.t[tray]
        } else {
            self.v[tray]
        }
    }

    /// Set one variable.
    fn set_variable(&mut self, tray: usize, k: usize, value: f64) {
        if k < self.c {
            self.liq[tray][k] = value;
        } else if k == self.c {
            self.t[tray] = value;
        } else {
            self.v[tray] = value;
        }
    }

    /// The residual vector, one block per tray.
    fn residual(&self) -> Vec<f64> {
        let mut f = vec![0.0; self.n * self.m];
        for j in 0..self.n {
            let block = self.residual_for_tray(j);
            f[j * self.m..(j + 1) * self.m].copy_from_slice(&block);
        }
        f
    }

    /// One tray's `C + 2` residuals: the component balances, the energy balance (or the
    /// temperature specification) and the summation equation.
    fn residual_for_tray(&self, j: usize) -> Vec<f64> {
        let mut f = vec![0.0; self.m];
        let lj = self.l[j];
        for i in 0..self.c {
            let mut balance = self.liq[j][i] + self.vap[j][i];
            if j + 1 < self.n {
                balance -= self.internal_liquid_fraction[j + 1] * self.liq[j + 1][i];
            }
            if j > 0 {
                balance -= self.internal_vapour_fraction[j - 1] * self.vap[j - 1][i];
            }
            balance -= self.feed_liq[j][i] + self.feed_vap[j][i];
            f[i] = balance / self.flow_scale;
        }

        if self.fixed_temperature[j].is_finite() {
            f[self.c] = (self.t[j] - self.fixed_temperature[j]) / self.temp_scale;
        } else {
            let mut h = lj * self.hl[j] + self.v[j] * self.hv[j];
            if j + 1 < self.n {
                h -= self.internal_liquid_fraction[j + 1] * self.l[j + 1] * self.hl[j + 1];
            }
            if j > 0 {
                h -= self.internal_vapour_fraction[j - 1] * self.v[j - 1] * self.hv[j - 1];
            }
            h -= self.feed_l_total[j] * self.feed_hl + self.feed_v_total[j] * self.feed_hv;
            let tray_flow = (lj + self.v[j]).max(self.flow_scale);
            let latent = (self.hv[j] - self.hl[j]).abs().max(1.0e3);
            f[self.c] = h / (tray_flow * latent);
        }

        let mut sum_kx = 0.0;
        for i in 0..self.c {
            let x = self.liq[j][i] / lj.max(1.0e-20);
            sum_kx += self.k[j][i] * x;
        }
        f[self.c + 1] = sum_kx - 1.0;
        f
    }

    /// The forced-phase fugacity coefficients' logarithms at `(T, P, z)`.
    ///
    /// The class's `computeSinglePhaseFugacityCoefficients`: one phase, its root **forced** to
    /// the side asked for rather than chosen by a stability test. That is the whole point of
    /// the K-value here - `K_i = phi_L_i(x) / phi_V_i(y)` at the tray's *own* compositions,
    /// where a flash would move `x` to the split's `x*` and leave a summation residual floor
    /// (the class's own javadoc records `~1e-3` in mole fraction).
    fn ln_phi(&self, z: &[f64], t: f64, p: f64, side: RootSide) -> Result<Vec<f64>> {
        let reduced = self.mixture.reduced_parameters(kelvins(t), pascals(p))?;
        Ok(self.mixture.phase_state(&reduced, z, side)?.ln_phi)
    }

    /// One forced phase's molar enthalpy at `(T, P, z)`, J/mol.
    fn single_phase_enthalpy(&self, z: &[f64], t: f64, p: f64, side: RootSide) -> Result<f64> {
        let reduced = self.mixture.reduced_parameters(kelvins(t), pascals(p))?;
        let root = self.mixture.phase_state(&reduced, z, side)?.z;
        let state = molar_enthalpy_entropy(
            &self.mixture,
            &self.ideal_gas,
            kelvins(t),
            pascals(p),
            z,
            root,
        )?;
        Ok(state.h.value)
    }

    /// The Wilson correlation for one component, which is the class's fallback and seed.
    fn wilson_k(&self, i: usize, t: f64, p_bar: f64) -> f64 {
        let (tc, pc, omega) = self.wilson[i];
        (pc / p_bar) * (5.37 * (1.0 + omega) * (1.0 - tc / t)).exp()
    }

    /// Evaluate every tray's thermodynamics.
    fn evaluate_thermo(&mut self) -> Result<()> {
        for j in 0..self.n {
            self.evaluate_thermo_for_tray(j)?;
        }
        Ok(())
    }

    /// Evaluate one tray: its K-values, its enthalpies and its derived vapour flows.
    ///
    /// `V[j]` is **not** recomputed - it is a free variable of the solver - and the vapour
    /// component flows are derived from it and the K-values at the end.
    fn evaluate_thermo_for_tray(&mut self, j: usize) -> Result<()> {
        let c = self.c;
        self.l[j] = self.liq[j].iter().sum::<f64>().max(1.0e-20);
        let x: Vec<f64> = (0..c).map(|i| self.liq[j][i] / self.l[j]).collect();
        let p = self.p[j];
        let t = self.t[j];

        let phi_l = self.ln_phi(&x, t, p, RootSide::Liquid)?;

        // The vapour guess from the previous K-values, or Wilson where there is none.
        let mut y = vec![0.0; c];
        let mut sum_guess = 0.0;
        for i in 0..c {
            let guess = self.k[j][i];
            let k_guess = if (1.0e-20..1.0e15).contains(&guess) {
                guess
            } else {
                self.wilson_k(i, t, p / 1.0e5)
            };
            y[i] = k_guess * x[i];
            sum_guess += y[i];
        }
        for value in &mut y {
            *value = if sum_guess > 1.0e-20 {
                *value / sum_guess
            } else {
                0.0
            };
        }
        if sum_guess <= 1.0e-20 {
            y.copy_from_slice(&x);
        }

        // The self-consistency sweeps: K = phi_L / phi_V, then y = K x / sum(K x).
        for _ in 0..THERMO_K_VALUE_ITERATIONS {
            let phi_v = self.ln_phi(&y, t, p, RootSide::Vapour)?;
            let mut sum = 0.0;
            for i in 0..c {
                let mut k_new = (phi_l[i] - phi_v[i]).exp();
                k_new = k_new.clamp(1.0e-15, 1.0e15);
                self.k[j][i] = k_new;
                y[i] = k_new * x[i];
                sum += y[i];
            }
            let max_update = (0..c)
                .map(|i| (self.k[j][i]).ln().abs())
                .fold(0.0_f64, f64::max);
            for value in &mut y {
                *value = if sum > 1.0e-20 { *value / sum } else { 0.0 };
            }
            if max_update <= THERMO_K_VALUE_TOLERANCE && sum > 1.0e-20 {
                break;
            }
            if sum <= 1.0e-20 {
                y.copy_from_slice(&x);
            }
        }

        // The efficiency proxy is the identity at the ends, so it is applied where it is not.
        let eta = self.tray_eta(j);
        if eta < 1.0 - 1.0e-10 {
            for i in 0..c {
                if self.k[j][i] > 0.0 {
                    self.k[j][i] = self.k[j][i].powf(eta);
                }
            }
        }

        let mut sum_kx = 0.0;
        for i in 0..c {
            sum_kx += self.k[j][i] * x[i];
        }
        for i in 0..c {
            y[i] = if sum_kx > 1.0e-20 {
                self.k[j][i] * x[i] / sum_kx
            } else {
                x[i]
            };
        }

        self.hl[j] = self.single_phase_enthalpy(&x, t, p, RootSide::Liquid)?;
        self.hv[j] = self.single_phase_enthalpy(&y, t, p, RootSide::Vapour)?;
        for i in 0..c {
            self.vap[j][i] = self.k[j][i] * x[i] * self.v[j];
        }
        Ok(())
    }

    /// The Murphree proxy's efficiency for one tray.
    ///
    /// Both ends are rigorous equilibrium - the class forces `trayEta` to one there - and this
    /// port refuses the parameter at its boundary, so every tray is one.
    fn tray_eta(&self, _tray: usize) -> f64 {
        1.0
    }

    /// The numerical Jacobian, one column per variable.
    ///
    /// **Only the perturbed tray's thermodynamics is re-evaluated**, which is what makes the
    /// matrix block-tridiagonal: a tray's K-values and enthalpies depend on its own
    /// temperature and composition and on nothing else, so a perturbed variable can reach the
    /// equations of its own tray and of its two neighbours and no further.
    fn jacobian(&mut self, f0: &[f64]) -> Result<BlockTridiagonal> {
        let mut jac = BlockTridiagonal::new(self.n, self.m)?;
        // **The base is frozen, and restored without re-evaluating.** Every column must be a
        // derivative of the *same* residual vector: `evaluateThermoForTray` advances its own
        // K-value fixed point, so a restore that re-evaluated would leave the next column
        // starting from a different base - the class's own note on
        // `restoreDerivedThermodynamicStateForTray`.
        let base = self.save_derived();
        for tray in 0..self.n {
            for k in 0..self.m {
                let original = self.variable(tray, k);
                let h = (original.abs() * PERTURBATION).max(MIN_PERTURBATION);
                self.set_variable(tray, k, original + h);
                self.evaluate_thermo_for_tray(tray)?;

                let start = tray.saturating_sub(1);
                let end = (tray + 1).min(self.n - 1);
                for j in start..=end {
                    let perturbed = self.residual_for_tray(j);
                    let block = jac.block_mut(j, tray)?;
                    for (equation, row) in block.iter_mut().enumerate() {
                        row[k] = (perturbed[equation] - f0[j * self.m + equation]) / h;
                    }
                }

                self.set_variable(tray, k, original);
                self.restore_derived_for_tray(tray, &base);
            }
        }
        Ok(jac)
    }

    /// The derived state a Jacobian build must hold still: the K-values, the vapour component
    /// flows, the liquid totals and the two enthalpies.
    fn save_derived(&self) -> DerivedState {
        (
            self.k.clone(),
            self.vap.clone(),
            self.l.clone(),
            self.hl.clone(),
            self.hv.clone(),
        )
    }

    /// One tray's derived state, put back as it was - **without** re-evaluating it.
    fn restore_derived_for_tray(&mut self, tray: usize, base: &DerivedState) {
        self.k[tray].clone_from(&base.0[tray]);
        self.vap[tray].clone_from(&base.1[tray]);
        self.l[tray] = base.2[tray];
        self.hl[tray] = base.3[tray];
        self.hv[tray] = base.4[tray];
    }

    /// `computeMassBalanceError`: the products against the feed, relative.
    fn mass_balance_error(&self) -> f64 {
        let mut total_feed = 0.0;
        for j in 0..self.n {
            for i in 0..self.c {
                total_feed += self.feed_liq[j][i] + self.feed_vap[j][i];
            }
        }
        let mut product = self.internal_vapour_fraction[self.n - 1] * self.v[self.n - 1]
            + self.internal_liquid_fraction[0] * self.l[0];
        for j in 0..self.n {
            product += (1.0 - self.internal_vapour_fraction[j]) * self.v[j];
            product += (1.0 - self.internal_liquid_fraction[j]) * self.l[j];
        }
        (product - total_feed).abs() / total_feed.max(1.0e-20)
    }

    /// `computeMaxComponentImbalance`: the worst per-tray, per-component leak, relative to the
    /// feed.
    fn max_component_imbalance(&self) -> f64 {
        let mut total_feed = 0.0;
        for j in 0..self.n {
            for i in 0..self.c {
                total_feed += self.feed_liq[j][i] + self.feed_vap[j][i];
            }
        }
        let denom = total_feed.max(1.0e-20);
        let mut worst: f64 = 0.0;
        for j in 0..self.n {
            for i in 0..self.c {
                let mut balance = self.liq[j][i] + self.vap[j][i];
                if j + 1 < self.n {
                    balance -= self.internal_liquid_fraction[j + 1] * self.liq[j + 1][i];
                }
                if j > 0 {
                    balance -= self.internal_vapour_fraction[j - 1] * self.vap[j - 1][i];
                }
                balance -= self.feed_liq[j][i] + self.feed_vap[j][i];
                worst = worst.max(balance.abs() / denom);
            }
        }
        worst
    }

    /// `meshClosureAcceptable`: the component balances close, relative to the feed.
    fn mesh_closure_acceptable(&self) -> bool {
        self.max_component_imbalance() <= COMPONENT_IMBALANCE_TOLERANCE
    }

    /// Save the primary state.
    fn save(&self) -> (Vec<Vec<f64>>, Vec<f64>, Vec<f64>) {
        (self.liq.clone(), self.t.clone(), self.v.clone())
    }

    /// Restore the primary state.
    fn restore(&mut self, saved: &(Vec<Vec<f64>>, Vec<f64>, Vec<f64>)) {
        self.liq.clone_from(&saved.0);
        self.t.clone_from(&saved.1);
        self.v.clone_from(&saved.2);
    }
}

/// The class's `solveNaphtaliSandholm`: initialize, seed, and iterate the Newton correction.
///
/// Returns the column's outcome, or a refusal naming what could not be reached.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a configuration this port does not carry - see the module
/// documentation - and [`AzothError::SolverNotConverged`] when the residual misses its gate.
pub fn solve(setup: &ColumnSetup, pressures: &[Pressure]) -> Result<ColumnOutcome> {
    if setup.reboiler_temperature.is_none() {
        return Err(AzothError::invalid_input(
            "solver_type",
            "naphtali_sandholm without a reboiler temperature would run \
             `solveBubblePointMethod`, whose V-cascade this port does not carry: pin the \
             reboiler's temperature, which is the class's own `useOverallMBClosure` condition",
        ));
    }
    if setup.top_specification.is_some() || setup.bottom_specification.is_some() {
        return Err(AzothError::invalid_input(
            "solver_type",
            "a product specification re-runs the whole column inside \
             `solveWithSpecificationTargets`, which is a second integration this port does not \
             carry: state the ends' temperatures instead",
        ));
    }

    let mut mesh = Mesh::initialize(setup, pressures)?;
    let mut residual = mesh.residual();
    let mut norm = vector_norm(&residual);

    let mut best = mesh.save();
    let mut best_norm = norm;
    let mut failed_steps = 0_usize;
    let mut stagnation_count = 0_usize;
    let mut non_descent = 0_usize;
    let mut completed = 0_usize;
    let mut previous_best = best_norm;

    for iteration in 1..=MAX_ITERATIONS {
        completed = iteration;
        if norm < TOLERANCE {
            return mesh.outcome(iteration - 1, norm);
        }

        let jacobian = mesh.jacobian(&residual)?;
        let rhs: Vec<Vec<f64>> = (0..mesh.n)
            .map(|j| residual[j * mesh.m..(j + 1) * mesh.m].to_vec())
            .collect();
        // **The class's own fallback**: where the block route reports the system singular it
        // eliminates the same matrix densely, and only a dense failure counts as a failed step.
        let step = match jacobian.solve(&rhs)? {
            Some(solution) => Some(solution),
            None => dense_solve(&jacobian.to_dense(), &rhs.concat()).map(|flat| {
                (0..mesh.n)
                    .map(|j| flat[j * mesh.m..(j + 1) * mesh.m].to_vec())
                    .collect::<Vec<Vec<f64>>>()
            }),
        };
        let Some(step) = step else {
            mesh.restore(&best);
            mesh.evaluate_thermo()?;
            residual = mesh.residual();
            norm = vector_norm(&residual);
            failed_steps += 1;
            if failed_steps > 5 {
                return Err(not_converged(iteration, norm));
            }
            continue;
        };
        let mut dx: Vec<f64> = step.concat();
        apply_trust_region(&mesh, &mut dx);
        let Some((_alpha, trial_residual, trial_norm)) = line_search(&mut mesh, &dx, norm)? else {
            mesh.restore(&best);
            mesh.evaluate_thermo()?;
            break;
        };

        if !trial_norm.is_finite() {
            mesh.restore(&best);
            mesh.evaluate_thermo()?;
            residual = mesh.residual();
            norm = vector_norm(&residual);
            failed_steps += 1;
            if failed_steps > 5 {
                return Err(not_converged(iteration, norm));
            }
            continue;
        }

        // A step that increases the residual by more than half is rejected outright.
        if trial_norm > 1.5 * norm && trial_norm > TOLERANCE {
            mesh.restore(&best);
            mesh.evaluate_thermo()?;
            residual = mesh.residual();
            norm = vector_norm(&residual);
            failed_steps += 1;
            if failed_steps > 10 {
                break;
            }
            continue;
        }

        if trial_norm >= norm {
            non_descent += 1;
            if non_descent >= MAX_NON_DESCENT_LINE_SEARCH_STEPS {
                mesh.restore(&best);
                mesh.evaluate_thermo()?;
                break;
            }
        }

        residual = trial_residual;
        norm = trial_norm;
        failed_steps = 0;

        if norm < best_norm {
            best_norm = norm;
            best = mesh.save();
        }

        if iteration % 10 == 0 {
            if best_norm > 0.95 * previous_best {
                stagnation_count += 1;
                if stagnation_count >= 2 {
                    mesh.restore(&best);
                    mesh.evaluate_thermo()?;
                    if mesh.mass_balance_error() < 0.005 && mesh.mesh_closure_acceptable() {
                        return mesh.outcome(iteration, best_norm);
                    }
                    break;
                }
            } else {
                stagnation_count = 0;
            }
            previous_best = best_norm;
        }
    }

    mesh.restore(&best);
    mesh.evaluate_thermo()?;
    let norm = vector_norm(&mesh.residual());
    // The class's own partial-convergence test, and the same refusal when it fails.
    if norm < TOLERANCE * 100.0 && mesh.mesh_closure_acceptable() {
        return mesh.outcome(completed, norm);
    }
    Err(not_converged(completed, norm))
}

/// The class's trust-region clamp, which limits each variable's step and scales the whole
/// vector by the smallest ratio it had to apply.
fn apply_trust_region(mesh: &Mesh, dx: &mut [f64]) {
    const MAX_DT: f64 = 10.0;
    let mut scale = 1.0_f64;
    for j in 0..mesh.n {
        let base = j * mesh.m;
        for i in 0..mesh.c {
            let cap = 0.5 * mesh.liq[j][i] + 1.0e-3 * mesh.flow_scale;
            let step = dx[base + i].abs();
            if step > cap && cap > 0.0 {
                scale = scale.min(cap / step);
            }
        }
        if !mesh.fixed_temperature[j].is_finite() {
            let step = dx[base + mesh.c].abs();
            if step > MAX_DT {
                scale = scale.min(MAX_DT / step);
            }
        }
        let cap = 0.5 * mesh.v[j] + 0.05 * mesh.flow_scale;
        let step = dx[base + mesh.c + 1].abs();
        if step > cap && cap > 0.0 {
            scale = scale.min(cap / step);
        }
    }
    if scale < 1.0 {
        for value in dx.iter_mut() {
            *value *= scale;
        }
    }
}

/// `applyUpdate`: the step, with the physical bounds the class enforces.
fn apply_update(mesh: &mut Mesh, dx: &[f64], alpha: f64) {
    for j in 0..mesh.n {
        let base = j * mesh.m;
        for i in 0..mesh.c {
            mesh.liq[j][i] = (mesh.liq[j][i] - alpha * dx[base + i]).max(1.0e-20);
        }
        if !mesh.fixed_temperature[j].is_finite() {
            mesh.t[j] = (mesh.t[j] - alpha * dx[base + mesh.c]).clamp(100.0, 1000.0);
        }
        mesh.v[j] = (mesh.v[j] - alpha * dx[base + mesh.c + 1]).max(0.0);
    }
}

/// The class's Armijo backtracking on `||F||`, which tries `alpha = 1, 1/2, 1/4, ...` and
/// accepts the **best** trial it saw - not necessarily the first that descended.
fn line_search(
    mesh: &mut Mesh,
    dx: &[f64],
    current_norm: f64,
) -> Result<Option<(f64, Vec<f64>, f64)>> {
    const ARMIJO: f64 = 1.0e-4;
    const BACKTRACK: f64 = 0.5;
    const MAX_BACKTRACK: usize = 15;

    let saved = mesh.save();
    let mut alpha = 1.0;
    let mut best: Option<(f64, Vec<f64>, f64)> = None;
    let mut last = f64::NAN;
    for _ in 0..MAX_BACKTRACK {
        mesh.restore(&saved);
        apply_update(mesh, dx, alpha);
        mesh.evaluate_thermo()?;
        let residual = mesh.residual();
        let norm = vector_norm(&residual);
        last = alpha;
        if norm.is_finite()
            && best
                .as_ref()
                .is_none_or(|(_, _, best_norm)| norm < *best_norm)
        {
            best = Some((alpha, residual.clone(), norm));
        }
        if norm.is_finite() && norm < (1.0 - ARMIJO * alpha) * current_norm {
            return Ok(Some((alpha, residual, norm)));
        }
        if alpha < 0.01 {
            break;
        }
        alpha *= BACKTRACK;
    }

    let Some((best_alpha, best_residual, best_norm)) = best else {
        mesh.restore(&saved);
        return Ok(None);
    };
    if best_alpha != last {
        mesh.restore(&saved);
        apply_update(mesh, dx, best_alpha);
        mesh.evaluate_thermo()?;
        let residual = mesh.residual();
        let norm = vector_norm(&residual);
        return Ok(Some((best_alpha, residual, norm)));
    }
    Ok(Some((best_alpha, best_residual, best_norm)))
}

/// A dense Gaussian elimination of the same system, for a diagnostic.
fn dense_solve(matrix: &[Vec<f64>], rhs: &[f64]) -> Option<Vec<f64>> {
    let size = rhs.len();
    let mut aug = vec![vec![0.0; size + 1]; size];
    for row in 0..size {
        aug[row][..size].copy_from_slice(&matrix[row][..size]);
        aug[row][size] = rhs[row];
    }
    for column in 0..size {
        let mut pivot = column;
        let mut largest = aug[column][column].abs();
        for row in column + 1..size {
            if aug[row][column].abs() > largest {
                largest = aug[row][column].abs();
                pivot = row;
            }
        }
        if largest < 1.0e-300 {
            return None;
        }
        aug.swap(column, pivot);
        for row in column + 1..size {
            let factor = aug[row][column] / aug[column][column];
            for k in column..=size {
                aug[row][k] -= factor * aug[column][k];
            }
        }
    }
    let mut dx = vec![0.0; size];
    for i in (0..size).rev() {
        let mut sum = aug[i][size];
        for j in i + 1..size {
            sum -= aug[i][j] * dx[j];
        }
        dx[i] = sum / aug[i][i];
    }
    Some(dx)
}

/// The Euclidean norm of a vector, which is the class's `vectorNorm`.
fn vector_norm(v: &[f64]) -> f64 {
    v.iter().map(|value| value * value).sum::<f64>().sqrt()
}

/// A refusal naming the residual the solve stopped at.
///
/// **The residual is the message.** `AzothError::SolverNotConverged` carries the iteration
/// count, the residual and the tolerance, and a caller reading them can tell a solve that
/// stopped one order short from one that never descended - which is what the class's own
/// `logger.warn("Naphtali-Sandholm did not converge ...")` prints.
fn not_converged(iterations: usize, residual: f64) -> AzothError {
    AzothError::SolverNotConverged {
        iterations: iterations as u32,
        residual,
        tolerance: TOLERANCE,
    }
}

impl Mesh {
    /// `initialize`: the feed's split, the pins, and the seed.
    ///
    /// The class reads the feed from the column's own external streams and flashes each at its
    /// state, splitting it into a vapour and a liquid contribution on the tray it enters. This
    /// port's column has one feed, at `feed_stage`, and its flash is the same one the
    /// substitution solve's feed-stage tray runs.
    fn initialize(setup: &ColumnSetup, pressures: &[Pressure]) -> Result<Self> {
        let n = pressures.len();
        let feed = &setup.feed;
        let names: Vec<&str> = feed.components.iter().map(String::as_str).collect();
        let (mixture, ideal_gas) =
            azoth_eos::databank::mixture_of(&names, azoth_eos::Cubic::Pr, None)?;
        let c = feed.components.len();
        let m = c + 2;
        let p: Vec<f64> = pressures.iter().map(|value| value.value).collect();
        let wilson: Vec<(f64, f64, f64)> = mixture
            .components()
            .iter()
            .map(|component| {
                (
                    component.tc.value,
                    component.pc.value / 1.0e5,
                    component.omega,
                )
            })
            .collect();

        let mut fixed_temperature = vec![f64::NAN; n];
        if setup.has_reboiler {
            fixed_temperature[0] = setup
                .reboiler_temperature
                .map_or(f64::NAN, |value| value.value);
        }
        if setup.has_condenser {
            fixed_temperature[n - 1] = setup
                .condenser_temperature
                .map_or(f64::NAN, |value| value.value);
        }

        // The feed's own flash, split into the vapour and liquid it carries onto its tray.
        let flash = azoth_eos::pt_flash::pt_flash(&mixture, feed.t, feed.p, &feed.z)?;
        let (n_vap, z_vap, n_liq, z_liq) = match flash.phase {
            azoth_eos::Phase::TwoPhase => {
                let beta = flash.beta.ok_or_else(|| {
                    AzothError::invalid_input(
                        "feed",
                        "the feed's flash is two-phase with no vapour fraction",
                    )
                })?;
                (
                    feed.n * beta,
                    flash.y.clone(),
                    feed.n * (1.0 - beta),
                    flash.x.clone(),
                )
            }
            azoth_eos::Phase::AllVapour => (feed.n, feed.z.clone(), 0.0, feed.z.clone()),
            azoth_eos::Phase::AllLiquid => (0.0, feed.z.clone(), feed.n, feed.z.clone()),
            azoth_eos::Phase::Trivial => {
                return Err(AzothError::invalid_input(
                    "feed",
                    "the feed's flash is a trivial solution, so nothing proves which phase it \
                     is and the solver has no split to seed from",
                ));
            }
        };

        let mut mesh = Self {
            n,
            c,
            m,
            t: vec![0.0; n],
            p,
            v: vec![0.0; n],
            l: vec![0.0; n],
            liq: vec![vec![0.0; c]; n],
            vap: vec![vec![0.0; c]; n],
            k: vec![vec![0.0; c]; n],
            hl: vec![0.0; n],
            hv: vec![0.0; n],
            feed_liq: vec![vec![0.0; c]; n],
            feed_vap: vec![vec![0.0; c]; n],
            feed_hl: 0.0,
            feed_hv: 0.0,
            feed_l_total: vec![0.0; n],
            feed_v_total: vec![0.0; n],
            internal_vapour_fraction: vec![1.0; n],
            internal_liquid_fraction: vec![1.0; n],
            fixed_temperature,
            flow_scale: 1.0,
            temp_scale: 100.0,
            mixture,
            ideal_gas,
            wilson,
            components: feed.components.clone(),
        };
        for (i, &value) in z_liq.iter().enumerate() {
            mesh.feed_liq[setup.feed_stage][i] = n_liq * value;
        }
        for (i, &value) in z_vap.iter().enumerate() {
            mesh.feed_vap[setup.feed_stage][i] = n_vap * value;
        }
        if n_liq > 0.0 {
            mesh.feed_hl =
                mesh.single_phase_enthalpy(&z_liq, feed.t.value, feed.p.value, RootSide::Liquid)?;
        }
        if n_vap > 0.0 {
            mesh.feed_hv =
                mesh.single_phase_enthalpy(&z_vap, feed.t.value, feed.p.value, RootSide::Vapour)?;
        }

        // ---- mol/s to mol/hr, which is the unit every constant of this solver is written in.
        for j in 0..n {
            for i in 0..c {
                mesh.feed_liq[j][i] *= MOL_PER_HOUR;
                mesh.feed_vap[j][i] *= MOL_PER_HOUR;
            }
        }
        mesh.feed_l_total[setup.feed_stage] = n_liq * MOL_PER_HOUR;
        mesh.feed_v_total[setup.feed_stage] = n_vap * MOL_PER_HOUR;

        let total_feed: f64 = mesh
            .feed_liq
            .iter()
            .chain(mesh.feed_vap.iter())
            .flat_map(|tray| tray.iter())
            .sum();
        mesh.flow_scale = (total_feed / n as f64).max(1.0);

        // The class pins a tray's temperature before it seeds, which is what anchors the
        // bottom of the profile.
        for j in 0..n {
            if mesh.fixed_temperature[j].is_finite() {
                mesh.t[j] = mesh.fixed_temperature[j];
            }
        }

        mesh.initialize_from_column(setup, total_feed)?;
        mesh.evaluate_thermo()?;
        Ok(mesh)
    }

    /// `initializeTrayStateFromColumn`: the MESH variables from the column's own converged
    /// trays.
    ///
    /// **This is the class's warm start, and this port always takes it.** NeqSim reaches it
    /// when the column has been solved before and its thermodynamics have not moved; a port
    /// that owns both solves can always take it, so the mesh solve here starts from the
    /// substitution core's answer.
    ///
    /// **The class's cold seed is not ported, and the reason is measured.**
    /// `initializeTrayState` + `seedSolutionForDirectNewton` build a profile whose per-component
    /// shoot ends with the top tray's light components projected onto the overhead, and the
    /// step that repairs it is `runBostonSullivanRefinement` - 464 lines of inside-out
    /// refinement from a local Antoine fit. Without it the seed's top tray carries a liquid
    /// whose light component is `1e-15` of its feed, and the first Newton step from that state
    /// moves the variables by `1e11` times their own size: measured on the captured binary
    /// column, the solver does not descend at all. With the warm start the same column
    /// converges in five iterations to the state NeqSim's own mesh solve reaches in four.
    ///
    /// # Errors
    /// [`AzothError::SolverNotConverged`] from the substitution core, and
    /// [`AzothError::InvalidInput`] where the seeded profile has a tray the mesh solve cannot
    /// start from.
    fn initialize_from_column(&mut self, setup: &ColumnSetup, total_feed: f64) -> Result<()> {
        let mut substitution = setup.clone();
        substitution.solver_type = SolverType::DirectSubstitution;
        let seed = crate::kernels::distillation_column::distillation_column(&substitution)?;

        // The class scales the previous flows to the new total feed, which is what makes the
        // seed usable after an input change as well as after a solve.
        let previous_total = seed.trays.last().map_or(0.0, |tray| tray.gas_n)
            + seed.trays.first().map_or(0.0, |tray| tray.liquid_n);
        let factor = if previous_total > 1.0e-12 {
            total_feed / (previous_total * MOL_PER_HOUR)
        } else {
            0.0
        };
        if !factor.is_finite() || factor <= 0.0 {
            return Err(AzothError::invalid_input(
                "column",
                format!(
                    "the substitution core's products carry {previous_total} mol/s, so there is \
                     no scale to seed the mesh solve from"
                ),
            ));
        }

        for (j, tray) in seed.trays.iter().enumerate() {
            // `!(x > 0.0)`, which is the class's own test and which a `NaN` fails - written
            // as a positive comparison because clippy reads the negation of a partial order as
            // ambiguous, and this one is not: `NaN` is neither greater nor less.
            if tray.gas_n <= 0.0
                || tray.gas_n.is_nan()
                || tray.liquid_n <= 0.0
                || tray.liquid_n.is_nan()
                || !tray.temperature.value.is_finite()
            {
                return Err(AzothError::invalid_input(
                    "column",
                    format!(
                        "tray {j} left the substitution core with {} mol/s of vapour, {} of \
                         liquid and {} K, and the mesh solve cannot start from that",
                        tray.gas_n, tray.liquid_n, tray.temperature.value
                    ),
                ));
            }
            self.t[j] = if self.fixed_temperature[j].is_finite() {
                self.fixed_temperature[j]
            } else {
                tray.temperature.value
            };
            self.v[j] = tray.gas_n * MOL_PER_HOUR * factor;
            self.l[j] = tray.liquid_n * MOL_PER_HOUR * factor;
            for i in 0..self.c {
                self.liq[j][i] = (tray.liquid_z[i] * self.l[j]).max(1.0e-20);
            }
        }
        Ok(())
    }

    /// The solution, as the kernel's outcome: the trays, the two products and the closure.
    fn outcome(&self, iterations: usize, residual: f64) -> Result<ColumnOutcome> {
        let components = self.components.clone();
        let mut trays = Vec::with_capacity(self.n);
        for j in 0..self.n {
            // ---- back to mol/s, which is what a stream carries.
            let liquid_n = self.l[j] / MOL_PER_HOUR;
            let gas_n = self.v[j] / MOL_PER_HOUR;
            // The compositions are ratios, so they come from the mol/hr quantities and not
            // from the converted flows - dividing a mol/hr component flow by a mol/s total is
            // a factor of 3600, and the library refuses a composition that does not sum to one.
            let liquid_z: Vec<f64> = self.liq[j].iter().map(|value| value / self.l[j]).collect();
            let gas_z: Vec<f64> = self.vap[j].iter().map(|value| value / self.v[j]).collect();
            trays.push(TrayProfile {
                temperature: kelvins(self.t[j]),
                pressure: pascals(self.p[j]),
                gas_n,
                liquid_n,
                gas_z,
                liquid_z,
            });
        }
        let top = &trays[self.n - 1];
        let bottom = &trays[0];
        let distillate = Stream::from_pt(
            components.clone(),
            top.gas_z.clone(),
            top.gas_n,
            top.pressure,
            top.temperature,
        )?;
        let bottoms = Stream::from_pt(
            components,
            bottom.liquid_z.clone(),
            bottom.liquid_n,
            bottom.pressure,
            bottom.temperature,
        )?;
        let (condenser_duty, reboiler_duty) = self.duties()?;
        // J/hr from the mesh, W on the products' side: the closure is a rate, so the feed
        // comes across the same conversion the duties do.
        let feed_enthalpy: f64 = (0..self.n)
            .map(|j| self.feed_l_total[j] * self.feed_hl + self.feed_v_total[j] * self.feed_hv)
            .sum::<f64>()
            / MOL_PER_HOUR;
        let products = distillate.n * distillate.h.value + bottoms.n * bottoms.h.value;
        let energy_residual = if feed_enthalpy.abs() > 0.0 {
            (feed_enthalpy + condenser_duty.value + reboiler_duty.value - products).abs()
                / feed_enthalpy.abs()
        } else {
            0.0
        };
        let mut mass_residual: f64 = 0.0;
        for i in 0..self.c {
            let supplied: f64 = (0..self.n)
                .map(|j| self.feed_liq[j][i] + self.feed_vap[j][i])
                .sum::<f64>()
                / MOL_PER_HOUR;
            if supplied.abs() > 1.0e-12 {
                let delivered = distillate.n * distillate.z[i] + bottoms.n * bottoms.z[i];
                mass_residual = mass_residual.max((supplied - delivered).abs() / supplied.abs());
            }
        }
        let tray_count = trays.len();
        Ok(ColumnOutcome {
            trays,
            // **This solve has no reactive route, and the model refuses the pair rather than
            // ignoring the flag**: its MESH equations take their fugacities from the mesh's own
            // mixture, where the class's trays take them from the tray's own flash - so a
            // reactive section here could only be silently non-reactive.
            warnings: Vec::new(),
            // **The mesh solve carries no side draws**, and the model refuses the pair rather
            // than ignoring the fractions: these equations take their fugacities from the mesh's
            // own mixture, so a draw has no tray outlet here to split.
            gas_side_draws: vec![None; tray_count],
            liquid_side_draws: vec![None; tray_count],
            pumparounds: vec![None; tray_count],
            distillate,
            bottoms,
            condenser_duty,
            reboiler_duty,
            iterations: iterations as u32,
            // **The temperature residual is the scaled MESH norm here**, which is this
            // solver's own convergence measure: NeqSim's column reports `NaN` under
            // `NAPHTALI_SANDHOLM` because the class's mean tray-temperature change is a
            // substitution solve's diagnostic and this solve does not form one.
            temperature_residual: residual,
            mass_residual,
            energy_residual,
        })
    }

    /// The two duties, from the ends' energy balance against their traffic.
    fn duties(&self) -> Result<(azoth_core::units::Power, azoth_core::units::Power)> {
        use azoth_core::units::watts;
        let top = self.n - 1;
        let mut condenser_in = if top > 0 {
            self.internal_vapour_fraction[top - 1] * self.v[top - 1] * self.hv[top - 1]
        } else {
            0.0
        };
        let mut reboiler_in = if self.n > 1 {
            self.internal_liquid_fraction[1] * self.l[1] * self.hl[1]
        } else {
            0.0
        };
        // The feed is the duty's *inlet* only on the tray it enters, which is the same
        // arithmetic the substitution kernel's `inlet_enthalpy` does over its streams.
        reboiler_in += self.feed_l_total[0] * self.feed_hl + self.feed_v_total[0] * self.feed_hv;
        condenser_in +=
            self.feed_l_total[top] * self.feed_hl + self.feed_v_total[top] * self.feed_hv;
        let condenser_out = self.v[top] * self.hv[top] + self.l[top] * self.hl[top];
        let reboiler_out = self.v[0] * self.hv[0] + self.l[0] * self.hl[0];
        // J/hr to W: the mesh carries the class's mol/hr, and a duty is a rate in watts.
        Ok((
            watts((condenser_out - condenser_in) / MOL_PER_HOUR),
            watts((reboiler_out - reboiler_in) / MOL_PER_HOUR),
        ))
    }
}
