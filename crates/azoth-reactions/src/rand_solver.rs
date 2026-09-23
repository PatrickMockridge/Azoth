//! The modified-RAND reactive equilibrium solve, from
//! `flashops/reactiveflash/ModifiedRANDSolver.java`.
//!
//! Simultaneous chemical and phase equilibrium by minimising the Gibbs energy subject to the
//! element balances: the unknowns are the moles `n[j][i]` and one Lagrange multiplier `λ_k`
//! per element row, and each pass corrects them together:
//!
//! ```text
//! e[j][i] = g0_i + ln x_j,i + ln phi_j,i - sum_k lambda_k A_ki     the potential error
//! C_kl    = sum_j sum_i A_ki A_li n_j,i                            the RAND matrix
//! rhs_k   = (b_k - sum A n) + sum A n e                            the inventory plus the error
//! n_j,i  *= exp(alpha (-e_j,i + sum_k A_ki dlambda_k))             the damped step
//! lambda += alpha dlambda
//! ```
//!
//! # What is reproduced and what is not, yet
//!
//! Reproduced: the potential error, the matrix and its right-hand side, the Tikhonov term,
//! the symmetric diagonal scaling, the damped step with its backtracking line search, the two
//! convergence tests, the single-phase damping rule and the *sliding-window* rule the
//! multiphase branch uses, the `NR = 0` short circuit, and the DIIS extrapolation the
//! multipliers are driven with.
//!
//! # The ionic branch, as three rules the caller states
//!
//! An ionic fluid changes three things, and each is a property of the *model* rather than of
//! the solve - so each is part of [`IonicPhase`] rather than something this module works out:
//! an ion's potential is its aqueous formation Gibbs energy less its reference-state log
//! fugacity coefficient, an ion's mole number in a gas phase is the floor, and the Lagrange
//! multipliers start from a non-gas phase because an ion exists only in solution.
//!
//! **What that costs is the reference itself.** `lnPhiRef` is
//! `ph.getLogInfiniteDiluteFugacity(i, solvent)`, and NeqSim computes it only when
//! `system instanceof SystemFurstElectrolyteEos` - on every other system the correction is
//! zero and an ion's potential is its formation Gibbs energy alone. The measurement is
//! `captures/rand_solver_ionic_probe.tsv`: a `SystemFurstElectrolyteEos` brine reaches
//! `ln_phi_ref = [0, -127.4, -171.5, -173.4, -160.8]` while a neutral control reports zeros
//! and `is_electrolyte_eos = false` on the same solver.
//!
//! # The split is not always determined, and the port does not pretend otherwise
//!
//! Where two phases converge to the *same* composition the Gibbs energy is flat along the
//! direction that trades moles between them, so every split satisfies the equilibrium
//! conditions and the one a run reports is decided by its path. The captured 300 K water-gas
//! shift is such a state: NeqSim ends at `(0.0141, 0.9859)` and this port at `(0.016, 0.984)`,
//! while both agree on the composition - the moles summed over the phases - to `1.2e-5`. The
//! same fluid at 600 K agrees to `1e-6`. Neither code converges a split tighter than its own
//! relaxed multiphase tolerance of `1e-4`, so that is the precision the comparison lives at.
//!
//! # `ln phi` is a closure, not a dependency
//!
//! The solver needs a fugacity coefficient at a *trial* composition, and NeqSim gets one by
//! writing the moles into its own phase and calling `init(1)`. Here the caller supplies it, per
//! phase, because the class's phases each carry their own root - which keeps this crate free of
//! an equation of state, and keeps the question of *which* phase model answers a separate one.

use azoth_core::{AzothError, Result};

use crate::diis::DiisAccelerator;
use crate::formula_matrix::matrix_rank;

/// The residual both the potential error and the element balance must come in under, from
/// `ModifiedRANDSolver.TOL`.
pub const TOL: f64 = 1.0e-9;

/// The pass cap, from `MAX_ITER`.
pub const MAX_ITERATIONS: usize = 500;

/// The floor a mole number is kept above, from `EPS`.
pub const EPS: f64 = 1.0e-30;

/// The gas constant NeqSim's solver carries, from `R_GAS`. **Not `azoth-core`'s**: this is a
/// literal in the class, and the standard potentials are computed with it.
pub const R_GAS: f64 = 8.314462;

/// The reference pressure the `ln(P/P_ref)` term is taken against, from `P_REF`.
pub const P_REF: f64 = 1.0;

/// The temperature the Cp polynomials are integrated from, from `computeG0`'s `T0`.
pub const T0: f64 = 298.15;

/// The singular-pivot floor `solveLinear` refuses below, from `mx < 1e-30`.
pub const LINEAR_PIVOT_FLOOR: f64 = 1.0e-30;

/// The DIIS history's length, from `DIIS_DEPTH`.
pub const DIIS_DEPTH: usize = 6;

/// The pass before DIIS may extrapolate, from `DIIS_START`.
pub const DIIS_START: usize = 5;

/// The sliding window the multiphase damping rule averages over, from the class's `WINDOW`.
pub const DAMPING_WINDOW: u32 = 10;

/// The potential error a multiphase solve accepts, from `maxE < 1.0e-4` - a looser test than
/// [`TOL`], and the reason the same fluid can report `1e-14` on one path and `1e-5` on another.
pub const MULTIPHASE_ERROR_TOLERANCE: f64 = 1.0e-4;

/// The element residual a *neutral* multiphase solve accepts, from `elementTolerance`, which is
/// `1e-8` when a charge row is present and `1e-4` when it is not.
pub const MULTIPHASE_ELEMENT_TOLERANCE: f64 = 1.0e-4;

/// A component's ideal-gas heat-capacity polynomial and its two formation properties.
///
/// The five coefficients are `Cp = A + B T + C T^2 + D T^3 + E T^4`, read from the component
/// databank by `eos.ideal_gas_cp`'s reader.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermoData {
    /// Ideal-gas enthalpy of formation at [`T0`], in J/mol.
    pub enthalpy_of_formation: f64,
    /// Ideal-gas absolute entropy at [`T0`], in J/(mol*K).
    pub absolute_entropy: f64,
    /// Gibbs energy of formation, in J/mol - the fallback when neither of the other two is
    /// stated.
    pub gibbs_energy_of_formation: f64,
    /// The heat-capacity polynomial's five coefficients, `A` through `E`.
    pub cp: [f64; 5],
}

/// The standard chemical potentials `g0`, from `computeG0`'s neutral branch.
///
/// `g0_i = (hT - T sT) / (R T) + ln(P/P_ref)`, where `hT` and `sT` are the formation
/// properties corrected from [`T0`] to `T` by integrating the Cp polynomial:
///
/// ```text
/// hT = dHf + integral(T0,T) Cp dT
/// sT = S0  + integral(T0,T) Cp/T dT
/// ```
///
/// The three fallbacks are the class's own, in its order: with `dHf` and `S0` but no Cp data
/// the polynomial is dropped, with only `dGf298` the Gibbs energy is used directly, and with
/// nothing at all the potential is `ln(P/P_ref)` alone - which is the ideal-gas potential of a
/// substance whose formation data the databank does not carry.
/// **An ion takes a different branch**, and it is NeqSim's own: `g0_i = dGf_aq / (R T) -
/// lnPhiRef_i`, with no Cp integration and no entropy term, because an ion's standard state is
/// the aqueous one its formation Gibbs energy is already measured against. `is_ion` empty
/// means a fluid with no ion and `ln_phi_ref` empty means every reference is zero - the
/// non-electrolyte case, where the class uses `dGf_aq / (R T)` alone.
#[must_use]
pub fn standard_potentials(
    data: &[ThermoData],
    temperature: f64,
    pressure: f64,
    is_ion: &[bool],
    ln_phi_ref: &[f64],
) -> Vec<f64> {
    let rt = R_GAS * temperature;
    let ln_p = (pressure / P_REF).ln();
    data.iter()
        .enumerate()
        .map(|(index, entry)| {
            if is_ion.get(index).copied().unwrap_or(false) {
                let reference = ln_phi_ref.get(index).copied().unwrap_or(0.0);
                return entry.gibbs_energy_of_formation / rt - reference;
            }
            let [cp_a, cp_b, cp_c, cp_d, cp_e] = entry.cp;
            let has_cp_data = cp_a.abs() > 1.0e-10 || cp_b.abs() > 1.0e-10;
            let has_thermo = entry.enthalpy_of_formation.abs() > 1.0e-10
                || entry.absolute_entropy.abs() > 1.0e-10;

            if has_thermo && has_cp_data {
                let dt = temperature - T0;
                let dt2 = temperature.powi(2) - T0.powi(2);
                let dt3 = temperature.powi(3) - T0.powi(3);
                let dt4 = temperature.powi(4) - T0.powi(4);
                let dt5 = temperature.powi(5) - T0.powi(5);
                let ln_ratio = (temperature / T0).ln();
                let delta_h = cp_a * dt
                    + cp_b / 2.0 * dt2
                    + cp_c / 3.0 * dt3
                    + cp_d / 4.0 * dt4
                    + cp_e / 5.0 * dt5;
                let delta_s = cp_a * ln_ratio
                    + cp_b * dt
                    + cp_c / 2.0 * dt2
                    + cp_d / 3.0 * dt3
                    + cp_e / 4.0 * dt4;
                let h_t = entry.enthalpy_of_formation + delta_h;
                let s_t = entry.absolute_entropy + delta_s;
                (h_t - temperature * s_t) / rt + ln_p
            } else if has_thermo {
                (entry.enthalpy_of_formation - temperature * entry.absolute_entropy) / rt + ln_p
            } else if entry.gibbs_energy_of_formation.abs() > 1.0e-10 {
                entry.gibbs_energy_of_formation / rt + ln_p
            } else {
                ln_p
            }
        })
        .collect()
}

/// A phase model: a phase's index and a mole-fraction composition in, each component's
/// `ln(phi_i)` there out. The index is what lets a caller answer with the phase's own root.
pub type PhaseLogPhi<'a> = &'a mut dyn FnMut(usize, &[f64]) -> Result<Vec<f64>>;

/// What only an ionic fluid's caller knows, from `hasIonicSpecies`' three uses.
///
/// Three facts, and each is a property of the *model* rather than of the solve: which
/// components are ions, which phases are gas, and each ion's reference-state log fugacity
/// coefficient. The last is `computeLnPhiRef`'s `lnPhiRef`, which NeqSim computes only when
/// `isElectrolyteEOS` - `system instanceof SystemFurstElectrolyteEos` - and leaves at zero
/// everywhere else, because for a plain cubic the infinite-dilution coefficient of an ion is
/// not a number its parameters can produce.
#[derive(Debug, Clone, PartialEq)]
pub struct IonicPhase {
    /// `isIon[i]`: whether a component carries a charge.
    pub is_ion: Vec<bool>,
    /// `isGasPhase[j]`: whether a phase is a gas, which is where an ion may not be.
    pub is_gas: Vec<bool>,
    /// `lnPhiRef[i]`: an ion's `ln(phi_i)` at infinite dilution in the solvent, zero wherever
    /// the phase model is not an electrolyte one.
    pub ln_phi_ref: Vec<f64>,
}

/// One phase's starting state, which `initialize` reads off the phase object: the composition
/// and the fraction the phase carries.
#[derive(Debug, Clone, PartialEq)]
pub struct PhaseFeed {
    /// `phase.getComponent(i).getx()`, one per component.
    pub fractions: Vec<f64>,
    /// `phase.getBeta()`, floored at [`EPS`].
    pub beta: f64,
}

/// What the solve answers with.
#[derive(Debug, Clone, PartialEq)]
pub struct RandSolution {
    /// The overall moles, summed over the phases - the composition the fluid ends at.
    pub moles: Vec<f64>,
    /// The moles in each phase, `n[j][i]`.
    pub phase_moles: Vec<Vec<f64>>,
    /// Each phase's share, `nPhase[j] / totalMoles` - **the solver's own split**, which is not
    /// the phase objects' `beta` after the driver writes them back.
    pub phase_amounts: Vec<f64>,
    /// `totalMoles`: `nPhase` summed.
    pub total_moles: f64,
    /// The element Lagrange multipliers at the answer, one per element row.
    pub lambda: Vec<f64>,
    /// Passes taken.
    pub iterations: u32,
    /// The largest absolute potential error over the species holding more than `1e-10` of the
    /// total.
    pub max_error: f64,
    /// `computeElementResidual`: the scaled root-mean-square element deviation, which is a
    /// *normalised* measure and not a bare difference.
    pub element_residual: f64,
    /// `getFinalResidual`: `max(max_error, element_residual)`, lowered by any extrapolated step
    /// DIIS kept - which is why it can sit below the pair it is the maximum of.
    pub final_residual: f64,
    /// Whether the residual passed the convergence test that applied - the strict one, or the
    /// relaxed multiphase one.
    pub converged: bool,
    /// `getDiisStepsAccepted`: extrapolated steps the solve kept.
    pub diis_steps: u32,
}

/// The single-phase neutral RAND solve: [`solve`] with one phase at `beta = 1`.
///
/// The class has one `solve()`, and a one-phase system is the case where its sums over the
/// phases collapse; the wrapper is that case named, so the callers that have a single phase do
/// not have to build a phase list to say so.
///
/// # Errors
/// [`AzothError::InvalidInput`] where the shapes disagree, and whatever `ln_phi` raises.
pub fn solve_single_phase(
    a_matrix: &[Vec<f64>],
    g0: &[f64],
    b: &[f64],
    feed_moles: &[f64],
    ln_phi: &mut dyn FnMut(&[f64]) -> Result<Vec<f64>>,
) -> Result<RandSolution> {
    let total: f64 = feed_moles.iter().sum();
    if total <= 0.0 {
        return Err(AzothError::InvalidInput {
            field: "feed_moles".to_string(),
            reason: "the feed holds no moles".to_string(),
        });
    }
    let fractions: Vec<f64> = feed_moles.iter().map(|moles| moles / total).collect();
    let mut one_phase = |_phase: usize, x: &[f64]| ln_phi(x);
    solve(
        a_matrix,
        g0,
        b,
        total,
        &[PhaseFeed {
            fractions,
            beta: 1.0,
        }],
        &mut one_phase,
        None,
    )
}

/// The neutral RAND solve over a phase list, `ModifiedRANDSolver.solve`.
///
/// `ln_phi` is handed a phase's index and its *mole-fraction* composition and answers each
/// component's `ln(phi_i)` there; the index is what a caller needs to pick the phase's own
/// root. NeqSim writes the moles into its own phase objects and reads
/// `getFugacityCoefficient` back, and this is the same question asked of whatever phase model
/// the caller has.
///
/// `total_moles` is the system's `getTotalNumberOfMoles`, which scales the phase moles at
/// `initialize`; `b` is the frozen element inventory, or the feed's own `A n` where the caller
/// has none.
///
/// # The rules that only apply above one phase
///
/// `np > 1` changes four things and each is the class's own: the step starts at a tenth rather
/// than the whole of the Newton direction, a *sliding window* over ten iterations replaces the
/// iteration-to-iteration damping rule, a residual of `1e-4` counts as converged rather than
/// [`TOL`], and DIIS is allowed to extrapolate the multipliers. The single-phase branch keeps
/// the strict rule, which is why the same fluid can report `1.4e-14` on one path and `2.0e-5`
/// on the other.
///
/// # Errors
/// [`AzothError::InvalidInput`] where the shapes disagree, and whatever `ln_phi` raises.
// Indexed rather than iterated, deliberately: every loop here is NeqSim's own loop, over a
// matrix, a phase and its components, and the correspondence is what makes the port checkable
// against the class line by line.
#[allow(clippy::needless_range_loop, clippy::too_many_arguments)]
pub fn solve(
    a_matrix: &[Vec<f64>],
    g0: &[f64],
    b: &[f64],
    total_moles: f64,
    phases: &[PhaseFeed],
    ln_phi: PhaseLogPhi<'_>,
    ions: Option<&IonicPhase>,
) -> Result<RandSolution> {
    let ne = b.len();
    let nc = g0.len();
    let np = phases.len();
    if a_matrix.len() != ne || np == 0 {
        return Err(AzothError::InvalidInput {
            field: "a_matrix".to_string(),
            reason: format!(
                "{ne} element row(s), {} potential(s) and {np} phase(s)",
                g0.len()
            ),
        });
    }
    for row in a_matrix {
        if row.len() != nc {
            return Err(AzothError::InvalidInput {
                field: "a_matrix".to_string(),
                reason: format!("a row has {} entries against {nc} component(s)", row.len()),
            });
        }
    }
    for phase in phases {
        if phase.fractions.len() != nc {
            return Err(AzothError::InvalidInput {
                field: "phases".to_string(),
                reason: format!(
                    "a phase has {} entries against {nc} component(s)",
                    phase.fractions.len()
                ),
            });
        }
    }

    // `initialize`: the moles and fractions the phase objects hold, floored.
    let mut n = vec![vec![0.0_f64; nc]; np];
    let mut fractions = vec![vec![0.0_f64; nc]; np];
    let mut n_phase = vec![0.0_f64; np];
    let mut beta = vec![0.0_f64; np];
    for j in 0..np {
        beta[j] = phases[j].beta.max(EPS);
        for i in 0..nc {
            fractions[j][i] = phases[j].fractions[i];
            n[j][i] = (fractions[j][i] * beta[j] * total_moles).max(EPS);
        }
        n_phase[j] = beta[j] * total_moles;
    }
    // `initialize`'s own pin, and `recalcTotals`' before every pass: **an ion is not in a gas
    // phase**. NeqSim floors such a mole number at `EPS` rather than removing it, so the
    // component keeps its column and the constraint is a value rather than a shape.
    pin_ions_in_gas(&mut n, ions);
    let mut total: f64 = n_phase.iter().sum();

    let mut ln_phi_here = vec![vec![0.0_f64; nc]; np];
    for j in 0..np {
        ln_phi_here[j] = ln_phi(j, &fractions[j])?;
    }

    // `NR = 0` is not an iteration: the element balance alone fixes the composition, and the
    // class returns the feed at once rather than driving a matrix that is singular where
    // `NE > rank(A)`.
    let independent_reactions = nc.saturating_sub(matrix_rank(a_matrix, nc));
    if independent_reactions == 0 {
        let element = element_residual(a_matrix, b, &n, total);
        return Ok(RandSolution {
            moles: (0..nc).map(|i| (0..np).map(|j| n[j][i]).sum()).collect(),
            phase_moles: n,
            phase_amounts: (0..np).map(|j| n_phase[j] / total).collect(),
            total_moles: total,
            lambda: vec![0.0; ne],
            iterations: 0,
            max_error: element,
            element_residual: element,
            final_residual: element,
            converged: true,
            diis_steps: 0,
        });
    }

    // `initializeLambda` starts from a **non-gas** phase when ions are present: an ion exists
    // only in solution, and a gas phase's ionic mole fraction is the floor rather than a state.
    let reference_phase = ions
        .and_then(|constraints| constraints.is_gas.iter().position(|gas| !gas))
        .unwrap_or(0);
    let mut lambda = initial_lambda(
        a_matrix,
        g0,
        &fractions[reference_phase],
        &ln_phi_here[reference_phase],
    );

    // The accelerator runs on the multipliers, with the element deviation as its error.
    let mut diis = DiisAccelerator::new(ne, DIIS_DEPTH);
    let mut diis_steps = 0_u32;

    // The class starts a multiphase solve at a tenth of the direction.
    let mut damping = if np > 1 { 0.1 } else { 1.0 };
    let mut previous_residual = f64::MAX;
    let mut stagnation = 0_u32;
    let mut window_residual = f64::MAX;
    let mut iterations_since_window = 0_u32;
    let element_tolerance = MULTIPHASE_ELEMENT_TOLERANCE;

    let mut final_error = f64::MAX;
    let mut final_element = f64::MAX;
    let mut final_residual = f64::MAX;
    let mut converged = false;
    let mut iterations = 0_u32;

    for iteration in 0..MAX_ITERATIONS {
        iterations = iteration as u32 + 1;

        // The potential error, per phase and component.
        let mut error = vec![vec![0.0_f64; nc]; np];
        for j in 0..np {
            for i in 0..nc {
                let x_i = fractions[j][i].max(EPS);
                let mut row_sum = 0.0;
                for k in 0..ne {
                    row_sum += lambda[k] * a_matrix[k][i];
                }
                let value = g0[i] + x_i.ln() + ln_phi_here[j][i] - row_sum;
                error[j][i] = if value.is_finite() { value } else { 0.0 };
            }
        }

        // The RAND matrix and its right-hand side, summed over the phases.
        let mut c = vec![vec![0.0_f64; ne]; ne];
        let mut rhs = vec![0.0_f64; ne];
        for k in 0..ne {
            let mut element_sum = 0.0;
            let mut element_error = 0.0;
            for j in 0..np {
                for i in 0..nc {
                    element_sum += a_matrix[k][i] * n[j][i];
                    element_error += a_matrix[k][i] * n[j][i] * error[j][i];
                }
            }
            rhs[k] = (b[k] - element_sum) + element_error;
            for l in 0..ne {
                let mut value = 0.0;
                for j in 0..np {
                    for i in 0..nc {
                        value += a_matrix[k][i] * a_matrix[l][i] * n[j][i];
                    }
                }
                c[k][l] = value;
            }
            // Tikhonov regularization, as the class does it.
            c[k][k] += (1.0e-10 * c[k][k].abs()).max(1.0e-14);
        }

        // Symmetric diagonal scaling, for the conditioning the ions cause.
        let scale: Vec<f64> = (0..ne)
            .map(|k| {
                let diagonal = c[k][k].abs();
                if diagonal > 1.0e-30 {
                    1.0 / diagonal.sqrt()
                } else {
                    1.0
                }
            })
            .collect();
        for k in 0..ne {
            for l in 0..ne {
                c[k][l] *= scale[k] * scale[l];
            }
            rhs[k] *= scale[k];
        }

        let Some(mut delta) = solve_linear(&c, &rhs) else {
            break;
        };
        for k in 0..ne {
            delta[k] *= scale[k];
        }

        let n_old = n.clone();
        let lambda_old = lambda.clone();

        // The damped step, with up to five halvings while the element residual does not fall.
        let mut alpha = damping;
        let mut accepted = false;
        for _ in 0..5 {
            n.clone_from_slice(&n_old);
            lambda.clone_from_slice(&lambda_old);
            for j in 0..np {
                for i in 0..nc {
                    let mut correction = alpha * -error[j][i];
                    for k in 0..ne {
                        correction += alpha * a_matrix[k][i] * delta[k];
                    }
                    if !correction.is_finite() {
                        correction = 0.0;
                    }
                    let correction = correction.clamp(-3.0, 3.0);
                    let stepped = n[j][i] * correction.exp();
                    n[j][i] = if !stepped.is_finite() || stepped < EPS {
                        EPS
                    } else {
                        stepped
                    };
                }
            }
            for k in 0..ne {
                lambda[k] += alpha * delta[k];
            }
            recalc_totals(&mut n_phase, &mut total, &mut n, ions);
            let element_residual = element_residual(a_matrix, b, &n, total);
            if element_residual < previous_residual * 1.5 || alpha < 0.05 {
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }

        if !accepted {
            // The class's own fallback: restore and take a tenth of the step.
            n.clone_from_slice(&n_old);
            lambda.clone_from_slice(&lambda_old);
            let alpha = 0.1_f64;
            for j in 0..np {
                for i in 0..nc {
                    let mut correction = alpha * -error[j][i];
                    for k in 0..ne {
                        correction += alpha * a_matrix[k][i] * delta[k];
                    }
                    if !correction.is_finite() {
                        correction = 0.0;
                    }
                    let correction = correction.clamp(-3.0, 3.0);
                    let stepped = n[j][i] * correction.exp();
                    n[j][i] = if !stepped.is_finite() || stepped < EPS {
                        EPS
                    } else {
                        stepped
                    };
                }
            }
            for k in 0..ne {
                lambda[k] += alpha * delta[k];
            }
            recalc_totals(&mut n_phase, &mut total, &mut n, ions);
        }

        for j in 0..np {
            fractions[j] = n[j].iter().map(|moles| moles / n_phase[j]).collect();
        }
        for j in 0..np {
            ln_phi_here[j] = ln_phi(j, &fractions[j])?;
        }

        // The convergence test: the worst potential error among the species that matter, and
        // the element balance.
        let mut max_error = 0.0_f64;
        for j in 0..np {
            for i in 0..nc {
                if n[j][i] > 1.0e-10 * total {
                    max_error = max_error.max(error[j][i].abs());
                }
            }
        }
        let element = element_residual(a_matrix, b, &n, total);
        final_error = max_error;
        final_element = element;
        let residual = max_error.max(element);
        final_residual = residual;

        if max_error < TOL && element < TOL {
            converged = true;
            break;
        }
        if np > 1 && max_error < MULTIPHASE_ERROR_TOLERANCE && element < element_tolerance {
            converged = true;
            break;
        }

        if np == 1 {
            // Single-phase: the classic per-iteration damping.
            if residual < previous_residual * 0.9 {
                damping = (damping * 1.5).min(1.0);
                stagnation = 0;
            } else if residual > previous_residual * 1.1 {
                damping = (damping * 0.5).max(0.01);
                stagnation += 1;
            } else {
                stagnation += 1;
            }
        } else {
            // Multi-phase: an immediate penalty, plus a ten-iteration window, because a small
            // damping makes under one percent of progress per pass and the per-iteration rule
            // would then hold the damping down forever.
            if residual > previous_residual * 1.5 {
                damping = (damping * 0.5).max(0.01);
                stagnation += 1;
            } else if residual > previous_residual * 1.05 {
                stagnation += 1;
            }
            iterations_since_window += 1;
            if iterations_since_window >= DAMPING_WINDOW {
                if residual < window_residual * 0.5 {
                    let ceiling = if residual > 10.0 {
                        0.3
                    } else if residual > 1.0 {
                        0.5
                    } else {
                        1.0
                    };
                    damping = (damping * 2.0).min(ceiling);
                    stagnation = 0;
                } else if residual > window_residual * 2.0 {
                    damping = (damping * 0.25).max(0.01);
                    stagnation += DAMPING_WINDOW;
                }
                window_residual = residual;
                iterations_since_window = 0;
            }
        }

        // A long stagnation at a small residual is accepted, as the class accepts it.
        if stagnation > 50 && max_error < 1.0e-3 && element < element_tolerance {
            converged = true;
            break;
        }

        // DIIS, on the multipliers, with the element deviation as its error. **Every** pass is
        // recorded - the class adds the entry outside the `DIIS_START` gate - and the gate only
        // decides whether an extrapolation is tried.
        let residual_vector = element_residual_vector(a_matrix, b, &n, total);
        diis.add_entry(&lambda, &residual_vector)?;
        if iteration >= DIIS_START && diis.can_extrapolate() {
            if let Some(extrapolated) = diis.extrapolate() {
                let saved_moles = n.clone();
                let saved_lambda = lambda.clone();
                for j in 0..np {
                    for i in 0..nc {
                        let mut correction = 0.0;
                        for k in 0..ne {
                            correction += (extrapolated[k] - lambda[k]) * a_matrix[k][i];
                        }
                        let correction = correction.clamp(-3.0, 3.0);
                        let stepped = n[j][i] * correction.exp();
                        n[j][i] = if !stepped.is_finite() || stepped < EPS {
                            EPS
                        } else {
                            stepped
                        };
                    }
                }
                lambda.clone_from_slice(&extrapolated);
                recalc_totals(&mut n_phase, &mut total, &mut n, ions);
                for j in 0..np {
                    fractions[j] = n[j].iter().map(|moles| moles / n_phase[j]).collect();
                }
                for j in 0..np {
                    ln_phi_here[j] = ln_phi(j, &fractions[j])?;
                }

                // The extrapolated state is judged on the whole residual - recomputed at the
                // extrapolated composition, because the potential error is not the
                // extrapolated one - and rolled back when it does not at least hold it.
                let mut diis_error = 0.0_f64;
                for j in 0..np {
                    for i in 0..nc {
                        if n[j][i] > 1.0e-10 * total {
                            let x_i = fractions[j][i].max(EPS);
                            let mut row_sum = 0.0;
                            for k in 0..ne {
                                row_sum += lambda[k] * a_matrix[k][i];
                            }
                            let value = g0[i] + x_i.ln() + ln_phi_here[j][i] - row_sum;
                            diis_error = diis_error.max(value.abs());
                        }
                    }
                }
                let diis_element = element_residual(a_matrix, b, &n, total);
                let diis_residual = diis_error.max(diis_element);
                if diis_residual < residual * 1.1 {
                    final_residual = final_residual.min(diis_residual);
                    diis_steps += 1;
                } else {
                    n = saved_moles;
                    lambda = saved_lambda;
                    recalc_totals(&mut n_phase, &mut total, &mut n, ions);
                    for j in 0..np {
                        fractions[j] = n[j].iter().map(|moles| moles / n_phase[j]).collect();
                    }
                    for j in 0..np {
                        ln_phi_here[j] = ln_phi(j, &fractions[j])?;
                    }
                }
            }
        }

        // `prevResidual` is the *last* word on the pass, after any extrapolation DIIS kept -
        // which is what makes a kept extrapolation carry into the next pass's damping and line
        // search rather than being compared against the value it replaced.
        previous_residual = final_residual;
    }

    Ok(RandSolution {
        moles: (0..nc).map(|i| (0..np).map(|j| n[j][i]).sum()).collect(),
        phase_moles: n,
        phase_amounts: (0..np).map(|j| n_phase[j] / total).collect(),
        total_moles: total,
        lambda,
        iterations,
        max_error: final_error,
        element_residual: final_element,
        final_residual,
        converged,
        diis_steps,
    })
}

/// `recalcTotals`: the phase amounts from the moles, floored, and the total. The mole
/// fractions follow from the two and are not stored here, because the class recomputes them
/// where it needs them.
fn recalc_totals(
    n_phase: &mut [f64],
    total: &mut f64,
    n: &mut [Vec<f64>],
    ions: Option<&IonicPhase>,
) {
    pin_ions_in_gas(n, ions);
    *total = 0.0;
    for (j, phase_moles) in n.iter().enumerate() {
        n_phase[j] = phase_moles.iter().sum::<f64>().max(EPS);
        *total += n_phase[j];
    }
    *total = total.max(EPS);
}

/// `enforceIonPhaseConstraints`: every ion's mole number in a gas phase is the floor.
fn pin_ions_in_gas(n: &mut [Vec<f64>], ions: Option<&IonicPhase>) {
    let Some(constraints) = ions else {
        return;
    };
    for (j, phase) in n.iter_mut().enumerate() {
        if !constraints.is_gas.get(j).copied().unwrap_or(false) {
            continue;
        }
        for (i, moles) in phase.iter_mut().enumerate() {
            if constraints.is_ion.get(i).copied().unwrap_or(false) {
                *moles = EPS;
            }
        }
    }
}

/// `computeElementResidual`: the scaled root-mean-square deviation of `A n` from `b`.
///
/// **The scaling is the class's and it is not cosmetic**: each element's deviation is divided
/// by `max(|b_k|, max(totalMoles 1e-6, 1e-10))`, so a charge row whose inventory is near zero
/// reports a bounded number rather than amplifying its own round-off.
#[allow(clippy::needless_range_loop)] // the class's own loops over the rows and phases
fn element_residual(a_matrix: &[Vec<f64>], b: &[f64], n: &[Vec<f64>], total: f64) -> f64 {
    let scale_floor = (total * 1.0e-6).max(1.0e-10);
    let mut sum = 0.0_f64;
    for k in 0..b.len() {
        let mut element_sum = 0.0;
        for phase in n {
            for i in 0..phase.len() {
                element_sum += a_matrix[k][i] * phase[i];
            }
        }
        let scale = b[k].abs().max(scale_floor);
        let deviation = (element_sum - b[k]) / scale;
        sum += deviation * deviation;
    }
    sum.sqrt()
}

/// `computeElementResidualVector`: the same deviations, one per element row, which is what the
/// DIIS accelerator takes as its error.
#[allow(clippy::needless_range_loop)] // the class's own loops over the rows and phases
fn element_residual_vector(
    a_matrix: &[Vec<f64>],
    b: &[f64],
    n: &[Vec<f64>],
    total: f64,
) -> Vec<f64> {
    let scale_floor = (total * 1.0e-6).max(1.0e-10);
    (0..b.len())
        .map(|k| {
            let mut element_sum = 0.0;
            for phase in n {
                for i in 0..phase.len() {
                    element_sum += a_matrix[k][i] * phase[i];
                }
            }
            (element_sum - b[k]) / b[k].abs().max(scale_floor)
        })
        .collect()
}

/// `initializeLambda`: the least-squares multipliers for the initial potentials,/// `initializeLambda`: the least-squares multipliers for the initial potentials,
/// `initializeLambda`: the least-squares multipliers for the initial potentials,
/// `lambda = (A A^T)^-1 A h` with `h_i = g0_i + ln x_i + ln phi_i`.
#[allow(clippy::needless_range_loop)] // the same index loops the class writes
fn initial_lambda(
    a_matrix: &[Vec<f64>],
    g0: &[f64],
    fractions: &[f64],
    ln_phi: &[f64],
) -> Vec<f64> {
    let ne = a_matrix.len();
    let nc = g0.len();
    let h: Vec<f64> = (0..nc)
        .map(|i| g0[i] + fractions[i].max(EPS).ln() + ln_phi[i])
        .collect();
    let mut ata = vec![vec![0.0_f64; ne]; ne];
    let mut ath = vec![0.0_f64; ne];
    for k in 0..ne {
        for l in 0..ne {
            for i in 0..nc {
                ata[k][l] += a_matrix[k][i] * a_matrix[l][i];
            }
        }
        for i in 0..nc {
            ath[k] += a_matrix[k][i] * h[i];
        }
    }
    solve_linear(&ata, &ath).unwrap_or_else(|| vec![0.0; ne])
}

/// `solveLinear`: Gaussian elimination with partial pivoting, refusing a pivot under
/// `1e-30`.
#[allow(clippy::needless_range_loop)] // the elimination's own loops, as the class writes them
fn solve_linear(matrix: &[Vec<f64>], rhs: &[f64]) -> Option<Vec<f64>> {
    let dim = rhs.len();
    let mut augment = vec![vec![0.0_f64; dim + 1]; dim];
    for i in 0..dim {
        augment[i][..dim].copy_from_slice(&matrix[i][..dim]);
        augment[i][dim] = rhs[i];
    }

    for column in 0..dim {
        let mut pivot = column;
        let mut largest = augment[column][column].abs();
        for row in (column + 1)..dim {
            if augment[row][column].abs() > largest {
                largest = augment[row][column].abs();
                pivot = row;
            }
        }
        if largest < LINEAR_PIVOT_FLOOR {
            return None;
        }
        augment.swap(column, pivot);
        for row in (column + 1)..dim {
            let factor = augment[row][column] / augment[column][column];
            for k in column..=dim {
                augment[row][k] -= factor * augment[column][k];
            }
        }
    }

    let mut solution = vec![0.0_f64; dim];
    for i in (0..dim).rev() {
        let mut sum = augment[i][dim];
        for k in (i + 1)..dim {
            sum -= augment[i][k] * solution[k];
        }
        solution[i] = sum / augment[i][i];
    }
    Some(solution)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_linear_solve_refuses_a_singular_pivot() {
        let matrix = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        assert_eq!(solve_linear(&matrix, &[1.0, 2.0]), None);
        let good = vec![vec![2.0, 0.0], vec![0.0, 4.0]];
        assert_eq!(solve_linear(&good, &[2.0, 8.0]), Some(vec![1.0, 2.0]));
    }

    /// The three fallbacks, in the class's own order.
    #[test]
    fn the_standard_potentials_fall_back_the_way_the_class_does() {
        let with_cp = ThermoData {
            enthalpy_of_formation: -110_525.0,
            absolute_entropy: 197.66,
            gibbs_energy_of_formation: -137_168.0,
            cp: [25.56759, 6.096130, 4.054656, -2.671301, 0.131021],
        };
        let without_cp = ThermoData {
            cp: [0.0; 5],
            ..with_cp
        };
        let only_gibbs = ThermoData {
            enthalpy_of_formation: 0.0,
            absolute_entropy: 0.0,
            ..with_cp
        };
        let nothing = ThermoData {
            enthalpy_of_formation: 0.0,
            absolute_entropy: 0.0,
            gibbs_energy_of_formation: 0.0,
            cp: [0.0; 5],
        };

        let potentials = standard_potentials(
            &[with_cp, without_cp, only_gibbs, nothing],
            600.0,
            1.0,
            &[],
            &[],
        );
        // The Cp polynomial and the constant-Cp fallback are different numbers, which is the
        // point of integrating it; the last two are the `dGf` and `ln(P)` branches.
        assert!(
            (potentials[0] - potentials[1]).abs() > 1.0,
            "{potentials:?}"
        );
        assert!((potentials[2] - (-137_168.0 / (R_GAS * 600.0))).abs() < 1e-12);
        assert_eq!(potentials[3], 0.0, "ln(P/P_ref) at the reference pressure");
    }

    /// **An ion takes the aqueous branch and no other.**
    ///
    /// `computeG0`'s two: a neutral is `(h_T - T s_T)/(R T) + ln(P/P_ref)` from its Cp
    /// polynomial, an ion is `dGf_aq/(R T) - lnPhiRef` with no integration at all - its
    /// standard state is the aqueous one its formation Gibbs energy is already measured
    /// against. The last assertion is the class's own degenerate case: with no electrolyte
    /// phase model `lnPhiRef` is zero, and the correction is the bare formation term.
    #[test]
    fn an_ion_takes_the_aqueous_standard_state() {
        let neutral = ThermoData {
            enthalpy_of_formation: -110_500.0,
            absolute_entropy: 197.7,
            gibbs_energy_of_formation: -137_168.0,
            cp: [30.87, -1.29e-2, 2.79e-5, -1.27e-8, 0.0],
        };
        let ion = ThermoData {
            enthalpy_of_formation: -240_100.0,
            absolute_entropy: 111.0,
            gibbs_energy_of_formation: -261_905.0,
            cp: [0.0; 5],
        };
        let data = [neutral, ion];
        let reference = -127.398_508_666_942_85;
        let rt = R_GAS * 298.15;

        let with_reference =
            standard_potentials(&data, 298.15, 1.01325, &[false, true], &[0.0, reference]);
        assert!(
            (with_reference[1] - (data[1].gibbs_energy_of_formation / rt - reference)).abs()
                < 1.0e-12,
            "{}",
            with_reference[1]
        );

        // Without an electrolyte phase model the reference is zero, and the ion's potential is
        // its formation Gibbs energy alone.
        let without_reference = standard_potentials(&data, 298.15, 1.01325, &[false, true], &[]);
        assert!((without_reference[1] - data[1].gibbs_energy_of_formation / rt).abs() < 1.0e-12);

        // And a mask that names no ion leaves every potential on the neutral branch.
        let all_neutral = standard_potentials(&data, 298.15, 1.01325, &[], &[]);
        assert!((all_neutral[1] - without_reference[1]).abs() > 1.0e-6);
    }

    /// **An ion is not in a gas phase**, which is `enforceIonPhaseConstraints`: its mole number
    /// there is the floor rather than a computed amount.
    ///
    /// The constraint is a *value* and not a shape - the component keeps its column - so it is
    /// asserted on the phase's moles.
    #[test]
    fn an_ion_is_pinned_out_of_a_gas_phase() {
        let a_matrix = vec![vec![1.0, 1.0]];
        let g0 = vec![0.0, 0.0];
        let b = vec![1.0];
        let phases = vec![
            PhaseFeed {
                fractions: vec![0.5, 0.5],
                beta: 0.5,
            },
            PhaseFeed {
                fractions: vec![0.5, 0.5],
                beta: 0.5,
            },
        ];
        let mut ln_phi = |_phase: usize, x: &[f64]| -> Result<Vec<f64>> { Ok(vec![0.0; x.len()]) };

        let with_gas = IonicPhase {
            is_ion: vec![true, false],
            is_gas: vec![true, false],
            ln_phi_ref: vec![0.0, 0.0],
        };
        let solution = solve(
            &a_matrix,
            &g0,
            &b,
            1.0,
            &phases,
            &mut ln_phi,
            Some(&with_gas),
        )
        .expect("the ionic solve runs");
        assert_eq!(
            solution.phase_moles[0][0], EPS,
            "the gas phase holds no ion"
        );
        assert!(
            solution.phase_moles[1][0] > 1.0e-9,
            "and the condensed phase holds it: {}",
            solution.phase_moles[1][0]
        );

        // The same fluid with no gas phase keeps the ion where the solve puts it.
        let condensed = IonicPhase {
            is_ion: vec![true, false],
            is_gas: vec![false, false],
            ln_phi_ref: vec![0.0, 0.0],
        };
        let free = solve(
            &a_matrix,
            &g0,
            &b,
            1.0,
            &phases,
            &mut ln_phi,
            Some(&condensed),
        )
        .expect("the ionic solve runs");
        assert!(
            free.phase_moles[0][0] > 1.0e-9,
            "nothing pins it where no phase is a gas"
        );
    }
}
