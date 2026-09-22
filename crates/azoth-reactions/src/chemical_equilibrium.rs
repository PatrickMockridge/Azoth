//! `reactions.chemical_equilibrium` - the Smith-Missen reactive solve.
//!
//! Spec: `specs/models/reactions/chemical_equilibrium.toml`, which carries the boundary
//! this takes, why the three optional modes are not in it, and what the capture shows
//! about a fluid NeqSim cannot solve.
//!
//! # The iteration
//!
//! Each pass forms the `M` matrix, the reduced chemical potentials, and the two
//! Lagrange systems, then takes a damped step:
//!
//! ```text
//! M[i][k]     = delta_ik / n_i                        (the ideal form)
//! mu[i]       = mu_ref[i] + ln(n_i) - ln(n_t) + ln(gamma_i)
//! AMA         = A M^-1 A^T ,   AMU = A M^-1 mu
//! [AMA  c^T] [lambda]   [AMU + correction]
//! [c     0 ] [ tau  ] = [ sum(n_i mu_i)   ]
//! dn          = M^-1 (A^T lambda - mu) + n tau
//! ```
//!
//! `c` is a conservation coupling row and `tau` the total-moles multiplier. The
//! conservation row is `b` unless the phase's current element amounts disagree with it,
//! in which case `c` becomes those amounts and `b - A n` is carried as a correction.
//!
//! # What this is not
//!
//! **`M` is the ideal form only.** NeqSim's `useFugacityDerivatives` adds
//! `d(ln phi_i)/d(n_j)` to it, `useFullMMatrix` subtracts `1/n_t` from every entry, and
//! `useAdaptiveDerivatives` switches the first on after five iterations. The first two
//! have **no caller in `src/main`**; the third is reached, but only from `solveChemEq`'s
//! second refinement onward, and the derivative it needs is a phase-model surface this
//! library does not have. See the spec's assumptions.
//!
//! **The mole-fraction standard state only.** `getLogReactionActivity` also has a
//! solute-molality branch, which is `SystemPitzer`'s, and a Deshmukh-Mather branch.

use azoth_core::{AzothError, CalcResult, Result, Warning, apply_checks};

use crate::equilibrium_constant::GAS_CONSTANT;
use crate::linalg::solve_lu;
use crate::model_gen;

/// The floor every mole number is raised to before a logarithm, from
/// `ChemicalEquilibrium.MIN_MOLES`.
pub const MIN_MOLES: f64 = 1e-60;

/// The relative tolerance the conservation-correction check uses, from
/// `CONSERVATION_CORRECTION_TOLERANCE`.
pub const CONSERVATION_CORRECTION_TOLERANCE: f64 = 1e-8;

/// Consecutive non-improving iterations before the solve gives up, from
/// `STAGNATION_LIMIT`.
pub const STAGNATION_LIMIT: u32 = 10;

/// Whether the conservation correction runs, which is NeqSim's
/// `system.getNumberOfPhases() == 1`.
///
/// A phase that is the whole system conserves the element amounts it was handed; a phase
/// that is one of several does not, and NeqSim then takes `b` as given. **Which branch is
/// taken changes the answer and not only its digits**: on the aqueous phase of a flashed
/// CO2-water fluid the correction drives the solve away from the fixed point it would
/// otherwise reach, and it does not converge at all.
pub const CONSERVATION_CORRECTION_PHASES: usize = 1;

/// Which standard state the activity term is taken on, NeqSim's
/// `ChemicalReactionConcentrationBasis`.
///
/// `SystemPitzer` is the only system that selects [`Self::SoluteMolality`], and the branch
/// changes the term for every component whose reference state is `solute` from
/// `ln n_i - ln n_t` to `ln n_i - ln w_solvent`, where `w_solvent` is the solvent's own
/// mass. **The Jacobian does not follow it**: `M[i][k]` stays `delta_ik / n_i`, which is
/// the derivative of the *mole-fraction* form, so a molality-basis solve iterates a
/// residual its own matrix does not describe. NeqSim does the same
/// (`ChemicalEquilibrium.java:210-235` builds `M` with no reference to the basis), and
/// reproducing it rather than correcting it is the tie-breaker the port follows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConcentrationBasis {
    /// `ln n_i - ln n_t` for every component, and the only form the captured fluids take.
    MoleFraction,
    /// `ln n_i - ln w_solvent` for a solute and the mole-fraction form for a solvent.
    SoluteMolality,
}

impl std::str::FromStr for ConcentrationBasis {
    type Err = AzothError;

    /// A basis this library does not carry is refused rather than defaulted, for the reason
    /// the reaction sources are: the two answer different questions.
    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        match value.trim().to_ascii_lowercase().as_str() {
            "mole_fraction" | "mole-fraction" => Ok(Self::MoleFraction),
            "solute_molality" | "solute-molality" => Ok(Self::SoluteMolality),
            other => Err(AzothError::InvalidInput {
                field: "concentration_basis".to_string(),
                reason: format!("`{other}` is not one of `mole_fraction`, `solute_molality`"),
            }),
        }
    }
}

/// The concentration basis together with the two facts the component vectors cannot state.
///
/// The reference-state split and the solvent's mass are databank properties of the
/// substances, and this id is handed vectors; the caller supplies them as the live phase
/// has them, which is what `calculateSolventWeight` computes there.
#[derive(Debug, Clone, Copy)]
struct ActivityBasis<'a> {
    /// Which standard state the term is taken on.
    basis: ConcentrationBasis,
    /// The solvent's own mass in kg, floored at [`MIN_MOLES`] as NeqSim floors it.
    solvent_weight: f64,
    /// One entry per component: 1 where its reference state is `solvent`.
    solvent_mask: &'a [f64],
}

impl ActivityBasis<'_> {
    /// The `ln` of each component's denominator, one per component: `n_t`, or the solvent's
    /// mass for a solute on the molality basis.
    ///
    /// `total` is NeqSim's `n_t`, which is **the phase's total moles read once when the
    /// solver is constructed and never updated** (`ChemicalEquilibrium.java:199`). It is
    /// not the running sum of the reactive set: `carbonate` splits one species into two, so
    /// a sum that followed the trial would drift away from the value NeqSim divides by for
    /// the whole solve. Where the phase holds only its reactive substances the two agree at
    /// the first pass and part afterwards.
    fn log_denominators(&self, total: f64, species: usize) -> Vec<f64> {
        let total = total.max(MIN_MOLES);
        let solvent = self.solvent_weight.max(MIN_MOLES);
        (0..species)
            .map(|i| {
                let is_solvent = self.solvent_mask.get(i).copied().unwrap_or(0.0) != 0.0;
                match self.basis {
                    ConcentrationBasis::MoleFraction => total.ln(),
                    ConcentrationBasis::SoluteMolality if !is_solvent => solvent.ln(),
                    ConcentrationBasis::SoluteMolality => total.ln(),
                }
            })
            .collect()
    }
}

/// Result of `reactions.chemical_equilibrium`.
#[derive(Debug, Clone, PartialEq)]
pub struct ChemicalEquilibriumResult {
    /// The moles of each species at the answer, in the caller's component order.
    ///
    /// **Present whether or not the solve converged**, because NeqSim returns the last
    /// iterate either way and its callers read it: a fluid the solver cannot solve still
    /// leaves a composition, and `converged` is what says whether to believe it.
    pub moles: Vec<f64>,
    /// Passes taken.
    pub iterations: u32,
    /// The final error: the sum of `|dn_i| / n_i` over the species that moved.
    pub error: f64,
    /// Whether the error came in under the tolerance. A fluid NeqSim cannot solve
    /// reports `false` and this port reproduces that rather than converging further.
    pub converged: bool,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ChemicalEquilibriumResult {
    const CALC_ID: &'static str = "reactions.chemical_equilibrium";
    const FIELDS: &'static [&'static str] =
        &["moles", "iterations", "error", "converged", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// The reactive equilibrium composition of a phase.
///
/// `a_matrix` is the element matrix, `n_elements` rows by one column per component, with
/// the electroneutrality row among them and `b`'s matching entry zero. `moles` is the
/// starting composition, `chem_ref` the components' reduced standard-state potentials
/// (`mu/RT`) and `log_activity` their `ln(gamma)`.
///
/// # Errors
/// [`AzothError::InvalidInput`] where the shapes disagree, and [`AzothError::OutOfRange`]
/// for a tolerance or iteration cap the spec bounds.
#[allow(clippy::too_many_arguments)] // the element matrix, its two vectors, the feed, the potentials, T
pub fn chemical_equilibrium(
    a_matrix: &[Vec<f64>],
    b: &[f64],
    whole_system: bool,
    moles: &[f64],
    chem_ref: &[f64],
    log_activity: &[f64],
    temperature: f64,
    max_iterations: u32,
    tolerance: f64,
    concentration_basis: ConcentrationBasis,
    solvent_weight: f64,
    solvent_mask: &[f64],
    phase_moles: f64,
) -> Result<ChemicalEquilibriumResult> {
    let spec = &model_gen::CHEMICAL_EQUILIBRIUM_SPEC;
    let mut warnings = Vec::new();

    let n_elements = b.len();
    let species = moles.len();

    // **One call over every declared input.** Calling this once per name with a
    // name-filtering closure emits a warning for every check the closure does not feed,
    // which is how this first read four warnings on a state that violates none.
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(temperature),
            "tolerance" => Some(tolerance),
            "max_iterations" => Some(f64::from(max_iterations)),
            _ => None,
        },
        &mut warnings,
    )?;

    if a_matrix.len() != n_elements {
        return Err(AzothError::InvalidInput {
            field: "a_matrix".to_string(),
            reason: format!(
                "the element matrix has {} row(s) against {} element(s)",
                a_matrix.len(),
                n_elements
            ),
        });
    }
    for (row, values) in a_matrix.iter().enumerate() {
        if values.len() != species {
            return Err(AzothError::InvalidInput {
                field: "a_matrix".to_string(),
                reason: format!(
                    "row {row} has {} entries against {species} component(s)",
                    values.len()
                ),
            });
        }
    }
    if chem_ref.len() != species || log_activity.len() != species {
        return Err(AzothError::InvalidInput {
            field: "chem_ref".to_string(),
            reason: format!(
                "{} reference potential(s) and {} activity coefficient(s) against {species} \
                 component(s)",
                chem_ref.len(),
                log_activity.len()
            ),
        });
    }

    // **The mask is checked where it is read.** On the mole-fraction basis the denominator
    // is `n_t` for every component whatever the mask says, and the operation's own call -
    // which cannot state the molality basis - passes the split it never reads.
    if concentration_basis == ConcentrationBasis::SoluteMolality && solvent_mask.len() != species {
        return Err(AzothError::InvalidInput {
            field: "solvent_mask".to_string(),
            reason: format!(
                "{} mask entr(ies) against {species} component(s), and the solute-molality \
                 basis reads one per component",
                solvent_mask.len()
            ),
        });
    }
    // **`n_t`, read once.** NeqSim takes the phase's total moles when the solver is
    // constructed and divides by that for the whole solve - `ChemicalEquilibrium.java:199`
    // - so it is not the running sum of the reactive set, and it is not the sum of `moles`
    // either unless the phase holds nothing else.
    let n_t = phase_moles.max(MIN_MOLES);
    let activity_basis = ActivityBasis {
        basis: concentration_basis,
        solvent_weight,
        solvent_mask,
    };

    // **Two compositions, because NeqSim has two.** `committed` is what the phase holds
    // and `n_mol` is the trial this pass computed; a pass whose error did not improve is
    // *not* written back, so the next pass re-derives from the same committed state
    // rather than from the trial it just rejected. Collapsing the two into one array
    // makes the solve wander instead of converge, which is what the first version of this
    // did.
    let mut committed = moles.to_vec();
    let mut n_mol = moles.to_vec();
    // `error` starts above any threshold so the first pass always commits; `err_old`
    // is assigned by the loop's first statement and has no initial value to read.
    let mut error = 1.0e10;
    let mut err_old: f64;
    let mut max_error = tolerance;
    let mut iterations: u32 = 0;
    let mut stagnation: u32 = 0;
    let mut best_error = f64::MAX;

    loop {
        iterations += 1;
        err_old = error;
        error = 0.0;

        // The trial starts from what the phase holds, not from the last trial.
        n_mol.copy_from_slice(&committed);

        let (dn, a_lambda) = chem_solve(
            a_matrix,
            b,
            whole_system,
            &committed,
            chem_ref,
            log_activity,
            &activity_basis,
            n_t,
        )?;
        let step = step_of(
            &committed,
            chem_ref,
            log_activity,
            &dn,
            &a_lambda,
            temperature,
            &activity_basis,
            n_t,
        );

        // The error and the trial in one pass, as NeqSim takes them: a species whose step
        // is under `1e-15` of its own moles contributes nothing and does not move.
        for i in 0..species {
            if committed[i] < MIN_MOLES {
                continue;
            }
            if !dn[i].is_finite() {
                error = f64::NAN;
                break;
            }
            if (dn[i] / committed[i]).abs() > 1e-15 {
                error += (dn[i] / committed[i]).abs();
                n_mol[i] = dn[i] * step + committed[i];
            }
        }

        if !error.is_finite() {
            break;
        }

        if error < best_error {
            best_error = error;
            stagnation = 0;
        } else {
            stagnation += 1;
        }
        if stagnation >= STAGNATION_LIMIT {
            break;
        }

        // `updateMoles()`: the trial is written back to the phase only when the error did
        // not get worse. A rejected pass leaves the phase where it was.
        if error <= err_old {
            committed.copy_from_slice(&n_mol);
        }

        // The loop condition, as NeqSim writes it: at least two passes, and then while
        // both the previous and the current error exceed the tolerance.
        let continuing =
            (err_old > max_error && error.abs() > max_error && iterations < max_iterations)
                || iterations < 2;
        if !continuing {
            break;
        }

        // The tolerance is relaxed on a solve that is taking a long time, which is
        // NeqSim's own behaviour and the reason a long solve can report convergence.
        if iterations > 15 {
            max_error *= 1.5;
        }
    }

    let converged = error.is_finite() && error < max_error;

    Ok(ChemicalEquilibriumResult {
        moles: committed,
        iterations,
        error,
        converged,
        warnings,
    })
}

/// One `chemSolve`: the two Lagrange systems and the Newton direction.
#[allow(clippy::too_many_arguments)] // the trial composition, the potentials and the activity basis
fn chem_solve(
    a_matrix: &[Vec<f64>],
    b: &[f64],
    whole_system: bool,
    n_mol: &[f64],
    chem_ref: &[f64],
    log_activity: &[f64],
    activity_basis: &ActivityBasis<'_>,
    n_t: f64,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let n_elements = b.len();
    let species = n_mol.len();
    let log_denominator = activity_basis.log_denominators(n_t, species);

    let mut m = vec![vec![0.0; species]; species];
    let mut chem_pot = vec![0.0; species];
    for i in 0..species {
        let n_i = n_mol[i].max(MIN_MOLES);
        m[i][i] = 1.0 / n_i;
        chem_pot[i] = chem_ref[i] + n_i.ln() - log_denominator[i] + log_activity[i];
    }

    // `M^-1 A^T`, then `A M^-1 A^T`.
    let mut a_transpose = vec![vec![0.0; n_elements]; species];
    for e in 0..n_elements {
        for i in 0..species {
            a_transpose[i][e] = a_matrix[e][i];
        }
    }
    let m_inv_at = solve_columns(&m, &a_transpose)?;

    let mut ama = vec![vec![0.0; n_elements]; n_elements];
    for e in 0..n_elements {
        for f in 0..n_elements {
            let mut sum = 0.0;
            for i in 0..species {
                sum += a_matrix[e][i] * m_inv_at[i][f];
            }
            ama[e][f] = sum;
        }
    }

    // `M^-1 mu`, then `A M^-1 mu`.
    let mu_column: Vec<Vec<f64>> = chem_pot.iter().map(|value| vec![*value]).collect();
    let m_inv_mu = solve_columns(&m, &mu_column)?;
    let mut amu = vec![0.0; n_elements];
    for e in 0..n_elements {
        let mut sum = 0.0;
        for i in 0..species {
            sum += a_matrix[e][i] * m_inv_mu[i][0];
        }
        amu[e] = sum;
    }

    // `sum(n_i mu_i)`, the total-moles row's right-hand side.
    let mut nmu = 0.0;
    for i in 0..species {
        nmu += n_mol[i] * chem_pot[i];
    }

    // The conservation coupling: `b` unless the phase's own element amounts disagree
    // with it, in which case the amounts are what is conserved and `b - A n` is carried.
    //
    // **Gated on the phase being the whole system**, which is NeqSim's
    // `system.getNumberOfPhases() == 1`. The gate is not a detail: on the aqueous phase of
    // a flashed CO2-water fluid the correction fires on the charge row and drives the
    // solve away from the fixed point it otherwise reaches in the same nineteen passes.
    let mut coupling = b.to_vec();
    let mut correction = vec![0.0; n_elements];
    if whole_system {
        for e in 0..n_elements {
            let mut current = 0.0;
            for i in 0..species {
                current += a_matrix[e][i] * n_mol[i];
            }
            if (b[e] - current).abs() > CONSERVATION_CORRECTION_TOLERANCE * b[e].abs().max(1.0) {
                coupling = (0..n_elements)
                    .map(|f| {
                        let mut sum = 0.0;
                        for i in 0..species {
                            sum += a_matrix[f][i] * n_mol[i];
                        }
                        sum
                    })
                    .collect();
                correction = (0..n_elements).map(|f| b[f] - coupling[f]).collect();
                break;
            }
        }
    }

    // The bordered system: `[[AMA, c^T], [c, 0]] [lambda, tau] = [AMU + correction, nmu]`.
    let size = n_elements + 1;
    let mut larger = vec![vec![0.0; size]; size];
    for e in 0..n_elements {
        for f in 0..n_elements {
            larger[e][f] = ama[e][f];
        }
        larger[e][n_elements] = coupling[e];
        larger[n_elements][e] = coupling[e];
    }
    let mut rhs = vec![0.0; size];
    for e in 0..n_elements {
        rhs[e] = amu[e] + correction[e];
    }
    rhs[n_elements] = nmu;

    let solved = solve_lu(&larger, &rhs)?;
    let tau = solved[n_elements];

    let mut dn = vec![0.0; species];
    // **`A^T lambda` is carried out because the step needs it**: NeqSim's `step()` reads
    // the same product back off `x_solve` to weigh the two Gibbs measures.
    let mut a_lambda = vec![0.0; species];
    let mut lambda_rhs = vec![vec![0.0]; species];
    for i in 0..species {
        let mut sum = 0.0;
        for e in 0..n_elements {
            sum += a_matrix[e][i] * solved[e];
        }
        a_lambda[i] = sum;
        lambda_rhs[i][0] = sum - chem_pot[i];
    }
    let direction = solve_columns(&m, &lambda_rhs)?;
    for i in 0..species {
        dn[i] = direction[i][0] + n_mol[i] * tau;
    }
    Ok((dn, a_lambda))
}

/// `step()`: the damped step length, from the two Gibbs measures.
///
/// **The negative-moles branch is not reproduced.** NeqSim sends a `n_omega` that has
/// gone negative to `innerStep`, which bisects the step back toward the current state;
/// that path is a refusal here, because it is reached only by a step the damping below
/// was supposed to prevent and no captured state takes it.
#[allow(clippy::too_many_arguments)] // the trial composition, the potentials and the activity basis
fn step_of(
    n_mol: &[f64],
    chem_ref: &[f64],
    log_activity: &[f64],
    dn: &[f64],
    a_lambda: &[f64],
    temperature: f64,
    activity_basis: &ActivityBasis<'_>,
    n_t: f64,
) -> f64 {
    let species = n_mol.len();
    let log_denominator = activity_basis.log_denominators(n_t, species);

    // **`R T` is not decoration.** `step()` recomputes its own `chem_pot` as
    // `R T (mu_ref + ln(n_i/n_t) + ln gamma_i)` - the *dimensional* potential - while the
    // `A^T lambda` it subtracts from it came out of `chemSolve`'s **reduced** solve. The
    // two arrays are in different units, so the factor does not cancel between `G_0` and
    // `G_1`, and dropping it drives the step to zero. That is an inconsistency in NeqSim
    // rather than a convention, and it is reproduced rather than corrected.
    let r_t = GAS_CONSTANT * temperature;
    let n_omega: Vec<f64> = (0..species).map(|i| n_mol[i] + dn[i]).collect();

    // **The negative-moles branch, which the first draft of this assumed was
    // unreachable and is reached on the very first pass.** A species whose trial goes
    // below zero cannot be logged, and NeqSim does not floor it - it takes the largest
    // step that keeps every species non-negative, three percent short of the boundary,
    // and returns it without computing either Gibbs measure. Flooring instead puts
    // `1/MIN_MOLES = 1e60` into `G_1` and drives the step to zero, which is how this was
    // found.
    if let Some(first) = n_omega.iter().position(|value| *value < 0.0) {
        let mut shortest = (-n_mol[first] / dn[first]) * 0.97;
        for i in first..species {
            if n_mol[i] + dn[i] < 0.0 {
                let candidate = (-n_mol[i] / dn[i]) * 0.97;
                if candidate < shortest {
                    shortest = candidate;
                }
            }
        }
        return if shortest > 1.0 { 1.0 } else { shortest };
    }

    let potential = |value: f64, index: usize| {
        r_t * (chem_ref[index] + value.max(MIN_MOLES).ln() - log_denominator[index]
            + log_activity[index])
    };

    // **The derivative does not follow the residual onto the molality basis.** These two
    // `1/n` terms are the derivative of the mole-fraction form, and NeqSim's `M` matrix is
    // the same expression whatever basis the system selects - so a molality-basis solve
    // damps with a gradient its own residual is not the integral of. Reproduced, because
    // the alternative is a step rule NeqSim does not have.
    let mut g_1 = 0.0;
    for i in 0..species {
        g_1 += (potential(n_omega[i], i) - a_lambda[i])
            * dn[i]
            * (1.0 / n_omega[i].max(MIN_MOLES) - 1.0 / n_t);
    }

    let mut step = 1.0;
    if g_1 > 0.0 {
        let mut g_0 = 0.0;
        for i in 0..species {
            g_0 += (potential(n_mol[i], i) - a_lambda[i])
                * dn[i]
                * (1.0 / n_mol[i].max(MIN_MOLES) - 1.0 / n_t);
        }
        let denominator = g_0 - g_1;
        if denominator.abs() > 1e-30 {
            step = g_0 / denominator;
        }
    }

    // Clamped to [0, 1]: a negative or overshooting step diverges and breaks element
    // conservation, so it is replaced by a full step.
    if !(0.0..=1.0).contains(&step) || !step.is_finite() {
        step = 1.0;
    }
    step
}

/// Solve `m x = rhs` for a square `m` and several right-hand sides.
fn solve_columns(m: &[Vec<f64>], rhs: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
    let columns = rhs.first().map_or(0, Vec::len);
    let mut out = vec![vec![0.0; columns]; rhs.len()];
    for c in 0..columns {
        let b: Vec<f64> = rhs.iter().map(|row| row[c]).collect();
        let solved = solve_lu(m, &b)?;
        for (i, value) in solved.iter().enumerate() {
            out[i][c] = *value;
        }
    }
    Ok(out)
}
