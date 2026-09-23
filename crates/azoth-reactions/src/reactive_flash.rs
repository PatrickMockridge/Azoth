//! The reactive flash driver, from
//! `flashops/reactiveflash/ReactiveMultiphaseTPflash.java`.
//!
//! `run()` is a sequence of decisions around the two solvers this crate already carries: it
//! builds the formula matrix, **short-circuits at `NR = 0`** (no independent reaction means
//! the element balance fully determines the composition), initialises a phase split with a
//! conventional VLE flash when more than one phase is allowed, asks the stability analysis
//! whether a second phase forms, and then drives the modified-RAND solver in an outer loop
//! that removes the phases that turn out negligible.
//!
//! # The phase list the driver starts from
//!
//! **`SystemSrkEos` constructs with two phase objects**, each holding the whole feed at
//! `beta = 1.0`, and `init(0)` restores that count - so `setNumberOfPhases(1)` reaches only
//! as far as the next init. Since the captured systems are built, initialised and handed to
//! the driver, the driver sees `np = 2` and `maxPhases = 2` and its `skipStability` rule
//! (`np >= 2 && maxPhases >= 2`) takes it into the outer loop **without any trial phase ever
//! being added**. The two phases the capture ends with are that pair, not a stability result.
//!
//! # The weight `computeGibbsEnergy` gives each phase is stale
//!
//! `updateSystem` writes the solver's converged fractions into each phase
//! (`ph.setBeta(beta[j])`) and then calls `system.init(1)` **inside the same method**, which
//! re-initialises every phase from the *system's* beta array - and that array is refreshed
//! only on the ionic branch (`system.setBeta(j, beta[j])` sits behind `hasIonicSpecies`). So
//! on a neutral fluid the converged split never reaches the system, and `computeGibbsEnergy`,
//! which weighs each phase by `phase.getBeta()`, uses whatever the last system-level write
//! left: `(1 - V, V)` from the VLE initialisation, or `(1.0, 1.0)` from the constructor.
//!
//! The capture measures both. The fluids the class's own tests use are never VLE-initialised,
//! so their two phases carry `beta = 1.0` each and the reported Gibbs energy is **twice** the
//! thermodynamic value - `-2.2595355` against this crate's single-phase `-1.1297676`. Forced
//! to one phase, the same shift at 300 K takes the other path: the VLE initialisation adds a
//! phase, the betas come out normalised `(0.2106, 0.7894)`, and the reported energy is
//! `-0.7162112` - one phase's worth. The driver's phase split, meanwhile, is
//! `(0.014138, 0.985862)` by the solver's own count, and nothing reads it.
//!
//! # One defect, three consequences
//!
//! The stale fraction is not only the Gibbs measure's problem. `removeNegligiblePhases` tests
//! `system.getBeta(j)` - the same stale array - so a phase the *solve* has driven to nothing is
//! **not removed**, and one the stale array calls small would be. And the next outer iteration's
//! `initialize` re-reads `phase.getBeta()`, so the fraction a solve starts from never moves
//! either, while the composition does. A port that read the solve's own `phase_amounts` in the
//! removal step would be a different algorithm from the class's.
//!
//! `removePhaseKeepTotalComposition` is the removal itself, and its name describes an intent it
//! does not carry out: it shifts the phase *index* array and decrements the count, touching no
//! moles at all.
//!
//! # What is ported here and what is not, yet
//!
//! Ported: the driver's **single-phase branch** (`solveSinglePhaseChemicalEquilibrium`, in
//! [`single_phase_equilibrium`]), its **VLE initialisation** ([`vle_initialization`]), its
//! **Gibbs measure** over a phase list ([`total_gibbs_energy`]), its **outer loop**
//! ([`outer_loop`]) with its acceptance rules, [`remove_negligible_phases`] and
//! [`add_trial_phase`]. The multiphase solve the loop drives is
//! [`crate::rand_solver::solve`].
//!
//! **Not ported**: the trace-ion short circuit (it needs `hasIonicSpecies`, which this crate's
//! callers refuse before reaching here), `runNonReactiveFlash` - the `NR = 0` multiphase
//! fallback to a conventional flash - and the reactive-stability *orchestrator* the loop's
//! recheck calls, which is [`crate::reactive_stability`] driven from outside.
//!
//! **No capture reaches the add.** The pair `SystemSrkEos` constructs means the stability
//! analysis is skipped on every captured state, so no trial is ever produced; reaching it needs
//! a fluid that is unstable *at one phase*, with the phase list cut back to one after the last
//! `init`. The removal *is* reached, and returns false.

use azoth_core::{AzothError, Result};

use crate::rand_solver::{PhaseFeed, PhaseLogPhi, RandSolution, solve_single_phase};
use crate::reactive_stability::CriticalConstants;

/// The floor `initializeWithVLEFlash` keeps a trial mole fraction above, which is its own
/// `1e-30` and not [`MIN_MOLES`]'s `1e-30` in the stability analysis - the same number, a
/// different constant in the source.
pub const VLE_MIN_MOLES: f64 = 1.0e-30;

/// The bound `initializeWithVLEFlash` rejects a Rachford-Rice answer at, from
/// `V < 1e-10 || V > 1.0 - 1e-10`.
pub const VLE_SINGLE_PHASE_TOLERANCE: f64 = 1.0e-10;

/// The `|ln K|` above which a component counts as volatile, from `hasVolatile`.
pub const VLE_VOLATILE_LOG_K: f64 = 0.1;

/// What `initializeWithVLEFlash` leaves behind: the vapour fraction it solved for and the two
/// normalised compositions.
#[derive(Debug, Clone, PartialEq)]
pub struct VleInitialisation {
    /// The Rachford-Rice root, in `(0, 1)` - the driver's `V` and the beta of the vapour.
    pub vapour_fraction: f64,
    /// The liquid composition `z_i / (1 + V (K_i - 1))`, floored and normalised.
    pub liquid: Vec<f64>,
    /// The vapour composition `K_i x_i`, floored and normalised.
    pub vapour: Vec<f64>,
}

/// `initializeWithVLEFlash`: Wilson K-values and a Rachford-Rice solve on the **overall**
/// composition.
///
/// `None` is the class's early return, which means one of two things and reports both the
/// same way: no component is volatile (`|ln K| <= 0.1` throughout), or the Rachford-Rice root
/// came back at a bound - `V <= 1e-10` or `V >= 1 - 1e-10` - so there is no split to
/// initialise. The driver then adds no phase, and whether the system holds one phase or two
/// afterwards is the phase list's business and not this function's.
///
/// The composition is `getz()`, the *overall* one, not the phase's `x` - the driver reads it
/// per component at `system.getPhase(0).getComponent(i).getz()`.
///
/// The Rachford-Rice loop is the class's own: 100 passes, a Newton step clamped to `[0, 1]`
/// after each, a `1e-30` floor on the denominator that *skips* the component rather than
/// guarding it, and a `1e-10` bound on the residual.
#[must_use]
#[allow(clippy::manual_range_contains)] // the class's own `V < 1e-10 || V > 1.0 - 1e-10`
pub fn vle_initialization(
    fractions: &[f64],
    constants: &[CriticalConstants],
    temperature: f64,
    pressure: f64,
) -> Option<VleInitialisation> {
    let nc = fractions.len();
    let mut log_k = vec![0.0_f64; nc];
    let mut has_volatile = false;
    for i in 0..nc {
        let entry = constants[i];
        if entry.pc > 0.0 && entry.tc > 0.0 {
            log_k[i] = (entry.pc / pressure).ln()
                + 5.373 * (1.0 + entry.omega) * (1.0 - entry.tc / temperature);
            if log_k[i].abs() > VLE_VOLATILE_LOG_K {
                has_volatile = true;
            }
        }
    }
    if !has_volatile {
        return None;
    }

    let mut vapour = 0.5_f64;
    for _ in 0..100 {
        let mut f = 0.0_f64;
        let mut derivative = 0.0_f64;
        for i in 0..nc {
            let k = log_k[i].exp();
            let denominator = 1.0 + vapour * (k - 1.0);
            if denominator.abs() < 1e-30 {
                continue;
            }
            f += fractions[i] * (k - 1.0) / denominator;
            derivative -= fractions[i] * (k - 1.0) * (k - 1.0) / (denominator * denominator);
        }
        if f.abs() < 1e-10 {
            break;
        }
        if derivative.abs() > 1e-30 {
            vapour -= f / derivative;
        }
        vapour = vapour.clamp(0.0, 1.0);
    }

    if vapour < VLE_SINGLE_PHASE_TOLERANCE || vapour > 1.0 - VLE_SINGLE_PHASE_TOLERANCE {
        return None;
    }

    let unfloored: Vec<f64> = (0..nc)
        .map(|i| fractions[i] / (1.0 + vapour * (log_k[i].exp() - 1.0)))
        .collect();
    let liquid: Vec<f64> = unfloored.iter().map(|x| x.max(VLE_MIN_MOLES)).collect();
    let vapour_composition: Vec<f64> = unfloored
        .iter()
        .enumerate()
        .map(|(i, x)| (log_k[i].exp() * x).max(VLE_MIN_MOLES))
        .collect();

    Some(VleInitialisation {
        vapour_fraction: vapour,
        liquid: normalise(&liquid),
        vapour: normalise(&vapour_composition),
    })
}

/// A stability analysis: a phase list in, the **unstable trial compositions** out - empty where
/// the analysis found nothing. That is `ReactiveStabilityAnalysis.run` and
/// `getUnstableTrialCompositions` in one call, because the driver's two branches on them both
/// end its loop.
pub type StabilityCheck<'a> = &'a mut dyn FnMut(&[PhaseFeed]) -> Result<Vec<Vec<f64>>>;

/// The outer loop's pass cap, from `MAX_OUTER_ITER`.
pub const MAX_OUTER_ITERATIONS: usize = 20;

/// The phase fraction below which `removeNegligiblePhases` drops a phase, from
/// `MIN_PHASE_FRACTION`.
pub const MIN_PHASE_FRACTION: f64 = 1.0e-12;

/// The fraction `addTrialPhases` gives a new phase, from `initialBeta`.
pub const TRIAL_PHASE_BETA: f64 = 0.01;

/// The residual `run` accepts a *non-converged* multiphase solve at, from `5.0e-3`.
pub const NEAR_CONVERGED_RESIDUAL: f64 = 5.0e-3;

fn normalise(values: &[f64]) -> Vec<f64> {
    let total: f64 = values.iter().sum();
    if total <= 0.0 {
        return values.to_vec();
    }
    values.iter().map(|value| value / total).collect()
}

/// What the driver reports for a fluid it finds stable in one phase.
#[derive(Debug, Clone, PartialEq)]
pub struct SinglePhaseOutcome {
    /// The equilibrium moles, one per component.
    pub moles: Vec<f64>,
    /// Their sum, which is the driver's `getTotalMoles` - **and not one**, because the shift
    /// changes the number of moles: `CO + H2O = CO2 + H2` leaves it alone, but a fluid where
    /// one species splits into two does not.
    pub total_moles: f64,
    /// `computeGibbsEnergy`: `sum_i x_i (ln x_i + ln phi_i)` at the answer, which is `G/RT`
    /// for one mole of mixture rather than a quantity in joules.
    pub gibbs_energy: f64,
    /// Passes the solve took.
    pub iterations: u32,
    /// Whether it converged.
    pub converged: bool,
    /// The solve itself, for the multipliers and the residuals.
    pub solution: RandSolution,
}

/// `solveSinglePhaseChemicalEquilibrium`: the driver's answer for a fluid that is stable in
/// one phase.
///
/// The element inventory `b` is the *frozen* one the driver captured from the feed
/// (`setElementBalance`), not one recomputed as the composition moves - which is what keeps
/// the balance a constraint on the whole fluid rather than on the current iterate.
///
/// **This branch never computes the Gibbs energy.** It returns before the driver's
/// `computeGibbsEnergy` call, so `getFinalGibbsEnergy()` reports its field's initial `0.0` for
/// every fluid that takes it - which the capture's forced-one-phase 600 K block shows.
///
/// # Errors
/// Whatever [`solve`] raises, and whatever `ln_phi` raises.
pub fn single_phase_equilibrium(
    a_matrix: &[Vec<f64>],
    g0: &[f64],
    b: &[f64],
    feed_moles: &[f64],
    ln_phi: &mut dyn FnMut(&[f64]) -> Result<Vec<f64>>,
) -> Result<SinglePhaseOutcome> {
    let solution = solve_single_phase(a_matrix, g0, b, feed_moles, ln_phi)?;
    let total_moles: f64 = solution.moles.iter().sum();
    let fractions: Vec<f64> = solution
        .moles
        .iter()
        .map(|moles| moles / total_moles)
        .collect();
    let ln_phi_here = ln_phi(&fractions)?;
    let gibbs_energy = gibbs_energy(&fractions, &ln_phi_here);

    Ok(SinglePhaseOutcome {
        total_moles,
        gibbs_energy,
        iterations: solution.iterations,
        converged: solution.converged,
        moles: solution.moles.clone(),
        solution,
    })
}

/// `computeGibbsEnergy`: `sum_i x_i (ln x_i + ln phi_i)` over one phase.
///
/// NeqSim sums this over the phases and weights each by its `beta`, so the sum is `G/RT` per
/// mole of mixture - **dimensionless**, which is why the capture's number is `-2.26` and not
/// a number of joules - and only when the weights sum to one. The constructor's two phases
/// both carry `beta = 1.0`, which is what doubles the captured number; see
/// [`total_gibbs_energy`].
#[must_use]
pub fn gibbs_energy(fractions: &[f64], ln_phi: &[f64]) -> f64 {
    fractions
        .iter()
        .zip(ln_phi)
        .map(|(x, phi)| x * (x.ln() + phi))
        .sum()
}

/// One phase as `computeGibbsEnergy` reads it: the weight it is given, its composition and
/// its fugacity coefficients.
#[derive(Debug, Clone, PartialEq)]
pub struct PhaseGibbs {
    /// `phase.getBeta()` - the phase object's own number, **not** the solver's fraction.
    pub beta: f64,
    /// The phase's mole fractions.
    pub fractions: Vec<f64>,
    /// Each component's `ln phi_i` in that phase.
    pub ln_phi: Vec<f64>,
}

/// `computeGibbsEnergy`: the beta-weighted sum over the driver's phase list.
///
/// **The weights are the phase objects' own betas, and on a neutral fluid those are stale**:
/// the solver's converged fractions are written into the phases and then overwritten by
/// `system.init(1)` from a system array the RAND solver never refreshes. What the caller
/// passes here is what the driver reads, which is the point of the function - a port that
/// substituted the solver's fractions would reproduce a number NeqSim does not report.
#[must_use]
pub fn total_gibbs_energy(phases: &[PhaseGibbs]) -> f64 {
    phases
        .iter()
        .map(|phase| phase.beta * gibbs_energy(&phase.fractions, &phase.ln_phi))
        .sum()
}

/// `normalizeBeta`: divide every phase's fraction by their sum.
///
/// # Errors
/// [`AzothError::InvalidInput`] where the fractions sum to zero, which the class does not
/// guard against and would divide by.
pub fn normalise_betas(phases: &mut [PhaseFeed]) -> Result<()> {
    let total: f64 = phases.iter().map(|phase| phase.beta).sum();
    if total <= 0.0 {
        return Err(AzothError::InvalidInput {
            field: "phases".to_string(),
            reason: "the phase fractions sum to zero".to_string(),
        });
    }
    for phase in phases.iter_mut() {
        phase.beta /= total;
    }
    Ok(())
}

/// `removeNegligiblePhases`: drop every phase whose fraction is under
/// [`MIN_PHASE_FRACTION`], from the last one backwards, never dropping the last phase there is.
///
/// **The fraction it tests is the one the *system* holds, and on a neutral fluid that is the
/// stale array** - `ModifiedRANDSolver.updateSystem` writes the solver's fractions into the
/// phase objects and then `system.init(1)` restores them from a system array only refreshed on
/// the ionic branch. So a phase the *solve* has driven to nothing is not removed by this test,
/// and one the stale array calls small is. Reading the solve's own `phase_amounts` here would
/// be a different algorithm from the class's.
///
/// `removePhaseKeepTotalComposition` is the class's removal and it **touches no moles**: it
/// shifts the phase *index* array and decrements the count, so the removed phase's moles simply
/// leave the list and the next solve's element balance redistributes what is left. The name
/// describes an intent the method does not carry out.
#[must_use]
pub fn remove_negligible_phases(phases: &mut Vec<PhaseFeed>) -> bool {
    let mut removed = false;
    let mut j = phases.len();
    while j > 0 {
        j -= 1;
        if phases.len() <= 1 {
            break;
        }
        if phases[j].beta < MIN_PHASE_FRACTION {
            phases.remove(j);
            removed = true;
        }
    }
    if removed {
        // The class calls `normalizeBeta` and `init(1)` here; the fractions are already
        // normalised ones and the betas are scaled, so only the scaling is left to do - and it
        // cannot fail on a list that still holds a phase.
        let _ = normalise_betas(phases);
    }
    removed
}

/// `addTrialPhases`: append **one** phase - the first of the trials, which the stability
/// analysis sorted most-unstable-first - at [`TRIAL_PHASE_BETA`], then renormalise every beta.
///
/// The class's guard is on `system.getMaxNumberOfPhases()`, which `init` may have reset, and
/// not on the effective maximum the driver resolved at entry. Returns whether a phase was
/// added.
///
/// # Errors
/// [`AzothError::InvalidInput`] where a trial's length disagrees with the phases'.
pub fn add_trial_phase(
    phases: &mut Vec<PhaseFeed>,
    trials: &[Vec<f64>],
    max_phases: usize,
) -> Result<bool> {
    let Some(trial) = trials.first() else {
        return Ok(false);
    };
    if phases.len() >= max_phases {
        return Ok(false);
    }
    if let Some(existing) = phases.first() {
        if trial.len() != existing.fractions.len() {
            return Err(AzothError::InvalidInput {
                field: "trials".to_string(),
                reason: format!(
                    "a trial has {} entries against {} component(s)",
                    trial.len(),
                    existing.fractions.len()
                ),
            });
        }
    }
    phases.push(PhaseFeed {
        fractions: normalise(
            &trial
                .iter()
                .map(|x| x.max(VLE_MIN_MOLES))
                .collect::<Vec<f64>>(),
        ),
        beta: TRIAL_PHASE_BETA,
    });
    normalise_betas(phases)?;
    Ok(true)
}

/// What the driver's outer loop ends with.
#[derive(Debug, Clone, PartialEq)]
pub struct OuterLoopOutcome {
    /// The phase list the loop stopped on, with each phase's composition as the last solve
    /// left it and its fraction as the *system* holds it.
    pub phases: Vec<PhaseFeed>,
    /// Whether the loop reached one of the class's convergence branches.
    pub converged: bool,
    /// `getTotalIterations`: every solve's passes, summed over the outer iterations.
    pub total_iterations: u32,
    /// `getEquilibriumTotalMoles`: the last solve's total.
    pub equilibrium_total_moles: f64,
    /// The last solve, for its multipliers, residuals and DIIS count.
    pub solution: RandSolution,
}

/// `run`'s outer loop, from the point the driver has a phase list to iterate on.
///
/// Each pass solves, removes what the removal step calls negligible, and then decides: a phase
/// list already at the maximum is done; otherwise the stability analysis is asked again, and
/// either it finds nothing, or its trials cannot be added, or a phase is added and the loop
/// goes round. A solve that did not converge is accepted only at
/// [`NEAR_CONVERGED_RESIDUAL`] - and the class's ionic guard on that acceptance is not ported,
/// because the callers of this crate refuse an ionic fluid before reaching here.
///
/// `stability` answers the *unstable trial compositions*, empty where the analysis found
/// nothing - which is `ReactiveStabilityAnalysis.run` and `getUnstableTrialCompositions` in
/// one call, because the class's two branches on them both end the loop.
///
/// # Errors
/// Whatever the solve raises, whatever `ln_phi` raises, and whatever `stability` raises.
#[allow(clippy::too_many_arguments)]
pub fn outer_loop(
    a_matrix: &[Vec<f64>],
    g0: &[f64],
    b: &[f64],
    total_moles: f64,
    initial: Vec<PhaseFeed>,
    max_phases: usize,
    ln_phi: PhaseLogPhi<'_>,
    stability: StabilityCheck<'_>,
) -> Result<OuterLoopOutcome> {
    let mut phases = initial;
    let mut converged = false;
    let mut total_iterations = 0_u32;
    let mut equilibrium_total_moles = total_moles;
    let mut outcome: Option<RandSolution> = None;

    for _ in 0..MAX_OUTER_ITERATIONS {
        let solution =
            crate::rand_solver::solve(a_matrix, g0, b, total_moles, &phases, &mut *ln_phi)?;
        total_iterations += solution.iterations;
        equilibrium_total_moles = solution.total_moles;
        let rand_converged = solution.converged;
        outcome = Some(solution);

        let solved = outcome.as_ref().expect("just assigned");
        if !rand_converged && phases.len() > 1 && solved.final_residual < NEAR_CONVERGED_RESIDUAL {
            converged = true;
            break;
        }

        let phase_removed = remove_negligible_phases(&mut phases);

        if rand_converged && !phase_removed {
            // A phase list already at the ceiling is accepted without asking the stability
            // analysis, which is what the class does to keep an immiscible fluid from
            // re-solving a single-phase state that is ill-posed for it.
            if phases.len() >= max_phases {
                converged = true;
                break;
            }
            let trials = stability(&phases)?;
            let phases_before = phases.len();
            // The class's two branches here - no trials returned, and a trial the phase
            // ceiling refused - both end the loop with the current equilibrium accepted, so
            // the count is the test rather than the return value.
            let _ = add_trial_phase(&mut phases, &trials, max_phases)?;
            if phases.len() == phases_before {
                converged = true;
                break;
            }
        }

        if rand_converged && phase_removed {
            converged = true;
            break;
        }
    }

    Ok(OuterLoopOutcome {
        phases,
        converged,
        total_iterations,
        equilibrium_total_moles,
        solution: outcome.expect("MAX_OUTER_ITERATIONS is 20, so the loop runs"),
    })
}
