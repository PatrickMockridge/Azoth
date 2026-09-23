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
//! # What is ported here and what is not, yet
//!
//! Ported: the driver's **single-phase branch** (`solveSinglePhaseChemicalEquilibrium`, in
//! [`single_phase_equilibrium`]), its **VLE initialisation** ([`vle_initialization`]) and its
//! **Gibbs measure** over a phase list ([`total_gibbs_energy`]).
//!
//! **Not ported**: the multiphase modified-RAND solve the outer loop drives, and with it
//! `addTrialPhases`, `removeNegligiblePhases`, the outer loop's acceptance rules and the
//! trace-ion short circuit. No capture reaches the add or the remove: no fluid this probe
//! drives has an unstable trial seed, and every beta stays far above the `1e-12` removal
//! floor.

use azoth_core::Result;

use crate::rand_solver::{RandSolution, solve};
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

/// `Phase.normalize()`: divide by the sum, leaving a vector that is already zero alone.
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
    let solution = solve(a_matrix, g0, b, feed_moles, ln_phi)?;
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
