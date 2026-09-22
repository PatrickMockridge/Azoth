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
//! # What is ported here and what is not, yet
//!
//! This is the driver's **single-phase branch**: the path a fluid takes when the stability
//! analysis finds no second phase, which is `solveSinglePhaseChemicalEquilibrium` - drop any
//! phase but the first, solve, and report. The captured water-gas shift is exactly that case,
//! and the driver's own numbers for it are pinned: the equilibrium moles, the total moles and
//! the Gibbs measure the driver computes.
//!
//! **The driver's own Gibbs measure counts a phase that is not there.** `computeGibbsEnergy`
//! sums over the phase *list*, and the captured system still holds two phase objects with the
//! same composition and `beta = 1.0` each - the trial phase the removal step never dropped - so
//! the number it reports is **twice** the thermodynamic value. That is a measurement, not a
//! rounding difference: the capture's `final_gibbs_energy` is `-2.2595355` against this port's
//! single-phase `-1.1297676`, and the test asserts the ratio.
//!
//! **Not ported**: the VLE initialisation, `addTrialPhases`, `removeNegligiblePhases`, the
//! outer loop's acceptance rules, and the trace-ion short circuit. Those are the *multiphase*
//! machinery, and the capture says this fluid does not reach them - the stability analysis
//! returns stable, so no trial phase is added. What the captured flash does leave behind is
//! two phase objects with identical compositions and `beta = 1.0`, which is that machinery's
//! bookkeeping rather than a chemical or phase-split result.

use azoth_core::Result;

use crate::rand_solver::{RandSolution, solve};

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
/// NeqSim loops over the phases and weights each by its `beta`; with one phase that weight is
/// one, and the sum is `G/RT` per mole of mixture - **dimensionless**, which is why the
/// capture's number is `-2.26` and not a number of joules.
#[must_use]
pub fn gibbs_energy(fractions: &[f64], ln_phi: &[f64]) -> f64 {
    fractions
        .iter()
        .zip(ln_phi)
        .map(|(x, phi)| x * (x.ln() + phi))
        .sum()
}
