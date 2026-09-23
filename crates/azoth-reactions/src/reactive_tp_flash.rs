//! The reactive flash as a model: `reactions.reactive_tp_flash`.
//!
//! Everything the driver does is in [`crate::reactive_flash`], which takes its phase model and
//! its reactive solve as closures. This module is what wires those closures to the databank and
//! a cubic, which is why this crate depends on `azoth-eos` after all: the driver needs a
//! fugacity coefficient at a trial composition, and the model is where the equation of state is
//! chosen.
//!
//! # The cubic, and the root every phase takes
//!
//! **SRK**, because the oracle's systems are `SystemSrkEos` and every captured number is an SRK
//! number. And **every phase takes the cubic's vapour root**, which is not a shortcut: NeqSim's
//! reactive flash clones phase 0 - a `gas` phase - for the new phase it adds, so both phases of
//! a VLE-initialised split are `gas`-typed there. The measurement is in the capture: rooting the
//! water-rich phase liquid collapses the 300 K split to a single phase, where the vapour root
//! reproduces it.
//!
//! # What the phases are called
//!
//! Nothing here names them. NeqSim's phase *types* come from `init(1)`'s assignment - `gas`,
//! `liquid`, `oil`, `aqueous` - and the capture records that the VLE-initialised liquid of the
//! 300 K state is `oil` while the constructor's pair are both `gas`. A type is a property of
//! NeqSim's system bookkeeping and not of the state this model computes, so the phases come
//! back in the driver's own order with their mole numbers and fractions and nothing else.
//!
//! # The one closure that is not the class's
//!
//! The stability analysis brings its reference and its trials to equilibrium through the
//! caller's solve, and NeqSim's own does that on a **clone of phase 0**, whose element
//! inventory is that phase's current `A n` rather than the driver's frozen one. A composition
//! is all this closure is handed, so it uses the driver's frozen inventory - and the two agree
//! on the first call, where phase 0 holds the feed, and differ where a recheck runs.

use azoth_core::units::{kelvins, pascals};
use azoth_core::warning::Warning;
use azoth_core::{AzothError, CalcResult, Result};
use azoth_eos::databank::mixture_of;
use azoth_eos::{Cubic, RootSide};

use crate::databank::{formation_properties, ionic_charge};
use crate::formula_matrix::FormulaMatrix;
use crate::rand_solver::{PhaseFeed, ThermoData, solve_single_phase, standard_potentials};
use crate::reactive_flash::{DriverState, run};
use crate::reactive_stability::CriticalConstants;

/// The cubic this model flashes with, from the oracle's `SystemSrkEos`.
pub const CUBIC: Cubic = Cubic::Srk;

/// What the flash answers with.
///
/// **Flat, because the transport layer's field names are**: a matrix and a vector rather than a
/// list of phase objects, so that `phase_moles[j][i]` crosses to Python as the spec's matrix
/// output and the three-way name agreement can be asserted.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactiveTpFlashResult {
    /// How many phases the driver stopped on.
    pub phase_count: usize,
    /// Each phase's moles, one row per phase and one column per component, in the driver's own
    /// phase order - **and no phase type is promised**, because a type is NeqSim's system
    /// bookkeeping and not a state this model computes.
    pub phase_moles: Vec<Vec<f64>>,
    /// Each phase's share of the fluid, the fraction the driver weighs it by.
    pub phase_fraction: Vec<f64>,
    /// `isConverged`. True on the branches that accept an answer without a converged solve,
    /// which is the class's own overwriting.
    pub converged: bool,
    /// Every solve's passes, summed. **Zero on the `NR = 0` fallback**, which the driver never
    /// counts.
    pub total_iterations: u32,
    /// The last solve's total moles, which one mole of feed does not fix: a reaction that
    /// splits one species into two moves it.
    pub equilibrium_total_moles: f64,
    /// `computeGibbsEnergy`, dimensionless - and weighed by the fractions the driver's phase
    /// list carries, which on a neutral fluid are the ones its own bookkeeping left.
    pub gibbs_energy: f64,
    /// `getFinalResidual`: the worst potential error against the element balance, which is what
    /// the solve stopped on. **The answer is only as tight as this**, and on a relaxed
    /// multiphase stop it is `1e-4` rather than `1e-9`.
    pub residual: f64,
    /// `getFinalElementResidual`, the scaled root-mean-square element deviation on its own.
    pub element_residual: f64,
    /// Caveats. This model has no range check of its own, so it is empty wherever the
    /// caller's state was accepted.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ReactiveTpFlashResult {
    const CALC_ID: &'static str = "reactions.reactive_tp_flash";
    const FIELDS: &'static [&'static str] = &[
        "phase_count",
        "phase_moles",
        "phase_fraction",
        "converged",
        "total_iterations",
        "equilibrium_total_moles",
        "gibbs_energy",
        "residual",
        "element_residual",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// `reactions.reactive_tp_flash`: simultaneous chemical and phase equilibrium at fixed
/// temperature and pressure.
///
/// `moles` is the **overall component amounts and not a composition**, which is what the
/// driver's `getOverallMoles` reads and what its frozen element inventory is built from.
///
/// # Errors
/// [`AzothError::InvalidInput`] for an ionic component (the RAND solver's ionic branches are
/// not ported, and this refuses rather than answering a neutral fluid's number), a length
/// disagreement, or a component the databank does not carry.
pub fn reactive_tp_flash(
    components: &[String],
    temperature: f64,
    pressure: f64,
    moles: &[f64],
    max_phases: usize,
) -> Result<ReactiveTpFlashResult> {
    if components.len() != moles.len() {
        return Err(AzothError::InvalidInput {
            field: "moles".to_string(),
            reason: format!(
                "{} entry(ies) against {} component(s)",
                moles.len(),
                components.len()
            ),
        });
    }
    if components.is_empty() || max_phases == 0 {
        return Err(AzothError::InvalidInput {
            field: "components".to_string(),
            reason: "a fluid with no component, or a phase ceiling of zero".to_string(),
        });
    }
    let total_moles: f64 = moles.iter().sum();
    if total_moles <= 0.0 {
        return Err(AzothError::InvalidInput {
            field: "moles".to_string(),
            reason: "the feed holds no moles".to_string(),
        });
    }

    // The ionic branch is refused rather than approximated: NeqSim pins a gas-phase ion to
    // `EPS` and corrects an electrolyte phase's reference state through
    // `getLogInfiniteDiluteFugacity`, and neither is ported.
    let mut charges = Vec::with_capacity(components.len());
    for component in components {
        let charge = ionic_charge(component)?.ok_or_else(|| AzothError::InvalidInput {
            field: "components".to_string(),
            reason: format!("the component databank has no row for `{component}`"),
        })?;
        charges.push(charge);
    }
    if let Some(index) = charges.iter().position(|charge| *charge != 0.0) {
        return Err(AzothError::InvalidInput {
            field: "components".to_string(),
            reason: format!(
                "`{}` is charged, and this model's phase model is a cubic: an ion's critical \
                 constants are the component table's filler rather than a measurement, and its \
                 reference state needs an electrolyte phase. The solver carries the ionic \
                 branch - see `rand_solver::IonicPhase` - and this model supplies none of it",
                components[index]
            ),
        });
    }

    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, ideal) = mixture_of(&names, CUBIC, None)?;

    let constants: Vec<CriticalConstants> = mixture
        .components()
        .iter()
        .map(|component| CriticalConstants {
            tc: component.tc.value,
            pc: component.pc.value / 1.0e5,
            omega: component.omega,
        })
        .collect();

    let matrix = FormulaMatrix::build(components)?;
    let b: Vec<f64> = matrix
        .matrix
        .iter()
        .map(|row| row.iter().zip(moles).map(|(a, n)| a * n).sum())
        .collect();

    let mut data = Vec::with_capacity(components.len());
    for (index, component) in components.iter().enumerate() {
        let formation =
            formation_properties(component)?.ok_or_else(|| AzothError::InvalidInput {
                field: "components".to_string(),
                reason: format!(
                    "the component databank carries no formation row for `{component}`"
                ),
            })?;
        data.push(ThermoData {
            enthalpy_of_formation: formation.enthalpy_of_formation,
            absolute_entropy: formation.absolute_entropy,
            gibbs_energy_of_formation: formation.gibbs_energy_of_formation,
            cp: [
                ideal.cp_a[index],
                ideal.cp_b[index],
                ideal.cp_c[index],
                ideal.cp_d[index],
                ideal.cp_e[index],
            ],
        });
    }

    // NeqSim's `getPressure()` is in bara, which is the unit its `ln(P/P_ref)` term is taken
    // against; the model's pressure is in Pa and is converted here and nowhere else.
    let pressure_bar = pressure / 1.0e5;
    // The ionic branch's three facts are not supplied: this model refuses a charged component
    // before it reaches here, so every potential takes the neutral branch. The solver carries
    // the branch - see `rand_solver::IonicPhase` - and what keeps it out of *this* model is the
    // cubic's own parameters, which are the databank's filler for an ion rather than a
    // measurement.
    let g0 = standard_potentials(&data, temperature, pressure_bar, &[], &[]);
    let reduced = mixture.reduced_parameters(kelvins(temperature), pascals(pressure))?;

    let feed_fractions: Vec<f64> = moles.iter().map(|moles| moles / total_moles).collect();
    let initial = vec![
        PhaseFeed {
            fractions: feed_fractions.clone(),
            beta: 1.0,
        },
        PhaseFeed {
            fractions: feed_fractions,
            beta: 1.0,
        },
    ];

    let mut phase_ln_phi = |_phase: usize, x: &[f64]| -> Result<Vec<f64>> {
        Ok(mixture
            .phase_state(&reduced, x, RootSide::Vapour)?
            .ln_phi
            .clone())
    };
    let mut single_ln_phi = |x: &[f64]| -> Result<Vec<f64>> {
        Ok(mixture
            .phase_state(&reduced, x, RootSide::Vapour)?
            .ln_phi
            .clone())
    };
    let mut ce = |x: &[f64]| -> Result<Vec<f64>> {
        let mut model = |y: &[f64]| -> Result<Vec<f64>> {
            Ok(mixture
                .phase_state(&reduced, y, RootSide::Vapour)?
                .ln_phi
                .clone())
        };
        let solved = solve_single_phase(&matrix.matrix, &g0, &b, x, &mut model)?;
        if !solved.converged {
            // The class's own behaviour with a failed solve: the composition it was given.
            return Ok(x.to_vec());
        }
        let total: f64 = solved.moles.iter().sum();
        Ok(solved.moles.iter().map(|moles| moles / total).collect())
    };

    let outcome = run(
        DriverState {
            feed_moles: moles,
            a_matrix: &matrix.matrix,
            g0: &g0,
            b: &b,
            total_moles,
            constants: &constants,
            charges: &charges,
            temperature,
            pressure: pressure_bar,
            max_phases,
            phases: initial,
        },
        &mut phase_ln_phi,
        &mut single_ln_phi,
        &mut ce,
        None,
    )?;

    // The moles per phase are the solve's own `n[j][i]`; the phases the *single-phase* branch
    // returns have a solve too, so the two branches are read the same way.
    let phases: Vec<Vec<f64>> = outcome
        .phases
        .iter()
        .enumerate()
        .map(|(index, phase)| {
            outcome
                .solution
                .as_ref()
                .and_then(|solution| solution.phase_moles.get(index).cloned())
                .unwrap_or_else(|| {
                    phase
                        .fractions
                        .iter()
                        .map(|x| x * total_moles * phase.beta)
                        .collect()
                })
        })
        .collect();

    let (residual, element_residual) = outcome
        .solution
        .as_ref()
        .map(|solution| (solution.final_residual, solution.element_residual))
        .unwrap_or((0.0, 0.0));

    Ok(ReactiveTpFlashResult {
        phase_count: phases.len(),
        phase_moles: phases,
        phase_fraction: outcome.phases.iter().map(|phase| phase.beta).collect(),
        converged: outcome.converged,
        total_iterations: outcome.total_iterations,
        equilibrium_total_moles: outcome.equilibrium_total_moles,
        gibbs_energy: outcome.gibbs_energy,
        residual,
        element_residual,
        warnings: Vec::new(),
    })
}
