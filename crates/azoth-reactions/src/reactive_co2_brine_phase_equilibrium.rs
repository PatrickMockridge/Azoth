//! `reactions.reactive_co2_brine_phase_equilibrium` - the phase/reaction coupling.
//!
//! Spec: `specs/models/reactions/reactive_co2_brine_phase_equilibrium.toml`.
//!
//! NeqSim's `ReactiveCO2BrinePhaseEquilibrium.run`: a phase selection with the species
//! inventory fixed, then chemistry on an isolated copy of the aqueous phase whose
//! conservative changes are the only ones written back. See the spec's assumptions for what
//! the loop does and for the two things this port states rather than hides - the activity
//! coefficients it holds fixed and the pass count its path decides.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, CalcResult, Result, Warning, apply_checks};
use azoth_eos::co2_brine_phase_equilibrium::co2_brine_phase_equilibrium;
use azoth_eos::databank as eos_databank;
use azoth_eos::electrolyte_cpa_statoil_phase::electrolyte_cpa_statoil_phase;

use crate::chemical_equilibrium::{ConcentrationBasis, chemical_equilibrium};
use crate::databank::{ReactionDataSource, ionic_charge};
use crate::equilibrium_constant::GAS_CONSTANT;
use crate::formula_matrix::FormulaMatrix;
use crate::lp_seed;
use crate::model_gen;
use crate::reactive_phase_equilibrium::{certify, solvent_mask, solvent_weight};
use crate::reference_potentials::reference_potentials;

/// The class's `MAXIMUM_ITERATIONS`, and the two gates its accept test applies.
const MAXIMUM_ITERATIONS: u32 = 30;
/// The reaction log residual the accept test holds.
const REACTION_TOLERANCE: f64 = 2.0e-6;
/// The net charge the accept test holds, in moles of elementary charge.
const CHARGE_TOLERANCE_MOLES: f64 = 1.0e-8;

/// The element-conservation tolerance `verifyElements` uses.
const ELEMENT_TOLERANCE_FLOOR: f64 = 1.0e-9;
const ELEMENT_TOLERANCE_RELATIVE: f64 = 1.0e-10;

/// The chemistry solve's own parameters, from `solveChemEq`'s defaults.
const CHEMISTRY_MAX_ITERATIONS: u32 = 100;
const CHEMISTRY_TOLERANCE: f64 = 1.0e-8;

/// `Phase.initRefPhases`' binary: the component at `1e-10` mol against `10` mol of water.
const REFERENCE_SOLUTE_MOLES: f64 = 1.0e-10;
const REFERENCE_WATER_MOLES: f64 = 10.0;

/// Result of `reactions.reactive_co2_brine_phase_equilibrium`.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactiveCo2BrinePhaseEquilibriumResult {
    /// How many phases the accepted state has.
    pub phase_count: u32,
    /// The label of each row: the root the CO2-rich phase was solved on, then `aqueous`.
    pub phase_type: Vec<String>,
    /// The mole fraction of the feed in each phase, summing to one.
    pub phase_fractions: Vec<f64>,
    /// Each phase's composition, one row per phase and one column per component.
    pub phase_composition: Vec<Vec<f64>>,
    /// Which root the CO2-rich phase was solved on.
    pub accepted_root: String,
    /// The tangent-plane distance the phase selection decided on.
    pub minimum_trial_distance: f64,
    /// The coupling's passes. **Path-decided**: see the spec's assumptions.
    pub iterations: u32,
    /// The worst `|ln Q - ln K|` over the reaction set at the accepted state.
    pub reaction_log_residual: f64,
    /// The accepted inventory's net charge, in moles of elementary charge.
    pub charge_moles: f64,
    /// The worst `|A n - b|` over the element rows at the accepted state, in mol.
    pub element_balance_residual: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ReactiveCo2BrinePhaseEquilibriumResult {
    const CALC_ID: &'static str = "reactions.reactive_co2_brine_phase_equilibrium";
    const FIELDS: &'static [&'static str] = &[
        "phase_count",
        "phase_type",
        "phase_fractions",
        "phase_composition",
        "accepted_root",
        "minimum_trial_distance",
        "iterations",
        "reaction_log_residual",
        "charge_moles",
        "element_balance_residual",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// The simultaneous phase and chemical equilibrium of a water-rich CO2/brine feed.
///
/// # Errors
/// * [`AzothError::InvalidInput`] where the shapes disagree, a component carries no element
///   row, or the loop produces an inventory that is not one.
/// * [`AzothError::SolverNotConverged`] if a chemistry solve or the coupling does not settle,
///   or the element balance moves.
pub fn reactive_co2_brine_phase_equilibrium(
    components: &[String],
    moles: &[f64],
    temperature: ThermodynamicTemperature,
    pressure: Pressure,
) -> Result<ReactiveCo2BrinePhaseEquilibriumResult> {
    let spec = &model_gen::REACTIVE_CO2_BRINE_PHASE_EQUILIBRIUM_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(temperature.value),
            "P" => Some(pressure.value),
            _ => None,
        },
        &mut warnings,
    )?;

    if components.len() != moles.len() {
        return Err(AzothError::invalid_input(
            "moles",
            format!(
                "{} component(s) and {} mole number(s); the two must match",
                components.len(),
                moles.len()
            ),
        ));
    }
    if moles.iter().any(|value| !value.is_finite() || *value < 0.0) {
        return Err(AzothError::invalid_input(
            "moles",
            "every mole number must be finite and not negative",
        ));
    }

    let matrix = FormulaMatrix::build(components)?;
    let a_matrix = matrix.matrix.clone();
    let expected = element_amounts(&a_matrix, moles);
    let mask = solvent_mask(components)?;
    let charges = components
        .iter()
        .map(|name| {
            ionic_charge(name)?.ok_or_else(|| {
                AzothError::invalid_input(
                    "components",
                    format!(
                        "`{name}` has no row in the component databank, so its charge is unknown"
                    ),
                )
            })
        })
        .collect::<Result<Vec<f64>>>()?;

    let mut work = moles.to_vec();
    let mut residual = f64::NAN;
    for iteration in 1..=MAXIMUM_ITERATIONS {
        let flash = co2_brine_phase_equilibrium(components, temperature, pressure, &work)?;
        warnings.extend(flash.warnings.iter().cloned());

        let aqueous = if flash.phase_count == 2 { 1 } else { 0 };
        let beta = flash.phase_fractions[aqueous];
        let total: f64 = work.iter().sum();
        let before: Vec<f64> = flash.x[aqueous].iter().map(|x| x * beta * total).collect();
        let mut b = element_amounts(&a_matrix, &before);
        let charge_index = b.len() - 1;
        b[charge_index] = 0.0;
        let log_activity =
            log_activity_coefficients(components, &flash.ln_phi[aqueous], temperature, pressure)?;
        let phase_charge: f64 = charges
            .iter()
            .zip(&before)
            .map(|(charge, amount)| charge * amount)
            .sum();
        let (reaction, _, _) = certify(
            components,
            ReactionDataSource::Standard,
            temperature.value,
            &a_matrix,
            &b,
            &before,
            &before,
            &log_activity,
            phase_charge,
            before.iter().sum(),
        )?;
        residual = reaction;

        if reaction.is_finite() && reaction <= REACTION_TOLERANCE {
            let charge: f64 = charges
                .iter()
                .zip(&work)
                .map(|(charge, amount)| charge * amount)
                .sum();
            if !charge.is_finite() || charge.abs() > CHARGE_TOLERANCE_MOLES {
                return Err(AzothError::SolverNotConverged {
                    iterations: iteration,
                    residual: charge,
                    tolerance: CHARGE_TOLERANCE_MOLES,
                });
            }
            let element_residual = element_balance(&a_matrix, &expected, &work);
            let phase_count = flash.phase_count;
            let phase_type = if phase_count == 2 {
                vec![flash.accepted_root.clone(), "aqueous".to_string()]
            } else {
                vec!["aqueous".to_string()]
            };
            return Ok(ReactiveCo2BrinePhaseEquilibriumResult {
                phase_count,
                phase_type,
                phase_fractions: flash.phase_fractions,
                phase_composition: flash.x,
                accepted_root: flash.accepted_root,
                minimum_trial_distance: flash.minimum_trial_distance,
                iterations: iteration,
                reaction_log_residual: reaction,
                charge_moles: charge,
                element_balance_residual: element_residual,
                warnings,
            });
        }

        // The chemistry runs on an isolated copy of the aqueous phase; only its conservative
        // species changes are written back.
        let potentials =
            reference_potentials(components, ReactionDataSource::Standard, temperature)?;
        warnings.extend(potentials.warnings.iter().cloned());
        let reduced: Vec<f64> = potentials
            .potentials
            .iter()
            .map(|potential| potential / (GAS_CONSTANT * temperature.value))
            .collect();
        let start = match lp_seed::initial_estimate(&a_matrix, &b, &reduced) {
            Some(estimate) => estimate
                .iter()
                .map(|value| value.max(lp_seed::MIN_MOLES))
                .collect(),
            None => before.clone(),
        };
        let solved = chemical_equilibrium(
            &a_matrix,
            &b,
            true,
            &start,
            &reduced,
            &log_activity,
            temperature.value,
            CHEMISTRY_MAX_ITERATIONS,
            CHEMISTRY_TOLERANCE,
            ConcentrationBasis::MoleFraction,
            solvent_weight(components, &before, &mask),
            &mask,
            before.iter().sum(),
        )?;
        warnings.extend(solved.warnings.iter().cloned());
        // A non-converged chemistry solve contributes its last iterate; the loop's element
        // check and the flashed-state accept test are what decide the answer. The spec
        // records why this differs from `solveChemEq`'s throw.
        let mut updated = Vec::with_capacity(work.len());
        for index in 0..work.len() {
            let value = work[index] + (solved.moles[index] - before[index]);
            if !value.is_finite() || value < -1.0e-12 {
                return Err(AzothError::InvalidInput {
                    field: "moles".to_string(),
                    reason: format!(
                        "the coupling produced {value} mol of component {index}, which is not an inventory"
                    ),
                });
            }
            updated.push(value.max(0.0));
        }
        let actual = element_amounts(&a_matrix, &updated);
        for (wanted, got) in expected.iter().zip(&actual) {
            let tolerance = ELEMENT_TOLERANCE_FLOOR.max(wanted.abs() * ELEMENT_TOLERANCE_RELATIVE);
            if !got.is_finite() || (got - wanted).abs() > tolerance {
                return Err(AzothError::SolverNotConverged {
                    iterations: iteration,
                    residual: got - wanted,
                    tolerance,
                });
            }
        }
        work = updated;
    }

    Err(AzothError::SolverNotConverged {
        iterations: MAXIMUM_ITERATIONS,
        residual,
        tolerance: REACTION_TOLERANCE,
    })
}

/// `ln gamma_k`: the phase's own `ln phi_k` less its unsymmetric reference.
fn log_activity_coefficients(
    components: &[String],
    actual: &[f64],
    temperature: ThermodynamicTemperature,
    pressure: Pressure,
) -> Result<Vec<f64>> {
    let mut out = Vec::with_capacity(components.len());
    for (index, name) in components.iter().enumerate() {
        let reference = if eos_databank::entry(name, None)?.reference_state == eos_databank::SOLVENT
        {
            electrolyte_cpa_statoil_phase(
                std::slice::from_ref(name),
                temperature,
                pressure,
                &[1.0],
                "liquid",
            )?
            .ln_phi[0]
        } else {
            let fraction =
                REFERENCE_SOLUTE_MOLES / (REFERENCE_SOLUTE_MOLES + REFERENCE_WATER_MOLES);
            electrolyte_cpa_statoil_phase(
                &[name.clone(), "water".to_string()],
                temperature,
                pressure,
                &[fraction, 1.0 - fraction],
                "liquid",
            )?
            .ln_phi[0]
        };
        out.push(actual[index] - reference);
    }
    Ok(out)
}

/// `A n`: the element amounts a composition carries.
fn element_amounts(a_matrix: &[Vec<f64>], moles: &[f64]) -> Vec<f64> {
    a_matrix
        .iter()
        .map(|row| {
            row.iter()
                .zip(moles)
                .map(|(coefficient, amount)| coefficient * amount)
                .sum()
        })
        .collect()
}

/// The worst `|A n - b|` over the element rows, excluding the charge row.
fn element_balance(a_matrix: &[Vec<f64>], expected: &[f64], moles: &[f64]) -> f64 {
    let actual = element_amounts(a_matrix, moles);
    let mut worst = 0.0_f64;
    for index in 0..expected.len().saturating_sub(1) {
        worst = worst.max((actual[index] - expected[index]).abs());
    }
    worst
}
