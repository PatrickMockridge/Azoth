//! `reactions.reactive_phase_equilibrium` - the solve as a phase operation.
//!
//! Spec: `specs/models/reactions/reactive_phase_equilibrium.toml`, which carries the
//! boundary this takes and why azoth's flashes do not reach it.
//!
//! `ChemicalReactionOperations` is NeqSim's facade: it chooses the phase, builds the
//! element matrix and the element amounts from that phase, solves, and writes the
//! composition back. This is that, as an operation over a phase it is handed.
//!
//! # What is built unconditionally and what is not
//!
//! The phase's reaction state - the matrix, the element amounts, the reference
//! potentials - is built whether or not the phase is reactive, because that is what
//! NeqSim's own constructor does and because the matrix is a property of the substances
//! rather than of the decision to solve. **Only the solve is conditional**, and
//! `skipped` is the flag for it:
//!
//! ```text
//! A            elements the components carry, sorted, with the charge row last
//! b            A n , with the charge row replaced by the correction below
//! chem_ref     from the independent reaction basis, in J/mol
//! solve        only when the phase is aqueous, liquid or oil
//! ```
//!
//! # The charge row is a correction and not a zero
//!
//! NeqSim's `calcBVector` reads the phase's own charge, subtracts the reactive set's
//! share of it, and sets the row to the negation of the remainder:
//!
//! ```text
//! inert = sum(z_i n_i over the phase) - (A n)[charge]
//! b[charge] = |inert| <= 1e-10 max(1, n_phase) ? 0 : -inert
//! ```
//!
//! So the row is zero only when every ion in the phase is in the reaction set. A brine
//! whose chloride is not modelled carries its sodium's charge into the constraint, and
//! the solve then offsets it - which is what makes the *phase* electroneutral rather
//! than the reactive set. The captured bicarbonate brine is the fluid where the row is
//! not zero.
//!
//! # The element row order is this library's own
//!
//! NeqSim collects the elements into a `HashSet` and iterates it, so its row order is
//! Java's hash order for the element symbols - a property of a hash table, not of the
//! chemistry. Here the rows are sorted by symbol. The two matrices are the same up to a
//! permutation of rows, and the capture's `elements=` line is what makes its own order
//! readable.

use azoth_core::{AzothError, CalcResult, Result, Warning, apply_checks};

use crate::chemical_equilibrium::{ChemicalEquilibriumResult, chemical_equilibrium};
use crate::databank::{ReactionDataSource, element_composition, ionic_charge};
use crate::equilibrium_constant::GAS_CONSTANT;
use crate::model_gen;
use crate::reactive_phase::is_reactive_phase;
use crate::reference_potentials::reference_potentials;

/// The charge, in moles of elementary charge, below which the phase is treated as
/// neutral. From `calcBVector`'s `1e-10 * Math.max(1.0, phase moles)`.
pub const CHARGE_NOISE_MOLE_FRACTION: f64 = 1e-10;

/// The floor `updateMoles` raises every mole number to before writing it back, from
/// `Math.max(newMoles[i], 1e-45)`.
pub const MIN_WRITTEN_MOLES: f64 = 1e-45;

/// Result of `reactions.reactive_phase_equilibrium`.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactivePhaseEquilibriumResult {
    /// **Whether the solve was skipped**, which is NeqSim's `getReactivePhaseIndex`
    /// returning `-1`: the phase is not aqueous, liquid or oil, so there is no phase for a
    /// water-based equilibrium to be solved in. **A skip is an answer and not a failure**:
    /// the composition comes back as it went in, and a caller reading it as a failed solve
    /// would be reading two different results as one.
    pub skipped: bool,
    /// The element matrix: one row per element the reactive components carry, sorted by
    /// symbol, with the electroneutrality row last. Columns are the caller's components.
    pub a_matrix: Vec<Vec<f64>>,
    /// The element amounts the solve conserves, one per row of [`Self::a_matrix`]. **The
    /// last entry is the charge correction** and is zero only when every ion in the phase
    /// is in the reaction set.
    pub b: Vec<f64>,
    /// Each component's standard-state reference potential, in J/mol, from the
    /// independent reaction basis. Present whether or not the solve ran.
    pub chem_ref: Vec<f64>,
    /// The composition after the solve, or the caller's own on a skip.
    pub moles: Vec<f64>,
    /// Passes the solve took. Zero when skipped.
    pub iterations: u32,
    /// The solve's final error. Zero when skipped.
    pub error: f64,
    /// Whether the solve converged. False when skipped, which is why `skipped` is
    /// checked first.
    pub converged: bool,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ReactivePhaseEquilibriumResult {
    const CALC_ID: &'static str = "reactions.reactive_phase_equilibrium";
    const FIELDS: &'static [&'static str] = &[
        "skipped",
        "a_matrix",
        "b",
        "chem_ref",
        "moles",
        "iterations",
        "error",
        "converged",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// The reactive equilibrium composition of one phase.
///
/// `components` and `moles` are the phase's reactive substances and their amounts. `phase`
/// is its type name as NeqSim spells it. `phase_charge` and `phase_moles` are the whole
/// phase's own, which is what the charge correction and its tolerance are read from, and
/// `whole_system` is NeqSim's `getNumberOfPhases() == 1`, which decides whether the solve
/// corrects its conservation coupling at all.
///
/// # Errors
/// [`AzothError::InvalidInput`] where the shapes disagree or a component has no element
/// row or no charge, [`AzothError::OutOfRange`] for a bounded input, and whatever
/// [`reference_potentials`] refuses.
#[allow(clippy::too_many_arguments)] // the reaction set, the phase it is solved in, and the phase's two totals
pub fn reactive_phase_equilibrium(
    components: &[String],
    source: ReactionDataSource,
    phase: &str,
    moles: &[f64],
    phase_charge: f64,
    phase_moles: f64,
    whole_system: bool,
    log_activity: &[f64],
    temperature: f64,
    max_iterations: u32,
    tolerance: f64,
) -> Result<ReactivePhaseEquilibriumResult> {
    let spec = &model_gen::REACTIVE_PHASE_EQUILIBRIUM_SPEC;
    let mut warnings = Vec::new();

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

    if components.len() != moles.len() || components.len() != log_activity.len() {
        return Err(AzothError::InvalidInput {
            field: "components".to_string(),
            reason: format!(
                "{} component(s), {} mole(s) and {} activity coefficient(s)",
                components.len(),
                moles.len(),
                log_activity.len()
            ),
        });
    }

    let a_matrix = element_matrix(components)?;
    let mut b = element_amounts(&a_matrix, moles);

    // The charge row, as `calcBVector` computes it. `phase_charge` is the phase's own
    // `sum(z_i n_i)` over *every* component it holds, which is why it is an input: the
    // spectator ions are not in `components` and the operation cannot see them.
    let charge_row = b.len() - 1;
    let inert_charge = phase_charge - b[charge_row];
    let noise = CHARGE_NOISE_MOLE_FRACTION * phase_moles.max(1.0);
    b[charge_row] = if inert_charge.abs() <= noise {
        0.0
    } else {
        -inert_charge
    };

    let potentials =
        reference_potentials(components, source, azoth_core::units::kelvins(temperature))?;
    warnings.extend(potentials.warnings.iter().cloned());
    let chem_ref = potentials.potentials;

    if !is_reactive_phase(phase) {
        return Ok(ReactivePhaseEquilibriumResult {
            skipped: true,
            a_matrix,
            b,
            chem_ref,
            moles: moles.to_vec(),
            iterations: 0,
            error: 0.0,
            converged: false,
            warnings,
        });
    }

    // The solver takes the potentials reduced, which is the form its own iteration adds a
    // logarithm to. NeqSim divides by `R T` at the same boundary.
    let reduced: Vec<f64> = chem_ref
        .iter()
        .map(|potential| potential / (GAS_CONSTANT * temperature))
        .collect();

    let solved: ChemicalEquilibriumResult = chemical_equilibrium(
        &a_matrix,
        &b,
        whole_system,
        moles,
        &reduced,
        log_activity,
        temperature,
        max_iterations,
        tolerance,
    )?;
    warnings.extend(solved.warnings.iter().cloned());

    Ok(ReactivePhaseEquilibriumResult {
        skipped: false,
        a_matrix,
        b,
        chem_ref,
        // `updateMoles`'s floor, applied to the value the phase is left holding.
        moles: solved
            .moles
            .iter()
            .map(|moles| moles.max(MIN_WRITTEN_MOLES))
            .collect(),
        iterations: solved.iterations,
        error: solved.error,
        converged: solved.converged,
        warnings,
    })
}

/// The element matrix of a reaction set: the elements its substances carry, sorted, with
/// the electroneutrality row last.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a component with no element row or no charge, which
/// are the two absences that would otherwise become a silent zero.
fn element_matrix(components: &[String]) -> Result<Vec<Vec<f64>>> {
    let mut elements: Vec<String> = Vec::new();
    let mut composition = Vec::with_capacity(components.len());
    let mut charge = Vec::with_capacity(components.len());

    for component in components {
        let Some(counts) = element_composition(component)? else {
            return Err(AzothError::InvalidInput {
                field: "components".to_string(),
                reason: format!(
                    "`{component}` has no row in the element table, so it has no formula to \
                     put in the balance. NeqSim's `Element` is empty for it and the row is \
                     built from whatever components do have one; two fluids can therefore \
                     carry different matrices under the same name"
                ),
            });
        };
        for (element, _) in &counts {
            if !elements.contains(element) {
                elements.push(element.clone());
            }
        }
        composition.push(counts);
        let Some(z) = ionic_charge(component)? else {
            return Err(AzothError::InvalidInput {
                field: "components".to_string(),
                reason: format!(
                    "`{component}` has no row in the component databank, so its \
                                 charge is unknown and cannot be defaulted to zero"
                ),
            });
        };
        charge.push(z);
    }

    elements.sort();

    let mut matrix = vec![vec![0.0; components.len()]; elements.len() + 1];
    for (column, counts) in composition.iter().enumerate() {
        for (element, count) in counts {
            let row = elements
                .iter()
                .position(|candidate| candidate == element)
                .expect("every element was collected above");
            matrix[row][column] = *count;
        }
    }
    let charge_row = elements.len();
    for (column, z) in charge.iter().enumerate() {
        matrix[charge_row][column] = *z;
    }
    Ok(matrix)
}

/// `A n`: the element amounts the composition carries.
fn element_amounts(a_matrix: &[Vec<f64>], moles: &[f64]) -> Vec<f64> {
    a_matrix
        .iter()
        .map(|row| {
            row.iter()
                .zip(moles)
                .map(|(coefficient, moles)| coefficient * moles)
                .sum()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_element_rows_are_sorted_and_the_charge_row_is_last() {
        let components = vec!["CO2".to_string(), "H2S".to_string(), "water".to_string()];
        let matrix = element_matrix(&components).expect("every row exists");
        // C H O S sorted, then the charge row.
        assert_eq!(matrix.len(), 5);
        assert_eq!(matrix[0], vec![1.0, 0.0, 0.0], "carbon");
        assert_eq!(matrix[1], vec![0.0, 2.0, 2.0], "hydrogen");
        assert_eq!(matrix[2], vec![2.0, 0.0, 1.0], "oxygen");
        assert_eq!(matrix[3], vec![0.0, 1.0, 0.0], "sulphur");
        assert_eq!(matrix[4], vec![0.0, 0.0, 0.0], "charge");
    }

    #[test]
    fn the_charge_row_reads_the_component_databank() {
        let components = vec!["CO2".to_string(), "OH-".to_string(), "H3O+".to_string()];
        let matrix = element_matrix(&components).expect("every row exists");
        assert_eq!(matrix.last().expect("a charge row"), &vec![0.0, -1.0, 1.0]);
    }

    #[test]
    fn a_component_with_no_element_row_is_refused_rather_than_dropped() {
        // MEG is a P9 inhibitor staple and has no element row: a fluid carrying it cannot
        // enter the matrix, and a port that dropped it would answer a smaller question.
        let components = vec!["CO2".to_string(), "MEG".to_string()];
        let error = element_matrix(&components).expect_err("MEG has no formula row");
        assert!(error.to_string().contains("MEG"), "{error}");
    }

    #[test]
    fn the_charge_correction_is_zero_only_when_the_reactive_set_holds_every_ion() {
        // A brine whose sodium is not in the reaction set: the constraint carries it.
        let components = vec!["water".to_string(), "HCO3-".to_string()];
        let matrix = element_matrix(&components).expect("every row exists");
        let moles = vec![10.0, 0.02];
        let b = element_amounts(&matrix, &moles);
        let charge_row = b.len() - 1;
        assert!(
            (b[charge_row] + 0.02).abs() < 1e-12,
            "the reactive set is -0.02"
        );

        // Every ion in the phase is in the set, so the phase's own charge is that -0.02
        // and there is nothing left over.
        let inert = -0.02 - b[charge_row];
        let noise = CHARGE_NOISE_MOLE_FRACTION * 10.02_f64.max(1.0);
        assert!(inert.abs() <= noise, "{inert} against {noise}");
    }
}
