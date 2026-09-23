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

use crate::chemical_equilibrium::{
    ChemicalEquilibriumResult, MIN_MOLES, chemical_equilibrium,
};
use crate::databank::{
    ReactionDataSource, element_composition, ionic_charge, reactions, stoichiometry,
};
use crate::equilibrium_constant::{GAS_CONSTANT, equilibrium_constant};
use crate::model_gen;
use crate::reactive_phase::is_reactive_phase;
use crate::reference_potentials::reference_potentials;

/// The charge, in moles of elementary charge, below which the phase is treated as
/// neutral. From `calcBVector`'s `1e-10 * Math.max(1.0, phase moles)`.
pub const CHARGE_NOISE_MOLE_FRACTION: f64 = 1e-10;

/// The floor `updateMoles` raises every mole number to before writing it back, from
/// `Math.max(newMoles[i], 1e-45)`.
pub const MIN_WRITTEN_MOLES: f64 = 1e-45;

/// The reaction log residual `solveChemEq` requires, from `REACTION_LOG_RESIDUAL_TOLERANCE`.
pub const REACTION_LOG_RESIDUAL_TOLERANCE: f64 = 2.0e-6;

/// The net charge it requires, in moles of elementary charge, from
/// `REACTIVE_PHASE_CHARGE_TOLERANCE_MOLES`.
pub const REACTIVE_PHASE_CHARGE_TOLERANCE_MOLES: f64 = 1.0e-8;

/// The element-balance residual it requires, in mol, from
/// `ELEMENT_BALANCE_RESIDUAL_TOLERANCE_MOLES`.
pub const ELEMENT_BALANCE_RESIDUAL_TOLERANCE_MOLES: f64 = 1.0e-8;

/// Which starting composition the Newton solve is given.
///
/// NeqSim's production path always runs the linear program: `initCalc` is built in
/// `ChemicalReactionOperations`' constructor and is only ever cleared inside
/// `reinitializeForReactivePhase`, which re-arms it. `None` is the direct path the captured
/// cases model, where the potentials are solved from the caller's own composition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactionSeed {
    /// Start from the caller's composition.
    None,
    /// Start from `LinearProgrammingChemicalEquilibrium`'s estimate, where the program has
    /// a solution - see [`crate::lp_seed`] for what happens where it does not.
    LinearProgramming,
}

impl std::str::FromStr for ReactionSeed {
    type Err = AzothError;

    /// A seed this library does not carry is refused rather than defaulted: the two differ
    /// in where the solve starts, and a fallback would answer a question not asked.
    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        match value.trim().to_ascii_lowercase().as_str() {
            "none" => Ok(Self::None),
            "linear_programming" | "linear-programming" => Ok(Self::LinearProgramming),
            other => Err(AzothError::InvalidInput {
                field: "seed".to_string(),
                reason: format!("`{other}` is not one of `none`, `linear_programming`"),
            }),
        }
    }
}

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
    /// **`solveChemEq`'s return**: the solve converged *and* the three residuals came in
    /// under their tolerances. False when skipped, false for a solve that ran and did not
    /// converge, and **false for one that converged and was not certified** - which is the
    /// distinction the raw solver flag cannot make.
    pub converged: bool,
    /// NeqSim's `MAXIMUM_EQUILIBRIUM_REFINEMENTS` loop: how many refinements it ran. **The
    /// port takes the first and not the second**, which needs `useAdaptiveDerivatives` and a
    /// live phase - so this is 1 wherever a solve ran and 0 on a skip.
    pub refinements: u32,
    /// Whether all three residuals came in under their tolerances: `2e-6` on the reaction
    /// log residual, `1e-8` mol on the net charge, `1e-8` mol on the element balance.
    pub certified: bool,
    /// `max |ln Q - ln K|` over the surviving reactions, at the answer. **The residual the
    /// certificate exists for**: on all five captured fluids it sits between `13.9` and
    /// `29.6`, so NeqSim's own `solveChemEq` returns false on every one of them.
    pub max_reaction_log_residual: f64,
    /// The phase's net charge `sum(z_i n_i)`, in moles of elementary charge, over **every**
    /// component it holds. `NaN` on a skip, which is what NeqSim returns there.
    pub net_charge_moles: f64,
    /// `max |A n - b|` over the element rows, in mol, at the answer. **The charge row is not
    /// among them** - it is checked separately, as the net charge.
    pub max_element_residual: f64,
    /// **Whether the linear program's estimate became the starting composition.** False is
    /// not a failure: the estimate is a state the program either has or does not, and
    /// NeqSim's own entry point reads its absence as "keep the composition the phase has".
    pub seed_applied: bool,
    /// The composition the solve started from: the estimate where one exists and the
    /// caller's own `moles` where it does not. **Floored at `updateMoles`'s `1e-45`** when
    /// the estimate is used, which is what NeqSim writes back.
    pub seed_moles: Vec<f64>,
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
        "refinements",
        "certified",
        "max_reaction_log_residual",
        "net_charge_moles",
        "max_element_residual",
        "seed_applied",
        "seed_moles",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// NeqSim's three-residual certificate, evaluated at the answer.
///
/// The three are `solveChemEq`'s own gate (`ChemicalReactionOperations.java:652-654`), and
/// one of them needs the reaction set back: `max |ln Q - ln K|` over the reactions this
/// fluid can run, with `Q` from the same activity term the solve used. The other two are
/// `A n - b` over the element rows and the phase's net charge.
///
/// **A species a surviving reaction names but the caller did not supply is refused**, not
/// skipped: NeqSim reads it off the live phase, which holds every product the reaction
/// machinery added, and a caller's shorter list would silently drop the reaction and report
/// a smaller residual than the one that exists.
#[allow(clippy::too_many_arguments)] // the reaction set, the phase, and the answer
fn certify(
    components: &[String],
    source: ReactionDataSource,
    temperature: f64,
    a_matrix: &[Vec<f64>],
    b: &[f64],
    answer: &[f64],
    input: &[f64],
    log_activity: &[f64],
    phase_charge: f64,
    phase_moles: f64,
) -> Result<(f64, f64, f64)> {
    let mut worst_reaction = 0.0_f64;
    for row in reactions(source)? {
        if !row.use_reaction {
            continue;
        }
        let species = stoichiometry(&row.name)?;
        if species.is_empty() {
            continue;
        }
        // `removeJunkReactions`: **all reactants present, or all products present** - the
        // same rule `reference_potentials` applies, corrected there first. The class falls
        // through to the products when a reactant is missing, so a reaction whose products
        // the fluid holds is kept even though it cannot run forwards; it is then in the
        // list this residual certifies, and its residual is what the gate reads.
        let side_present = |negative: bool| {
            let side: Vec<&(String, f64)> = species
                .iter()
                .filter(|(_, coefficient)| (*coefficient < 0.0) == negative)
                .collect();
            !side.is_empty()
                && side
                    .iter()
                    .all(|(name, _)| components.iter().any(|held| held == name))
        };
        if !side_present(true) && !side_present(false) {
            continue;
        }
        let ln_k =
            equilibrium_constant(source, &row.name, azoth_core::units::kelvins(temperature))?.ln_k;
        let mut quotient = 0.0_f64;
        for (name, coefficient) in &species {
            let index = components
                .iter()
                .position(|held| held == name)
                .ok_or_else(|| AzothError::InvalidInput {
                    field: "components".to_string(),
                    reason: format!(
                        "the reaction `{}` names {name}, which this fluid does not carry, so \
                         its log residual cannot be evaluated here",
                        row.name
                    ),
                })?;
            let x = (answer[index] / phase_moles).max(MIN_MOLES);
            quotient += coefficient * (x.ln() + log_activity[index]);
        }
        worst_reaction = worst_reaction.max((quotient - ln_k).abs());
    }

    // The charge row is the phase's net charge, and the correction `b` carries is applied
    // to the reactive set only - so the phase's own charge moves by what the set gained.
    let charge_row = a_matrix.len() - 1;
    let mut net_charge = phase_charge;
    for i in 0..answer.len() {
        net_charge += a_matrix[charge_row][i] * (answer[i] - input[i]);
    }

    let mut worst_element = 0.0_f64;
    for row in 0..charge_row {
        let amount: f64 = (0..answer.len())
            .map(|i| a_matrix[row][i] * answer[i])
            .sum();
        worst_element = worst_element.max((amount - b[row]).abs());
    }

    Ok((worst_reaction, net_charge, worst_element))
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
    seed: ReactionSeed,
    concentration_basis: crate::chemical_equilibrium::ConcentrationBasis,
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
            // Nothing was solved, so nothing is certified: NeqSim's
            // `getReactivePhaseChargeMoles` is NaN here and its two other residuals have no
            // phase to be computed on.
            refinements: 0,
            certified: false,
            max_reaction_log_residual: f64::NAN,
            net_charge_moles: f64::NAN,
            max_element_residual: f64::NAN,
            // Nothing was solved, so no estimate was applied: the caller's own composition
            // comes back as both the seed and the answer.
            seed_applied: false,
            seed_moles: moles.to_vec(),
            warnings,
        });
    }

    // The solver takes the potentials reduced, which is the form its own iteration adds a
    // logarithm to. NeqSim divides by `R T` at the same boundary.
    let reduced: Vec<f64> = chem_ref
        .iter()
        .map(|potential| potential / (GAS_CONSTANT * temperature))
        .collect();

    // NeqSim's seed: the estimate is written back into the phase through `updateMoles`
    // before the solve starts, so it is the starting composition where one exists.
    let estimate = match seed {
        ReactionSeed::None => None,
        ReactionSeed::LinearProgramming => {
            crate::lp_seed::initial_estimate(&a_matrix, &b, &reduced)
        }
    };
    let seed_moles: Vec<f64> = estimate.clone().unwrap_or_else(|| moles.to_vec());
    let start: Vec<f64> = if estimate.is_some() {
        seed_moles
            .iter()
            .map(|value| value.max(MIN_WRITTEN_MOLES))
            .collect()
    } else {
        moles.to_vec()
    };

    // **The basis is the caller's, and its other two facts are derived here.** NeqSim reads
    // all three off its system; this operation is handed vectors, so the solvent mask is the
    // components whose reference state is `solvent` and the solvent's mass is `sum(n_j M_j)`
    // over them - both of them properties of `moles` and the component databank, and neither
    // of them something a caller should have to restate. On the mole-fraction basis the mass
    // is computed and never read, which is what the solver's own spec says of it.
    let solvent_mask = solvent_mask(components)?;
    let solvent_weight = solvent_weight(components, moles, &solvent_mask);
    let solved: ChemicalEquilibriumResult = chemical_equilibrium(
        &a_matrix,
        &b,
        whole_system,
        &start,
        &reduced,
        log_activity,
        temperature,
        max_iterations,
        tolerance,
        concentration_basis,
        solvent_weight,
        &solvent_mask,
        phase_moles,
    )?;
    warnings.extend(solved.warnings.iter().cloned());

    // NeqSim's certificate, on the composition the solve left. **`converged` is this and
    // not the solver's own flag**: `solveChemEq` returns false wherever one of the three
    // residuals is over its tolerance, and on all five captured fluids it is.
    let (max_reaction_log_residual, net_charge_moles, max_element_residual) = certify(
        components,
        source,
        temperature,
        &a_matrix,
        &b,
        &solved.moles,
        moles,
        log_activity,
        phase_charge,
        phase_moles,
    )?;
    let certified = solved.converged
        && max_reaction_log_residual <= REACTION_LOG_RESIDUAL_TOLERANCE
        && net_charge_moles.abs() <= REACTIVE_PHASE_CHARGE_TOLERANCE_MOLES
        && max_element_residual <= ELEMENT_BALANCE_RESIDUAL_TOLERANCE_MOLES;

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
        seed_applied: estimate.is_some(),
        seed_moles,
        converged: certified,
        refinements: 1,
        certified,
        max_reaction_log_residual,
        net_charge_moles,
        max_element_residual,
        warnings,
    })
}

/// The element matrix of a reaction set: the elements its substances carry, sorted, with
/// the electroneutrality row last.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a component with no element row or no charge, which
/// are the two absences that would otherwise become a silent zero.
/// Which components keep the mole-fraction form on the solute-molality basis.
///
/// The reference state is the component databank's, for the reason the ionic charge is: a
/// caller that had to state the mask could state a different one from the table's, and the
/// branch it selects is a property of the substance.
///
/// # Errors
/// [`AzothError::PropertyUnavailable`] if a component has no databank row.
pub(crate) fn solvent_mask(components: &[String]) -> Result<Vec<f64>> {
    components
        .iter()
        .map(|name| {
            let entry = azoth_eos::databank::entry(name, None)?;
            Ok(f64::from(
                u8::from(entry.reference_state == azoth_eos::databank::SOLVENT),
            ))
        })
        .collect()
}

/// `sum(n_j M_j)` over the components the mask marks, in kg.
///
/// NeqSim's `getPhase().getTotalVolume()`-side quantity is the phase's own; this is the
/// reactive subset's, which is the only part of the phase this operation can see.
#[must_use]
pub(crate) fn solvent_weight(components: &[String], moles: &[f64], mask: &[f64]) -> f64 {
    let mut weight = 0.0;
    for ((name, amount), marked) in components.iter().zip(moles).zip(mask) {
        if *marked < 0.5 {
            continue;
        }
        if let Ok(entry) = azoth_eos::databank::entry(name, None) {
            if let Some(mass) = entry.molar_mass {
                weight += amount.max(0.0) * mass;
            }
        }
    }
    weight
}

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
