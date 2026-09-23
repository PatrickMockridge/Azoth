//! `reactions.reference_potentials` - the independent basis, and the potentials from it.
//!
//! Spec: `specs/models/reactions/reference_potentials.toml`, which carries why the three
//! steps are in this order and what the caller's component order means.
//!
//! NeqSim builds this in three passes, and the port keeps them:
//!
//! 1. **`removeJunkReactions`** - a reaction is kept when *every* reactant it names is a
//!    component the caller supplied **or** every product is. Reactants are the negative
//!    coefficients; the second half is `reactantsContains`' fall-through, and the port
//!    had only the first until `PitzerStrictnessProbe` measured a fluid it changes.
//! 2. **`removeDependentReactions`** - a greedy rank test over the reactions, dropping
//!    each one that does not raise the accumulated matrix's rank.
//! 3. **`requireValidatedEvidenceForActiveReactions`** - the `pitzer` source refuses when
//!    any *survivor* is not `VALIDATED`. It runs after both removals, so it judges the
//!    active set and not the table, and it is the only source that asks.
//! 4. **`calcReferencePotentials`** - a greedy rank test over the *components*, a solve
//!    for the independent ones from `sum(nu_i mu_i) = -R T ln K`, then a propagation to
//!    the rest in dependency order.
//!
//! # The two places a naive port goes wrong
//!
//! **The rank tests see only the stoichiometry.** `reacGMatrix` is
//! `nRows x (nComponents + 1)`, the extra column holding `-R T ln K`, and every rank
//! call is on a matrix one column narrower - so the numbers being ranked are the
//! coefficients, and their integrality is what makes the test exact. See [`crate::linalg`].
//!
//! **The order of the caller's components is the order of the answer.** NeqSim's own
//! order comes from `ChemicalReactionList.getAllComponents`, which builds an array by
//! iterating a `HashSet` - unspecified, and a property of Java's hash table rather than
//! of the chemistry. `potentials` here is in the caller's order, and the capture records
//! what NeqSim's order happened to be for the fluid it was run on.
//!
//! # The two branches that answer from the databank
//!
//! `calcReferencePotentials` has two ways to answer from a component's Gibbs energy of
//! formation rather than from the reaction set, and **both are reproduced**:
//!
//! * when the component rank falls below the reaction count it returns null, and
//!   `ChemicalReactionOperations.calcChemRefPot` reads that null as *every* component's
//!   Gibbs energy of formation (`ChemicalReactionOperations.java:513-519`);
//! * when the propagation deadlocks it seeds the first uncomputed dependent component with
//!   that same number and carries on, and a final sweep seeds whatever is left
//!   (`ChemicalReactionList.java:526-557`).
//!
//! The number is [`crate::databank::formation_properties`]'s, out of
//! `GIBBSENERGYOFFORMATION`, and it is taken **raw and not negated** - NeqSim's line is
//! `result[depCol] = gf` against a solve that used `-RT ln K` - so a seeded potential is
//! not on the same footing as a solved one. That is what makes this a fallback rather
//! than an equation, and it is why `independent` is reported: a caller can tell which
//! potentials came from the basis and which were filled in.
//!
//! **The fallback is not an edge case.** Of the 14,333 subsets of the species the three
//! tables' loaded reactions name, 11,837 answer - 3,320 of them seeding at least one
//! component - and 2,496 are refused by the pitzer source's evidence gate.
//! `every_reaction_set_the_vendored_tables_admit_answers_and_the_seed_is_exercised`
//! asserts all three counts. **The product fallback moved the seed count**: it was 4,888
//! before, on the reactant-only rule, because the reactions the products keep are ones the
//! propagation no longer deadlocks on. **A component the databank has no row for is refused rather than
//! seeded with zero**: NeqSim's `Component` field defaults to zero when the row is absent
//! and such a component cannot be built there at all, while here the name is a caller's
//! string, and eight of the element table's eighty components are in that state.

use std::collections::HashMap;

use azoth_core::units::ThermodynamicTemperature;
use azoth_core::{AzothError, CalcResult, Result, Warning, apply_checks};

use crate::databank::{ReactionDataSource, formation_properties, stoichiometry};
use crate::equilibrium_constant::GAS_CONSTANT;
use crate::linalg::{rank_of_integer_matrix, solve_lu};
use crate::model_gen;

/// The magnitude below which a stoichiometric coefficient counts as absent, from
/// `calcReferencePotentials`.
const COEFFICIENT_FLOOR: f64 = 1e-10;

/// One side of `reactantsContains`: every component on this side is one the caller supplied,
/// and the side is not empty.
///
/// **All reactants present, or all products present** - not "every reactant", which is what
/// this port had and what its spec said. The class falls through to the *products* when a
/// reactant is missing, so a reaction whose products the fluid can hold is kept even though it
/// cannot run forwards. Measured on `SystemPitzer` for a `MDEA`/`water`/`CO2` fluid: `MDEAprot`
/// names `MDEA+` as a reactant, which the fluid does not carry, and `MDEA` and `H3O+` as
/// products, which it does - so NeqSim keeps it, and the pitzer source's evidence gate then
/// refuses the fluid over it.
///
/// Reactants are the negative coefficients and products everything else, including a zero,
/// which is how the class splits them. A side with no names does not satisfy it, because the
/// class's `test` starts false and only a completed scan sets it.
///
/// Shared with [`crate::reactive_hybrid_eos_ge_flash`], which asks the same question of the
/// same rows to find out which components a fluid drives: a second copy of this rule is how the
/// two would come to different answers about one fluid.
pub(crate) fn side_is_present(
    coefficients: &[(String, f64)],
    present: &HashMap<&str, usize>,
    negative: bool,
) -> bool {
    let mut seen = false;
    for (component, nu) in coefficients {
        if (*nu < 0.0) != negative {
            continue;
        }
        seen = true;
        if !present.contains_key(component.as_str()) {
            return false;
        }
    }
    seen
}

/// Result of `reactions.reference_potentials`.
#[derive(Debug, Clone, PartialEq)]
pub struct ReferencePotentialsResult {
    /// One standard-state reference potential per component, in J/mol, **in the order the
    /// caller gave them**.
    pub potentials: Vec<f64>,
    /// A mask over the components: 1 where the basis solved for the potential directly,
    /// 0 where it was propagated from the stoichiometry.
    pub independent: Vec<f64>,
    /// A mask over the source's **loaded** reactions - `use_reaction` rows in table order
    /// - 1 where the reaction survived the dependency test and 0 where it was dropped.
    ///
    /// **The order is the table's**, and it is part of the answer rather than a
    /// presentation choice: the reducer is greedy over the rows in the order they are
    /// read, so a reaction that survives in one order can be the one dropped in another.
    pub survivors: Vec<f64>,
    /// The rank the reaction basis reached, which is what the dependency test maximised.
    pub rank: usize,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ReferencePotentialsResult {
    const CALC_ID: &'static str = "reactions.reference_potentials";
    const FIELDS: &'static [&'static str] =
        &["potentials", "independent", "survivors", "rank", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// The standard-state reference potentials of a fluid's reactive components.
///
/// `components` is the reactive set **in the caller's order**, which is the order of the
/// returned potentials. The reaction set is not a parameter: it is what the source
/// carries, filtered to the ones this fluid can run, which is how NeqSim reaches it -
/// `readReactions` reads the whole table and `removeJunkReactions` drops what the fluid
/// cannot support - all reactants present, or all products present. So the order the greedy reducer sees is the table's physical order,
/// and that order is part of the answer.
///
/// # Errors
/// [`AzothError::PropertyUnavailable`] for a reaction name the source does not carry,
/// [`AzothError::InvalidInput`] where a rank test meets a non-integral coefficient, where
/// the propagation cannot reach a component, or where the `pitzer` source has an active row
/// that is not `VALIDATED` - which is `requireValidatedEvidenceForActiveReactions`, and it
/// names the rows it rejected.
pub fn reference_potentials(
    components: &[String],
    source: ReactionDataSource,
    temperature: ThermodynamicTemperature,
) -> Result<ReferencePotentialsResult> {
    let spec = &model_gen::REFERENCE_POTENTIALS_SPEC;
    let mut warnings = Vec::new();

    let t = temperature.value;
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t),
            _ => None,
        },
        &mut warnings,
    )?;

    let index_of: HashMap<&str, usize> = components
        .iter()
        .enumerate()
        .map(|(i, name)| (name.as_str(), i))
        .collect();
    let width = components.len();

    // 1. The candidates: the source's loaded rows, in table order, and
    //    `removeJunkReactions`.
    //
    //    **NeqSim filters with the fluid's own names before its ions are added**, and
    //    this filters with the set the caller supplies - which for a fluid that already
    //    carries its ions is the same set. The capture's three fluids are checked against
    //    that claim by the crate's tests rather than assumed here, and the augmentation is
    //    why the MDEA/CO2 case supplies the *augmented* names: the gate fires on the
    //    constructor's second pass.
    //    The `-R T ln K` each row contributes is taken **here**, from the row that
    //    survived, rather than looked up again by name: `MDEAprot` is duplicated in the
    //    table, and a second lookup would answer with the first row rather than the one
    //    the reducer kept.
    let mut rows: Vec<Vec<f64>> = Vec::new();
    let mut kept: Vec<String> = Vec::new();
    let mut rhs: Vec<f64> = Vec::new();
    //: The survivor mask is over the *loaded* rows, so the position of each candidate in
    //: that shorter list is carried with it.
    let mut loaded_positions: Vec<usize> = Vec::new();
    let mut statuses: Vec<Option<String>> = Vec::new();
    let mut loaded = 0usize;
    for row in crate::databank::reactions(source)? {
        if !row.use_reaction {
            continue;
        }
        let position = loaded;
        loaded += 1;

        let coefficients = stoichiometry(&row.name)?;
        // `reactantsContains`: **all reactants present, or all products present.**
        //
        // Not "every reactant", which is what this port had and what the spec said: the
        // class falls through to the *products* when a reactant is missing, so a reaction
        // whose products the fluid can hold is kept even though it cannot run forwards.
        // Measured on `SystemPitzer` for a `MDEA`/`water`/`CO2` fluid: `MDEAprot` names
        // `MDEA+` as a reactant, which the fluid does not carry, and `MDEA` and `H3O+` as
        // products, which it does - so NeqSim keeps it and the pitzer source's evidence
        // gate then refuses the fluid over it.
        //
        // Reactants are the negative coefficients and products everything else, including
        // a zero, which is how the class splits them. A reaction with no names on a side
        // does not satisfy that side, because the class's `test` starts false and only a
        // completed scan sets it.
        if !side_is_present(&coefficients, &index_of, true)
            && !side_is_present(&coefficients, &index_of, false)
        {
            continue;
        }

        let mut coefficients_row = vec![0.0; width];
        for (component, coefficient) in &coefficients {
            if let Some(&column) = index_of.get(component.as_str()) {
                coefficients_row[column] = *coefficient;
            }
        }

        let [k1, k2, k3, k4] = row.coefficients;
        let ln_k = k1 + k2 / t + k3 * t.ln() + k4 * t;

        rows.push(coefficients_row);
        kept.push(row.name.clone());
        statuses.push(row.validation_status.clone());
        rhs.push(-GAS_CONSTANT * t * ln_k);
        loaded_positions.push(position);
    }

    // 2. `removeDependentReactions`: greedy, and order-dependent.
    let mut independent_reactions: Vec<Vec<f64>> = Vec::new();
    let mut survivors: Vec<String> = Vec::new();
    let mut survivor_statuses: Vec<Option<String>> = Vec::new();
    let mut independent_rhs: Vec<f64> = Vec::new();
    let mut survivor_mask = vec![0.0f64; loaded];
    for ((((row, name), status), value), position) in rows
        .into_iter()
        .zip(kept)
        .zip(statuses)
        .zip(rhs)
        .zip(loaded_positions)
    {
        let mut candidate = independent_reactions.clone();
        candidate.push(row.clone());
        if rank_of_integer_matrix(&candidate)? > independent_reactions.len() {
            independent_reactions = candidate;
            survivors.push(name);
            survivor_statuses.push(status);
            independent_rhs.push(value);
            survivor_mask[position] = 1.0;
        }
    }

    // 2b. `requireValidatedEvidenceForActiveReactions`: **the pitzer source is the one that
    //     requires evidence, and it requires it of the *survivors*** - the check runs after
    //     both removals upstream, so a row the fluid cannot run is not judged. Measured on
    //     `SystemPitzer`: 43 of the table's 47 rows are not `VALIDATED`, and a plain
    //     CO2/water fluid survives on the four that are, while `MDEA`/`CO2`/water trips it
    //     on `MDEAprot` alone - reached on the constructor's *second* pass, after the loop
    //     has added the `MDEA+` that reaction needs.
    if source == ReactionDataSource::Pitzer {
        let mut unvalidated: Vec<&String> = survivors
            .iter()
            .zip(&survivor_statuses)
            .filter(|(_, status)| {
                !status
                    .as_deref()
                    .is_some_and(|value| value.trim().eq_ignore_ascii_case("VALIDATED"))
            })
            .map(|(name, _)| name)
            .collect();
        if !unvalidated.is_empty() {
            // Sorted, as the class sorts its own list before naming it.
            unvalidated.sort();
            let named: Vec<String> = unvalidated.iter().map(|name| (*name).clone()).collect();
            return Err(AzothError::InvalidInput {
                field: "components".to_string(),
                reason: format!(
                    "the `{}` source rejects unvalidated active rows, and the reactions \
                     without validated evidence are {}",
                    source.identifier(),
                    named.join(", ")
                ),
            });
        }
    }

    if survivors.is_empty() {
        return Ok(ReferencePotentialsResult {
            potentials: vec![0.0; width],
            independent: vec![0.0; width],
            survivors: survivor_mask,
            rank: 0,
            warnings,
        });
    }

    // 3. The greedy column selection, which is the reaction basis transposed: a column is
    //    independent when adding it raises the rank of the rows so far.
    let n_rows = survivors.len();
    let mut current: Vec<Vec<f64>> = Vec::new();
    let mut independent_columns: Vec<usize> = Vec::new();
    let mut dependent_columns: Vec<usize> = Vec::new();

    for column in 0..width {
        // The candidate matrix: every column kept so far, plus this one. With nothing
        // kept yet that is each row's single coefficient.
        let next: Vec<Vec<f64>> = independent_reactions
            .iter()
            .enumerate()
            .map(|(i, row)| {
                let mut values = if current.is_empty() {
                    Vec::new()
                } else {
                    current[i].clone()
                };
                values.push(row[column]);
                values
            })
            .collect();

        let current_rank = if current.is_empty() {
            0
        } else {
            rank_of_integer_matrix(&current)?
        };
        if rank_of_integer_matrix(&next)? > current_rank {
            current = next;
            independent_columns.push(column);
            if independent_columns.len() == n_rows {
                for k in (column + 1)..width {
                    dependent_columns.push(k);
                }
                break;
            }
        } else {
            dependent_columns.push(column);
        }
    }

    // **A rank-deficient basis is an answer and not an error.** `calcReferencePotentials`
    // returns null when the component rank falls below the reaction count, and
    // `ChemicalReactionOperations.calcChemRefPot` reads that null as every component's
    // Gibbs energy of formation. Nothing is marked independent there, because nothing was
    // solved for; `rank` still reports the rank it fell short of.
    let rank = independent_columns.len();
    let mut potentials = vec![0.0; width];
    let mut computed = vec![false; width];
    let mut independent_mask = vec![0.0; width];

    if rank < n_rows {
        for column in 0..width {
            potentials[column] = formation_seed(components, column)?;
            computed[column] = true;
        }
    } else if n_rows > 0 {
        // The solve: `A_indep x = -B`, with `A_indep` square by the test above.
        let mut square = vec![vec![0.0; n_rows]; n_rows];
        for (i, row) in current.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                square[i][j] = *value;
            }
        }
        let negated: Vec<f64> = independent_rhs.iter().map(|value| -value).collect();
        let solved = solve_lu(&square, &negated)?;
        for (i, &column) in independent_columns.iter().enumerate() {
            potentials[column] = solved[i];
            computed[column] = true;
            independent_mask[column] = 1.0;
        }

        // The propagation: a dependent component is computable once a surviving reaction
        // exists in which every *other* component it names is known.
        let max_iterations = dependent_columns.len() * 2 + 1;
        for _ in 0..max_iterations {
            let mut progress = false;
            for &column in &dependent_columns {
                if computed[column] {
                    continue;
                }
                for (r, row) in independent_reactions.iter().enumerate() {
                    let nu = row[column];
                    if nu.abs() < COEFFICIENT_FLOOR {
                        continue;
                    }
                    let all_others_known = (0..width)
                        .all(|j| j == column || row[j].abs() <= COEFFICIENT_FLOOR || computed[j]);
                    if !all_others_known {
                        continue;
                    }
                    let mut sum_others = 0.0;
                    for j in 0..width {
                        if j != column {
                            sum_others += row[j] * potentials[j];
                        }
                    }
                    potentials[column] = (independent_rhs[r] - sum_others) / nu;
                    computed[column] = true;
                    progress = true;
                    break;
                }
            }
            if progress {
                continue;
            }
            // **The deadlock fallback, reproduced.** NeqSim seeds the first uncomputed
            // dependent component with its Gibbs energy of formation and carries on, one
            // per round that made no progress (`ChemicalReactionList.java:526-541`). The
            // propagation then continues from that value, which is what makes the seed a
            // step and not a resting place.
            let next = dependent_columns
                .iter()
                .copied()
                .find(|&column| !computed[column]);
            match next {
                Some(column) => {
                    potentials[column] = formation_seed(components, column)?;
                    computed[column] = true;
                }
                None => break,
            }
        }

        // The final sweep (`:546-557`): whatever the rounds above left uncomputed takes
        // the same seed. It is reached only when a round had nothing left to seed.
        for &column in &dependent_columns {
            if !computed[column] {
                potentials[column] = formation_seed(components, column)?;
            }
        }
    }

    apply_checks(
        spec.derived_checks(),
        |name| (name == "rank").then_some(rank as f64),
        &mut warnings,
    )?;

    Ok(ReferencePotentialsResult {
        potentials,
        independent: independent_mask,
        survivors: survivor_mask,
        rank,
        warnings,
    })
}

/// NeqSim's seed for a component the propagation cannot reach: **the component's Gibbs
/// energy of formation, taken raw and not negated** (`ChemicalReactionList.java:531`).
///
/// **A component with no databank row is refused rather than seeded with zero.** NeqSim's
/// `Component` field defaults to zero when the row is absent, and a component with no row
/// cannot be built in NeqSim at all; here the name is a caller's string, and eight of the
/// element table's eighty components are in that state - `DEAH+`, `N2O`, `NO`, `SO3`,
/// `acetaldehyde`, `dimethyl ether`, `formaldehyde` and `propylene` - so the seed has no
/// value to take and the refusal names the component.
fn formation_seed(components: &[String], column: usize) -> Result<f64> {
    let name = &components[column];
    match formation_properties(name)? {
        Some(row) => Ok(row.gibbs_energy_of_formation),
        None => Err(AzothError::InvalidInput {
            field: "reactions".to_string(),
            reason: format!(
                "the reference potential of {name} cannot be reached by propagation and \
                 the component databank has no row for it, so there is no Gibbs energy of \
                 formation to seed it with as NeqSim's fallback does"
            ),
        }),
    }
}
