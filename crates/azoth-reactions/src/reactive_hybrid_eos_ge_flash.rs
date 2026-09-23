//! The reactive coupling on the fixed-role hybrid flash, from `TPHybridEosGeFlash.run`.
//!
//! `eos.hybrid_eos_ge_flash` solves one gas-oil-brine split at a *fixed* species inventory.
//! NeqSim's `run` wraps it in a second loop when `system.isChemicalSystem()`: each pass
//! re-equilibrates the brine's chemistry, projects the resulting mole changes onto the
//! stoichiometric conservation null space, writes the adjusted inventory back, and only then
//! solves the fractions again. This module is that loop, and it lives here rather than beside
//! the flash because the reaction tables do: `azoth-eos` cannot see them.
//!
//! # The loop, in the class's own order
//!
//! ```text
//! prepare the roles
//! solveFixedTopologyPhaseEquilibrium          the feed-balanced split, before chemistry
//! repeat at most 100 passes:
//!     solveAqueousChemicalEquilibrium         chemistry, projection, write-back
//!     solveFixedTopologyPhaseEquilibrium      fractions at the adjusted inventory
//!     stop when pass >= 3 and both tests pass
//! ```
//!
//! **The first solve is not part of the loop.** Chemistry would otherwise meet the role seeds
//! rather than a split, and NeqSim establishes a feed-balanced one before it does.
//!
//! # What the passes share
//!
//! The split, which is NeqSim's `system`: each pass re-enters the fraction solve from the
//! fractions the last one converged on. `solve_fixed_topology_from` is that re-entry point, and
//! re-seeding each pass would throw away the only thing the passes have in common.
//!
//! **One divergence, and it is measured.** NeqSim's write-back also runs
//! `synchronizeHybridEosGeOverallComposition`'s `initBeta`/`normalizeBeta`, which recompute each
//! role's fraction as `phase moles / coupled total` from the phase objects' own stored counts -
//! quantities that move in the eighth digit between passes. This carries the fractions instead.
//! The answer is unaffected: the fraction solve converges from either starting point to the
//! captured endpoint, and the pass count is unaffected because the earliest a pass may converge
//! is the class's own minimum of three.
//!
//! # The chemistry's input is the phase's activity vector
//!
//! `solveChemEq` reads `phase.getLogActivityCoefficient(i, water)` once, before its solve, and
//! that is the `log_activity` this loop hands
//! [`crate::reactive_phase_equilibrium::reactive_phase_equilibrium`] - which is what lets that
//! operation state its concentration basis at all. See [`log_activities`] for how the vector is
//! recovered from `eos.pitzer_phase`, and the spec for what the recovery does not claim.

use std::collections::HashMap;

use azoth_core::spec::ModelAlgorithm;
use azoth_core::{AzothError, CalcResult, Result, Warning};
use azoth_eos::Cubic;
use azoth_eos::databank::mixture_of_with_ions;
use azoth_eos::hybrid_eos_ge_flash::{seed_phases, solve_fixed_topology_from};
use azoth_eos::multiphase::MultiphasePhase;
use azoth_eos::pitzer_phase::pitzer_phase;

use crate::chemical_equilibrium::ConcentrationBasis;
use crate::databank::{ReactionDataSource, ionic_charge, reactions, stoichiometry};
use crate::linalg::{project_onto_null_space, row_space_basis};
use crate::reactive_phase_equilibrium::{ReactionSeed, element_matrix, reactive_phase_equilibrium};
use crate::reference_potentials::side_is_present;

/// The reaction source every fluid of this model runs.
///
/// `SystemPitzer.getChemicalReactionDataSource()` is `PITZER`, and the hybrid route is reachable
/// only from an EoS/GE system - so the source is a property of the model and not an input. It is
/// also the source that carries the validation-status column, whose gate is inside
/// [`reactive_phase_equilibrium`]'s own potentials.
pub const SOURCE: ReactionDataSource = ReactionDataSource::Pitzer;

/// The outer loop's cap, `MAXIMUM_REACTIVE_ITERATIONS`.
pub const MAXIMUM_REACTIVE_PASSES: u32 = 100;

/// The fewest passes that may certify a coupled state, `MINIMUM_REACTIVE_ITERATIONS`.
///
/// The chemistry and the split are each other's input, and two passes are the least that lets
/// either see the other's answer. Measured on the captured fluids, the second pass fails both
/// clauses anyway - `2 >= 3` and a deviation of `1.011e-10` against `1e-10` - so the floor is
/// not what stops the loop there, but it is what makes a two-pass answer unreportable in
/// general.
pub const MINIMUM_REACTIVE_PASSES: u32 = 3;

/// The sum of absolute aqueous mole-fraction changes below which the chemistry has settled,
/// `REACTIVE_COMPOSITION_TOLERANCE`.
pub const REACTIVE_COMPOSITION_TOLERANCE: f64 = 1.0e-10;

/// The fixed-topology solve's own tolerance, `HYBRID_SOLVER_TOLERANCE`.
pub const HYBRID_SOLVER_TOLERANCE: f64 = 1.0e-10;

/// The floor the coupled inventory is raised to, `MINIMUM_COUPLED_COMPONENT_MOLES`.
pub const MINIMUM_COUPLED_MOLES: f64 = 1.0e-45;

/// The raw reaction delta below which the conservation projection is **not** applied,
/// `REACTION_DELTA_CONSERVATION_TOLERANCE_MOLES`.
pub const REACTION_DELTA_CONSERVATION_TOLERANCE_MOLES: f64 = 1.0e-8;

/// The negative overall amount that is a failure rather than rounding, `-1.0e-9`.
pub const NEGATIVE_INVENTORY_TOLERANCE_MOLES: f64 = 1.0e-9;

/// The brine Newton's pass cap, `ChemicalEquilibrium`'s own default - `solveChemEq` is called
/// without one, so the class's is what runs.
pub const CHEMISTRY_MAX_ITERATIONS: u32 = 100;

/// The error that Newton is trying to reach, the class's own default for the same reason.
pub const CHEMISTRY_TOLERANCE: f64 = 1.0e-8;

/// The water the class's two-component reference phase holds, in mol, from `initRefPhases`.
pub const REFERENCE_SOLVENT_MOLES: f64 = 10.0;

/// The solute it holds, in mol, from the same place.
pub const REFERENCE_SOLUTE_MOLES: f64 = 1.0e-10;

/// The aqueous role's index, `azoth_eos::hybrid_eos_ge_flash`'s role order.
const AQUEOUS_ROLE: usize = 2;

/// The reactive hybrid flash's answer.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactiveHybridEosGeFlashResult {
    /// Each role's mole fraction of the feed, in `[gas, oil, aqueous]` order, summing to one.
    pub beta: Vec<f64>,
    /// Each role's composition, one row per role and one column per component.
    pub x: Vec<Vec<f64>>,
    /// The reaction-adjusted overall inventory at the answer, one entry per component. **Not the
    /// feed**: chemistry changes species amounts and the total moves with them.
    pub coupled_moles: Vec<f64>,
    /// The brine's species amounts at the answer, in the reactive set's own order.
    pub aqueous_moles: Vec<f64>,
    /// Coupled passes taken, the class's own count and not the fraction solve's.
    ///
    /// **A loop that does not certify is an error and not an answer**, which is NeqSim's own
    /// `IllegalStateException`: a returned result is one where all three conditions held, so
    /// there is no flag to report beside it.
    pub passes: u32,
    /// The last pass's composition deviation - the sum of `|x_old - x_new|` over the brine.
    pub chemical_deviation: f64,
    /// The last fraction solve's own residual: the larger of its step norm and its gradient norm.
    pub residual: f64,
    /// The worst `|z_i - sum_k beta_k x_ik|` at the answer, over the coupled inventory.
    pub max_material_balance_residual: f64,
    /// The worst spread of `ln(x_i phi_i P)` over the roles that are there.
    pub max_log_fugacity_residual: f64,
    /// The worst element residual of the coupled inventory against the feed, `|A n - A n_feed|`.
    pub element_residual: f64,
    /// The net charge the coupled inventory carries against the feed, in moles of elementary
    /// charge.
    pub charge_residual: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ReactiveHybridEosGeFlashResult {
    const CALC_ID: &'static str = "reactions.reactive_hybrid_eos_ge_flash";
    const FIELDS: &'static [&'static str] = &[
        "beta",
        "x",
        "coupled_moles",
        "aqueous_moles",
        "passes",
        "chemical_deviation",
        "residual",
        "max_material_balance_residual",
        "max_log_fugacity_residual",
        "element_residual",
        "charge_residual",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// What the coupled loop needs from the model's spec, beside the fluid.
///
/// `algorithm` is the fixed-topology solve's own, which `azoth-eos` owns; the chemistry's two
/// are `ChemicalEquilibrium`'s, and NeqSim's caller passes its defaults straight through.
#[derive(Debug, Clone, Copy)]
pub struct CoupledAlgorithm<'a> {
    /// The fraction solve's cap and tolerance.
    pub algorithm: &'a ModelAlgorithm,
    /// `solveChemEq`'s pass cap for the brine's Newton.
    pub chemistry_max_iterations: u32,
    /// The error that Newton is trying to reach.
    pub chemistry_tolerance: f64,
}

/// The coupled reactive flash, as the model `reactions.reactive_hybrid_eos_ge_flash`.
///
/// The boundary the spec declares: names rather than resolved records, and the model's own
/// checks applied before anything is solved. Everything below it is [`solve_coupled`].
///
/// # Errors
/// [`AzothError::OutOfRange`] if `T` or `P` is outside the declared range, and every error
/// [`solve_coupled`] raises.
#[allow(non_snake_case)] // `T` and `P` are the symbols in the model's own name
pub fn reactive_hybrid_eos_ge_flash(
    components: &[String],
    cubic: Cubic,
    T: f64,
    P: f64,
    moles: &[f64],
) -> Result<ReactiveHybridEosGeFlashResult> {
    let spec = &crate::model_gen::REACTIVE_HYBRID_EOS_GE_FLASH_SPEC;
    let mut warnings = Vec::new();
    azoth_core::apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T),
            "P" => Some(P),
            _ => None,
        },
        &mut warnings,
    )?;

    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let mut result = solve_coupled(
        &names,
        cubic,
        T,
        P,
        moles,
        CoupledAlgorithm {
            // The fraction solve's own algorithm, which is the underlying model's: this loop
            // re-enters it rather than running a second one.
            algorithm: azoth_eos::algorithm_of(&azoth_eos::model_gen::HYBRID_EOS_GE_FLASH_SPEC)?,
            chemistry_max_iterations: CHEMISTRY_MAX_ITERATIONS,
            chemistry_tolerance: CHEMISTRY_TOLERANCE,
        },
    )?;
    warnings.extend(result.warnings);
    result.warnings = warnings;
    Ok(result)
}

/// The coupled loop itself, over names already resolved against the databank.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the feed is not a mole vector, a component has no formula
///   or no charge row, or the chemistry drives an amount past [`NEGATIVE_INVENTORY_TOLERANCE_MOLES`].
/// * [`AzothError::PropertyUnavailable`] if a phase model cannot answer at a trial composition,
///   or the source's evidence gate refuses an active reaction row.
/// * [`AzothError::SolverNotConverged`] if a Newton correction cannot be solved, if the brine
///   cannot hold its ions, or if the coupled loop reaches its cap without converging.
pub fn solve_coupled(
    names: &[&str],
    cubic: Cubic,
    t: f64,
    p: f64,
    moles: &[f64],
    coupled_algorithm: CoupledAlgorithm<'_>,
) -> Result<ReactiveHybridEosGeFlashResult> {
    if names.len() != moles.len() {
        return Err(AzothError::invalid_input(
            "moles",
            format!(
                "a feed for {} component(s) has {} mole(s)",
                names.len(),
                moles.len()
            ),
        ));
    }

    let (mixture, _) = mixture_of_with_ions(names, cubic, None)?;
    let mut warnings: Vec<Warning> = Vec::new();
    let reduced =
        mixture.reduced_parameters(azoth_core::units::kelvins(t), azoth_core::units::pascals(p))?;
    let algorithm = coupled_algorithm.algorithm;

    let ion: Vec<bool> = names
        .iter()
        .map(|name| {
            azoth_eos::databank::entry(name, None)
                .map(|entry| entry.class == azoth_eos::databank::ION)
        })
        .collect::<Result<_>>()?;

    // The reactive set the fluid drives, and the matrix the projection is built on. Both are
    // properties of the substances rather than of the decision to solve.
    let reactive = reactive_components(names)?;
    let a_matrix = element_matrix(&reactive)?;
    let basis = row_space_basis(&a_matrix)?;
    let reactive_at: Vec<usize> = reactive
        .iter()
        .map(|name| {
            names
                .iter()
                .position(|candidate| candidate == name)
                .ok_or_else(|| {
                    AzothError::invalid_input(
                        "components",
                        format!("`{name}` is reactive and not in the fluid"),
                    )
                })
        })
        .collect::<Result<_>>()?;

    let total: f64 = moles.iter().sum();
    if !total.is_finite() || total <= 0.0 {
        return Err(AzothError::invalid_input(
            "moles",
            format!("the feed sums to {total}, which is not a mole count"),
        ));
    }
    let z: Vec<f64> = moles.iter().map(|value| value / total).collect();
    let mut phases = seed_phases(&mixture, &z, &ion)?;

    // `run`'s first solve, outside the loop.
    let seeded =
        solve_fixed_topology_from(&mixture, &reduced, moles, &ion, algorithm, &mut phases)?;

    let mut coupled = moles.to_vec();
    let mut passes = 0;
    let mut coupled_converged = false;
    let mut chemical_deviation = f64::INFINITY;
    let mut residual = seeded.residual;
    let mut last = seeded;

    for pass in 0..MAXIMUM_REACTIVE_PASSES {
        passes = pass + 1;
        chemical_deviation = chemical_step(
            names,
            &reactive,
            &reactive_at,
            &phases,
            &mut coupled,
            &a_matrix,
            &basis,
            t,
            p,
            pass == 0,
            coupled_algorithm,
        )?;

        // The adjusted inventory is what the next fraction solve balances, which is what
        // `synchronizeHybridEosGeOverallComposition` does on NeqSim's own system.
        last =
            solve_fixed_topology_from(&mixture, &reduced, &coupled, &ion, algorithm, &mut phases)?;
        residual = last.residual;

        coupled_converged = passes >= MINIMUM_REACTIVE_PASSES
            && chemical_deviation <= REACTIVE_COMPOSITION_TOLERANCE
            && residual.is_finite()
            && residual <= HYBRID_SOLVER_TOLERANCE;
        if coupled_converged {
            break;
        }
    }

    if !coupled_converged {
        return Err(AzothError::SolverNotConverged {
            iterations: passes,
            residual: residual.max(chemical_deviation),
            tolerance: HYBRID_SOLVER_TOLERANCE,
        });
    }

    let (element_residual, charge_residual) =
        conservation_residuals(names, &a_matrix, &reactive_at, &phases, &coupled)?;
    // `getNumberOfMolesInPhase` for a species of the brine: its mole fraction there times the
    // phase's own mole count, which is the adjusted inventory's total.
    let coupled_total: f64 = coupled.iter().sum();
    let aqueous_moles: Vec<f64> = reactive_at
        .iter()
        .map(|&index| {
            phases[AQUEOUS_ROLE].composition[index] * phases[AQUEOUS_ROLE].fraction * coupled_total
        })
        .collect();
    warnings.extend(last.warnings.iter().cloned());

    Ok(ReactiveHybridEosGeFlashResult {
        beta: last.beta,
        x: last.x,
        coupled_moles: coupled,
        aqueous_moles,
        passes,
        chemical_deviation,
        residual,
        max_material_balance_residual: last.max_material_balance_residual,
        max_log_fugacity_residual: last.max_log_fugacity_residual,
        element_residual,
        charge_residual,
        warnings,
    })
}

/// One coupled pass: solve the brine, project the delta, write the inventory back.
///
/// Returns the composition deviation the convergence test reads - the sum of `|x_old - x_new|`
/// over the brine's components, **taken after the write-back and before the fraction solve**,
/// which is where `solveAqueousChemicalEquilibrium` reads it.
#[allow(clippy::too_many_arguments)] // the loop's own state, and each of these is a part of it
fn chemical_step(
    names: &[&str],
    reactive: &[String],
    reactive_at: &[usize],
    phases: &[MultiphasePhase],
    coupled: &mut [f64],
    a_matrix: &[Vec<f64>],
    basis: &[Vec<f64>],
    t: f64,
    p: f64,
    initialise: bool,
    coupled_algorithm: CoupledAlgorithm<'_>,
) -> Result<f64> {
    let aqueous = &phases[AQUEOUS_ROLE];
    let total: f64 = coupled.iter().sum();
    let phase_moles = aqueous.fraction * total;

    // The phase's whole mole vector, which is what `solveChemEq` writes into and reads back.
    let mut phase_amounts: Vec<f64> = aqueous
        .composition
        .iter()
        .map(|fraction| (fraction * phase_moles).max(0.0))
        .collect();
    let old_fractions: Vec<f64> = reactive_at
        .iter()
        .map(|&index| aqueous.composition[index])
        .collect();
    let old_moles: Vec<f64> = reactive_at
        .iter()
        .map(|&index| phase_amounts[index])
        .collect();

    // The phase's own charge over *every* component it holds, which is the quantity
    // `calcBVector` subtracts the reactive set's share from.
    let mut phase_charge = 0.0;
    for (index, name) in names.iter().enumerate() {
        if let Some(charge) = ionic_charge(name)? {
            phase_charge += charge * phase_amounts[index];
        }
    }

    let activity = log_activities(names, t, p, &aqueous.composition)?;
    let log_activity: Vec<f64> = reactive_at.iter().map(|&index| activity[index]).collect();

    let solved = reactive_phase_equilibrium(
        reactive,
        SOURCE,
        "aqueous",
        &old_moles,
        phase_charge,
        phase_moles,
        false,
        &log_activity,
        t,
        coupled_algorithm.chemistry_max_iterations,
        coupled_algorithm.chemistry_tolerance,
        if initialise {
            ReactionSeed::LinearProgramming
        } else {
            ReactionSeed::None
        },
        // `SystemPitzer.getChemicalReactionConcentrationBasis()` is `SOLUTE_MOLALITY`, and the
        // captured numbers are on that basis.
        ConcentrationBasis::SoluteMolality,
    )?;

    let raw: Vec<f64> = solved
        .moles
        .iter()
        .zip(&old_moles)
        .map(|(new, old)| new - old)
        .collect();
    let delta = conservative_delta(a_matrix, basis, coupled, reactive_at, &raw);
    apply_delta(coupled, reactive_at, &delta)?;

    // `updateMoles` writes the answer into the phase and the deviation is read from there, so
    // the phase's own total moves with the reactive components and its other components do not.
    for (position, &index) in reactive_at.iter().enumerate() {
        phase_amounts[index] = (old_moles[position] + raw[position]).max(MINIMUM_COUPLED_MOLES);
    }
    let new_total: f64 = phase_amounts.iter().sum();
    let mut deviation = 0.0;
    if new_total > 0.0 {
        for (position, &index) in reactive_at.iter().enumerate() {
            deviation += (old_fractions[position] - phase_amounts[index] / new_total).abs();
        }
    }
    Ok(deviation)
}

/// The conservative part of a reaction delta: `delta - A+ A delta`, with the two shortcuts in
/// front of it.
///
/// `getConservativeReactionDeltas`. A delta that already satisfies the element balance to
/// [`REACTION_DELTA_CONSERVATION_TOLERANCE_MOLES`] and leaves no inventory negative is used as it
/// stands; otherwise it is projected onto the conservation null space, which
/// [`project_onto_null_space`] computes as the component outside the matrix's row space.
/// **Measured: the shortcut is the branch the captured fluids take on every pass** - the worst
/// raw residual is `2.0e-10` - so the reactive oracle does not exercise the projection, and the
/// projection's arithmetic is pinned by its own test against the captured `A+`.
///
/// The last step scales the whole delta uniformly when an inventory would go negative, which is
/// `feasibleStep`: scaling together keeps the delta a reaction direction, where clipping one
/// component would not.
fn conservative_delta(
    a_matrix: &[Vec<f64>],
    basis: &[Vec<f64>],
    coupled: &[f64],
    reactive_at: &[usize],
    raw: &[f64],
) -> Vec<f64> {
    let worst = a_matrix
        .iter()
        .map(|row| {
            row.iter()
                .zip(raw)
                .map(|(coefficient, delta)| coefficient * delta)
                .sum::<f64>()
                .abs()
        })
        .fold(0.0_f64, f64::max);
    let non_negative = raw
        .iter()
        .zip(reactive_at)
        .all(|(delta, &index)| coupled[index] + delta >= -NEGATIVE_INVENTORY_TOLERANCE_MOLES);

    let mut delta = if worst <= REACTION_DELTA_CONSERVATION_TOLERANCE_MOLES && non_negative {
        raw.to_vec()
    } else {
        project_onto_null_space(basis, raw)
    };

    let mut step = 1.0_f64;
    for (position, value) in delta.iter().enumerate() {
        if *value >= 0.0 {
            continue;
        }
        let available = (coupled[reactive_at[position]] - MINIMUM_COUPLED_MOLES).max(0.0);
        step = step.min(available / -value);
    }
    if step < 1.0 {
        let scale = step.max(0.0);
        for value in delta.iter_mut() {
            *value *= scale;
        }
    }
    delta
}

/// Add a reaction delta to the coupled inventory, under `updateCoupledOverallComposition`'s own
/// rules: a negative result past the tolerance is a failure and not rounding, and every entry is
/// floored at [`MINIMUM_COUPLED_MOLES`].
fn apply_delta(coupled: &mut [f64], reactive_at: &[usize], delta: &[f64]) -> Result<()> {
    for (position, &index) in reactive_at.iter().enumerate() {
        let updated = coupled[index] + delta[position];
        if !updated.is_finite() || updated < -NEGATIVE_INVENTORY_TOLERANCE_MOLES {
            return Err(AzothError::InvalidInput {
                field: "moles".to_string(),
                reason: format!(
                    "aqueous chemical equilibrium produced an invalid overall amount ({updated}) \
                     for component {index}"
                ),
            });
        }
        coupled[index] = updated.max(MINIMUM_COUPLED_MOLES);
    }
    Ok(())
}

/// `ChemicalEquilibrium.calcRefPot`'s vector, recovered from `eos.pitzer_phase`.
///
/// `log_activity_i = ln phi_i - ln phi_i^reference`, where the reference is the pure component
/// for a solvent and a **two-component reference phase** - the solute at
/// [`REFERENCE_SOLUTE_MOLES`], water at [`REFERENCE_SOLVENT_MOLES`] - for everything else. That is
/// `getLogActivityCoefficient`'s own definition, and `pitzer_phase` reports both halves of it:
/// `ln_gamma` is the live arm's `ln phi` up to the branch's own factors, and `gamma_inf` is the
/// infinite-dilution coefficient the ionic arm divides by.
///
/// **Measured against the captured vector**, on two passes and every component: water matches to
/// `3e-16` as `ln_gamma` alone, every ion to `3e-15` as `ln_gamma - ln gamma_inf`, and a neutral
/// that is not water takes the same rule plus the reference phase's own `ln(m/x)` shift,
/// `-ln x_water + ln x_water^reference` - exact to the last digit of NeqSim's `1e-11` reference
/// geometry, because `m/x` is `1/(M_w x_water)` in both phases.
///
/// # What this does not claim
///
/// The neutral arm's `gamma_inf` is reported as one where `pitzer_phase` does not divide by it, so
/// a strongly non-ideal neutral keeps its full `ln gamma` where NeqSim's reference phase would
/// divide it out. No captured fluid exercises that - the brines here hold no such neutral - and
/// the spec states it as a divergence rather than an agreement.
///
/// # Errors
/// [`AzothError::InvalidInput`] if the fluid carries no water, which is the component the
/// reference is taken against.
pub fn log_activities(names: &[&str], t: f64, p: f64, composition: &[f64]) -> Result<Vec<f64>> {
    let activity = pitzer_phase(names, t, p, composition)?;
    let water = names
        .iter()
        .position(|name| name.eq_ignore_ascii_case("water"))
        .ok_or_else(|| {
            AzothError::invalid_input(
                "components",
                "the activity vector is taken against water and this fluid carries none",
            )
        })?;

    let reference_water =
        REFERENCE_SOLVENT_MOLES / (REFERENCE_SOLVENT_MOLES + REFERENCE_SOLUTE_MOLES);
    let live_water = composition[water];

    let mut out = Vec::with_capacity(names.len());
    for (index, name) in names.iter().enumerate() {
        let mut value = activity.ln_gamma[index] - activity.gamma_inf[index].ln();
        let charged = ionic_charge(name)?.is_some_and(|charge| charge != 0.0);
        if !charged && index != water {
            value += -live_water.ln() + reference_water.ln();
        }
        out.push(value);
    }
    Ok(out)
}

/// The components a fluid drives out of [`SOURCE`]'s table, **in the fluid's own order**.
///
/// NeqSim reaches this through `readReactions` -> `removeJunkReactions(componentNames)` ->
/// `getAllComponents()`: the reactions the fluid can support are kept, and the set is the union of
/// *their* components. A component the fluid carries but no surviving reaction names is not
/// reactive and must not be solved for - which is how `Ca++` and `Cl-`, both of which the element
/// table covers and neither of which the pitzer source's reactions name for a carbonate brine,
/// stay out of its chemistry.
///
/// **The order is this library's own.** NeqSim's is a `HashSet` iteration, so it is a property of
/// Java's string hash and not of the chemistry; every vector that depends on this order moves with
/// it, so the order is a presentation choice, stated rather than reproduced. On the captured
/// fluids the fluid's order is the capture's.
pub fn reactive_components(names: &[&str]) -> Result<Vec<String>> {
    let present: HashMap<&str, usize> = names
        .iter()
        .enumerate()
        .map(|(index, name)| (*name, index))
        .collect();
    let mut named: Vec<String> = Vec::new();
    for row in reactions(SOURCE)? {
        if !row.use_reaction {
            continue;
        }
        let coefficients = stoichiometry(&row.name)?;
        if !side_is_present(&coefficients, &present, true)
            && !side_is_present(&coefficients, &present, false)
        {
            continue;
        }
        for (component, _) in &coefficients {
            if present.contains_key(component.as_str()) && !named.contains(component) {
                named.push(component.clone());
            }
        }
    }
    Ok(names
        .iter()
        .filter(|name| named.iter().any(|kept| kept == *name))
        .map(|name| (*name).to_string())
        .collect())
}

/// The split's element and charge residuals against the coupled inventory it was solved at.
///
/// The quantity the capture's `conserved` prints and `SystemHybridEosGeFlashTest` asserts: the
/// coupled inventory `A n` against the same rows applied to the amounts the phases actually
/// hold, with the charge row taken over **every** component's charge rather than the reactive
/// set's - which is how that test builds it. An element row is skipped where the inventory is
/// an ion's, and the charge row is the one that covers them.
///
/// A role the solver drove to nothing contributes its own `1e-50` floor here and not an exact
/// zero, which is the same floor NeqSim's vanished phase holds.
fn conservation_residuals(
    names: &[&str],
    a_matrix: &[Vec<f64>],
    reactive_at: &[usize],
    phases: &[MultiphasePhase],
    coupled: &[f64],
) -> Result<(f64, f64)> {
    let total: f64 = coupled.iter().sum();
    let summed: Vec<f64> = (0..names.len())
        .map(|index| {
            phases
                .iter()
                .map(|phase| phase.fraction * phase.composition[index])
                .sum::<f64>()
                * total
        })
        .collect();

    let rows = a_matrix.len();
    let mut worst = 0.0_f64;
    for (row, coefficients) in a_matrix.iter().enumerate() {
        if row + 1 == rows {
            continue;
        }
        let difference: f64 = coefficients
            .iter()
            .enumerate()
            .map(|(column, coefficient)| {
                let index = reactive_at[column];
                coefficient * (coupled[index] - summed[index])
            })
            .sum();
        worst = worst.max(difference.abs());
    }

    let mut charge = 0.0;
    for (index, name) in names.iter().enumerate() {
        if let Some(value) = ionic_charge(name)? {
            charge += value * (coupled[index] - summed[index]);
        }
    }
    Ok((worst, charge.abs()))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The captured conservation matrix, whose oxygen row is dependent on the other three.
    fn captured_matrix() -> Vec<Vec<f64>> {
        vec![
            vec![1.0, 0.0, 1.0, 1.0, 0.0, 0.0],
            vec![0.0, 2.0, 1.0, 0.0, 1.0, 3.0],
            vec![2.0, 1.0, 3.0, 3.0, 1.0, 1.0],
            vec![0.0, 0.0, -1.0, -2.0, -1.0, 1.0],
        ]
    }

    /// The captured coupled inventory, in the fluid's order.
    fn captured_inventory() -> Vec<f64> {
        vec![
            5.0,
            0.049997206913475076,
            55.499994364714574,
            6.0e-4,
            2.0e-4,
            0.0010027670978232208,
            2.619878633189711e-8,
            1.1443550870656359e-8,
            2.830738529430478e-6,
        ]
    }

    /// The reactive components' positions in that inventory: `CO2, water, OH-, H3O+, HCO3-,
    /// CO3--` - the fluid's order, spectators dropped.
    const REACTIVE_AT: [usize; 6] = [1, 2, 5, 6, 7, 8];

    fn conservation(matrix: &[Vec<f64>], delta: &[f64]) -> f64 {
        matrix
            .iter()
            .map(|row| {
                row.iter()
                    .zip(delta)
                    .map(|(coefficient, value)| coefficient * value)
                    .sum::<f64>()
                    .abs()
            })
            .fold(0.0_f64, f64::max)
    }

    /// **The projection runs where the raw delta is not conservative.**
    ///
    /// No captured fluid exercises this: on all six passes of both, `max|A delta|` is
    /// `2.0e-10` against a `1e-8` shortcut, so the unprojected delta is used and the
    /// pseudo-inverse's stand-in never runs. This is the branch that does, on a delta with a
    /// component normal to the conservation space: the answer satisfies `A delta = 0` and is
    /// not the delta that went in.
    #[test]
    fn a_non_conservative_delta_is_projected_onto_the_null_space() {
        let matrix = captured_matrix();
        let basis = row_space_basis(&matrix).expect("a basis");
        let coupled = captured_inventory();
        // A unit step on water: the element rows it violates are carbon's and oxygen's.
        let raw = vec![0.0, 0.0, 1.0e-3, 0.0, 0.0, 0.0];
        assert!(
            conservation(&matrix, &raw) > REACTION_DELTA_CONSERVATION_TOLERANCE_MOLES,
            "this delta has to be the one the shortcut refuses"
        );

        let projected = conservative_delta(&matrix, &basis, &coupled, &REACTIVE_AT, &raw);
        assert!(
            conservation(&matrix, &projected) < 1.0e-16,
            "the projected delta conserves to rounding: {}",
            conservation(&matrix, &projected)
        );
        assert!(
            projected
                .iter()
                .zip(&raw)
                .any(|(value, original)| (value - original).abs() > 1.0e-6),
            "and it is not the delta that went in"
        );
    }

    /// **The shortcut is taken where the raw delta already conserves**, which is the branch
    /// the captured fluids take on every pass.
    #[test]
    fn a_conservative_delta_is_used_as_it_stands() {
        let matrix = captured_matrix();
        let basis = row_space_basis(&matrix).expect("a basis");
        let coupled = captured_inventory();
        // `A`'s own row space: its fourth row is `2 C + 0.5 H - 0.5 charge`, so a combination
        // of the rows is a direction the elements already allow.
        let mut raw = vec![0.0; 6];
        for (column, coefficient) in matrix[2].iter().enumerate() {
            raw[column] = coefficient * 1.0e-10;
        }
        assert!(conservation(&matrix, &raw) < REACTION_DELTA_CONSERVATION_TOLERANCE_MOLES);

        let delta = conservative_delta(&matrix, &basis, &coupled, &REACTIVE_AT, &raw);
        assert_eq!(delta, raw, "the delta is used as it stands");
    }

    /// **A delta that would take an inventory negative is scaled whole.**
    ///
    /// `feasibleStep`: the carbonate sits at `2.6e-8` and a step of `-1e-6` on it is not
    /// admissible. Scaling one component would take the delta out of the null space, so the
    /// whole step is scaled and the direction survives.
    #[test]
    fn a_negative_inventory_scales_the_whole_delta() {
        let matrix = captured_matrix();
        let basis = row_space_basis(&matrix).expect("a basis");
        let coupled = captured_inventory();
        let mut raw = vec![0.0; 6];
        raw[5] = -1.0e-6;

        let delta = conservative_delta(&matrix, &basis, &coupled, &REACTIVE_AT, &raw);
        for (position, value) in delta.iter().enumerate() {
            let index = REACTIVE_AT[position];
            assert!(
                coupled[index] + value >= -NEGATIVE_INVENTORY_TOLERANCE_MOLES,
                "component {index} would go negative: {} + {value}",
                coupled[index]
            );
        }
        // The step a feasible delta may take is bounded by the carbonate's own inventory over
        // the size of its negative step, so the scaled delta cannot be larger than that.
        let available = coupled[REACTIVE_AT[5]] - MINIMUM_COUPLED_MOLES;
        assert!(delta[5] >= -available * (1.0 + 1.0e-12));
        assert!(delta[5] < 0.0, "the direction survives the scaling");
        assert!(
            conservation(&matrix, &delta) < 1.0e-16,
            "and the scaled delta still conserves"
        );
    }
}
