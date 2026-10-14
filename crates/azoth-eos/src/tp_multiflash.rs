//! `eos.tp_multiflash` - how many phases a feed splits into, and how much of each.
//!
//! Spec: `specs/models/eos/tp_multiflash.toml`. NeqSim's `TPmultiflash`, reached through
//! `TPflash.runInternal()` when `system.doMultiPhaseCheck()` is true - so the plain
//! [`crate::pt_flash`] is this model with the stability seeding left out, and the two agree
//! everywhere the seeding has nothing to add.
//!
//! **The seeding is the whole of the difference, and it is measured rather than assumed.**
//! `validation/neqsim/TpMultiFlashSweep.java` grids sixteen mixtures over 207 states each -
//! 3,312 states, 6,624 flashes - and the flag changes the answer at exactly three kinds of
//! place: a feed the two-phase flash returned as single-phase or trivial that a tangent-plane
//! trial finds unstable (`C1/nC4 50/50` at 370 K / 40 bar, one state in 207); a state that
//! splits **three** ways where the plain flash splits two (`CO2/C1/nc10 40/30/30`, 12 states;
//! `N2/CO2/n-octane 10/40/50`, 9 states, both CO2-rich, cold and at low pressure); and a
//! water-bearing feed, where the seeding adds a hydrocarbon liquid to an aqueous one (139 of
//! 207 for water/n-hexane). Every other mixture is identical to the last bit or differs by
//! under `3e-9`.
//!
//! So this is three steps, in the order upstream runs them:
//!
//! 1. **The two-phase flash**, [`crate::pt_flash`], for the phases and a starting split.
//! 2. **The stability trial**, [`crate::stability_test`]: the stationary points of the
//!    tangent-plane distance. A trial below the plane whose composition is *not* one of the
//!    existing phases is a phase the flash has not found, and it is added at upstream's
//!    starting fraction - the feed's mole fraction of the trial's largest component.
//! 3. **The fraction solve**, [`crate::multiphase::solve_phase_fractions`], then the merge.
//!
//! **What is not here.** Upstream's 3,816 lines are mostly rescues for systems this crate
//! does not have: the ionic rules (P8), the CPA rules (P7) and the hydrate coupling (P9). What
//! it also has, and what a water-bearing cubic feed needs, is **four** further ways a phase can
//! enter the set. This paragraph used to name two of them as "the aqueous seeds" and to say
//! they were "gated on a water component and on nothing else"; both halves were wrong, so each
//! is written out here:
//!
//! 1. **`stabilityAnalysisEnhanced()`** - upstream's **second public stability analysis**,
//!    which `run()` calls where the basic one found nothing. It adds Wilson-K, `1/K` and an
//!    acentric-perturbation "LLE-specific" trial, does not skip non-hydrocarbon components, and
//!    tests against every existing phase rather than phase 0. **Not ported, and named nowhere in
//!    this tree until now.**
//! 2. **`seedAdditionalPhaseFromFeed()`** - seeds a **GAS** phase. It fires only where an
//!    **AQUEOUS** phase already exists, no gas phase sits above `beta = 1e-6`, water's feed
//!    fraction is at least `1e-4` and some hydrocarbon's is at least `1e-4`.
//! 3. **`seedHydrocarbonLiquidFromFeed()`** - seeds an **OIL** phase, and it too requires an
//!    existing **AQUEOUS** one, along with `waterZ >= 1e-6` and heavy hydrocarbons - molar mass
//!    above `0.045` kg/mol - totalling at least `5e-3`.
//! 4. The **aqueous** seed, which is inline in `run()` rather than a method: it fires where no
//!    aqueous phase exists at all and `waterZ > 1e-6`, prescribing one from the feed - `1e-16`
//!    by default, water at `max(waterZ, 1e-12)`, other non-hydrocarbon non-inert components at
//!    `min(z * 1e-2, 1e-8)` - at `beta = 1e-5`.
//!
//! **2 and 3 cannot fire in the state that matters here** - a gas carrying ppm water, with one
//! phase - because both require an aqueous phase to be there already. 4 can, but only once
//! `multiPhaseTest` is true, and only 1 or the basic analysis can set that. So a feed carrying
//! water comes back with the phases the tangent-plane trial finds, and the measured 139-of-207
//! divergence for water/n-hexane is upstream's seeding reaching states this does not; the spec
//! says so.
//!
//! **And that is what the multiphase water dew point was waiting on - which is not one of the
//! four above.** Measured on the gas NeqSim's own `WaterDewPointTest` builds, at 70 bara, this
//! flash was **single-phase at `260.545` K**, the very temperature NeqSim's dew point puts an
//! aqueous phase at, and the phase appeared here only at `217` K - forty-three kelvin colder.
//! `validation/neqsim/WaterDewPointSeedProbe.java` attributed it: neither the enhanced
//! analysis nor either named seed, but that **`TPmultiflash` overrides `stabilityAnalysis`**
//! with a second method of the same name - the same Wilson pair plus a **pure-component
//! fallback** - and this tree had ported the *base*. That override is
//! [`crate::multiflash_stability_test`], and with it this flash answers at the same temperature
//! NeqSim does and `eos.water_dew_point_multiphase` is an id.
//!
//! **A stalled two-phase flash no longer stops the seeding either.** Upstream's flash is where
//! this model's phase set *starts* and not what it answers, and there are states - a water
//! dew point's incipient aqueous phase among them - where it reaches its cap and stops short:
//! measured, `300` steps to `5.3e-5` against the `1e-10` it declares. It is now a warning and
//! the set falls back to the feed, which is what makes `eos.water_dew_point_line`'s sweep
//! possible at all.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result, Warning, WarningCode, apply_checks};

use crate::algorithm_of;
use crate::mixture::{Mixture, PhaseState, ReducedParameters, RootSide};
use crate::model_gen;
use crate::multiflash_stability_test::multiflash_stability_test;
use crate::multiphase::{MultiphasePhase, MultiphaseSplit, PhaseKind, solve_phase_fractions};
use crate::pt_flash::pt_flash;
use crate::results::{MultiflashSeedPath, Phase, TpMultiflashResult, TpMultiflashSeed};
use crate::stability_test::feed_state;

/// A phase below this is dropped. Upstream's `1.1 * phaseFractionMinimumLimit`.
///
/// The solve floors a fraction at `1e-12` - [`crate::multiphase`]'s `FRACTION_FLOOR`, which is
/// upstream's `phaseFractionMinimumLimit` - so a phase *pinned* there is one the Newton could
/// not give an amount to. This sits a hair above that floor so a pinned phase is dropped while
/// one the solve settled just above it is kept.
const FRACTION_DROP: f64 = 1.1e-12;

/// Two phases at the same composition are one phase, **whatever roots they were labelled
/// with**: a seeded phase can land on a composition that sits on one root only, where the
/// labels differ and the phases do not. Upstream reaches the same conclusion from a density
/// comparison instead.
const DUPLICATE_TOLERANCE: f64 = 1.0e-06;

/// The phases a two-phase flash reported, as a starting set for the fraction solve.
fn phases_of(
    flash: &crate::results::PtFlashResult,
    mixture: &Mixture,
    reduced: &ReducedParameters,
    z: &[f64],
) -> Result<Vec<MultiphasePhase>> {
    match flash.vapour_fraction {
        Some(beta) => Ok(vec![
            MultiphasePhase {
                fraction: 1.0 - beta,
                composition: flash.x.clone(),
                kind: PhaseKind::Cubic(RootSide::Liquid),
            },
            MultiphasePhase {
                fraction: beta,
                composition: flash.y.clone(),
                kind: PhaseKind::Cubic(RootSide::Vapour),
            },
        ]),
        // No fraction: either every K-value was on one side, or the iteration reached
        // `x = y = z`. The first names the phase; the second does not, so the feed is placed
        // on whichever root of the cubic has the lower Gibbs energy - the same rule
        // `stability_test` uses for the state the tangent plane is measured from.
        None => {
            let side = match flash.phase {
                Phase::AllLiquid => RootSide::Liquid,
                Phase::AllVapour => RootSide::Vapour,
                _ => feed_state(mixture, reduced, z)?.1,
            };
            Ok(vec![MultiphasePhase {
                fraction: 1.0,
                composition: z.to_vec(),
                kind: PhaseKind::Cubic(side),
            }])
        }
    }
}

/// `sum_i |a_i - b_i|`, the trivial-solution distance upstream tests against `1e-4`.
fn distance(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).sum()
}

/// The multiphase flash: the two-phase flash, the tangent-plane trial, and the fraction solve.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `z` is the wrong length, negative, or does not sum to
///   one.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
/// * [`AzothError::SolverNotConverged`] if the fraction solve's Hessian is singular, which
///   means two phases have met. Reaching the cap is **not** an error: the split comes back
///   with its residual and a warning, as upstream returns it.
///
/// # Example
/// ```
/// use azoth_eos::Cubic;
/// use azoth_eos::databank::mixture_of;
/// use azoth_eos::tp_multiflash;
/// use azoth_core::units::{kelvins, pascals};
///
/// let names = ["CO2", "methane", "nc10"];
/// let (mixture, _) = mixture_of(&names, Cubic::Pr, None)?;
/// let r = tp_multiflash(&mixture, kelvins(200.0), pascals(1.0e6), &[0.4, 0.3, 0.3])?;
/// assert_eq!(r.phase_count, 3);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T`, `P` and `z` are the symbols in the chemistry
pub fn tp_multiflash(
    mixture: &Mixture,
    T: ThermodynamicTemperature,
    P: Pressure,
    z: &[f64],
) -> Result<TpMultiflashResult> {
    let spec = &model_gen::TP_MULTIFLASH_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T.value),
            "P" => Some(P.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let n = mixture.len();
    if z.len() != n {
        return Err(AzothError::invalid_input(
            "z",
            format!("a feed for {n} components has {} entries", z.len()),
        ));
    }
    if let Some(bad) = z.iter().position(|&value| value < 0.0) {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "z[{bad}] is {} but a mole fraction cannot be negative",
                z[bad]
            ),
        ));
    }
    let sum: f64 = z.iter().sum();
    if (sum - 1.0).abs() > 1.0e-09 {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "the mole fractions sum to {sum}, not to one. Renormalising them here \
                 would make a composition error invisible in every number downstream, \
                 so it is refused instead"
            ),
        ));
    }

    let algorithm = algorithm_of(spec)?;
    let reduced = mixture.reduced_parameters(T, P)?;

    // **A stalled two-phase flash is warned about rather than propagated.** Upstream's flash is
    // the *seed* this model's seeding starts from and not its answer, and there are states -
    // a water dew point's incipient aqueous phase among them - where it reaches its cap and
    // stops short: measured, `300` steps to a residual of `5.3e-5` against the `1e-10` it
    // declares, at `261.876` K and `111.6` bara on the gas `eos.water_dew_point_line` sweeps.
    // Refusing there loses an answer the seeding and the fraction solve would have found, and
    // `solve_phase_fractions`' own cap is already a warning for exactly that reason. The set
    // falls back to the feed on its lower-Gibbs root, which is the same set `phases_of` builds
    // where the flash reports no fraction at all.
    let (mut phases, flash_min_t_over_tc) = match pt_flash(mixture, T, P, z) {
        Ok(flash) => {
            warnings.extend(flash.warnings.iter().cloned());
            let phases = phases_of(&flash, mixture, &reduced, z)?;
            (phases, flash.min_t_over_tc)
        }
        Err(AzothError::SolverNotConverged {
            iterations,
            residual,
            tolerance,
        }) => {
            warnings.push(Warning::new(
                WarningCode::SolverNotConverged,
                format!(
                    "the two-phase flash reached its cap after {iterations} steps with a \
                     residual of {residual:.3e} against the {tolerance:.3e} it declares, so the \
                     seeding below starts from the feed instead of from the split it nearly \
                     reached"
                ),
            ));
            let (_, side) = feed_state(mixture, &reduced, z)?;
            (
                vec![MultiphasePhase {
                    fraction: 1.0,
                    composition: z.to_vec(),
                    kind: PhaseKind::Cubic(side),
                }],
                f64::INFINITY,
            )
        }
        Err(other) => return Err(other),
    };

    // The seeding. **This is `TPmultiflash`'s own `stabilityAnalysis` and not the base
    // method `eos.stability_test` ports** - the two share a name and nothing else - and the
    // difference is the whole of what this model was missing: the override adds a
    // pure-component trial per eligible component, which is what finds a water-rich phase
    // in a gas carrying ppm water and a CO2-rich one in the N2/CO2/n-octane case.
    let seeding = multiflash_stability_test(mixture, T, P, z)?;
    warnings.extend(seeding.warnings.iter().cloned());

    // The phase count the two-phase flash reached, which is what `seeded` is measured against:
    // a phase added and then merged away left the answer the flash's own, and reporting it as
    // seeded would say the model found something it did not keep.
    let base_count = phases.len();
    let mut seeded = false;
    if seeding.path != MultiflashSeedPath::Stable {
        let (_, side) = feed_state(mixture, &reduced, &seeding.seed_composition)?;
        phases.push(MultiphasePhase {
            fraction: seeding.seed_fraction,
            composition: seeding.seed_composition.clone(),
            kind: PhaseKind::Cubic(side),
        });
        seeded = true;
    }

    // **Nothing seeded, nothing to solve.** When the trial added no phase the answer is the
    // two-phase flash's, and re-solving the fractions it already converged would replace a
    // converged split with an iterate of the same equations - measurably, it moved methane/
    // n-butane's 250 K betas in the tenth digit and took 41 steps to do it. Upstream's own
    // numbers say the same: with the flag on and off, every state the sweep found no seeding
    // at reports identical betas to all seventeen digits.
    let split = if seeded {
        let total: f64 = phases.iter().map(|phase| phase.fraction).sum();
        for phase in &mut phases {
            phase.fraction /= total;
        }
        solve_phase_fractions(mixture, &reduced, z, &mut phases, algorithm)?
    } else {
        MultiphaseSplit {
            fractions: phases.iter().map(|phase| phase.fraction).collect(),
            compositions: phases
                .iter()
                .map(|phase| phase.composition.clone())
                .collect(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        }
    };

    // The merge. A phase the solve pinned at the floor is a phase that is not there.
    let mut kept: Vec<MultiphasePhase> = Vec::with_capacity(phases.len());
    for phase in phases.iter() {
        if phase.fraction < FRACTION_DROP {
            continue;
        }
        // Two phases at the same composition are one phase however their roots were labelled.
        //
        // **The side is deliberately not part of the test.** A phase the seeding added is put
        // on the lower-Gibbs root at the *trial* composition, and the solve can then drive it
        // onto a composition that sits on one root only - where the two sides are the same
        // number and the labels differ while the phases do not. Requiring the same side there
        // reports one phase twice, which is what this merge exists to prevent; measured, it
        // did exactly that at N2/CO2/n-octane 180 K / 10 bar and methane/n-butane 250 K /
        // 50 bar before the test was dropped.
        //
        // Upstream reaches the same conclusion from the other end - its rule is a density
        // comparison, and it folds the pair when the compositions agree. The mass moves into
        // the phase already kept, which is upstream's `mergeAndRemoveDuplicatePhase`; dropping
        // it instead would lose the material.
        if let Some(existing) = kept
            .iter_mut()
            .find(|other| distance(&other.composition, &phase.composition) < DUPLICATE_TOLERANCE)
        {
            existing.fraction += phase.fraction;
            continue;
        }
        kept.push(phase.clone());
    }

    // A dropped phase's share returns to the others, which is upstream's
    // `removePhaseKeepTotalComposition`: removing a phase is a statement that it is not there,
    // not that a part of the feed has gone. Measured, without this the fractions sum to
    // `1 - 1e-12` - the floor of the phase that was dropped - and the answer is a feed with a
    // component missing.
    let total: f64 = kept.iter().map(|phase| phase.fraction).sum();
    if total > 0.0 {
        for phase in &mut kept {
            phase.fraction /= total;
        }
    }

    // Whether the solve contributed to the answer, which is the same predicate `seeded`
    // reports. When the merge folded the seeded phase away the answer is the two-phase flash's
    // and the solve produced nothing - so the steps it took and the residual it left are not
    // properties of the answer, and reporting them would report work that was discarded.
    // Measured: at N2/CO2/n-octane 190 K / 1 bar the seeding fires and is then merged away, and
    // the discarded solve took 32 steps in one kernel and 34 in the other.
    let grew = kept.len() > base_count;
    let iterations = if grew { split.iterations } else { 0 };
    let residual = if grew { split.residual } else { 0.0 };
    if grew && !split.converged {
        warnings.push(Warning::new(
            WarningCode::SolverNotConverged,
            format!(
                "the fraction solve left by its cap after {iterations} steps with a residual \
                 of {residual:.3e}, against the {} the spec declares. The fractions are the \
                 iterate it reached, not a state it settled at, and upstream accepts this \
                 rather than refusing the flash",
                algorithm.tolerance
            ),
        ));
    }
    let mut geometry = Vec::with_capacity(kept.len());
    for phase in &kept {
        geometry.push(state_of(mixture, &reduced, phase)?);
    }

    Ok(TpMultiflashResult {
        phase_count: kept.len() as u32,
        phase_fractions: kept.iter().map(|phase| phase.fraction).collect(),
        x: kept.iter().map(|phase| phase.composition.clone()).collect(),
        z_factor: geometry
            .iter()
            .map(|state| state.as_ref().map_or(f64::NAN, |state| state.z))
            .collect(),
        ln_phi: geometry
            .iter()
            .map(|state| {
                state
                    .as_ref()
                    .map_or_else(Vec::new, |state| state.ln_phi.clone())
            })
            .collect(),
        seeded: if grew {
            TpMultiflashSeed::StabilitySeeded
        } else {
            TpMultiflashSeed::TwoPhaseFlash
        },
        tm: seeding.tm.clone(),
        iterations,
        residual,
        min_t_over_tc: flash_min_t_over_tc.min(seeding.min_t_over_tc),
        warnings,
    })
}

/// The state of a phase, or `None` if its composition is degenerate.
///
/// A phase the merge left at a composition the cubic cannot evaluate - which means a
/// component it does not have, so `ln` of a zero - still has a fraction worth reporting, so
/// the geometry is optional rather than the phase being dropped.
fn state_of(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    phase: &MultiphasePhase,
) -> Result<Option<PhaseState>> {
    let PhaseKind::Cubic(side) = phase.kind else {
        // A wax phase is not a root of the cubic, so there is no `PhaseState` for it and
        // nothing here to report: this routine exists to hand a *cubic* state to the
        // comparison that follows, and a caller that put a solid in the set has already
        // said it is not one.
        return Ok(None);
    };
    match mixture.phase_state(reduced, &phase.composition, side) {
        Ok(state) => Ok(Some(state)),
        Err(AzothError::OutOfRange { .. }) => Ok(None),
        Err(other) => Err(other),
    }
}
