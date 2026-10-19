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
//! water-bearing feed, where the flag's difference is common (139 of 207 for water/n-hexane).
//! Every other mixture is identical to the last bit or differs by under `3e-9`.
//!
//! `validation/neqsim/TpMultiflashSeedProbe.java` attributes the water-bearing difference: at
//! water/n-hexane 350 K / 20 bar the **basic analysis** adds the aqueous phase, and neither
//! named feed seed survives the merge anywhere in the probe's grid. Their direct invocation is
//! compared against NeqSim's composition by this module's unit tests.
//!
//! So this is three steps, in the order upstream runs them:
//!
//! 1. **The two-phase flash**, [`crate::pt_flash`], for the phases and a starting split.
//! 2. **The stability trial**, [`crate::multiflash_stability_test`]: the stationary points of
//!    the tangent-plane distance, measured from **phase 0 of the flashed set** - the
//!    least-dense phase, as upstream orders phases by density - with a trial trivial when it
//!    reproduces **any** existing phase. A trial below the plane is a phase the flash has not
//!    found, and it is added at upstream's starting fraction - the feed's mole fraction of the
//!    trial's largest component.
//! 3. **The fraction solve**, [`crate::multiphase::solve_phase_fractions`], then the merge.
//!
//! **The other ways a phase can enter.** Upstream's 3,816 lines are mostly rescues for systems
//! this crate does not have: the ionic rules (P8), the CPA rules (P7) and the hydrate coupling
//! (P9). Beside those it has four further seeding constructs, all four ported here:
//!
//! 1. **`stabilityAnalysisEnhanced()`** - upstream's second public analysis, run **after** the
//!    fraction solve on a two-phase result. It tries Wilson `K`, `1/K` and an
//!    acentric-perturbation trial per existing phase and re-solves when one is admitted. Its
//!    opt-in pre-solve call needs `setEnhancedMultiPhaseCheck(true)`, a system switch this
//!    model has no input for.
//! 2. **`seedAdditionalPhaseFromFeed()`** - seeds a **GAS** where an aqueous phase exists, no
//!    gas sits above `beta = 1e-6`, water is at least `1e-4` of the feed and a hydrocarbon is.
//! 3. **`seedHydrocarbonLiquidFromFeed()`** - seeds an **OIL**, also requiring an aqueous
//!    phase, with `waterZ >= 1e-6` and heavy hydrocarbons (molar mass above `0.045` kg/mol)
//!    totalling at least `5e-3` and outweighing the water.
//! 4. The **aqueous** seed, inline in `run()` rather than a method: it is the only one that can
//!    fire on a feed with no aqueous phase already, and it does not re-solve.
//!
//! The phase-type gates are `hydrate_inhibitor_wt::label`'s, which is NeqSim's `PhaseEos.init`
//! rule. **Measured at every state the probe scans: neither named seed's addition survives the
//! merge**, so their end-to-end effect is the full-chain answers the cases already pin; their
//! arithmetic is compared to NeqSim's own compositions by the unit tests at the end of this
//! file.
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
use crate::databank::ION;
use crate::hydrate_inhibitor_wt::{PhaseLabel, label};
use crate::mixture::{Mixture, PhaseState, ReducedParameters, RootSide, wilson_k};
use crate::model_gen;
use crate::multiflash_stability_test::{enhanced_against_phases, seeding_against_reference};
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
    // **A two-phase flash's pair, and only its pair.** A *negative* flash reports an
    // extrapolated composition rather than a state - the amount of the absent phase that
    // would bring the feed to saturation - so seeding from it starts a one-phase answer at a
    // composition the feed does not have. Measured, methane/n-heptane at 293.15 K and 150
    // bara came back as one phase at `x_methane = 0.547` against a feed of `0.0265`.
    match (flash.phase, flash.vapour_fraction) {
        (Phase::TwoPhase, Some(beta)) => Ok(vec![
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
        // No two-phase split: either every K-value was on one side, the iteration reached
        // `x = y = z`, or the flash is negative. The first names the phase; the rest do not,
        // so the feed is placed on whichever root of the cubic has the lower Gibbs energy -
        // the same rule `stability_test` uses for the state the tangent plane is measured
        // from.
        _ => {
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
    //
    // **Measured from the least-dense existing phase, as upstream is.** `TPflash` orders
    // phases by density and the analysis reads its reference potentials from phase 0; the
    // two-phase set here is ordered liquid first, so the reference is chosen by the largest
    // compressibility. A trial that reproduces **any** existing phase is trivial, which is
    // upstream's `xTrivialCheck0`/`xTrivialCheck1` pair. On a single-phase feed the reference
    // is the feed and this is the model's own single-phase case.
    let mut reference: Option<(usize, RootSide)> = None;
    let mut reference_z = f64::NEG_INFINITY;
    for (index, phase) in phases.iter().enumerate() {
        if let PhaseKind::Cubic(side) = phase.kind
            && let Some(state) = state_of(mixture, &reduced, phase)?
            && state.z.is_finite()
            && state.z > reference_z
        {
            reference_z = state.z;
            reference = Some((index, side));
        }
    }
    let existing: Vec<Vec<f64>> = phases
        .iter()
        .map(|phase| phase.composition.clone())
        .collect();
    let seeding = if let Some((index, side)) = reference {
        seeding_against_reference(
            mixture,
            T,
            P,
            z,
            &phases[index].composition,
            side,
            &existing,
            Vec::new(),
        )?
    } else {
        let (_, side) = feed_state(mixture, &reduced, z)?;
        seeding_against_reference(mixture, T, P, z, z, side, &existing, Vec::new())?
    };
    warnings.extend(seeding.warnings.iter().cloned());

    // The phase count the two-phase flash reached, which is what `seeded` is measured against:
    // a phase added and then merged away left the answer the flash's own, and reporting it as
    // seeded would say the model found something it did not keep.
    let base_count = phases.len();
    let mut seeded = false;
    // Which construct added each phase. The report names the first surviving addition, so a
    // phase that was added and merged away does not claim the answer.
    let mut provenance: Vec<TpMultiflashSeed> = vec![TpMultiflashSeed::TwoPhaseFlash; phases.len()];
    if seeding.path != MultiflashSeedPath::Stable {
        let (_, side) = feed_state(mixture, &reduced, &seeding.seed_composition)?;
        phases.push(MultiphasePhase {
            fraction: seeding.seed_fraction,
            composition: seeding.seed_composition.clone(),
            kind: PhaseKind::Cubic(side),
        });
        provenance.push(TpMultiflashSeed::StabilitySeeded);
        seeded = true;
    }

    // **The two named feed seeds.** Both require an **aqueous** phase to exist already; the
    // gas seed runs only when the stability analysis added nothing (upstream's
    // `!multiPhaseTest`), the oil seed unconditionally. Measured by
    // `TpMultiflashSeedProbe.java`: the gas seed fires across water/n-hexane, the oil seed on
    // the eight-component water-bearing gas.
    //
    // The set before them is kept, because these are the seeds a failed solve reverts.
    let before_feed_seeds = phases.clone();
    let mut feed_seeded = false;
    if !seeded && let Some(phase) = seed_gas_from_feed(mixture, T, P, z, &phases, &reduced)? {
        phases.push(phase);
        provenance.push(TpMultiflashSeed::GasSeeded);
        seeded = true;
        feed_seeded = true;
    }
    if let Some(phase) = seed_oil_from_feed(mixture, T, P, z, &phases, &reduced)? {
        phases.push(phase);
        provenance.push(TpMultiflashSeed::OilSeeded);
        seeded = true;
        feed_seeded = true;
    }

    // **Nothing seeded, nothing to solve.** When the trial added no phase the answer is the
    // two-phase flash's, and re-solving the fractions it already converged would replace a
    // converged split with an iterate of the same equations - measurably, it moved methane/
    // n-butane's 250 K betas in the tenth digit and took 41 steps to do it. Upstream's own
    // numbers say the same: with the flag on and off, every state the sweep found no seeding
    // at reports identical betas to all seventeen digits.
    let mut split = if seeded {
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

    // **A feed seed whose solve leaves by its cap is not evidence of a phase.** The fraction
    // solve floors and renormalises each step, so a seed it cannot settle drives the whole set
    // to the floors rather than leaving the set before it alone - measured on the 1 bara water
    // dew point's marginal gas (224.2 K, 35 ppm water), where the oil seed's addition took the
    // aqueous phase to the floor with it and the answer became one phase where the two-phase
    // flash had two. The set *before* the feed seeds is restored and solved again; the basic
    // and enhanced additions are not reverted, because they are not what failed.
    if feed_seeded && !split.converged {
        phases = before_feed_seeds;
        provenance.truncate(phases.len());
        seeded = phases.len() > base_count;
        split = if seeded {
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
        warnings.push(Warning::new(
            WarningCode::FallbackToInput,
            "a feed seed's fraction solve left by its cap, so the seed was removed and the set \
             before it is the answer: a seed the solve cannot settle is not a phase",
        ));
    }

    // **Upstream's post-flash enhanced check.** `run()` runs `stabilityAnalysisEnhanced` once
    // on a two-phase multiphase result, after the fraction solve, and a phase it adds triggers
    // a bounded rerun of the whole flash. This is that call; the rerun is the re-solve below,
    // which is the observable part - the recursive `run()` also re-runs the basic analysis,
    // and on the captured state that finds nothing. Measured: C1/nC4 + 1e-4 water at 370 K /
    // 40 bar, where this check's addition is removed again by the merge, and its presence is
    // what decides whether the inline aqueous seed (not yet ported) can fire.
    if seeded && phases.len() == 2 {
        let existing: Vec<(Vec<f64>, RootSide)> = phases
            .iter()
            .filter_map(|phase| match phase.kind {
                PhaseKind::Cubic(side) => Some((phase.composition.clone(), side)),
                _ => None,
            })
            .collect();
        if let Some(seed) = enhanced_against_phases(mixture, T, P, z, &existing)? {
            phases.push(MultiphasePhase {
                fraction: seed.fraction,
                composition: seed.composition,
                kind: PhaseKind::Cubic(seed.side),
            });
            provenance.push(TpMultiflashSeed::EnhancedSeeded);
            let total: f64 = phases.iter().map(|phase| phase.fraction).sum();
            for phase in &mut phases {
                phase.fraction /= total;
            }
            split = solve_phase_fractions(mixture, &reduced, z, &mut phases, algorithm)?;
        }
    }

    // **The inline aqueous seed.** `run()` adds an AQUEOUS phase where none exists, the water
    // fraction is above `1e-6` and fewer than three phases are present, after the fraction
    // solve and the enhanced check - and **does not re-solve**: the other phases keep the
    // fractions they converged to, rescaled by `1 - beta`. Measured in the probe's ablation at
    // C1/nC4 + 1e-4 water, 370 K / 40 bar: aqueous beta `9.99990000099999e-6`, its
    // composition the feed's water over `1e-16` traces. In the default chain of that state the
    // enhanced check's addition is what suppresses this seed.
    if seeded && phases.len() < 3 {
        let labels = phase_labels(mixture, &reduced, &phases)?;
        if !labels.contains(&Some(PhaseLabel::Aqueous))
            && let Some(water) = water_index(mixture)
            && z[water] > 1.0e-6
        {
            let mut composition: Vec<f64> = Vec::with_capacity(phases.len());
            for (index, component) in mixture.components().iter().enumerate() {
                let value = if index == water {
                    z[index].max(1.0e-12)
                } else if matches!(component.class.as_str(), "hc" | "inert") {
                    1.0e-16
                } else {
                    (z[index] * 1.0e-2).min(1.0e-8)
                };
                composition.push(value);
            }
            let sum: f64 = composition.iter().sum();
            if sum > 0.0 {
                for value in &mut composition {
                    *value /= sum;
                }
            }
            phases.push(MultiphasePhase {
                fraction: 1.0e-5,
                composition,
                kind: PhaseKind::Cubic(RootSide::Liquid),
            });
            provenance.push(TpMultiflashSeed::AqueousSeeded);
            let total: f64 = phases.iter().map(|phase| phase.fraction).sum();
            for phase in &mut phases {
                phase.fraction /= total;
            }
        }
    }

    // The merge. A phase the solve pinned at the floor is a phase that is not there.
    let mut kept: Vec<MultiphasePhase> = Vec::with_capacity(phases.len());
    let mut kept_provenance: Vec<TpMultiflashSeed> = Vec::with_capacity(phases.len());
    for (index, phase) in phases.iter().enumerate() {
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
        kept_provenance.push(provenance[index]);
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
        // The **last** construct whose addition survived: when the basic analysis and the
        // inline aqueous seed both contributed, the answer's distinguishing phase is the later
        // one, and that is the construct a case is pinning.
        seeded: if grew {
            kept_provenance
                .iter()
                .rev()
                .copied()
                .find(|path| *path != TpMultiflashSeed::TwoPhaseFlash)
                .unwrap_or(TpMultiflashSeed::TwoPhaseFlash)
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

/// NeqSim's `PhaseEos.init` label for every phase, or `None` where the cubic cannot evaluate
/// one. The seeds gate on these labels; an unlabelled phase is neither a gas, an oil nor an
/// aqueous phase.
fn phase_labels(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    phases: &[MultiphasePhase],
) -> Result<Vec<Option<PhaseLabel>>> {
    let mut out = Vec::with_capacity(phases.len());
    for phase in phases {
        match state_of(mixture, reduced, phase)? {
            Some(state) => out.push(Some(label(mixture, reduced, &phase.composition, state.z)?)),
            None => out.push(None),
        }
    }
    Ok(out)
}

/// The feed index of `water`, by name.
fn water_index(mixture: &Mixture) -> Option<usize> {
    mixture
        .names()?
        .iter()
        .position(|name| name.eq_ignore_ascii_case("water"))
}

/// Molar mass above which upstream's oil seed counts a hydrocarbon as heavy, in kg/mol.
const HEAVY_HYDROCARBON: f64 = 0.045;

/// `seedAdditionalPhaseFromFeed`: a **GAS** added to a set that has an aqueous phase and no
/// gas above a trace.
///
/// The composition is the feed weighted by Wilson `K` in log space - `ln z + ln K`, clipped to
/// `[-50, 50]`, normalised by its largest entry and floored at `1e-16` - and the fraction is
/// `1e-3`. The gates are upstream's: an aqueous phase exists, no gas sits above `1e-6`, water
/// is at least `1e-4` of the feed, a hydrocarbon is at least `1e-4`, and fewer than three
/// phases are present.
fn seed_gas_from_feed(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    phases: &[MultiphasePhase],
    reduced: &ReducedParameters,
) -> Result<Option<MultiphasePhase>> {
    if phases.len() >= 3 {
        return Ok(None);
    }
    let labels = phase_labels(mixture, reduced, phases)?;
    if phases
        .iter()
        .zip(&labels)
        .any(|(phase, label)| *label == Some(PhaseLabel::Gas) && phase.fraction > 1.0e-6)
    {
        return Ok(None);
    }
    if !labels.contains(&Some(PhaseLabel::Aqueous)) {
        return Ok(None);
    }
    let Some(water) = water_index(mixture) else {
        return Ok(None);
    };
    if z[water] < 1.0e-4 {
        return Ok(None);
    }
    if !mixture
        .components()
        .iter()
        .enumerate()
        .any(|(index, component)| component.class == "hc" && z[index] > 1.0e-4)
    {
        return Ok(None);
    }

    let k = wilson_k(mixture, t, p);
    let mut logs: Vec<f64> = Vec::with_capacity(z.len());
    for (index, component) in mixture.components().iter().enumerate() {
        let mut value = z[index].max(1.0e-100).ln();
        if z[index] > 0.0 && component.tc.value > 0.0 && component.pc.value > 0.0 {
            let log_k = k[index].max(1.0e-20).ln();
            if log_k.is_finite() {
                value += log_k.clamp(-50.0, 50.0);
            }
        }
        logs.push(value);
    }
    let maximum = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut composition: Vec<f64> = logs
        .iter()
        .map(|value| (value - maximum).exp().max(1.0e-16))
        .collect();
    let total: f64 = composition.iter().sum();
    for fraction in &mut composition {
        *fraction /= total;
    }
    Ok(Some(MultiphasePhase {
        fraction: 1.0e-3,
        composition,
        kind: PhaseKind::Cubic(RootSide::Vapour),
    }))
}

/// `seedHydrocarbonLiquidFromFeed`: an **OIL** phase added to a water-bearing set.
///
/// Heavy hydrocarbons keep their feed fraction, light hydrocarbons and water are damped to
/// `min(z * 1e-2, 1e-8)`, ions and every other component sit at `1e-16`, and the fraction is
/// `1e-5`. The gates are upstream's: an aqueous phase exists, no oil does, water is at least
/// `1e-6` of the feed, and the heavy hydrocarbons total at least `5e-3` and outweigh the water.
fn seed_oil_from_feed(
    mixture: &Mixture,
    _t: ThermodynamicTemperature,
    _p: Pressure,
    z: &[f64],
    phases: &[MultiphasePhase],
    reduced: &ReducedParameters,
) -> Result<Option<MultiphasePhase>> {
    if phases.len() >= 3 {
        return Ok(None);
    }
    let labels = phase_labels(mixture, reduced, phases)?;
    if labels.contains(&Some(PhaseLabel::Oil)) {
        return Ok(None);
    }
    if !labels.contains(&Some(PhaseLabel::Aqueous)) {
        return Ok(None);
    }
    let Some(water) = water_index(mixture) else {
        return Ok(None);
    };
    if z[water] < 1.0e-6 {
        return Ok(None);
    }
    let heavy: f64 = mixture
        .components()
        .iter()
        .enumerate()
        .filter(|(index, component)| {
            component.class == "hc"
                && z[*index] > 1.0e-6
                && component
                    .molar_mass
                    .is_some_and(|mass| mass > HEAVY_HYDROCARBON)
        })
        .map(|(index, _)| z[index])
        .sum();
    if heavy < 5.0e-3 || heavy <= z[water] {
        return Ok(None);
    }

    let mut composition: Vec<f64> = Vec::with_capacity(z.len());
    for (index, component) in mixture.components().iter().enumerate() {
        let value = if component.class == ION {
            1.0e-16
        } else if component.class == "hc" {
            if component
                .molar_mass
                .is_some_and(|mass| mass > HEAVY_HYDROCARBON)
            {
                z[index].max(1.0e-12)
            } else {
                (z[index] * 1.0e-2).min(1.0e-8)
            }
        } else if mixture
            .names()
            .and_then(|names| names.get(index))
            .is_some_and(|name| name.eq_ignore_ascii_case("water"))
        {
            (z[index] * 1.0e-2).min(1.0e-8)
        } else {
            1.0e-16
        };
        composition.push(value);
    }
    let total: f64 = composition.iter().sum();
    if total > 0.0 {
        for fraction in &mut composition {
            *fraction /= total;
        }
    }
    Ok(Some(MultiphasePhase {
        fraction: 1.0e-5,
        composition,
        kind: PhaseKind::Cubic(RootSide::Liquid),
    }))
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

#[cfg(test)]
mod tests {
    use azoth_core::units::{kelvins, pascals};

    use super::*;
    use crate::Cubic;
    use crate::databank::mixture_of;

    /// The gas seed's composition on the state `TpMultiflashSeedProbe` measures it at:
    /// water/n-hexane 50/50, 350 K, 20 bar, after the basic analysis has added the aqueous
    /// phase (two phases, both aqueous-labelled, no gas).
    #[test]
    fn the_gas_seed_reproduces_neqsims_composition() {
        let (mixture, _) =
            mixture_of(&["water", "n-hexane"], Cubic::Pr, None).expect("the databank resolves");
        let reduced = mixture
            .reduced_parameters(kelvins(350.0), pascals(20.0e5))
            .expect("the state is in range");
        let phases = [
            MultiphasePhase {
                fraction: 2.0 / 3.0,
                composition: vec![0.5, 0.5],
                kind: PhaseKind::Cubic(RootSide::Liquid),
            },
            MultiphasePhase {
                fraction: 1.0 / 3.0,
                composition: vec![1.0, 0.0],
                kind: PhaseKind::Cubic(RootSide::Liquid),
            },
        ];
        let seed = seed_gas_from_feed(
            &mixture,
            kelvins(350.0),
            pascals(20.0e5),
            &[0.5, 0.5],
            &phases,
            &reduced,
        )
        .expect("the label arithmetic is available")
        .expect("the seed fires: an aqueous phase exists and no gas does");
        assert_eq!(seed.kind, PhaseKind::Cubic(RootSide::Vapour));
        assert!((seed.fraction - 1.0e-3).abs() < 1.0e-15);
        // NeqSim: water 0.26943831469243495, n-hexane 0.73056168530756500.
        assert!((seed.composition[0] - 0.269_438_314_692_434_95).abs() < 1.0e-9);
        assert!((seed.composition[1] - 0.730_561_685_307_565).abs() < 1.0e-9);
    }

    /// The oil seed on the state the probe attributes it at: the eight-component water-bearing
    /// gas at 230 K / 80 bar, with the aqueous phase the basic analysis added.
    #[test]
    fn the_oil_seed_reproduces_neqsims_composition() {
        let names = [
            "CO2", "nitrogen", "methane", "ethane", "propane", "i-butane", "n-butane", "water",
        ];
        let (mixture, _) = mixture_of(&names, Cubic::Srk, None).expect("the databank resolves");
        let x = 35.0e-6;
        let z = [0.02, 0.01, 0.9 - x, 0.05, 0.01, 0.005, 0.005, x];
        let reduced = mixture
            .reduced_parameters(kelvins(230.0), pascals(80.0e5))
            .expect("the state is in range");
        let phases = [
            MultiphasePhase {
                fraction: 1.0,
                composition: z.to_vec(),
                kind: PhaseKind::Cubic(RootSide::Vapour),
            },
            MultiphasePhase {
                fraction: 3.0e-5,
                composition: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
                kind: PhaseKind::Cubic(RootSide::Liquid),
            },
        ];
        let seed = seed_oil_from_feed(
            &mixture,
            kelvins(230.0),
            pascals(80.0e5),
            &z,
            &phases,
            &reduced,
        )
        .expect("the label arithmetic is available")
        .expect("the seed fires: an aqueous phase exists, no oil does, heavy ends are present");
        assert_eq!(seed.kind, PhaseKind::Cubic(RootSide::Liquid));
        assert!((seed.fraction - 1.0e-5).abs() < 1.0e-15);
        let expected = [
            9.999960000159797e-15,
            9.999960000159797e-15,
            9.999960000159798e-7,
            9.999960000159798e-7,
            9.999960000159798e-7,
            0.49999800000798994,
            0.49999800000798994,
            9.999960000159798e-7,
        ];
        for (got, want) in seed.composition.iter().zip(expected) {
            assert!((got - want).abs() < 1.0e-12, "{:?}", seed.composition);
        }
    }
}
