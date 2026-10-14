//! `eos.multiflash_stability_test` - the seeding `TPmultiflash` overrides.
//!
//! Spec: `specs/models/eos/multiflash_stability_test.toml`. Oracle:
//! `validation/neqsim/captures/water_dew_point_seed_probe.tsv`.
//!
//! # Two methods, one name
//!
//! `Flash.stabilityAnalysis` is a two-trial Michelsen test - the vapour-like and liquid-like
//! trials from Wilson K-values - and **it reports instability by setting K-values; it never
//! adds a phase**. [`crate::stability_test`] is that method, and its spec names it. But
//! `TPmultiflash` **overrides** it with a different method of the same name, and
//! `TPmultiflash.run()` calls the override. This is the override.
//!
//! The difference is not a variation on the base; it is a second procedure:
//!
//! 1. **the same two Wilson trials**, run first, and the method **returns as soon as one of
//!    them lies below the plane** - so the fallback below is only reached when Wilson finds
//!    nothing;
//! 2. **a pure-component fallback**: one trial per component, each seeded at that component
//!    alone (`1.0`, every other at `1e-12`), scanned from the last component down, and
//!    **skipping every hydrocarbon but the heaviest and the lightest**;
//! 3. **a selection** that takes the first trial below the plane in that same descending
//!    order, and starts the new phase at the *trial component's own* feed fraction.
//!
//! # That fallback is the whole of two recorded divergences, measured
//!
//! `eos.tp_multiflash` seeds from [`crate::stability_test`], the base method, so it misses
//! everything only the fallback finds. Two places in this tree record the loss and neither
//! knew what it was losing: the `n2/CO2/n-octane` case, whose own spec text says a third
//! phase is reached upstream "through a second seeding mechanism this model does not port",
//! and the multiphase water dew point, which is single-phase here at the temperature
//! NeqSim calls the water dew point.
//!
//! `validation/neqsim/WaterDewPointSeedProbe.java` attributes it: on NeqSim's own
//! `WaterDewPointTest` gas at 70 bara the **base** method is what adds an aqueous phase at
//! `260.545` K (beta `1.75e-5`, water at `x = 0.9999939`), and the class's *enhanced*
//! analysis adds nothing at any state probed. Iterating that gas's pure-water trial here
//! reaches `x = 0.999994` with `tm = -9.2e-7` - inside the `-1e-8` the selection needs, and
//! barely, which is why the dew point sits where it does and why a shallow instability is
//! the one that has to survive the arithmetic.
//!
//! # What is not carried
//!
//! The iteration upstream interleaves is **successive substitution alone** here: it
//! accelerates with a Wegstein step every seventh iteration in the Wilson pair and a
//! dominant-eigenvalue step every fifth in the fallback loop, and past `150` steps the
//! fallback switches to a Newton solve on composition derivatives. None of the three is
//! carried, which is the omission [`crate::stability_test`] already makes against the base
//! method's own acceleration; the spec says so rather than leaving it to be discovered. The
//! measured states converge in a single step, so nothing here rests on it.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result, apply_checks};

use crate::algorithm_of;
use crate::mixture::{Mixture, ReducedParameters, RootSide, normalise, wilson_k};
use crate::model_gen;
use crate::results::{MultiflashSeedPath, MultiflashStabilityTestResult};
use crate::stability_test::{feed_state, ln};

/// A component at or below this feed fraction is not in the fluid at all.
///
/// Upstream's `1e-100`, used for the reference potentials, the eligibility of a
/// pure-component trial, and the tangent-plane sum.
const PRESENT: f64 = 1.0e-100;

/// The logarithm a component that is not present carries.
///
/// `exp(SENTINEL)` underflows to zero, so an absent component contributes nothing to the
/// tangent-plane sum and stays at zero in every trial composition - which is what upstream's
/// own `-10000.0` does.
const SENTINEL: f64 = -10_000.0;

/// The Wilson K-value floor, upstream's.
const WILSON_FLOOR: f64 = 1.0e-20;

/// Below the largest `|ln K|` this, the Wilson pair is skipped.
///
/// Upstream's guard against a feed near or above its critical point, where the correlation
/// says nothing. It also skips them outright for two electrolyte systems this crate does not
/// carry, which is why the rule here is the condition alone.
const WILSON_NEUTRAL: f64 = 0.01;

/// The mole number every component but the trial's own is seeded at, in a pure-component
/// trial: upstream's `1e-12`.
const TRACE: f64 = 1.0e-12;

/// A trial this close to the feed, in sum-absolute composition, **is** the feed.
///
/// Upstream's `xTrivialCheck`, and a trivial trial is marked stable rather than unstable:
/// the plane passes through the feed, so a trial that reproduces it has a tangent-plane
/// distance of zero by construction.
const TRIVIAL: f64 = 1.0e-04;

/// Below this a trial is below the tangent plane, and the phase it found is real.
const TM_LIMIT: f64 = -1.0e-08;

/// The tangent-plane distance a trial that failed or went trivial is marked with.
///
/// Upstream's `10.0`, chosen to be far above any real distance so the selection cannot pick
/// it. Zero would not do: a skipped component's `tm` stays at Java's default `0.0`, and the
/// selection reads every entry.
const FAILED_TM: f64 = 10.0;

/// Steps the Wilson pair may take, upstream's `maxiter`.
const WILSON_CAP: u32 = 50;

/// Steps a pure-component trial may take, upstream's `maxiter`.
const PURE_CAP: u32 = 200;

/// The change in `ln W` a trial must fall below to have settled.
///
/// **A sum of absolute changes, not their root mean square** - upstream's `err` - which is
/// what `eos.stability_test`'s own iteration deliberately is not.
const TOLERANCE: f64 = 1.0e-09;

/// One trial's outcome: its composition, `exp` of its converged mole numbers summed, and
/// the steps it took.
struct Settled {
    composition: Vec<f64>,
    totals: f64,
    iterations: u32,
    failed: bool,
}

/// Iterate one trial to a stationary point of the tangent-plane distance.
///
/// Unlike [`crate::stability_test`]'s own iteration this **does not raise** when the cap is
/// reached: upstream computes the distance from whatever iterate it has and leaves the
/// decision to the threshold, so a trial that never settled is reported rather than refused.
fn trial(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    d: &[f64],
    valid: &[bool],
    seed: &[f64],
    side: RootSide,
    cap: u32,
) -> Result<Settled> {
    let mut w = seed.to_vec();
    normalise(&mut w);
    let mut ln_w: Vec<f64> = w.iter().map(|&value| ln(value)).collect();
    let mut previous = f64::INFINITY;
    let mut iterations = 0;
    let mut totals = 0.0;

    for step in 1..=cap {
        iterations = step;
        let state = match mixture.phase_state(reduced, &w, side) {
            Ok(state) => state,
            // Upstream catches an init failure, marks the trial and moves on. So does this:
            // a state the cubic cannot evaluate is a trial that says nothing, not a reason
            // to refuse the flash.
            Err(AzothError::OutOfRange { .. }) => {
                return Ok(Settled {
                    composition: w,
                    totals: 0.0,
                    iterations,
                    failed: true,
                });
            }
            Err(other) => return Err(other),
        };
        let new: Vec<f64> = ln_w
            .iter()
            .enumerate()
            .map(|(index, _)| {
                if valid[index] {
                    d[index] - state.ln_phi[index]
                } else {
                    SENTINEL
                }
            })
            .collect();
        let error: f64 = new
            .iter()
            .zip(&ln_w)
            .map(|(&fresh, &old)| (fresh - old).abs())
            .sum();
        let mole_numbers: Vec<f64> = new.iter().map(|&value| value.exp()).collect();
        totals = mole_numbers.iter().sum();
        w = mole_numbers;
        normalise(&mut w);
        ln_w = new;
        if error <= TOLERANCE && error <= previous {
            return Ok(Settled {
                composition: w,
                totals,
                iterations,
                failed: false,
            });
        }
        previous = error;
    }

    Ok(Settled {
        composition: w,
        totals,
        iterations,
        failed: false,
    })
}

/// Whether a component counts as a hydrocarbon for the fallback's skip rule.
///
/// NeqSim's `isHydrocarbon` is `HC` **or** a TBP or plus fraction, and this crate has no
/// pseudo-components, so the class is the whole of it. It is *not* the predicate
/// [`crate::hydrate_inhibitor_wt`]'s phase label uses, which counts `isInert()` on the same
/// side - the two are different questions and get different answers for carbon dioxide.
fn is_hydrocarbon(mixture: &Mixture, index: usize) -> bool {
    mixture.components()[index].class == "hc"
}

/// Which trials the pure-component fallback runs, and from which seed.
///
/// Upstream skips a component whose feed fraction is absent, one that carries an ion, and
/// **every hydrocarbon but the heaviest and the lightest** - the two a hydrocarbon mixture's
/// ends are named by, so a mixture of twenty alkanes still gets three trials rather than
/// twenty.
///
/// # Errors
/// [`AzothError::PropertyUnavailable`] if a hydrocarbon carries no molar mass. The ordering
/// this rule rests on is by molar mass, and a component with none has no place in it - so the
/// zero it would otherwise stand in as is refused rather than used, which is what the Python
/// reference and the boundary's `_molar_masses` guard both say too.
fn pure_trials(mixture: &Mixture, z: &[f64]) -> Result<Vec<usize>> {
    let count = mixture.len();
    let mut heaviest = 0;
    let mut lightest = 0;
    let mut heaviest_mass = 0.0_f64;
    let mut lightest_mass = f64::INFINITY;
    for index in 0..count {
        if !is_hydrocarbon(mixture, index) {
            continue;
        }
        let Some(mass) = mixture.components()[index].molar_mass else {
            return Err(AzothError::property_unavailable(
                "component",
                "molar mass",
                "a card-added component needs its own molar mass",
            ));
        };
        if mass > heaviest_mass {
            heaviest_mass = mass;
            heaviest = index;
        }
        if mass < lightest_mass {
            lightest_mass = mass;
            lightest = index;
        }
    }

    let mut out = Vec::new();
    for index in (0..count).rev() {
        // Upstream also skips an ion. A `Mixture` cannot carry one - `mixture_of` refuses to
        // build a cubic over an `ION`, because the table's ionic rows carry filler critical
        // constants - so the clause has nothing to skip here and is recorded rather than
        // restated.
        if z[index] < PRESENT {
            continue;
        }
        if is_hydrocarbon(mixture, index) && index != heaviest && index != lightest {
            continue;
        }
        out.push(index);
    }
    Ok(out)
}

/// The smallest `T / Tc` over the components, which is the near-critical warning's measure.
fn min_t_over_tc(mixture: &Mixture, t: ThermodynamicTemperature) -> f64 {
    mixture
        .components()
        .iter()
        .map(|component| t.value / component.tc.value)
        .fold(f64::INFINITY, f64::min)
}

/// The index of the largest entry, upstream's `dominantComp`.
/// The seeding `TPmultiflash` overrides `Flash.stabilityAnalysis` with.
///
/// Returns the phase to add - its composition and the fraction to start it at - or that the
/// feed is stable as it stands.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
/// * [`AzothError::InvalidInput`] if `z` is the wrong length, has a negative entry, or does
///   not sum to one.
/// * Propagates the kernels' range checks.
pub fn multiflash_stability_test(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<MultiflashStabilityTestResult> {
    let spec = &model_gen::MULTIFLASH_STABILITY_TEST_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "P" => Some(p.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let count = mixture.len();
    if z.len() != count {
        return Err(AzothError::invalid_input(
            "z",
            format!("a feed for {count} components has {} entries", z.len()),
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
            format!("the mole fractions sum to {sum}, not to one"),
        ));
    }

    let algorithm = algorithm_of(spec)?;
    let reduced = mixture.reduced_parameters(t, p)?;
    warnings.extend(reduced.warnings.iter().cloned());

    let valid: Vec<bool> = z.iter().map(|&value| value > PRESENT).collect();
    let present: Vec<f64> = z
        .iter()
        .map(|&value| if value > PRESENT { value } else { 1.0 })
        .collect();

    // The reference potentials, from the feed on its lower-Gibbs root.
    let (feed, _) = feed_state(mixture, &reduced, z)?;
    let d: Vec<f64> = z
        .iter()
        .zip(&feed.ln_phi)
        .zip(&present)
        .map(|((&zi, &lp), &held)| if zi > PRESENT { ln(held) + lp } else { 0.0 })
        .collect();

    let mut tm: Vec<f64> = Vec::new();
    let mut iterations = 0;

    // The Wilson pair, which upstream runs first and returns from.
    let k = wilson_k(mixture, t, p);
    let seeded: Vec<(f64, RootSide, Vec<f64>)> = {
        let liquid: Vec<f64> = (0..count)
            .map(|index| {
                if valid[index] {
                    z[index] / k[index].max(WILSON_FLOOR)
                } else {
                    0.0
                }
            })
            .collect();
        let vapour: Vec<f64> = (0..count)
            .map(|index| {
                if valid[index] {
                    z[index] * k[index].max(WILSON_FLOOR)
                } else {
                    0.0
                }
            })
            .collect();
        vec![
            (z[dominant(&liquid)], RootSide::Liquid, liquid),
            (z[dominant(&vapour)], RootSide::Vapour, vapour),
        ]
    };
    let neutral = k
        .iter()
        .zip(&valid)
        .filter(|pair| *pair.1)
        .map(|(&value, _)| value.max(WILSON_FLOOR).ln().abs())
        .fold(0.0_f64, f64::max);
    if neutral >= WILSON_NEUTRAL {
        for (fraction, side, composition) in seeded {
            let settled = trial(
                mixture,
                &reduced,
                &d,
                &valid,
                &composition,
                side,
                WILSON_CAP.min(algorithm.max_iterations),
            )?;
            iterations += settled.iterations;
            let distance = seed_of(&settled, z);
            tm.push(distance);
            if distance < TM_LIMIT {
                return Ok(MultiflashStabilityTestResult {
                    path: if side == RootSide::Liquid {
                        MultiflashSeedPath::WilsonLiquid
                    } else {
                        MultiflashSeedPath::WilsonVapour
                    },
                    tm,
                    seed_composition: settled.composition,
                    seed_fraction: fraction,
                    seed_component_index: -1,
                    iterations,
                    min_t_over_tc: min_t_over_tc(mixture, t),
                    warnings,
                });
            }
        }
    }

    // The pure-component fallback.
    for index in pure_trials(mixture, z)? {
        let composition: Vec<f64> = (0..count)
            .map(|other| {
                if !valid[other] {
                    0.0
                } else if other == index {
                    1.0
                } else {
                    TRACE
                }
            })
            .collect();
        let settled = trial(
            mixture,
            &reduced,
            &d,
            &valid,
            &composition,
            RootSide::Liquid,
            PURE_CAP,
        )?;
        iterations += settled.iterations;
        let distance = seed_of(&settled, z);
        tm.push(distance);
        if distance < TM_LIMIT {
            return Ok(MultiflashStabilityTestResult {
                path: MultiflashSeedPath::PureComponent,
                tm,
                seed_composition: settled.composition,
                seed_fraction: z[index],
                seed_component_index: index as i64,
                iterations,
                min_t_over_tc: min_t_over_tc(mixture, t),
                warnings,
            });
        }
    }

    Ok(MultiflashStabilityTestResult {
        path: MultiflashSeedPath::Stable,
        tm,
        seed_composition: Vec::new(),
        seed_fraction: 0.0,
        seed_component_index: -1,
        iterations,
        min_t_over_tc: min_t_over_tc(mixture, t),
        warnings,
    })
}

/// The distance a trial reports: `1 - sum(W)`, or `FAILED_TM` where it failed or went
/// trivial.
///
/// Upstream's `xTrivialCheck` on a single-phase fluid - the sum-absolute distance from the
/// feed - is the whole of the trivial test here.
fn seed_of(settled: &Settled, z: &[f64]) -> f64 {
    let distance: f64 = settled
        .composition
        .iter()
        .zip(z)
        .map(|(&a, &b)| (a - b).abs())
        .sum();
    if settled.failed || distance < TRIVIAL {
        FAILED_TM
    } else {
        1.0 - settled.totals
    }
}

/// The index of the largest entry, upstream's `dominantComp`.
fn dominant(values: &[f64]) -> usize {
    let mut best = 0;
    for (index, &value) in values.iter().enumerate() {
        if value > values[best] {
            best = index;
        }
    }
    best
}
