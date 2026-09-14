//! Spec-driven tests for the `eos.stability_test` model.
//!
//! The spec's cases pin the two implementations to each other. These tests are for
//! the things a case *cannot* say: the identities that hold at any answer, the
//! cross-model check against `eos.pt_flash` - which is the scientific content of
//! the pair - and the two decisions that were made rather than transcribed, each of
//! which is a wrong answer that looks like a right one.
//!
//! # What the cases do not pin here, and what pins it instead
//!
//! The generator that turns a model spec into [`model_gen`] carries a case's scalar
//! and vector expectations and has no field for a *matrix* one, so the `w` block a
//! case declares never reaches `TestCase`. `tm` and `iterations` do, and they are
//! asserted below; the compositions are covered two other ways - structurally here
//! (every row is a composition, and the trivial trial's row is the feed), and
//! exactly by `python/tests/models/test_stability_test.py`, which reads the raw
//! spec and compares both backends' `w` on every case. Rust's `w` matching Python's
//! on every case, and Python's matching the spec's, is the same claim transitively.

use azoth_core::units::{kelvins, pascals};
use azoth_core::{AzothError, CalcResult, WarningCode};
use azoth_eos::mixture::{Component, Mixture, RootSide};
use azoth_eos::{
    Phase, StabilityTestResult, StabilityVerdict, databank, model_gen, pt_flash, stability_test,
};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.stability_test";

/// The methane/n-butane pair the spec's cases use, resolved through the databank so
/// the pair the sweeps run and the pair the cases run are the same fluid, `kij`
/// included.
fn methane_butane() -> Mixture {
    databank::mixture_of(&["methane", "n-butane"], None)
        .expect("the pair resolves")
        .0
}

fn mixture_from_case(case: &azoth_core::spec::TestCase) -> Mixture {
    let names = case
        .list("components")
        .expect("the case declares components");
    databank::mixture_of(names, None)
        .expect("the case's components resolve")
        .0
}

fn verdict_case(case: &azoth_core::spec::TestCase) -> StabilityTestResult {
    let mixture = mixture_from_case(case);
    stability_test(
        &mixture,
        kelvins(common::input(case, "T")),
        pascals(common::input(case, "P")),
        case.vector("z").expect("the case declares z"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

/// Compare a trial distance, where an expectation of zero has no scale.
///
/// `common::assert_close` is *relative*, deliberately - the numbers in this registry
/// span nine orders of magnitude - and a relative comparison against an expected
/// zero divides by the smallest positive double and fails on a value that is zero to
/// within a few ulps. A trivial trial's distance is exactly that: `tm = 1 - sum(W)`
/// where `sum(W)` is one to rounding, so the distance is `0.0` or `-4.4e-16`
/// depending on which way the last `exp` rounded - and, as the spec's second case
/// records, it can land on either side. Comparing two such numbers relatively is
/// comparing round-off.
///
/// So an expectation at or below `1e-12` - which no non-trivial distance is - is
/// compared absolutely, at the case's own tolerance.
#[track_caller]
fn assert_tm(actual: f64, expected_value: f64, tolerance: f64, context: &str) {
    if expected_value.abs() <= 1e-12 {
        assert!(
            actual.abs() <= tolerance,
            "{context}: got {actual}, expected a distance of zero to within {tolerance:e}"
        );
    } else {
        common::assert_close(actual, expected_value, tolerance, context);
    }
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = verdict_case(case);
        let context = &format!("{}::{}", spec.id, case.id);

        let expected_tm = case.expected_vector("tm").expect("the case declares tm");
        assert_eq!(result.tm.len(), expected_tm.len(), "{context}: tm length");
        for (i, (&actual, &expected_value)) in result.tm.iter().zip(expected_tm).enumerate() {
            assert_tm(
                actual,
                expected_value,
                case.tolerance,
                &format!("{context} (tm[{i}])"),
            );
        }

        assert_eq!(
            result.iterations.len(),
            case.expected_vector("iterations")
                .expect("the case declares iterations")
                .len(),
            "{context}: iterations length"
        );
        for (i, (&actual, &expected_value)) in result
            .iterations
            .iter()
            .zip(case.expected_vector("iterations").expect("declared"))
            .enumerate()
        {
            assert_eq!(
                f64::from(actual),
                expected_value,
                "{context}: iteration count for trial {i}"
            );
        }

        // The compositions, as far as the generated table can carry them: two rows,
        // each a composition. See this file's header for what pins their values.
        assert_eq!(result.w.len(), 2, "{context}: two trials, so two rows");
        for (i, row) in result.w.iter().enumerate() {
            assert_eq!(
                row.len(),
                case.vector("z").expect("the case declares z").len(),
                "{context}: w[{i}] length"
            );
            let sum: f64 = row.iter().sum();
            assert!(
                (sum - 1.0).abs() <= 1e-12,
                "{context}: w[{i}] sums to {sum}, not to one"
            );
            assert!(
                row.iter().all(|&value| value >= 0.0),
                "{context}: w[{i}] = {row:?} has a negative mole fraction"
            );
        }

        // Both cases are two-trial states, so neither verdict is a foregone
        // conclusion: the vapourization trial finds the feed and says nothing, and
        // it is the *other* row that decides. Asserted because a model that swapped
        // the two would satisfy every value above.
        let z = case.vector("z").expect("the case declares z");
        for (i, (&row, &zi)) in result.w[0].iter().zip(z).enumerate() {
            assert!(
                (row - zi).abs() <= 1e-10,
                "{context}: the first row should be the vapour-like trial, and on these \
                 two states its stationary point is the feed itself; w[0][{i}] is {row} \
                 against z[{i}] = {zi}"
            );
        }
    }
}

/// The scientific content of the pair: this model's verdict against the flash's phase.
///
/// The two models answer different questions about one state, and their answers are
/// related. A feed a flash splits into two phases is unstable as a single phase -
/// that is what the split *means* - and a feed it reports as one phase is not. So
/// over a sweep of states the verdict must be `unstable` exactly where the phase is
/// `two_phase`, and `stable` everywhere else: `all_liquid`, `all_vapour` and
/// `trivial` are all single-phase readings.
///
/// `trivial` is the one worth stating. The flash reached `x = y = z` and cannot say
/// which single phase the feed is; this model *can*, and says it is one phase. The
/// two are consistent rather than redundant - which is why the pair exists.
///
/// Nothing about this is a value comparison: the two answers are enums, computed by
/// separate iterations over separate equations, and the claim is that they never
/// contradict each other.
#[test]
fn the_verdict_agrees_with_the_flash_phase_on_every_state() {
    let mixtures = [
        ("methane/butane", methane_butane()),
        (
            "propane/butane",
            databank::mixture_of(&["propane", "n-butane"], None)
                .expect("the pair resolves")
                .0,
        ),
    ];

    let mut splits = 0;
    let mut checked = 0;
    for (name, mixture) in &mixtures {
        for t in [250.0, 280.0, 300.0, 330.0, 350.0, 400.0] {
            for p in [
                100_000.0,
                1_000_000.0,
                3_000_000.0,
                10_000_000.0,
                20_000_000.0,
                50_000_000.0,
            ] {
                for z in [[0.1, 0.9], [0.5, 0.5], [0.9, 0.1]] {
                    let Ok(flash) = pt_flash(mixture, kelvins(t), pascals(p), &z) else {
                        continue;
                    };
                    let Ok(test) = stability_test(mixture, kelvins(t), pascals(p), &z) else {
                        continue;
                    };
                    checked += 1;

                    let expected = if flash.phase == Phase::TwoPhase {
                        splits += 1;
                        StabilityVerdict::Unstable
                    } else {
                        StabilityVerdict::Stable
                    };
                    assert_eq!(
                        test.verdict,
                        expected,
                        "{name}, T={t}, P={p}, z={z:?}: the flash says `{}` and this model \
                         says `{}`. A state that splits is unstable as a single phase, and \
                         one that does not is not - a disagreement here is a contradiction, \
                         not a tolerance",
                        flash.phase.as_str(),
                        test.verdict.as_str()
                    );
                }
            }
        }
    }

    assert!(
        splits > 20,
        "the sweep should find plenty of two-phase states, found {splits}"
    );
    assert!(
        checked > 150,
        "the sweep should cover plenty of states, covered {checked}"
    );
}

/// A feed the flash calls `trivial` is one this model calls `stable`.
///
/// The contrast that is the model's reason to exist, on the spec's second case. At
/// 430 K and 60 bar the flash converges to `x = y = z` and reports `trivial`, which
/// says the feed is single phase and *not which one* - the K-values straddled one
/// throughout, so nothing in that model ever proved what the feed was. This model
/// answers the question it could not ask.
///
/// Both trials converge to the feed here, so the verdict is the weakest kind of
/// `stable`: two trials that found nothing, not two trials that separated. That is
/// the honest answer to the question asked, and the test asserts the shape of it -
/// both distances at zero to rounding, and `tm[1]` *negative*, which is exactly why
/// the threshold is `-1e-8` and not zero.
#[test]
fn a_trivial_flash_is_not_a_stable_feed() {
    let mixture = methane_butane();
    let (t, p, z) = (430.0, 6_000_000.0, [0.6, 0.4]);

    let flash = pt_flash(&mixture, kelvins(t), pascals(p), &z).expect("the flash computes");
    assert_eq!(
        flash.phase,
        Phase::Trivial,
        "the contrast is only a contrast if the flash cannot answer here"
    );

    let test = stability_test(&mixture, kelvins(t), pascals(p), &z).expect("the test computes");
    assert_eq!(test.verdict, StabilityVerdict::Stable);
    for (i, &distance) in test.tm.iter().enumerate() {
        assert!(
            distance.abs() <= 1e-8,
            "trial {i} reached a stationary point at tm = {distance}, so this is not the \
             two-trivial-trials state the case describes"
        );
    }
    // The *sign* of that round-off is deliberately not asserted. It says which way the
    // last `exp` rounded, and this test used to pin it as negative - a measurement of
    // the old illustrative `kij` rather than a property of the model, and one the
    // databank's fitted parameter rounds the other way. What the threshold has to
    // tolerate is a distance below zero, and that a run ever produces one is a fact
    // about cancellation rather than about 430 K.
}

/// The feed is placed on whichever admissible root has the lower Gibbs energy.
///
/// The decision this model is most likely to get quietly wrong, and the one the
/// spec's `notes` record a measurement for. The comparison is
/// `A^R/RT - ln Z + Z` and **not** `A^R/RT + Z`: the ideal part of `G/RT` carries
/// `-ln Z`, and `V = Z R T / P` differs between the two roots, so only `-ln Z`
/// separates them at one `(T, P, n)`.
///
/// Pure methane at 150 K and 1 bar is superheated vapour - its saturation pressure
/// there is about 10 bar - so the vapour-like root is the state the feed is in. The
/// two comparisons disagree about which root that is, and the numbers below are the
/// spec's measurement of it. The second half of the test is the consequence rather
/// than the rule: with the right comparison the feed is a stable single phase, and
/// with the wrong one both trials measure their distance from a state the feed is
/// not in and the model calls a plain vapour unstable.
#[test]
fn the_feed_is_placed_on_its_lower_gibbs_root() {
    let mixture = databank::mixture_of(&["methane"], None)
        .expect("methane resolves")
        .0;
    let (t, p) = (kelvins(150.0), pascals(100_000.0));
    let z = [1.0];

    let reduced = mixture.reduced_parameters(t, p).unwrap();
    let (a_mix, b_mix) = mixture.mixture_parameters(&reduced, &z);
    let roots = azoth_eos::pr_z_factor(a_mix, b_mix).unwrap();

    let gibbs = |compressibility: f64| -> (f64, f64, f64) {
        let residual = mixture
            .helmholtz_energy(&reduced, &z, compressibility)
            .unwrap();
        (
            residual,
            residual + compressibility,
            residual - compressibility.ln() + compressibility,
        )
    };
    let (liquid_a, liquid_naive, liquid_right) = gibbs(roots.z_min);
    let (vapour_a, vapour_naive, vapour_right) = gibbs(roots.z_max);

    common::assert_close(liquid_a, -2.5578138051194212, 1e-9, "liquid A^R/RT");
    common::assert_close(vapour_a, -0.015548927641895367, 1e-9, "vapour A^R/RT");
    common::assert_close(
        liquid_right,
        3.1458414529963945,
        1e-9,
        "liquid A^R/RT - ln Z + Z",
    );
    common::assert_close(
        vapour_right,
        0.9845725795485407,
        1e-9,
        "vapour A^R/RT - ln Z + Z",
    );

    assert!(
        liquid_naive < vapour_naive,
        "the wrong comparison is only worth recording if it picks the other root: \
         A^R/RT + Z gives {liquid_naive} for the liquid and {vapour_naive} for the vapour"
    );
    assert!(
        vapour_right < liquid_right,
        "A^R/RT - ln Z + Z must pick the vapour root for a superheated vapour: it gives \
         {vapour_right} for the vapour against {liquid_right} for the liquid"
    );

    // The consequence. A superheated vapour is a stable single phase, and the trial
    // seeded on its own side converges to the feed itself - `tm` exactly zero and
    // `w` exactly `z`, which is the shape a trivial stationary point has.
    let result = stability_test(&mixture, t, p, &z).unwrap();
    assert_eq!(
        result.verdict,
        StabilityVerdict::Stable,
        "the wrong root calls this unstable: {:?}",
        result.tm
    );
    assert_eq!(
        result.tm[0], 0.0,
        "the vapour-like trial should land on the feed, so this is the vapour root; \
         with the wrong comparison it measures from the liquid root instead and comes \
         back strongly negative"
    );
    assert_eq!(
        result.w[0], z,
        "the feed's own stationary point is the feed"
    );
    assert_eq!(result.iterations[0], 1, "it starts at the answer");
    assert!(
        result.tm[1] > 0.0,
        "the liquid-like trial is measured against the vapour root and finds nothing \
         below it; got {}",
        result.tm[1]
    );
}

/// A pure component's trials both sit at the feed, and one of them is trivial.
///
/// With one component every seed normalises to `[1.0]`, so both trials evaluate the
/// same composition and differ only in *which root* they claim. The trial claiming
/// the root the feed is actually on has `ln W = ln phi(z) - ln phi(z) = 0` at its
/// first step: `tm = 1 - sum(W) = 0` exactly, `w == z`, and it converges in one
/// step. That is the trivial stationary point stated exactly rather than approached.
///
/// **Which trial that is, is the feed's decision and not the seed's.** A subcooled
/// liquid sits on the liquid root, so it is the *second* trial that is trivial and
/// the first that reports a positive distance; a superheated vapour is the other way
/// round; and a state above the critical temperature has one admissible root, so
/// both trials claim it and both are trivial. Asserting only "the vapour-like trial
/// is trivial" would be asserting something false for half the states a caller
/// meets, which is what the state table below pins.
#[test]
fn a_pure_components_trials_are_both_the_feed_and_one_is_trivial() {
    for name in ["methane", "n-butane", "propane"] {
        let tc = databank::entry(name, None).expect("a databank entry").tc;
        let mixture = databank::mixture_of(&[name], None).expect("resolves").0;
        for t in [200.0, 280.0, 330.0] {
            for p in [100_000.0, 1_000_000.0, 5_000_000.0] {
                let z = [1.0];
                let Ok(result) = stability_test(&mixture, kelvins(t), pascals(p), &z) else {
                    continue;
                };
                let context = &format!("Tc={tc}, T={t}, P={p}");

                for (i, row) in result.w.iter().enumerate() {
                    assert_eq!(row, &z, "{context}: trial {i} did not return the feed");
                }
                let trivial: Vec<usize> = result
                    .tm
                    .iter()
                    .enumerate()
                    .filter(|&(_, &distance)| distance == 0.0)
                    .map(|(i, _)| i)
                    .collect();
                assert!(
                    !trivial.is_empty(),
                    "{context}: the trial claiming the root the feed is on has nothing to \
                     iterate, so at least one distance is exactly zero; got {:?}",
                    result.tm
                );
                for &trial in &trivial {
                    assert_eq!(
                        result.iterations[trial], 1,
                        "{context}: a trivial trial starts at its stationary point"
                    );
                }
                for (i, &distance) in result.tm.iter().enumerate() {
                    assert!(
                        distance >= 0.0,
                        "{context}: trial {i} reached {distance}, below the tangent plane, \
                         so a pure component would be reported unstable"
                    );
                }
                assert_eq!(
                    result.verdict,
                    StabilityVerdict::Stable,
                    "{context}: a pure component has no split of its own to find"
                );
            }
        }
    }

    // And which trial is the trivial one, on two states that differ only in side:
    // butane's saturation pressure at 330 K is about 6 bar, so 1 bar is superheated
    // vapour and 10 bar is subcooled liquid.
    let butane = databank::mixture_of(&["n-butane"], None)
        .expect("n-butane resolves")
        .0;
    let vapour = stability_test(&butane, kelvins(330.0), pascals(100_000.0), &[1.0]).unwrap();
    assert_eq!(
        vapour.tm[0], 0.0,
        "a superheated vapour is on the vapour root"
    );
    assert!(vapour.tm[1] > 0.0);

    let liquid = stability_test(&butane, kelvins(330.0), pascals(1_000_000.0), &[1.0]).unwrap();
    assert_eq!(
        liquid.tm[1], 0.0,
        "a subcooled liquid is on the liquid root"
    );
    assert!(liquid.tm[0] > 0.0);
}

/// A trial that cannot reach its stationary point raises, and is not discarded.
///
/// The rule the spec's port notes record as a decision: a tangent-plane distance
/// bounds stability only *at a stationary point*, so a `tm` read from a partial
/// iteration is a fact about the iteration rather than about the mixture - and
/// discarding it turns "could not tell" into "stable", which is the failure this
/// library is organised against. NeqSim breaks out of its loop and carries on; this
/// raises, exactly as the flash does.
///
/// The state is a feed with an **absent** component. Its reference potential is
/// `-inf` and its trial mole number is zero and stays zero, so the rms change in
/// `ln W` between two iterations is `-inf - -inf` - a `NaN` - at every step. The
/// iteration therefore never meets a tolerance and runs to the cap: 2000 steps, by
/// design, and an error rather than a verdict. Measured on the Python reference at
/// 0.09 s, so this costs nothing to assert.
///
/// The cap is the spec's and deliberately far above every state tested - the
/// measured worst case on the binaries is 507 steps - which is why it takes a
/// degenerate composition to reach it.
#[test]
fn a_trial_that_cannot_converge_raises_rather_than_calling_the_feed_stable() {
    let mixture = methane_butane();
    let err = stability_test(&mixture, kelvins(330.0), pascals(2_500_000.0), &[1.0, 0.0])
        .expect_err("a trial that never converges must not produce a verdict");
    match err {
        AzothError::SolverNotConverged {
            iterations,
            residual,
            tolerance,
        } => {
            assert_eq!(iterations, 2000, "the cap is the spec's");
            assert_eq!(tolerance, 1.0e-10, "the tolerance is the spec's");
            assert!(
                residual.is_nan(),
                "the residual is a difference of two `-inf`s, so it is NaN rather than \
                 large; got {residual}"
            );
        }
        other => panic!("expected a solver failure, got {other:?}"),
    }
}

/// The two trials come back vapour-like first, and the order is the contract.
///
/// `tm`, `w` and `iterations` are positional and a caller reads index 0 as the
/// vapour-like trial. Swapping them would silently relabel the two, and every value
/// in the result would still be a number the model produced - so nothing but a test
/// can catch it.
///
/// The spec's first case is the one that separates them: the vapour-like trial
/// converges to the feed (16 steps, `tm = 0`) and the liquid-like one proves the
/// feed unstable (9 steps, `tm = -0.215`). Different distances *and* different
/// counts, in the order the model documents.
#[test]
fn the_trials_are_reported_vapour_like_first() {
    let mixture = methane_butane();
    let result = stability_test(&mixture, kelvins(330.0), pascals(2_500_000.0), &[0.6, 0.4])
        .expect("computes");

    assert!(
        result.tm[0] > result.tm[1],
        "the vapour-like trial is the one that finds the feed and the liquid-like one is \
         the one that proves it unstable; got {:?}",
        result.tm
    );
    assert_eq!(result.iterations, vec![16, 9]);

    // And the rows are the two Wilson seeds' stationary points, not two guesses: at
    // this state light methane's K-value is well above one and butane's well below,
    // so `z_i K_i` is the methane-rich seed and `z_i / K_i` the butane-rich one. A
    // model that ran the liquid-like trial first would put the butane-rich row at
    // index 0.
    assert!(
        result.w[0][0] > result.w[1][0],
        "row 0 should be the vapour-like trial, whose seed is z_i K_i: {:?}",
        result.w
    );
    assert!(
        result.w[1][1] > result.w[0][1],
        "row 1 should be the liquid-like trial, whose seed is z_i / K_i: {:?}",
        result.w
    );
}

/// `tm` is never NaN, which is the one value a caller cannot reason about.
///
/// A `NaN` distance would satisfy neither `tm < -1e-8` nor its negation, so a feed
/// would be called `stable` by a comparison that never ran. The threshold happens to
/// be safe against that - `NaN < x` is false in Rust, so the verdict falls to
/// `stable` - which is exactly the wrong answer to give silently. This asserts the
/// distances are numbers wherever a result is returned at all.
#[test]
fn no_distance_is_ever_nan() {
    let mixtures = [
        ("methane/butane", methane_butane()),
        (
            "propane/butane",
            databank::mixture_of(&["propane", "n-butane"], None)
                .expect("the pair resolves")
                .0,
        ),
    ];

    let mut checked = 0;
    for (name, mixture) in &mixtures {
        for t in [250.0, 280.0, 300.0, 330.0, 350.0, 400.0] {
            for p in [
                100_000.0,
                1_000_000.0,
                3_000_000.0,
                10_000_000.0,
                50_000_000.0,
            ] {
                for z in [[0.1, 0.9], [0.5, 0.5], [0.9, 0.1]] {
                    let Ok(result) = stability_test(mixture, kelvins(t), pascals(p), &z) else {
                        continue;
                    };
                    checked += 1;
                    for (i, &distance) in result.tm.iter().enumerate() {
                        assert!(
                            distance.is_finite(),
                            "{name}, T={t}, P={p}, z={z:?}: tm[{i}] is {distance}"
                        );
                    }
                    assert!(
                        result.min_t_over_tc.is_finite() && result.min_t_over_tc > 0.0,
                        "{name}, T={t}, P={p}, z={z:?}: min_t_over_tc is {}",
                        result.min_t_over_tc
                    );
                }
            }
        }
    }
    assert!(checked > 100, "the sweep only covered {checked} states");
}

/// The verdict is the threshold applied to the two distances, and nothing else.
///
/// The identity that ties the two reported things together: `unstable` if and only
/// if some `tm` is below `-1e-8`. Asserted on returned results rather than
/// recomputed from the inputs, so a result reporting a verdict its own distances do
/// not support fails.
#[test]
fn the_verdict_is_the_threshold_applied_to_the_distances() {
    let mixture = methane_butane();
    let mut unstable = 0;
    for t in [250.0, 280.0, 300.0, 330.0, 350.0, 400.0] {
        for p in [
            100_000.0,
            1_000_000.0,
            3_000_000.0,
            10_000_000.0,
            50_000_000.0,
        ] {
            for z in [[0.1, 0.9], [0.5, 0.5], [0.9, 0.1]] {
                let Ok(result) = stability_test(&mixture, kelvins(t), pascals(p), &z) else {
                    continue;
                };
                let expected = if result.tm.iter().any(|&distance| distance < -1.0e-08) {
                    unstable += 1;
                    StabilityVerdict::Unstable
                } else {
                    StabilityVerdict::Stable
                };
                assert_eq!(result.verdict, expected, "T={t}, P={p}, z={z:?}");
            }
        }
    }
    assert!(unstable > 0, "the sweep found no unstable feed at all");
}

#[test]
fn a_feed_that_does_not_sum_to_one_is_refused() {
    let mixture = methane_butane();
    for z in [vec![0.6, 0.5], vec![0.6, 0.3], vec![0.6, 0.4, 0.0]] {
        let err = stability_test(&mixture, kelvins(330.0), pascals(2_500_000.0), &z).unwrap_err();
        assert!(
            matches!(err, AzothError::InvalidInput { .. }),
            "for z = {z:?}"
        );
    }
}

#[test]
fn a_negative_mole_fraction_is_refused() {
    let mixture = methane_butane();
    let err = stability_test(&mixture, kelvins(330.0), pascals(2_500_000.0), &[1.6, -0.6])
        .expect_err("a negative mole fraction is not a feed");
    assert!(matches!(err, AzothError::InvalidInput { .. }));
}

#[test]
fn a_non_positive_temperature_or_pressure_is_refused() {
    let mixture = methane_butane();
    for (t, p) in [(0.0, 2_500_000.0), (-1.0, 2_500_000.0), (330.0, 0.0)] {
        assert!(
            stability_test(&mixture, kelvins(t), pascals(p), &[0.6, 0.4]).is_err(),
            "T={t}, P={p}"
        );
    }
}

/// The near-critical bound warns without refusing, as it does on the flash.
///
/// Above `T / Tc = 0.9` for the nearest component the cubic's roots are close to
/// coalescing, and this model rests entirely on the trials finding *distinct*
/// stationary points - so a `stable` verdict here is the one to distrust. The
/// arithmetic is still defined and the answer is still the criterion's, which is why
/// it is a warning and not a refusal.
#[test]
fn a_near_critical_component_warns_rather_than_failing() {
    let mixture = methane_butane();
    // Butane's Tc is 425.12 K, so 400 K puts it at 0.941.
    let result = stability_test(&mixture, kelvins(400.0), pascals(2_000_000.0), &[0.6, 0.4])
        .expect("a near-critical state still computes");
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code == WarningCode::OutOfValidRange),
        "expected an out-of-valid-range warning, got {:?}",
        result.warnings
    );
    assert!(
        (result.min_t_over_tc - 400.0 / 425.12).abs() < 1e-12,
        "the reported bound subject should be the smallest T/Tc"
    );
}

#[test]
fn the_result_is_clean_at_an_ordinary_state() {
    let result = stability_test(
        &methane_butane(),
        kelvins(330.0),
        pascals(2_500_000.0),
        &[0.6, 0.4],
    )
    .unwrap();
    assert!(
        result.is_clean(),
        "unexpected warnings: {:?}",
        result.warnings
    );
}

/// A malformed `kij` is refused at construction, as it is for every other model.
///
/// Asserted here too rather than left to the flash's test because the refusal is
/// what makes the *model* layer's shared arithmetic safe: this model reaches the
/// mixing rule through the same `Mixture`, and a matrix that got past construction
/// would give every trial a different mixture than the caller described.
#[test]
fn a_malformed_kij_matrix_is_refused() {
    let components = || {
        vec![
            Component::new(kelvins(190.56), pascals(4_599_200.0), 0.01142).unwrap(),
            Component::new(kelvins(425.12), pascals(3_796_000.0), 0.2002).unwrap(),
        ]
    };
    assert!(Mixture::new(components(), vec![0.0; 3]).is_err());
    assert!(Mixture::new(components(), vec![0.1, 0.05, 0.05, 0.0]).is_err());
    assert!(Mixture::new(components(), vec![0.0, 0.05, 0.07, 0.0]).is_err());
    assert!(Mixture::new(components(), vec![0.0, 0.05, 0.05, 0.0]).is_ok());
}

/// Each reported trial is a stationary point of the iteration, and `tm` is its value.
///
/// The strongest data-free check on this model, and the one that ties it to the
/// cubic rather than to its own arithmetic. Everything below is recomputed here from
/// the public API - the feed's lower-Gibbs root included, by the rule the module
/// documents - and then the *returned* `w` is fed back through the model's own
/// stationarity equation:
///
/// ```text
/// W_i  = exp(d_i - ln phi_i(w))        must reproduce w after normalisation
/// tm   = 1 - sum_i W_i                 must reproduce the returned distance
/// ```
///
/// A model that reported a partially-converged iterate, a composition from the wrong
/// root, or a distance measured somewhere other than at the `w` it returned fails
/// this without any expected value to compare against.
#[test]
fn every_reported_trial_is_a_stationary_point() {
    let spec = model_gen::model(MODEL_ID).expect("the model");
    for case in spec.cases {
        let mixture = mixture_from_case(case);
        let (t, p) = (
            kelvins(common::input(case, "T")),
            pascals(common::input(case, "P")),
        );
        let z = case.vector("z").expect("the case declares z");
        let result = verdict_case(case);
        let context = &format!("{}::{}", spec.id, case.id);

        let reduced = mixture.reduced_parameters(t, p).unwrap();
        let (a_mix, b_mix) = mixture.mixture_parameters(&reduced, z);
        let roots = azoth_eos::pr_z_factor(a_mix, b_mix).unwrap();
        // The feed's root: the lower Gibbs energy of the admissible ones, which is
        // the comparison the model documents - `A^R/RT - ln Z + Z`.
        let gibbs = |compressibility: f64| {
            mixture
                .helmholtz_energy(&reduced, z, compressibility)
                .unwrap()
                - compressibility.ln()
                + compressibility
        };
        let root = if roots.z_min == roots.z_max || gibbs(roots.z_min) < gibbs(roots.z_max) {
            roots.z_min
        } else {
            roots.z_max
        };
        let feed = mixture.phase_state_at(&reduced, z, root).unwrap();
        let d: Vec<f64> = z
            .iter()
            .zip(&feed.ln_phi)
            .map(|(&zi, &lp)| zi.ln() + lp)
            .collect();

        for (trial, side) in [(0, RootSide::Vapour), (1, RootSide::Liquid)] {
            let w = &result.w[trial];
            let state = mixture.phase_state(&reduced, w, side).unwrap();

            let mole_numbers: Vec<f64> = d
                .iter()
                .zip(&state.ln_phi)
                .map(|(&di, &lp)| (di - lp).exp())
                .collect();
            let totals: f64 = mole_numbers.iter().sum();

            for (i, (back, &reported)) in mole_numbers
                .iter()
                .map(|&value| value / totals)
                .zip(w)
                .enumerate()
            {
                assert!(
                    (back - reported).abs() <= 1e-09,
                    "{context}: w[{trial}][{i}] is {reported}, but the iteration at that \
                     composition gives {back} - it is not a stationary point"
                );
            }
            // Absolute rather than relative: at a trivial trial both sides of this
            // are of order `1e-16`, and a relative comparison between two numbers
            // that size is a comparison of round-off. Every meaningful distance here
            // is `O(1)`, so `1e-9` separates them without flattering either.
            assert!(
                (1.0 - totals - result.tm[trial]).abs() <= 1.0e-09,
                "{context}: tm[{trial}] is reported as {}, but the iteration at the \
                 returned w gives {}",
                result.tm[trial],
                1.0 - totals
            );
        }
    }
}
