//! Spec-driven tests for the `eos.pt_flash` model.
//!
//! The spec's cases pin the two implementations to each other. These tests are for
//! the things a case *cannot* say: the identities that hold at any answer, the
//! reductions that tie the model layer to the registered kernels below it, and the
//! failure modes that were found by measurement rather than by construction.

use azoth_core::units::{kelvins, pascals};
use azoth_core::{AzothError, CalcResult, WarningCode};
use azoth_eos::mixture::{Component, Mixture, RootSide};
use azoth_eos::{Phase, databank, model_gen, pt_flash};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.pt_flash";

fn mixture_of(tc: &[f64], pc: &[f64], omega: &[f64], kij: Vec<f64>) -> Mixture {
    let components = (0..tc.len())
        .map(|i| {
            Component::new(kelvins(tc[i]), pascals(pc[i]), omega[i]).expect("a valid component")
        })
        .collect();
    Mixture::new(components, kij).expect("a valid mixture")
}

/// The methane/n-butane pair the spec's first and third cases use.
///
/// Resolved through the databank, so the pair the sweeps run and the pair the cases
/// run are the same fluid. The critical constants used to be typed out here and had
/// drifted from `COMP.csv` - methane at 0.01142 and 4 599 200 Pa, where NeqSim says
/// 0.0115 and 4 599 000 - and the `kij` was an "illustrative" 0.05 rather than the
/// 0.01289789 `INTER.csv` fits.
fn methane_butane() -> Mixture {
    databank::mixture_of(&["methane", "n-butane"])
        .expect("the pair resolves")
        .0
}

fn mixture_from_case(case: &azoth_core::spec::TestCase) -> Mixture {
    let names = case
        .list("components")
        .expect("the case declares components");
    databank::mixture_of(names)
        .expect("the case's components resolve")
        .0
}

fn flash_case(case: &azoth_core::spec::TestCase) -> azoth_eos::PtFlashResult {
    let mixture = mixture_from_case(case);
    pt_flash(
        &mixture,
        kelvins(common::input(case, "T")),
        pascals(common::input(case, "P")),
        case.vector("z").expect("the case declares z"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = flash_case(case);

        for name in ["z_liquid", "z_vapour", "min_t_over_tc"] {
            common::assert_close(
                match name {
                    "z_liquid" => result.z_liquid,
                    "z_vapour" => result.z_vapour,
                    _ => result.min_t_over_tc,
                },
                common::expected(case, name),
                case.tolerance,
                &format!("{}::{} ({name})", spec.id, case.id),
            );
        }
        // `residual` is deliberately *not* pinned as a case value: it is an rms of
        // `ln`-built quantities, so the two implementations differ in its last few
        // ulps and a value near 4.7e-11 has no significant figures to spare. What
        // is reproducible - and what the case does assert - is the iteration count.
        // See the spec's correction 4.
        assert!(
            result.residual <= spec.algorithm.expect("a procedure").tolerance,
            "{}::{}: residual {:e} exceeds the declared tolerance {:e}",
            spec.id,
            case.id,
            result.residual,
            spec.algorithm.expect("a procedure").tolerance
        );
        assert_eq!(
            result.iterations as f64,
            common::expected(case, "iterations"),
            "{}::{}: iteration count",
            spec.id,
            case.id
        );

        // The vapour fraction is optional, and every spec case is a state where it
        // exists - a trivial solution is a *behaviour* and is tested as one below,
        // not pinned as a value.
        let beta = result
            .beta
            .unwrap_or_else(|| panic!("{}::{}: expected a vapour fraction", spec.id, case.id));
        common::assert_close(
            beta,
            common::expected(case, "beta"),
            case.tolerance,
            &format!("{}::{} (beta)", spec.id, case.id),
        );

        for name in ["x", "y", "k", "ln_phi_liquid", "ln_phi_vapour"] {
            let actual: &[f64] = match name {
                "x" => &result.x,
                "y" => &result.y,
                "k" => &result.k,
                "ln_phi_liquid" => &result.ln_phi_liquid,
                _ => &result.ln_phi_vapour,
            };
            let expected = case.expected_vector(name).expect("declared");
            assert_eq!(actual.len(), expected.len(), "{name}: length");
            for (i, (a, e)) in actual.iter().zip(expected).enumerate() {
                common::assert_close(
                    *a,
                    *e,
                    case.tolerance,
                    &format!("{}::{} ({name}[{i}])", spec.id, case.id),
                );
            }
        }
    }
}

/// The five identities that hold at *any* correct answer.
///
/// Stronger than any single expected value, because a plausible-looking wrong answer
/// cannot satisfy them: material balance, both normalisations, the K-definition, and
/// the Rachford-Rice residual. Recomputed here from the result rather than read from
/// it, so a result reporting an identity it had not achieved fails.
#[test]
fn the_identities_hold_at_every_returned_answer() {
    let spec = model_gen::model(MODEL_ID).expect("the model");
    for case in spec.cases {
        let result = flash_case(case);
        let z = case.vector("z").expect("the case declares z");
        let beta = result.beta.expect("the spec's cases are all splits");
        let context = &format!("{}::{}", spec.id, case.id);

        for (i, &zi) in z.iter().enumerate() {
            let material = (1.0 - beta) * result.x[i] + beta * result.y[i];
            assert!(
                (material - zi).abs() <= 1e-12,
                "{context}: material balance at component {i} is {material}, not {zi}"
            );
            assert!(
                (result.k[i] - result.y[i] / result.x[i]).abs() <= 1e-12,
                "{context}: K[{i}] is {} but y/x is {}",
                result.k[i],
                result.y[i] / result.x[i]
            );
        }

        let sum_x: f64 = result.x.iter().sum();
        let sum_y: f64 = result.y.iter().sum();
        assert!(
            (sum_x - 1.0).abs() <= 1e-12 && (sum_y - 1.0).abs() <= 1e-12,
            "{context}: compositions sum to {sum_x} and {sum_y}, not to one"
        );

        let residual: f64 = z
            .iter()
            .zip(&result.k)
            .map(|(&zi, &ki)| zi * (ki - 1.0) / (1.0 + beta * (ki - 1.0)))
            .sum();
        assert!(
            residual.abs() <= 1e-12,
            "{context}: the Rachford-Rice residual is {residual:e}"
        );
    }
}

/// The answer is a minimum of the total Gibbs energy, not merely a fixed point.
///
/// The strongest available correctness test for a flash and the one that needs no
/// external data: perturbing `beta` at fixed `K` must *raise* `G / RT`, on both
/// sides, and raise it quadratically - which a converged-but-wrong answer cannot do.
/// This is the embryo of the stability test the model does not have.
#[test]
fn the_answer_is_a_gibbs_minimum() {
    let mixture = methane_butane();
    let z = [0.1, 0.9];
    let (t, p) = (kelvins(300.0), pascals(3_000_000.0));
    let result = pt_flash(&mixture, t, p, &z).unwrap();
    let beta = result.beta.expect("a split");

    let reduced = mixture.reduced_parameters(t, p).unwrap();
    let total = |beta: f64| -> f64 {
        let x: Vec<f64> = z
            .iter()
            .zip(&result.k)
            .map(|(&zi, &ki)| zi / (1.0 + beta * (ki - 1.0)))
            .collect();
        let y: Vec<f64> = result.k.iter().zip(&x).map(|(&ki, &xi)| ki * xi).collect();
        let liquid = mixture.phase_state(&reduced, &x, RootSide::Liquid).unwrap();
        let vapour = mixture.phase_state(&reduced, &y, RootSide::Vapour).unwrap();
        let g_liquid: f64 = x
            .iter()
            .zip(&liquid.ln_phi)
            .map(|(&xi, &lp)| xi * (xi.ln() + lp))
            .sum();
        let g_vapour: f64 = y
            .iter()
            .zip(&vapour.ln_phi)
            .map(|(&yi, &lp)| yi * (yi.ln() + lp))
            .sum();
        (1.0 - beta) * g_liquid + beta * g_vapour
    };

    let base = total(beta);
    for step in [1e-3, 1e-4, 1e-5] {
        for offset in [-step, step] {
            let raised = total(beta + offset) - base;
            assert!(
                raised > 0.0,
                "perturbing beta by {offset:e} lowered G/RT by {raised:e} - the answer \
                 is not a minimum"
            );
        }
    }
    // Quadratic, not merely positive: shrinking the step by 10 shrinks the rise by
    // 100. A converged-but-wrong answer, or a wrong root, fails this even when the
    // rise happens to have the right sign.
    let coarse = total(beta + 1e-3) - base;
    let fine = total(beta + 1e-4) - base;
    assert!(
        (coarse / fine - 100.0).abs() < 1.0,
        "the rise is not quadratic: {coarse:e} at 1e-3 against {fine:e} at 1e-4"
    );
}

/// The mixture form reduces to the registered pure-component kernel at `N = 1`.
///
/// The model layer's fugacity coefficient is the one piece of arithmetic here that no
/// registered calculation covers, and this is the check that keeps it honest:
/// `eos.pr_departure` is the same expression with the cross-sum factor collapsed to
/// one, so an implementation that got the factor wrong is caught by a kernel with its
/// own spec, its own worked example and its own source.
#[test]
fn the_mixture_form_reduces_to_pr_departure_at_one_component() {
    for (tc, pc, omega, tr, pr) in [
        (369.83, 4_248_000.0, 0.1523, 0.8, 0.25),
        (425.12, 3_796_000.0, 0.2002, 0.7, 0.4),
        (190.56, 4_599_200.0, 0.01142, 1.2, 0.9),
    ] {
        let mixture = mixture_of(&[tc], &[pc], &[omega], vec![0.0]);
        let t = kelvins(tc * tr);
        let p = pascals(pc * pr);
        let reduced = mixture.reduced_parameters(t, p).unwrap();
        let state = mixture
            .phase_state(&reduced, &[1.0], RootSide::Vapour)
            .unwrap();

        let kappa = azoth_eos::pr_kappa(omega).unwrap().kappa;
        let pure =
            azoth_eos::pr_departure(reduced.a[0], reduced.b[0], state.z, kappa, t.value / tc)
                .unwrap();

        assert!(
            (state.ln_phi[0] - pure.ln_phi).abs() < 1e-12,
            "Tc={tc}: the mixture form gives ln phi = {} and pr_departure gives {}",
            state.ln_phi[0],
            pure.ln_phi
        );
        assert!(
            (state.a_mix - reduced.a[0]).abs() < 1e-15
                && (state.b_mix - reduced.b[0]).abs() < 1e-15,
            "Tc={tc}: a pure component's mixture parameters should be its own"
        );
    }
}

/// The mixture parameters reduce to the registered binary kernel at `N = 2`.
///
/// A tolerance rather than bit-equality, and the reason is written where it bites:
/// `eos.vdw1f_mix_binary` evaluates its three terms longhand and this sums `i` then
/// `j`, so the two associate differently. A bit-equality claim here would be a claim
/// about summation order rather than about the mixing rule.
#[test]
fn the_mixture_parameters_reduce_to_the_binary_kernel() {
    let mixture = methane_butane();
    for (t, p) in [
        (330.0, 2_500_000.0),
        (300.0, 3_000_000.0),
        (350.0, 5_000_000.0),
    ] {
        let reduced = mixture.reduced_parameters(kelvins(t), pascals(p)).unwrap();
        for x1 in [0.1, 0.3, 0.6, 0.9] {
            let x = [x1, 1.0 - x1];
            let (a_mix, b_mix) = mixture.mixture_parameters(&reduced, &x);
            // The kernel takes the pair's interaction parameter as its last argument.
            // Read off the mixture rather than written here, so the two terms of the
            // comparison cannot describe different fluids - which is exactly what
            // happened while the fixture was literal and this argument was not.
            let kernel = azoth_eos::vdw1f_mix_binary(
                x1,
                reduced.a[0],
                reduced.a[1],
                reduced.b[0],
                reduced.b[1],
                mixture.kij(0, 1),
            )
            .unwrap();

            assert!(
                (a_mix - kernel.a_mix).abs() <= 1e-12 * a_mix.abs(),
                "T={t}, P={p}, x1={x1}: a_mix {a_mix} against the kernel's {}",
                kernel.a_mix
            );
            assert!(
                (b_mix - kernel.b_mix).abs() <= 1e-12 * b_mix.abs(),
                "T={t}, P={p}, x1={x1}: b_mix {b_mix} against the kernel's {}",
                kernel.b_mix
            );
        }
    }
}

/// The flash's vapour fraction reduces to the registered binary kernel at `N = 2`.
///
/// A *cross-layer* check no single-language test can replace: the closed form is
/// algebra that follows from the Rachford-Rice equation, and this bisection is a
/// general procedure that happens to agree with it. A bracket error in the model
/// layer would still pass every case whose K-values straddle one, and this catches it
/// wherever the two implementations disagree about which root they found.
#[test]
fn the_vapour_fraction_reduces_to_the_binary_kernel() {
    let mixture = methane_butane();
    for (t, p, z) in [
        (300.0, 3_000_000.0, [0.1, 0.9]),
        (330.0, 2_500_000.0, [0.6, 0.4]),
        (330.0, 3_000_000.0, [0.5, 0.5]),
    ] {
        let result = pt_flash(&mixture, kelvins(t), pascals(p), &z).unwrap();
        let kernel = azoth_eos::rachford_rice_binary(z[0], result.k[0], result.k[1]).unwrap();
        assert!(
            (result.beta.unwrap() - kernel.beta).abs() <= 1e-12,
            "T={t}, P={p}, z={z:?}: the bisection gives {} and the closed form {}",
            result.beta.unwrap(),
            kernel.beta
        );
    }
}

/// A feed whose K-values all sit on one side of one is single phase, by proof.
///
/// When every `K_i > 1` the Rachford-Rice equation has no root at all: `sum_i y_i =
/// sum_i K_i x_i = 1` with `sum_i x_i = 1` requires `1 > 1`. When every `K_i < 1` it
/// requires `1 < 1`. There is nothing to bisect toward and nothing to report a
/// vapour fraction of, which is why `beta` is absent and why the loop settles on the
/// **first** iteration - a diagnosis, not an iteration that failed to converge.
///
/// The compositions are the feed because a single phase *is* the feed. What the
/// model does not do is name the phase it is looking at beyond what the K-values
/// say; here they say it unambiguously.
#[test]
fn a_feed_with_no_rachford_rice_root_is_single_phase() {
    let mixture = methane_butane();
    for (t, p, expected) in [
        // Above the dew line: every K above one.
        (280.0, 100_000.0, Phase::AllVapour),
        (330.0, 100_000.0, Phase::AllVapour),
        // Well inside the liquid region: every K below one.
        (280.0, 50_000_000.0, Phase::AllLiquid),
        (330.0, 50_000_000.0, Phase::AllLiquid),
    ] {
        let z = [0.6, 0.4];
        let result = pt_flash(&mixture, kelvins(t), pascals(p), &z).unwrap();

        assert_eq!(result.phase, expected, "T={t}, P={p}");
        assert!(
            result.beta.is_none(),
            "T={t}, P={p}: no root means no vapour fraction, but {} was reported",
            result.beta.unwrap()
        );
        assert_eq!(
            result.iterations, 1,
            "T={t}, P={p}: a no-root diagnosis settles on the first step"
        );
        assert!(
            result.residual.is_nan(),
            "T={t}, P={p}: no step completed, so there is no residual to report"
        );
        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.code == WarningCode::TrivialSolution),
            "T={t}, P={p}: the caller is not told the feed is single phase"
        );
        assert_eq!(result.x, z, "T={t}, P={p}: a single phase is the feed");
        assert_eq!(result.y, z, "T={t}, P={p}: a single phase is the feed");
        assert!(result.k.iter().all(|&v| v != 1.0), "T={t}, P={p}");
    }
}

/// A converged-but-out-of-range vapour fraction is the negative flash, and is
/// reported rather than suppressed.
///
/// `beta > 1` says the feed is superheated vapour and this is how much vapour would
/// have to be condensed out of it. It is a real reading of the same equation, not a
/// failure of it - so unlike the no-root case there *is* a number, and unlike the
/// two-phase case it carries a warning.
#[test]
fn a_negative_flash_reports_its_vapour_fraction() {
    let mixture = methane_butane();
    let (t, p, z) = (330.0, 1_000_000.0, [0.6, 0.4]);
    let result = pt_flash(&mixture, kelvins(t), pascals(p), &z).unwrap();

    assert_eq!(result.phase, Phase::AllVapour);
    let beta = result
        .beta
        .expect("a negative flash still has a vapour fraction");
    assert!(beta > 1.0, "expected beta above one, got {beta}");
    assert_eq!(
        result.iterations, 9,
        "the iteration ran, so this is not a proof"
    );
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code == WarningCode::OutOfValidRange),
        "an out-of-range vapour fraction should say so: {:?}",
        result.warnings
    );
    // The out-of-range reading is still a valid solution of the equations it names,
    // so the material balance holds on it exactly as it does on a split.
    for (i, &zi) in z.iter().enumerate() {
        let material = (1.0 - beta) * result.x[i] + beta * result.y[i];
        assert!(
            (material - zi).abs() <= 1e-12,
            "the negative flash does not satisfy the material balance at {i}"
        );
        assert!(
            result.x[i] > 0.0 && result.y[i] > 0.0,
            "composition {i} is negative"
        );
    }
}

/// A feed the iteration cannot split converges to the trivial solution, and says so.
///
/// This is the model's honest gap made testable. Where the K-values straddle one
/// throughout, the no-root proof above never applies, and for a single-phase feed
/// successive substitution converges to `x = y = z`. There, every K-value is 1,
/// `g(beta)` is identically zero, and the vapour fraction is **indeterminate** - not
/// outside `[0, 1]`, but undefined. Measured at `-7.7e10` for one feed and `-2.2e11`
/// for another, neither reproducible across implementations.
///
/// So the contract is that no number is reported at all, and that the caller is told
/// why. A model returning a large negative `beta` here would be returning something
/// readable as "all liquid" from a state that is a compressed liquid only by
/// accident of the pressure - and the same code path, on a lower-pressure feed,
/// produced a *negative* number for a feed that is unambiguously vapour.
///
/// The result states the trivial solution exactly - `x = y = z`, `K = 1` - rather
/// than the last iterate that approached it, which is a function of where the
/// bisection stopped on an identically-zero function.
#[test]
fn a_feed_that_converges_to_the_trivial_solution_says_so() {
    let mixture = methane_butane();
    for t in [280.0, 300.0, 330.0] {
        let z = [0.6, 0.4];
        let result = pt_flash(&mixture, kelvins(t), pascals(20_000_000.0), &z).unwrap();

        assert_eq!(result.phase, Phase::Trivial, "T={t}");
        assert!(
            result.beta.is_none(),
            "T={t}: a trivial solution has no vapour fraction, but {} was reported",
            result.beta.unwrap()
        );
        assert!(
            result.iterations > 1,
            "T={t}: a trivial solution is the result of an iteration, not of a proof"
        );
        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.code == WarningCode::TrivialSolution),
            "T={t}: the caller is not told which failure this is"
        );
        for (i, &zi) in z.iter().enumerate() {
            assert_eq!(result.x[i], zi, "T={t}: x should be the feed");
            assert_eq!(result.y[i], zi, "T={t}: y should be the feed");
            assert_eq!(result.k[i], 1.0, "T={t}: K should be one");
        }
    }
}

/// Every `Phase` value is reachable, which is what stops one being dead.
///
/// `eos.pr_z_factor`'s `root_structure` had a `two_roots` variant that could not
/// occur - a value a caller branches on and never sees - and it was removed rather
/// than kept for symmetry. This is the same check for this model's enum, and it is
/// why `all_liquid` is in the enum at all: it was found by sweeping 49 states, not
/// assumed to exist.
#[test]
fn every_phase_value_is_reachable() {
    let mixtures = [
        ("methane/butane", methane_butane()),
        (
            "propane/butane",
            mixture_of(
                &[369.83, 425.12],
                &[4_248_000.0, 3_796_000.0],
                &[0.1523, 0.2002],
                vec![0.0, 0.02, 0.02, 0.0],
            ),
        ),
    ];

    let mut seen = std::collections::HashSet::new();
    for (_, mixture) in &mixtures {
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
                    if let Ok(result) = pt_flash(mixture, kelvins(t), pascals(p), &z) {
                        seen.insert(result.phase);
                    }
                }
            }
        }
    }

    for phase in [
        Phase::TwoPhase,
        Phase::AllLiquid,
        Phase::AllVapour,
        Phase::Trivial,
    ] {
        assert!(
            seen.contains(&phase),
            "{phase:?} is in the enum but no state in the sweep reaches it, so a \
             caller branches on a value that never arrives"
        );
    }
}

/// The Rachford-Rice bracket, corrected - the regression that is not a spec case.
///
/// The textbook bracket `1/(1 - K_max) < beta < 1/(1 - K_min)` is correct only when
/// the K-values straddle one. When every one of them is on the same side the interval
/// is empty and inverted, and bisection on it converges to whatever the arithmetic
/// lands on. Measured with the wrong bracket: a negative molar composition for a
/// subcooled liquid feed (the cubic refused `b_reduced = -4.78`) and `beta = -0.129`
/// from a bracket of `(-0.0022, 0)`.
///
/// The property that catches it needs no expected values: on the interval the
/// bisection searches, `1 + beta (K_i - 1)` is positive for **every** `i`. That is
/// the same fact as the bracket's correctness - it is what keeps both compositions
/// positive - so asserting it on the returned state is asserting the bracket.
#[test]
fn both_compositions_stay_positive_at_every_returned_answer() {
    let mixtures = [
        ("methane/butane", methane_butane()),
        (
            "propane/butane",
            mixture_of(
                &[369.83, 425.12],
                &[4_248_000.0, 3_796_000.0],
                &[0.1523, 0.2002],
                vec![0.0, 0.02, 0.02, 0.0],
            ),
        ),
    ];

    let mut splits = 0;
    for (name, mixture) in &mixtures {
        for t in [250.0, 280.0, 300.0, 330.0, 350.0] {
            for p in [100_000.0, 1_000_000.0, 3_000_000.0, 10_000_000.0] {
                for z in [[0.1, 0.9], [0.5, 0.5], [0.9, 0.1]] {
                    let Ok(result) = pt_flash(mixture, kelvins(t), pascals(p), &z) else {
                        continue;
                    };
                    assert!(
                        result.x.iter().all(|&v| v > 0.0),
                        "{name}, T={t}, P={p}, z={z:?}: a negative liquid mole fraction, \
                         which means the Rachford-Rice bracket was wrong"
                    );
                    assert!(
                        result.y.iter().all(|&v| v > 0.0),
                        "{name}, T={t}, P={p}, z={z:?}: a negative vapour mole fraction"
                    );
                    assert!(
                        result.x.iter().all(|&v| v <= 1.0) && result.y.iter().all(|&v| v <= 1.0),
                        "{name}, T={t}, P={p}, z={z:?}: a mole fraction above one"
                    );
                    if result.phase == Phase::TwoPhase {
                        splits += 1;
                    }
                }
            }
        }
    }
    assert!(
        splits > 20,
        "the sweep should find plenty of two-phase states, found {splits}"
    );
}

#[test]
fn a_feed_that_does_not_sum_to_one_is_refused() {
    let mixture = methane_butane();
    for z in [vec![0.6, 0.5], vec![0.6, 0.3], vec![0.6, 0.4, 0.0]] {
        let err = pt_flash(&mixture, kelvins(330.0), pascals(2_500_000.0), &z).unwrap_err();
        assert!(
            matches!(err, AzothError::InvalidInput { .. }),
            "for z = {z:?}"
        );
    }
}

#[test]
fn a_negative_mole_fraction_is_refused() {
    let mixture = methane_butane();
    let err = pt_flash(&mixture, kelvins(330.0), pascals(2_500_000.0), &[1.6, -0.6]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }));
}

/// A `kij` matrix that is not `N x N`, not symmetric, or not zero-diagonal is
/// refused at construction rather than at evaluation.
///
/// All three are mistakes a caller cannot see in the answer: a non-zero diagonal
/// silently rescales that component's attraction, and an asymmetric matrix has no
/// meaning in a double sum that is symmetric by construction.
#[test]
fn a_malformed_kij_matrix_is_refused() {
    let components = || {
        vec![
            Component::new(kelvins(190.56), pascals(4_599_200.0), 0.01142).unwrap(),
            Component::new(kelvins(425.12), pascals(3_796_000.0), 0.2002).unwrap(),
        ]
    };

    assert!(
        Mixture::new(components(), vec![0.0; 3]).is_err(),
        "too short"
    );
    assert!(
        Mixture::new(components(), vec![0.1, 0.05, 0.05, 0.0]).is_err(),
        "non-zero diagonal"
    );
    assert!(
        Mixture::new(components(), vec![0.0, 0.05, 0.07, 0.0]).is_err(),
        "asymmetric"
    );
    assert!(Mixture::new(components(), vec![0.0, 0.05, 0.05, 0.0]).is_ok());
    assert!(Mixture::new(vec![], vec![]).is_err(), "no components");
}

#[test]
fn a_non_positive_temperature_or_pressure_is_refused() {
    let mixture = methane_butane();
    for (t, p) in [(0.0, 2_500_000.0), (-1.0, 2_500_000.0), (330.0, 0.0)] {
        assert!(
            pt_flash(&mixture, kelvins(t), pascals(p), &[0.6, 0.4]).is_err(),
            "T={t}, P={p}"
        );
    }
}

#[test]
fn a_component_with_a_non_positive_critical_point_is_refused() {
    for (tc, pc) in [(0.0, 4_599_200.0), (190.56, 0.0), (-1.0, 4_599_200.0)] {
        assert!(
            Component::new(kelvins(tc), pascals(pc), 0.01142).is_err(),
            "Tc={tc}, Pc={pc}"
        );
    }
}

/// The near-critical bound warns without refusing.
///
/// Above `T / Tc = 0.9` for the nearest component the roots are close to coalescing
/// and successive substitution converges slowly; the answer is still the equation's,
/// so it is a warning and not an error.
#[test]
fn a_near_critical_component_warns_rather_than_failing() {
    let mixture = methane_butane();
    // Butane's Tc is 425.12 K, so 400 K puts it at 0.941.
    let result = pt_flash(&mixture, kelvins(400.0), pascals(2_000_000.0), &[0.6, 0.4])
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
    let result = pt_flash(
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

/// The flash reports the same roots the cubic gives at its own compositions.
///
/// Ties the model layer to `eos.pr_z_factor` rather than to its own arithmetic: a
/// `z_liquid` that was not actually the smallest admissible root would make every
/// `ln phi` downstream wrong while `beta` stayed plausible.
#[test]
fn the_reported_roots_are_the_cubics_roots_at_the_reported_compositions() {
    for case in model_gen::model(MODEL_ID).unwrap().cases {
        let mixture = mixture_from_case(case);
        let result = flash_case(case);
        let reduced = mixture
            .reduced_parameters(
                kelvins(common::input(case, "T")),
                pascals(common::input(case, "P")),
            )
            .unwrap();

        for (x, reported_z, reported_ln_phi, side) in [
            (
                &result.x,
                result.z_liquid,
                &result.ln_phi_liquid,
                RootSide::Liquid,
            ),
            (
                &result.y,
                result.z_vapour,
                &result.ln_phi_vapour,
                RootSide::Vapour,
            ),
        ] {
            let state = mixture.phase_state(&reduced, x, side).unwrap();
            assert!(
                (state.z - reported_z).abs() < 1e-15,
                "{}: the reported root {reported_z} is not the one the cubic gives",
                case.id
            );
            for (i, (&from_state, &from_result)) in
                state.ln_phi.iter().zip(reported_ln_phi).enumerate()
            {
                assert!(
                    (from_state - from_result).abs() < 1e-15,
                    "{}: ln_phi[{i}] is {from_result} but the mixture gives {from_state}",
                    case.id
                );
            }
        }
    }
}
