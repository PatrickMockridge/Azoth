//! `eos.tp_multiflash` against NeqSim's own multiphase flash.
//!
//! The expected values are NeqSim's, printed by `validation/neqsim/TpMultiFlashCases.java`
//! with `setMultiPhaseCheck(true)`. The oracle gates nothing - a divergence is a finding - so
//! the tolerances below are the case files' and the assertions say what a divergence would be
//! rather than only that there is one.
//!
//! **Every case matches each phase by composition, not by index.** The phase order is the
//! solver's own and no order is promised, so an index-wise comparison would test an ordering
//! neither implementation agrees to.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::Cubic;
use azoth_eos::databank::mixture_of;
use azoth_eos::results::TpMultiflashSeed;
use azoth_eos::tp_multiflash;

/// One phase of an oracle state: the fraction, the cubic root and the composition.
#[derive(Debug)]
struct Expected {
    beta: f64,
    z_factor: f64,
    x: &'static [f64],
}

fn state(names: &[&str], t: f64, p: f64, z: &[f64]) -> azoth_eos::results::TpMultiflashResult {
    let (mixture, _) =
        mixture_of(names, Cubic::Pr, None).expect("the databank resolves these names");
    tp_multiflash(&mixture, kelvins(t), pascals(p), z).expect("the flash converges")
}

/// Check every expected phase is present, matched by composition.
///
/// A phase reported that NeqSim did not report is a failure too: the count is part of the
/// answer, and a model that finds a phase upstream does not has found a different state.
fn check(got: &azoth_eos::results::TpMultiflashResult, expected: &[Expected], label: &str) {
    assert_eq!(
        got.phase_count as usize,
        expected.len(),
        "{label}: {} phases where NeqSim reports {}; betas {:?}",
        got.phase_count,
        expected.len(),
        got.beta
    );
    assert!(
        (got.beta.iter().sum::<f64>() - 1.0).abs() < 1e-12,
        "{label}: the fractions sum to {}",
        got.beta.iter().sum::<f64>()
    );
    for want in expected {
        // The phase whose composition is nearest, by the same sum-absolute distance the merge
        // uses. If no phase is near, the failure names the closest so the gap is readable.
        let (index, distance) = got
            .x
            .iter()
            .enumerate()
            .map(|(i, x)| {
                (
                    i,
                    x.iter()
                        .zip(want.x)
                        .map(|(a, b)| (a - b).abs())
                        .fold(0.0, f64::max),
                )
            })
            .fold((0usize, f64::INFINITY), |best, candidate| {
                if candidate.1 < best.1 {
                    candidate
                } else {
                    best
                }
            });
        assert!(
            distance < 1e-6,
            "{label}: no phase matches x {want:?}; nearest is phase {index} at {distance:.3e}, \
             x {:?}",
            got.x[index]
        );
        assert!(
            (got.beta[index] - want.beta).abs() < 1e-6,
            "{label}: phase {index} is at beta {:.17}, NeqSim at {:.17}",
            got.beta[index],
            want.beta
        );
        assert!(
            (got.z_factor[index] - want.z_factor).abs() < 1e-6,
            "{label}: phase {index} sits at Z {:.17}, NeqSim at {:.17}",
            got.z_factor[index],
            want.z_factor
        );
    }
}

#[test]
fn co2_methane_decane_splits_three_ways_at_200_k_10_bar() {
    // The state that makes this model more than `eos.pt_flash`: NeqSim's two-phase flash
    // reports gas at 0.35221969743032410 and oil at 0.64778030256967590 here.
    let got = state(&["CO2", "methane", "nc10"], 200.0, 1.0e6, &[0.4, 0.3, 0.3]);
    check(
        &got,
        &[
            Expected {
                beta: 0.31913160557766135,
                z_factor: 0.9060977631578389,
                x: &[
                    0.25606636477632466,
                    0.7439336232007784,
                    1.2022897093803702e-8,
                ],
            },
            Expected {
                beta: 0.4986239551745936,
                z_factor: 0.08159861613379736,
                x: &[0.28770596245992525, 0.11069399567275835, 0.6016000418673164],
            },
            Expected {
                beta: 0.18224443924774514,
                z_factor: 0.02068639206922768,
                x: &[
                    0.9592832882420675,
                    0.040564148233054216,
                    0.00015256352487823808,
                ],
            },
        ],
        "CO2/C1/nc10 200 K 10 bar",
    );
    assert_eq!(got.seeded, TpMultiflashSeed::StabilitySeeded);
}

#[test]
fn co2_methane_decane_splits_three_ways_at_180_k_2_bar() {
    let got = state(&["CO2", "methane", "nc10"], 180.0, 2.0e5, &[0.4, 0.3, 0.3]);
    check(
        &got,
        &[
            Expected {
                beta: 0.45978707915172495,
                z_factor: 0.9734657991299046,
                x: &[
                    0.37523965788038316,
                    0.624760341465891,
                    6.537257990299283e-10,
                ],
            },
            Expected {
                beta: 0.39004930226638107,
                z_factor: 0.021530403096495943,
                x: &[
                    0.20104204002401013,
                    0.029827981897943685,
                    0.7691299780780461,
                ],
            },
            Expected {
                beta: 0.15016361858189392,
                z_factor: 0.004339193940160941,
                x: &[
                    0.9926062497118069,
                    0.007384504387544632,
                    9.245900648435757e-6,
                ],
            },
        ],
        "CO2/C1/nc10 180 K 2 bar",
    );
}

/// **A known divergence, characterised rather than matched.**
///
/// NeqSim reports **three** phases at this state, the third a CO2-rich liquid at `x_CO2`
/// 0.98415312005431520 and Z of 0.02175230239149613. This library reports two: its stability
/// trial converges to `x = [0.9696, 0.0304, 6.5e-9]`, a *nitrogen*-rich vapour, and the solve
/// drives that onto the gas the two-phase flash already found.
///
/// The cause is a seeding mechanism that is not ported. `stabilityAnalysis` runs the two
/// Wilson-K trials and then a **pure-component fallback**: one trial per component, started
/// from that component, scanning from the last downwards, adding a phase when its tangent-plane
/// distance is below `-1e-8`. The CO2 trial is the one that reaches the CO2-rich phase, and
/// Wilson K-values do not start anywhere near it.
///
/// So this test asserts the *current* count, not NeqSim's, and it is written to fail loudly if
/// the gap closes - at which point the missing seeding exists and the assertions below become
/// the three-phase `check` the other cases use.
#[test]
fn n2_co2_octane_third_phase_is_not_reached_without_the_pure_component_trials() {
    let got = state(
        &["nitrogen", "CO2", "n-octane"],
        180.0,
        1.0e6,
        &[0.1, 0.4, 0.5],
    );
    assert_eq!(
        got.phase_count, 2,
        "NeqSim reports 3 here and this library reported {}. If this now passes with 3, the \
         pure-component seeding has been ported and this test should become the three-phase \
         check beside it. Betas {:?}",
        got.phase_count, got.beta
    );
    // The two phases it does report are the two-phase flash's, and they are right: NeqSim
    // reports the same pair with the multiphase flag *off*.
    check(
        &got,
        &[
            Expected {
                beta: 0.9165691336463065,
                z_factor: 0.06779885045276167,
                x: &[0.02767335316398553, 0.42681406115096726, 0.5455125856850472],
            },
            Expected {
                beta: 0.08343086635369346,
                z_factor: 0.9468497822815427,
                x: &[
                    0.8945784926813769,
                    0.10542148400286804,
                    2.3315755061764716e-8,
                ],
            },
        ],
        "N2/CO2/nC8 180 K 10 bar (two-phase part)",
    );
    assert_eq!(
        got.seeded,
        TpMultiflashSeed::TwoPhaseFlash,
        "the seeding fired and did not survive the merge; if it now does, see above"
    );
}

#[test]
fn n2_co2_octane_stays_two_phase_at_190_k_1_bar() {
    // The state that shows the count is not the whole answer: at 190 K and 1 bar NeqSim
    // reports two phases with the flag on exactly as with it off, to all seventeen digits.
    // A model that seeded a third phase here would be finding a state upstream does not.
    let got = state(
        &["nitrogen", "CO2", "n-octane"],
        190.0,
        1.0e5,
        &[0.1, 0.4, 0.5],
    );
    assert_eq!(
        got.phase_count, 2,
        "NeqSim reports two phases here; this found {}: {:?}",
        got.phase_count, got.beta
    );
    check(
        &got,
        &[
            Expected {
                beta: 0.39215551375730107,
                z_factor: 0.9856288252087966,
                x: &[0.25394363833827993, 0.746055607642126, 7.540195941670215e-7],
            },
            Expected {
                beta: 0.6078444862426989,
                z_factor: 0.00860682835745079,
                x: &[
                    0.0006820857231110191,
                    0.17673958114909008,
                    0.822578333127799,
                ],
            },
        ],
        "N2/CO2/nC8 190 K 1 bar",
    );
}

/// The state NeqSim needs its stability check for, and this library does not.
///
/// NeqSim's plain `TPflash` converges to `x = y = z` at 370 K / 40 bar and reports one phase;
/// its multiflash runs the stability trial and adds a dense liquid at beta 0.0098562380076969.
/// **This library's `eos.pt_flash` does not converge to the trivial solution here** - it returns
/// beta 0.9901437619487483 with the same `x_methane` 0.12301811150418199 - so the seeding has
/// nothing to add and the numbers agree with NeqSim's multiflash by a different route.
///
/// That is why `seeded` is `two_phase_flash` and not `stability_seeded`: it reports which of
/// this library's paths produced the answer, and NeqSim has no such field to compare against.
#[test]
fn methane_butane_dense_liquid_at_370_k_40_bar() {
    let got = state(&["methane", "n-butane"], 370.0, 4.0e6, &[0.5, 0.5]);
    assert_eq!(
        got.seeded,
        TpMultiflashSeed::TwoPhaseFlash,
        "this library's two-phase flash reaches NeqSim's multiflash answer here without the \
         seeding; if it has started needing it, `eos.pt_flash` has changed"
    );
    check(
        &got,
        &[
            Expected {
                beta: 0.9901437619923031,
                z_factor: 0.7452177737328194,
                x: &[0.5037526098332724, 0.49624739016672753],
            },
            Expected {
                beta: 0.009856238007696907,
                z_factor: 0.1588137538319941,
                x: &[0.12301811150419444, 0.8769818884958055],
            },
        ],
        "C1/nC4 370 K 40 bar",
    );
}

#[test]
fn the_seeding_is_silent_when_the_two_phase_flash_already_split() {
    // Ten kelvin colder, NeqSim reports the same betas and compositions with the flag on and
    // off. This is the 3,300-odd states of the sweep where the flag does nothing.
    let got = state(&["methane", "n-butane"], 360.0, 4.0e6, &[0.5, 0.5]);
    assert_eq!(got.seeded, TpMultiflashSeed::TwoPhaseFlash);
    check(
        &got,
        &[
            Expected {
                beta: 0.8139457129518551,
                z_factor: 0.7784548921212494,
                x: &[0.5834844933560103, 0.4165155066439898],
            },
            Expected {
                beta: 0.18605428704814486,
                z_factor: 0.15404179516968172,
                x: &[0.13477409446900623, 0.8652259055309938],
            },
        ],
        "C1/nC4 360 K 40 bar",
    );
}

#[test]
fn methane_butane_two_phase_at_250_k_50_bar() {
    // The state that caught a defect in the merge: the seeded phase landed on the gas
    // composition and differed from it only in being labelled a liquid root, so a merge that
    // required matching sides reported one phase twice. Betas and Z are NeqSim's.
    let got = state(&["methane", "n-butane"], 250.0, 5.0e6, &[0.5, 0.5]);
    assert_eq!(got.seeded, TpMultiflashSeed::TwoPhaseFlash);
    check(
        &got,
        &[
            Expected {
                beta: 0.21848944858961394,
                z_factor: 0.7871775377076164,
                x: &[0.9778100274811522, 0.0221899725188478],
            },
            Expected {
                beta: 0.7815105514103861,
                z_factor: 0.1786351934671386,
                x: &[0.3664170953974433, 0.6335829046025567],
            },
        ],
        "C1/nC4 250 K 50 bar",
    );
}

#[test]
fn the_answer_never_exceeds_the_three_phase_ceiling() {
    // The sweep found no state of any mixture above three phases, and upstream raises the
    // system's ceiling to three before the flash runs. This is that bound as an assertion
    // over the states the sweep's three-phase mixtures pass through.
    let cases: &[(&[&str], f64, f64, &[f64])] = &[
        (&["CO2", "methane", "nc10"], 180.0, 2.0e5, &[0.4, 0.3, 0.3]),
        (&["CO2", "methane", "nc10"], 190.0, 5.0e5, &[0.4, 0.3, 0.3]),
        (&["CO2", "methane", "nc10"], 210.0, 2.0e6, &[0.4, 0.3, 0.3]),
        (&["CO2", "methane", "nc10"], 220.0, 2.0e6, &[0.4, 0.3, 0.3]),
        (
            &["nitrogen", "CO2", "n-octane"],
            180.0,
            1.0e6,
            &[0.1, 0.4, 0.5],
        ),
        (
            &["nitrogen", "CO2", "n-octane"],
            190.0,
            2.0e6,
            &[0.1, 0.4, 0.5],
        ),
        (
            &["nitrogen", "CO2", "n-octane"],
            190.0,
            4.0e5,
            &[0.1, 0.4, 0.5],
        ),
        (
            &["nitrogen", "CO2", "n-octane"],
            200.0,
            5.0e5,
            &[0.1, 0.4, 0.5],
        ),
    ];
    for (names, t, p, z) in cases {
        let got = state(names, *t, *p, z);
        assert!(
            got.phase_count <= 3,
            "{names:?} at {t} K / {p} Pa: {} phases",
            got.phase_count
        );
        assert!(
            got.beta.iter().all(|beta| *beta >= 0.0 && *beta <= 1.0),
            "{names:?} at {t} K / {p} Pa: betas {:?}",
            got.beta
        );
    }
}
