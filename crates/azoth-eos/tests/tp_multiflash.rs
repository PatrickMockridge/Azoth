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
    phase_fractions: f64,
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
        got.phase_fractions
    );
    assert!(
        (got.phase_fractions.iter().sum::<f64>() - 1.0).abs() < 1e-12,
        "{label}: the fractions sum to {}",
        got.phase_fractions.iter().sum::<f64>()
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
            (got.phase_fractions[index] - want.phase_fractions).abs() < 1e-6,
            "{label}: phase {index} is at beta {:.17}, NeqSim at {:.17}",
            got.phase_fractions[index],
            want.phase_fractions
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
                phase_fractions: 0.31913160557766135,
                z_factor: 0.9060977631578389,
                x: &[
                    0.25606636477632466,
                    0.7439336232007784,
                    1.2022897093803702e-8,
                ],
            },
            Expected {
                phase_fractions: 0.4986239551745936,
                z_factor: 0.08159861613379736,
                x: &[0.28770596245992525, 0.11069399567275835, 0.6016000418673164],
            },
            Expected {
                phase_fractions: 0.18224443924774514,
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
                phase_fractions: 0.45978707915172495,
                z_factor: 0.9734657991299046,
                x: &[
                    0.37523965788038316,
                    0.624760341465891,
                    6.537257990299283e-10,
                ],
            },
            Expected {
                phase_fractions: 0.39004930226638107,
                z_factor: 0.021530403096495943,
                x: &[
                    0.20104204002401013,
                    0.029827981897943685,
                    0.7691299780780461,
                ],
            },
            Expected {
                phase_fractions: 0.15016361858189392,
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

/// **The state that was a recorded divergence until the reference phase was fixed.**
///
/// NeqSim reports three phases here, the third a CO2-rich liquid at `x_CO2` 0.98415312005431520
/// and Z of 0.02175230239149613, and this library used to report the two-phase flash's own
/// answer. The pure-component fallback was already ported; what was wrong was where it was
/// measured from. `stabilityAnalysis` reads its reference potentials from **phase 0 of the
/// flashed system** - the least-dense phase, because `TPflash` orders by density - and marks a
/// trial trivial when it reproduces **any** existing phase. This kernel derived both from the
/// feed, which is phase 0 only for a single-phase feed. With the reference fixed the CO2 trial
/// reaches the CO2-rich phase and the answer is NeqSim's.
#[test]
fn n2_co2_octane_splits_three_ways_at_180_k_10_bar() {
    let got = state(
        &["nitrogen", "CO2", "n-octane"],
        180.0,
        1.0e6,
        &[0.1, 0.4, 0.5],
    );
    check(
        &got,
        &[
            Expected {
                phase_fractions: 0.08437940065407534,
                z_factor: 0.9495405616635599,
                x: &[0.9132770769209018, 0.08672289733209693, 2.574700136832e-8],
            },
            Expected {
                phase_fractions: 0.7563697015044588,
                z_factor: 0.07716944734261964,
                x: &[0.02712244085245802, 0.31195737408395324, 0.6609201850635887],
            },
            Expected {
                phase_fractions: 0.1592508978414658,
                z_factor: 0.02175230239149613,
                x: &[
                    0.015218973076522303,
                    0.9841531200543152,
                    0.000627906869162516,
                ],
            },
        ],
        "N2/CO2/nC8 180 K 10 bar",
    );
    assert_eq!(got.seeded, TpMultiflashSeed::StabilitySeeded);
}

/// The state where upstream's **enhanced** analysis is what decides the answer.
///
/// Ablating `stabilityAnalysisEnhanced` from `run()` leaves an aqueous phase at beta `1e-5`
/// (the inline seed); with it, the enhanced analysis adds a phase, the bounded rerun solves
/// three, and the merge removes the addition again - two phases. Both chains are in
/// `captures/tp_multiflash_seed_probe.tsv`; the expected values below are the default one.
#[test]
fn the_post_flash_enhanced_check_decides_a_water_bearing_state() {
    // NeqSim's `getz` is the mole numbers over their total; the probe's feed is
    // 0.5 / 0.5 / 1e-4 mol.
    let got = state(
        &["methane", "n-butane", "water"],
        370.0,
        40.0e5,
        &[
            0.49995000499950004,
            0.49995000499950004,
            9.999000099990002e-5,
        ],
    );
    check(
        &got,
        &[
            Expected {
                phase_fractions: 0.9902639834091616,
                z_factor: 0.7452175536344108,
                x: &[
                    0.5036561370177296,
                    0.49624322017919587,
                    0.00010064280307448483,
                ],
            },
            Expected {
                phase_fractions: 0.009736016590838425,
                z_factor: 0.15881276682734707,
                x: &[
                    0.12299408866860612,
                    0.8769723187503989,
                    3.359258099504507e-5,
                ],
            },
        ],
        "C1/nC4/water 370 K 40 bar",
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
        got.phase_count, got.phase_fractions
    );
    check(
        &got,
        &[
            Expected {
                phase_fractions: 0.39215551375730107,
                z_factor: 0.9856288252087966,
                x: &[0.25394363833827993, 0.746055607642126, 7.540195941670215e-7],
            },
            Expected {
                phase_fractions: 0.6078444862426989,
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
                phase_fractions: 0.9901437619923031,
                z_factor: 0.7452177737328194,
                x: &[0.5037526098332724, 0.49624739016672753],
            },
            Expected {
                phase_fractions: 0.009856238007696907,
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
                phase_fractions: 0.8139457129518551,
                z_factor: 0.7784548921212494,
                x: &[0.5834844933560103, 0.4165155066439898],
            },
            Expected {
                phase_fractions: 0.18605428704814486,
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
                phase_fractions: 0.21848944858961394,
                z_factor: 0.7871775377076164,
                x: &[0.9778100274811522, 0.0221899725188478],
            },
            Expected {
                phase_fractions: 0.7815105514103861,
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
            got.phase_fractions
                .iter()
                .all(|beta| *beta >= 0.0 && *beta <= 1.0),
            "{names:?} at {t} K / {p} Pa: betas {:?}",
            got.phase_fractions
        );
    }
}

/// **The aqueous seeds this port declines, measured on NeqSim's own gas.**
///
/// `WaterDewPointTemperatureMultiphaseFlash` brackets on where a *multiphase flash* reports
/// an aqueous phase, so a water dew point built on it is only as good as the seeding
/// underneath. This is that seeding, on the gas NeqSim's own `WaterDewPointTest` builds -
/// CO2 0.02, nitrogen 0.01, methane 0.9 - x, ethane 0.05, propane 0.01, i-butane 0.005,
/// n-butane 0.005 and water x - on a plain SRK at 70 bara.
///
/// **At the temperature NeqSim calls the water dew point this flash is single-phase.** Its
/// own capture puts the multiphase dew point at `260.545463562012` K for 35 ppm water and
/// `284.302406311035` K for 200 ppm, each with an aqueous phase beside the gas; here the
/// flash answers one gas phase at both. The aqueous phase does appear - at `217` K and
/// `260` K - so the trial is not broken; it simply does not find the water-rich stationary
/// point at the temperatures upstream calls the dew point.
///
/// **This test was written as a canary and it fired, on 2026-10-13.** When it was added the
/// assertions were the other way round: this flash answered one gas phase at NeqSim's own
/// dew-point temperature, and the test said so. `eos.multiflash_stability_test` landed -
/// `TPmultiflash`'s **own** `stabilityAnalysis`, the method its `run()` calls and the one the
/// base class's two-trial test shadows - and with it the pure-component trial that reaches
/// the aqueous phase. So the assertions are inverted here rather than deleted, and the two
/// classes downstream, `eos.water_dew_point_multiphase` and `eos.water_dew_point_line`, are
/// portable now.
#[test]
fn the_aqueous_phase_appears_far_colder_than_neqsim_calls_the_dew_point() {
    use azoth_eos::hydrate_inhibitor_wt::{PhaseLabel, label};

    const NAMES: [&str; 8] = [
        "co2", "nitrogen", "methane", "ethane", "propane", "i-butane", "n-butane", "water",
    ];

    fn gas(ppm: f64) -> Vec<f64> {
        let x = ppm * 1e-6;
        vec![0.02, 0.01, 0.9 - x, 0.05, 0.01, 0.005, 0.005, x]
    }

    /// Whether the flash at `t` reports an aqueous phase, or `None` if it does not converge.
    fn aqueous_at(mixture: &azoth_eos::Mixture, t: f64, z: &[f64]) -> Option<bool> {
        let result = tp_multiflash(mixture, kelvins(t), pascals(7.0e6), z).ok()?;
        let reduced = mixture
            .reduced_parameters(kelvins(t), pascals(7.0e6))
            .expect("the state's reduced parameters");
        let labels = result
            .x
            .iter()
            .zip(&result.z_factor)
            .map(|(x, z_factor)| label(mixture, &reduced, x, *z_factor))
            .collect::<azoth_core::Result<Vec<_>>>()
            .expect("every phase carries a class");
        Some(labels.contains(&PhaseLabel::Aqueous))
    }

    for (ppm, neqsim_dew_point) in [(35.0, 260.545463562012), (200.0, 284.302406311035)] {
        let z = gas(ppm);
        let (mixture, _) =
            mixture_of(&NAMES, Cubic::Srk, None).expect("the databank resolves these names");

        assert_eq!(
            aqueous_at(&mixture, neqsim_dew_point, &z),
            Some(true),
            "{ppm} ppm: NeqSim's own multiphase dew point is {neqsim_dew_point} K, and this \
             flash answers one gas phase there"
        );

        // And ten kelvin warmer it does not, so this is a crossing rather than a flash that
        // reports an aqueous phase wherever it is asked.
        assert_eq!(
            aqueous_at(&mixture, neqsim_dew_point + 10.0, &z),
            Some(false),
            "{ppm} ppm: an aqueous phase ten kelvin above the dew point is one the seeding \
             invented"
        );
    }
}

/// **A stalled two-phase flash is a seed, and the seeding runs anyway.**
///
/// Upstream's flash is where this model's phase set *starts*, not what it answers, and there
/// are states where it reaches its cap and stops short. This is one: on the gas
/// `eos.water_dew_point_line` sweeps, at `111.55555` bara and `261.876` K, `eos.pt_flash`
/// takes its full three hundred steps to a residual of `5.3e-5` against the `1e-10` it
/// declares. Refusing there, which is what this model did until 2026-10-13, loses the answer
/// the seeding and the fraction solve would have found - and the two-phase flash's own id
/// keeps its strictness, because that one *is* answering rather than seeding.
///
/// Pinned here because it is the state a water dew point brackets through: without the
/// fallback the sweep at that pressure has no point at all.
#[test]
fn a_stalled_two_phase_flash_still_seeds() {
    use azoth_core::CalcResult;
    use azoth_core::warning::WarningCode;

    let names = [
        "co2", "nitrogen", "methane", "ethane", "propane", "i-butane", "n-butane", "water",
    ];
    let x = 35e-6;
    let z = [0.02, 0.01, 0.9 - x, 0.05, 0.01, 0.005, 0.005, x];
    let (mixture, _) =
        mixture_of(&names, Cubic::Srk, None).expect("the databank resolves these names");

    let got = tp_multiflash(&mixture, kelvins(261.876), pascals(1.1155555e7), &z)
        .expect("a stalled seed is not a failure of this model");

    assert_eq!(got.phase_count, 2, "the aqueous phase the seeding found");
    assert!(
        got.has_warning(WarningCode::SolverNotConverged),
        "the stalled flash is reported rather than hidden"
    );
}
