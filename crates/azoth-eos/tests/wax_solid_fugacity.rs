//! Spec-driven tests for `eos.wax_solid_fugacity`.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::spec_gen;
use azoth_eos::wax_solid_fugacity;
use azoth_test_support as common;

const CALC_ID: &str = "eos.wax_solid_fugacity";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::results::WaxSolidFugacityResult {
    wax_solid_fugacity(
        common::input(case, "molar_mass"),
        kelvins(common::input(case, "tc")),
        pascals(common::input(case, "pc")),
        common::input(case, "omega"),
        common::input(case, "heat_of_fusion"),
        kelvins(common::input(case, "triple_point_temperature")),
        kelvins(common::input(case, "T")),
        pascals(common::input(case, "P")),
        common::input_str(case, "eos"),
    )
    .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    common::assert_skips_are_explained(spec);

    let mut executed = 0;
    for case in spec.all_tests() {
        if !case.is_active() {
            continue;
        }
        match case.kind {
            "worked_example" | "reference" => {
                let result = call(case);
                let context = format!("{}::{}", spec.id, case.id);
                common::assert_close(
                    result.fugacity_coefficient,
                    common::expected(case, "fugacity_coefficient"),
                    case.tolerance,
                    &format!("{context} (fugacity coefficient)"),
                );
                common::assert_consistent(&result, &context);
            }
            "property" => match case.property {
                Some("unit_round_trip") => {}
                other => panic!("{}::{}: unhandled property {other:?}", spec.id, case.id),
            },
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 2,
        "expected several active cases, ran {executed}"
    );
}

/// **All twenty-one of the capture's states**, not only the two the spec records.
///
/// `validation/neqsim/captures/wax_reference_probe.tsv` prints the seven wax cuts at 285, 275
/// and 261 K, each with the reference liquid's own coefficient and molar volume beside the
/// product NeqSim reports. The spec pins two of them as cases; this runs the rest, because
/// the pair of assertions the two cases make is one shape and the twenty-one are the range -
/// the coefficient spans twenty-four decades across the cuts, so a term dropped at either end
/// is invisible at the other.
///
/// The inputs are the capture's own `component[*]` rows; the constants are stated here rather
/// than parsed so the test reads as the measurement it is.
#[test]
fn the_whole_capture() {
    // `molar_mass` kg/mol, `tc` K, `pc` Pa, `omega`, `dh_fus` J/mol, `t_tp` K.
    const CUTS: [(f64, f64, f64, f64, f64, f64); 7] = [
        (
            0.170,
            632.271665693285,
            1526386.94925152,
            0.674258471831839,
            26418.6063505621,
            260.290076470588,
        ),
        (
            0.289092472581105,
            747.128826726606,
            1037555.15125897,
            0.962417531753728,
            53900.9280050369,
            312.288568891688,
        ),
        (
            0.33021895027422,
            776.619447657179,
            1004153.62726608,
            1.04695278637736,
            63494.4152292418,
            322.055087380289,
        ),
        (
            0.37853703499591,
            809.455923999465,
            979158.7377216,
            1.13531038666064,
            74832.9995039473,
            331.116944513871,
        ),
        (
            0.440544507719916,
            849.445499029544,
            960855.768001826,
            1.229888278603,
            89490.8975249053,
            340.240259652719,
        ),
        (
            0.521130807870325,
            898.791034295216,
            950425.515252361,
            1.31836543537214,
            108720.173205326,
            349.429861434454,
        ),
        (
            0.685054415944094,
            993.100301671903,
            951903.245094917,
            1.36697955720612,
            148461.408421741,
            362.982038040391,
        ),
    ];
    // NeqSim's `wax_phi` at each of the three temperatures, in the capture's cut order.
    const CAPTURED: [[f64; 7]; 3] = [
        [
            3.498_989_803_6e-05,
            1.533_287_137_6e-10,
            2.905_240_942_9e-12,
            2.791_885_350_2e-14,
            7.873_533_151_7e-17,
            5.293_708_355_4e-20,
            1.705_038_741_1e-25,
        ],
        [
            8.429_736_227_1e-06,
            1.586_496_868_9e-11,
            2.278_228_556_3e-13,
            1.589_088_924_8e-15,
            3.011_767_224_2e-18,
            1.251_432_383_1e-21,
            1.828_810_752_3e-27,
        ],
        [
            9.640_371_758_6e-07,
            5.119_135_082_3e-13,
            4.856_743_297_0e-15,
            2.099_126_473_6e-17,
            2.199_975_387_5e-20,
            4.473_599_548_9e-24,
            2.052_002_534_4e-30,
        ],
    ];

    for (row, &t) in [285.0, 275.0, 261.0].iter().enumerate() {
        for (index, &(molar_mass, tc, pc, omega, dh, ttp)) in CUTS.iter().enumerate() {
            let result = wax_solid_fugacity(
                molar_mass,
                kelvins(tc),
                pascals(pc),
                omega,
                dh,
                kelvins(ttp),
                kelvins(t),
                pascals(5.0e5),
                "srk",
            )
            .expect("computes");
            common::assert_close(
                result.fugacity_coefficient,
                CAPTURED[row][index],
                1.0e-9,
                &format!("cut {index} at {t} K"),
            );
        }
    }
}

/// **The pressure term is referred to one bar, and it is not the `1e5`-small one.**
///
/// It is the term this port got wrong first: evaluated with `v` in m³/mol against a pressure
/// in **bar** it comes out `1e5` too small, and the coefficient is then out by a per cent that
/// grows with the cut's molar mass - `0.6 %` on the lightest cut and `1.4 %` on the heaviest,
/// which is close enough to look like round-off. NeqSim's `refPressure` is `1.0` in its
/// bar-valued code and its molar volume is in `1e-5 m³/mol`, so in pascals the term is
/// `(v - 0.9v)(P - 1.0e5)/(R T)` and the `1e-5` is not a conversion to apply.
///
/// The two candidate values are built here from azoth's own reference state, so this is a
/// check of the units rather than a restatement of the answer: the wrong one is the right one
/// times `exp(pres_correct - pres_wrong)`, and the test is that the function reports the
/// first and not the second.
#[test]
fn the_volume_term_is_not_the_one_with_the_pressure_in_bar() {
    use azoth_eos::mixture::{Component, Mixture, RootSide};
    use azoth_eos::{Alpha, Cubic};

    let (molar_mass, tc, pc, omega, dh, ttp) = (
        0.170,
        632.271665693285,
        1526386.94925152,
        0.674258471831839,
        26418.6063505621,
        260.290076470588,
    );
    let (t, p) = (261.0, 5.0e5);
    let r = 8.314_462_1;

    let reported = wax_solid_fugacity(
        molar_mass,
        kelvins(tc),
        pascals(pc),
        omega,
        dh,
        kelvins(ttp),
        kelvins(t),
        pascals(p),
        "srk",
    )
    .expect("computes")
    .fugacity_coefficient;

    // The same expression, with the two candidate volume terms, from the same reference state.
    let reference = Mixture::new(
        vec![Component::new(kelvins(tc), pascals(pc), omega).expect("positive")],
        vec![0.0],
    )
    .expect("one component")
    .with_cubic(Cubic::Srk)
    .with_alpha(Alpha::Srk);
    let state = reference
        .phase_state(
            &reference
                .reduced_parameters(kelvins(t), pascals(p))
                .expect("reduces"),
            &[1.0],
            RootSide::Liquid,
        )
        .expect("liquid root");
    let v_liquid = state.z * r * t / p;
    let mw_grams = molar_mass * 1000.0;
    let delta_cp = (0.3033 * mw_grams - 4.635e-4 * mw_grams * t) * 4.184;
    let ratio = ttp / t;
    let base = -dh / (r * t) * (1.0 - t / ttp) + delta_cp / r * (ratio - 1.0 - ratio.ln());
    let correct = -(v_liquid - 0.9 * v_liquid) * (p - 1.0e5) / r / t;
    let in_bar = -(v_liquid - 0.9 * v_liquid) * (p / 1.0e5 - 1.0) / r / t;

    let expected = state.ln_phi[0].exp() * (base + correct).exp();
    let wrong = state.ln_phi[0].exp() * (base + in_bar).exp();
    assert!(
        (reported / expected - 1.0).abs() < 1.0e-9,
        "the coefficient is {reported} and the term referred to one bar gives {expected}"
    );
    assert!(
        (wrong / reported - 1.0).abs() > 1.0e-3,
        "the bar-valued term gives {wrong} against {reported}, so this test cannot tell them \
         apart on this cut"
    );
}
