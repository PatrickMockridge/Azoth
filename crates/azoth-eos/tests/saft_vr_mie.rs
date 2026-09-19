//! SAFT-VR-Mie's parameters and the effective diameter, against NeqSim's own.
//!
//! Every expected number comes from `validation/neqsim/SaftVrMieProbe.java`, which prints
//! `PhaseSAFTVRMie`'s intermediates with the method each came from. The diameter is the
//! first layer and the one everything else is built on, so it is pinned before anything
//! uses it.

use azoth_eos::saft_vr_mie::{MieComponent, barker_henderson, mie_prefactor};

fn methane() -> MieComponent {
    MieComponent {
        m: 1.0,
        lambda_r: 12.65,
        lambda_a: 6.0,
        sigma: 3.7412e-10,
        epsik: 153.36,
    }
}

fn n_butane() -> MieComponent {
    MieComponent {
        m: 1.8514,
        lambda_r: 13.457,
        lambda_a: 6.0,
        sigma: 4.0887e-10,
        epsik: 273.64,
    }
}

/// To the probe's printed digits, which is fifteen.
fn matches(actual: f64, expected: f64, context: &str) {
    let relative = (actual / expected - 1.0).abs();
    assert!(
        relative < 1.0e-13,
        "{context}: {actual} against the probe's {expected}, a relative {relative:e}"
    );
}

#[test]
fn the_effective_diameter_is_neqsims() {
    matches(
        methane().d(300.0).expect("a diameter"),
        3.58753830702527e-10,
        "methane at 300 K",
    );
    matches(
        methane().d(350.0).expect("a diameter"),
        3.57210695782575e-10,
        "methane at 350 K",
    );
    matches(
        n_butane().d(350.0).expect("a diameter"),
        3.96665785055817e-10,
        "n-butane at 350 K",
    );
}

/// The prefactor is the constant that puts the well's minimum at `-epsilon`, and it is
/// checked at the two exponents the table actually carries rather than only at its own
/// definition: `C(12.65, 6)` is the number the three states above are built on.
#[test]
fn the_mie_prefactor_is_the_published_constant() {
    // `C(12, 6) = 4` exactly: the 12-6 Lennard-Jones case, where the prefactor's algebra
    // collapses to `2 * 2^1`.
    matches(mie_prefactor(12.0, 6.0), 4.0, "the 12-6 case");
    matches(
        mie_prefactor(12.65, 6.0),
        3.728_592_467_848_367,
        "methane's",
    );
    matches(
        mie_prefactor(13.457, 6.0),
        3.456_527_193_236_619,
        "n-butane's",
    );
}

/// **The weak branch is a ten-point rule and the rule is coarse**, so the check is
/// against the rule rather than against the integral.
///
/// Measured against a converged integral, NeqSim's ten Gauss-Legendre points carry the
/// Barker-Henderson integral to only `1e-4` at `theta = 0.9` and `4e-3` at `0.05`: the
/// integrand is a step-like `1 - exp(-u*)` and ten points cannot resolve the knee. No
/// state in the probe is on this branch, so there is no oracle number to pin - and the
/// temptation is to "fix" the rule, which would be porting a different model. What is
/// checked instead is a second transcription of NeqSim's own sum, node for node, which is
/// the failure this can actually have: a swapped node with its weight, or a loop that
/// starts or ends one early.
#[test]
fn the_weak_branch_is_neqsims_ten_point_rule() {
    const NODES: [f64; 10] = [
        0.013_046_735_74,
        0.067_468_316_65,
        0.160_295_215_85,
        0.283_302_302_94,
        0.425_562_830_50,
        0.574_437_169_50,
        0.716_697_697_06,
        0.839_704_784_15,
        0.932_531_683_35,
        0.986_953_264_26,
    ];
    const WEIGHTS: [f64; 10] = [
        0.033_335_672_15,
        0.074_725_674_58,
        0.109_543_181_26,
        0.134_633_359_65,
        0.147_762_112_36,
        0.147_762_112_36,
        0.134_633_359_65,
        0.109_543_181_26,
        0.074_725_674_58,
        0.033_335_672_15,
    ];

    // Two states on the weak branch, one shallow and one close to the crossing.
    for theta in [0.2, 0.9] {
        let mut expected = 0.0;
        for (node, weight) in NODES.iter().zip(WEIGHTS.iter()) {
            expected += weight * (1.0 - (-theta * (node.powf(-12.0) - node.powf(-6.0))).exp());
        }
        let actual = barker_henderson(theta, 12.0, 6.0);
        assert!(
            (actual / expected - 1.0).abs() < 1.0e-14,
            "at theta = {theta} the model gives {actual} and the transcribed rule {expected}"
        );
    }

    // And the rule's own error, so that the next reader does not take it for an integral.
    // Simpson on 2e6 points agrees with a composite 40-point rule to all digits above.
    let converged = [0.795_114_723_922_264, 0.931_751_726_994_897];
    for (theta, reference) in [0.05, 0.9].into_iter().zip(converged) {
        let rule = barker_henderson(theta, 12.0, 6.0);
        let error = (rule / reference - 1.0).abs();
        assert!(
            error > 1.0e-5,
            "at theta = {theta} the rule is within {error:e} of the integral, which is \
             tighter than this rule manages - one of the two has changed"
        );
    }
}

/// A component with no set is refused, and so is a state that is not one.
#[test]
fn absent_or_impossible_inputs_are_refused() {
    let absent = MieComponent {
        m: 0.0,
        lambda_r: 12.0,
        lambda_a: 6.0,
        sigma: 0.0,
        epsik: 0.0,
    };
    assert!(!absent.has_parameters(), "a zero segment number is absent");

    let methane = methane();
    assert!(methane.d(0.0).is_err(), "zero kelvin is not a state");

    // A shallow well: the effective sphere is smaller than the potential's range, which is
    // what the integral says a Barker-Henderson diameter is.
    let shallow = MieComponent {
        epsik: 90.0,
        ..methane
    };
    let d = shallow.d(400.0).expect("a diameter");
    assert!(d > 0.0 && d < shallow.sigma, "d = {d} against sigma");
    let inverted = MieComponent {
        lambda_r: 6.0,
        ..methane
    };
    assert!(
        inverted.d(300.0).is_err(),
        "an attractive exponent that is not below the repulsive one is not a potential"
    );
}
