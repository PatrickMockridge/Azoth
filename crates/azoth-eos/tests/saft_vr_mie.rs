//! SAFT-VR-Mie's parameters and the effective diameter, against NeqSim's own.
//!
//! Every expected number comes from `validation/neqsim/SaftVrMieProbe.java`, which prints
//! `PhaseSAFTVRMie`'s intermediates with the method each came from. The diameter is the
//! first layer and the one everything else is built on, so it is pinned before anything
//! uses it.

use azoth_eos::saft_vr_mie::{
    MieComponent, a_s1_bare, a1_mie, a2_mie, a3_mie, b_bare, barker_henderson, chain_contact_value,
    chain_g1, chain_g2, contact_value_0, dispersion_pair_sum, eta_effective, g_hs, k_hs, mie_alpha,
    mie_prefactor, pressure_over_rt, state,
};

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

/// **`1e-8`, one order looser than the diameter's `1e-13`**, and the difference is what
/// the oracle can offer rather than what the port achieves. The diameter is a function of
/// the parameters alone and reproduces to fifteen digits; these layers are read at the
/// probe's own `eta`, which its `volInit` builds from a volume the phase does not report
/// elsewhere - the `1.1e-9` its bookkeeping costs, recorded in `tests/pcsaft.rs` for the
/// same reason.
fn via_eta(actual: f64, expected: f64, context: &str) {
    let relative = (actual / expected - 1.0).abs();
    assert!(
        relative < 1.0e-8,
        "{context}: {actual} against the probe's {expected}, a relative {relative:e}"
    );
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

/// The effective packing fraction, the bare `aS1`/`B` pair and `K_HS`, at the two states
/// the probe prints.
#[test]
fn the_bare_layers_are_neqsims() {
    // Methane at 300 K: eta = 0.0315889086860159, x0 = 1.04283207030119.
    let (eta, x0) = (0.031_588_908_686_015_9, 1.042_832_070_301_19);
    // **`eta_eff` has no oracle of its own**: the probe does not print it, and it is not
    // `eta` - the softness of the potential moves the hard-sphere reference off the real
    // packing fraction, which is the whole point of the correction. What pins it is the two
    // layers built on it, `aS1Bare` and `B`, below.
    assert!(
        (eta_effective(eta, 6.0) - eta).abs() > 0.01,
        "the effective packing fraction should differ from the packing fraction"
    );
    via_eta(a_s1_bare(eta, 6.0), -0.346_542_502_308_464, "aS1Bare(6)");
    via_eta(
        a_s1_bare(eta, 12.65),
        -0.110_203_246_791_488,
        "aS1Bare(12.65)",
    );
    via_eta(b_bare(eta, 6.0, x0), 0.042_574_771_360_466_2, "B(6)");
    via_eta(b_bare(eta, 12.65, x0), 0.037_267_259_194_146_6, "B(12.65)");
    via_eta(k_hs(eta), 0.778_171_407_751_256, "K_HS");
    via_eta(
        mie_alpha(12.65, 6.0),
        0.856_481_516_794_185,
        "Lafitte alpha",
    );
}

/// The contact value's quartic and the two chain perturbations, per component.
#[test]
fn the_chain_layers_are_neqsims() {
    let methane = methane();
    let butane = n_butane();
    // The binary at 350 K: eta = 0.0280863049064808, and each component's own x0.
    let eta = 0.028_086_304_906_480_8;
    for (component, x0, g_hs_0, g1, g2) in [
        (
            methane,
            1.047_337_060_219_82,
            1.067_442_307_213_12,
            -0.139_081_494_124_591,
            0.086_337_566_007_181_9,
        ),
        (
            butane,
            1.030_766_996_811_85,
            1.069_756_432_606_19,
            -0.137_475_062_806_657,
            0.106_318_540_526_763,
        ),
    ] {
        let t = 350.0;
        let c_mie = mie_prefactor(component.lambda_r, component.lambda_a);
        let zeta = eta * x0 * x0 * x0;
        via_eta(contact_value_0(eta, x0), g_hs_0, "gHS0");
        via_eta(
            chain_g1(eta, component.lambda_r, component.lambda_a, c_mie, x0),
            g1,
            "g1",
        );
        via_eta(
            chain_g2(
                eta,
                zeta,
                component.lambda_r,
                component.lambda_a,
                component.epsik / t,
                c_mie,
                x0,
            ),
            g2,
            "g2",
        );
    }
}

/// **The mixture's `gHS` is the blended one and a pure fluid's is the contact value.**
///
/// This is the layer a port is most likely to stop short of: methane's `gHS` above is the
/// Carnahan-Starling number to every printed digit, and the binary's is `1.02809527174773`,
/// which the CS form does not give.
#[test]
fn the_mixture_contact_value_is_the_blended_one() {
    let eta = 0.028_086_304_906_480_8;
    let methane_d = methane().d(350.0).expect("a diameter");
    let butane_d = n_butane().d(350.0).expect("a diameter");
    let diameters = [methane_d, butane_d];

    let blended = chain_contact_value(
        &[methane(), n_butane()],
        &[0.6, 0.4],
        350.0,
        eta,
        &diameters,
    );
    via_eta(blended, 1.028_095_271_747_73, "the binary's gHS");

    // And the pure fluid keeps the Carnahan-Starling value, because `w = x (m - 1)` is zero.
    let pure_eta = 0.031_588_908_686_015_9;
    let pure = chain_contact_value(
        &[methane()],
        &[1.0],
        300.0,
        pure_eta,
        &[methane().d(300.0).expect("a diameter")],
    );
    assert_eq!(
        pure,
        g_hs(pure_eta),
        "a one-segment fluid keeps the CS value"
    );
}

/// The three dispersion terms, at methane's state and at the binary's two components.
///
/// Methane 300 K: the probe's `A1Disp`, `A2Disp` and `A3Disp` are `-0.192886652170265`,
/// `-0.0175616577606500` and `-0.00145830264597508`. The binary's are the pair sums, which
/// this does not yet reproduce - so only the sub-layers are checked there.
#[test]
fn the_dispersion_terms_are_neqsims() {
    let methane = methane();
    let (eta, t) = (0.031_588_908_686_015_9, 300.0);
    let x0 = 1.042_832_070_301_19;
    let eps = methane.epsik / t;
    let c_mie = mie_prefactor(methane.lambda_r, methane.lambda_a);
    let zeta = eta * x0 * x0 * x0;
    let (lr, la) = (methane.lambda_r, methane.lambda_a);

    via_eta(
        a1_mie(eta, lr, la, eps, c_mie, x0),
        -0.192_886_652_170_265,
        "A1",
    );
    via_eta(
        a2_mie(eta, zeta, lr, la, eps, c_mie, x0),
        -0.017_561_657_760_650_0,
        "A2",
    );
    via_eta(a3_mie(zeta, lr, la, eps), -0.001_458_302_645_975_08, "A3");

    // **The binary's components have their own `A1`, and the mixture's is a pair sum over
    // them** - so a per-component `A1` is checked against the *pure* states the probe
    // prints, and the sum is left to the layer that does it.
    let butane = n_butane();
    let beta = mie_alpha(butane.lambda_r, butane.lambda_a);
    assert!(
        beta > 0.0 && beta < 1.0,
        "Lafitte's alpha is a softness: {beta}"
    );
}

/// The mixture's dispersion is a **pair sum with its own cross parameters**, and this is
/// the check that says so: methane/n-butane at 350 K gives `-0.197247388294105`,
/// `-0.0230307443986750` and `-0.00229661243223986`.
#[test]
fn the_mixtures_dispersion_is_a_pair_sum() {
    let (a1, a2, a3) = dispersion_pair_sum(
        &[methane(), n_butane()],
        &[0.6, 0.4],
        350.0,
        0.028_086_304_906_480_8,
    )
    .expect("a pair sum");
    via_eta(a1, -0.197_247_388_294_105, "A1");
    via_eta(a2, -0.023_030_744_398_675_0, "A2");
    via_eta(a3, -0.002_296_612_432_239_86, "A3");

    // A pure fluid's pair sum is its direct value: the weight is one and the cross
    // parameters are the pure ones.
    let methane = methane();
    let eta = 0.031_588_908_686_015_9;
    let x0 = 1.042_832_070_301_19;
    let (s1, _, _) = dispersion_pair_sum(&[methane], &[1.0], 300.0, eta).expect("a pair sum");
    via_eta(
        s1,
        a1_mie(
            eta,
            methane.lambda_r,
            methane.lambda_a,
            methane.epsik / 300.0,
            mie_prefactor(methane.lambda_r, methane.lambda_a),
            x0,
        ),
        "the one-component pair sum is the direct value",
    );
}

/// **The whole Helmholtz energy at both states**, which is the layer everything above was
/// built for. The probe prints `F_hc`, `F_disp` and `F`; the assembly is
/// `m_bar a_hs - m_minus_1 ln g_hs + a_1 + a_2 + a_3`.
#[test]
fn the_helmholtz_energy_is_neqsims() {
    let (eta, v_methane) = (0.031_588_908_686_015_9, 4.609_634_632_033_79e-4);
    // The volume NeqSim converged to, so the state is the one the probe describes.
    let pure = state(&[methane()], &[1.0], 300.0, v_methane).expect("a state");
    via_eta(pure.eta, eta, "eta");
    via_eta(pure.f_hc(), 0.131_541_289_151_816, "F_hc");
    via_eta(pure.f_disp(), -0.211_906_612_576_890, "F_disp");
    via_eta(pure.f(), -0.080_365_323_425_074_2, "F");

    let binary = state(
        &[methane(), n_butane()],
        &[0.6, 0.4],
        350.0,
        8.260_537_757_932_58e-4,
    )
    .expect("a state");
    via_eta(binary.eta, 0.028_086_304_906_480_8, "eta");
    via_eta(binary.g_hs, 1.028_095_271_747_73, "g_hs");
    via_eta(binary.f_hc(), 0.146_641_004_513_517, "F_hc");
    via_eta(binary.f_disp(), -0.298_374_800_324_796, "F_disp");
    // **And the segment number is what makes that different from the plain sum**: the
    // three terms add to -0.222574745379874, which is the pure-fluid assembly's answer and
    // a 25% error for this mixture.
    via_eta(
        binary.a1 + binary.a2 + binary.a3,
        -0.222_574_745_379_874,
        "a1 + a2 + a3",
    );
    via_eta(
        binary.f_disp() / (binary.a1 + binary.a2 + binary.a3),
        binary.m_bar,
        "m_bar",
    );
    via_eta(binary.f(), -0.151_733_795_811_279, "F");
}

/// The pressure the kernel predicts at NeqSim's own converged volume.
///
/// `P v/(RT)` must be the compressibility the probe reports on the same phase:
/// `0.924019412709511` for methane and `0.851583764555351` for the binary. That ties the
/// `eta` derivative - the hard-sphere term's in closed form, the other two by the central
/// differences NeqSim uses - to the *state* the oracle converged to, so a sign, a term or a
/// power wrong in the chain shows up as a pressure that is not the one the volume came from.
#[test]
fn the_pressure_at_neqsims_volume_is_neqsims_compressibility() {
    // **`1e-6`, because the derivative is noise-limited on both sides.** The port takes
    // NeqSim's own step - a relative `1e-5` in the packing fraction - and at that step the
    // difference of two nearly equal energies is in charge: the methane state lands within
    // `1e-8` of the probe's `Z` and the binary within `5.9e-7`, and *shrinking* the step
    // makes the binary worse rather than better. At a step ten times larger the same port
    // gives `0.851583766815854`, which is `2.7e-9` from the probe - so the model is right
    // and this tolerance is the arithmetic.
    let v_methane = 4.609_634_632_033_79e-4;
    let pure = pressure_over_rt(&[methane()], &[1.0], 300.0, v_methane).expect("a pressure");
    assert!(
        (pure * v_methane / 0.924_019_412_709_511 - 1.0).abs() < 1.0e-6,
        "methane Z = {}",
        pure * v_methane
    );

    let v_binary = 8.260_537_757_932_58e-4;
    let binary = pressure_over_rt(&[methane(), n_butane()], &[0.6, 0.4], 350.0, v_binary)
        .expect("a pressure");
    assert!(
        (binary * v_binary / 0.851_583_764_555_351 - 1.0).abs() < 1.0e-6,
        "binary Z = {}",
        binary * v_binary
    );
}
