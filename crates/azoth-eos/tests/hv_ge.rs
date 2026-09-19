//! The Huron-Vidal GE model, checked against NeqSim 3.20.0's GE rules.
//!
//! Three checks, in increasing order of what they can catch: the fixed water/ethanol
//! matrices against CLASSIC_HV, the same numbers resolved from the databank by name, and
//! `CO2`/`water` against all three of CLASSIC_HV, Huron-Vidal and Wong-Sandler - which is
//! the one that can tell the two rules' temperature coefficients apart.
//!
//! **NeqSim's activity coefficients have to be read after evaluating them.**
//! `PhaseGE.getActivityCoefficient` returns a cached field that only
//! `getExcessGibbsEnergy` fills, so a probe that skips that call reads a plausible-looking
//! number that is not the model's answer - measured, `gamma = 1.000002` for water where
//! the model says 19.3. `validation/neqsim/GeGamma.java` is the driver that does it
//! correctly, and these values came from it.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::Cubic;
use azoth_eos::hv_ge::{hv_d_ln_gamma_dn, hv_ln_gamma};

/// The water/ethanol NRTL parameters NeqSim's database carries for CLASSIC_HV.
fn water_ethanol() -> ([f64; 2], [f64; 2]) {
    (
        [0.015_684_626_494_084_11, 0.036_425_601_555_317_06],
        [0.000_846_308_042_272_862_2, 0.002_417_275_020_583_037_4],
    )
}

/// NeqSim's `SRKHuronVidal2` (CLASSIC_HV) for water/ethanol at T = 300 K, x = 0.5/0.5.
///
/// The fitted `Dij` in Kelvin, its temperature coefficient `DijT`, and the symmetric
/// non-randomness, for the pair marked HV. `lambda` is SRK's `ln(2)`.
#[test]
fn reproduces_neqsims_water_ethanol_gamma() {
    let (a, b) = water_ethanol();
    let x = [0.5, 0.5];
    let kij = [0.0, 0.0, 0.0, 0.0];
    let hv_gij = [0.0, -2612.51, 2207.03, 0.0];
    let hv_gij_t = [0.0, 7.3, -4.6, 0.0];
    let hv_alpha = [0.0, 0.2245, 0.2245, 0.0];
    let hv_pairs = [false, true, true, false];

    let ln_gamma = hv_ln_gamma(
        &x,
        300.0,
        &a,
        &b,
        &kij,
        &hv_gij,
        &hv_gij_t,
        &hv_alpha,
        &hv_pairs,
        std::f64::consts::LN_2,
    );

    assert!(
        (ln_gamma[0] - 0.703_759_864_912_056_9).abs() < 1e-12,
        "water: {}",
        ln_gamma[0]
    );
    assert!(
        (ln_gamma[1] - 0.509_936_195_107_294_5).abs() < 1e-12,
        "ethanol: {}",
        ln_gamma[1]
    );
}

/// The composition derivative agrees with a finite difference of the coefficients.
///
/// `d ln gamma_i / dn_p` is the mole-number derivative NeqSim's `dlngammadn` returns;
/// perturbing `n_p` and renormalising is the direct numerical check of it.
#[test]
fn d_ln_gamma_dn_matches_finite_differences() {
    let (a, b) = water_ethanol();
    let x = [0.5, 0.5];
    let kij = [0.0, 0.0, 0.0, 0.0];
    let hv_gij = [0.0, -2612.51, 2207.03, 0.0];
    let hv_gij_t = [0.0, 7.3, -4.6, 0.0];
    let hv_alpha = [0.0, 0.2245, 0.2245, 0.0];
    let hv_pairs = [false, true, true, false];

    let d = hv_d_ln_gamma_dn(
        &x,
        300.0,
        &a,
        &b,
        &kij,
        &hv_gij,
        &hv_gij_t,
        &hv_alpha,
        &hv_pairs,
        std::f64::consts::LN_2,
    );

    let eps = 1e-6;
    for p in 0..2 {
        let mut x_plus = [0.0; 2];
        let mut x_minus = [0.0; 2];
        for k in 0..2 {
            x_plus[k] = (x[k] + if k == p { eps } else { 0.0 }) / (1.0 + eps);
            x_minus[k] = (x[k] - if k == p { eps } else { 0.0 }) / (1.0 - eps);
        }
        let ln_plus = hv_ln_gamma(
            &x_plus,
            300.0,
            &a,
            &b,
            &kij,
            &hv_gij,
            &hv_gij_t,
            &hv_alpha,
            &hv_pairs,
            std::f64::consts::LN_2,
        );
        let ln_minus = hv_ln_gamma(
            &x_minus,
            300.0,
            &a,
            &b,
            &kij,
            &hv_gij,
            &hv_gij_t,
            &hv_alpha,
            &hv_pairs,
            std::f64::consts::LN_2,
        );
        for i in 0..2 {
            let approx = (ln_plus[i] - ln_minus[i]) / (2.0 * eps);
            assert!(
                (approx - d[i * 2 + p]).abs() < 1e-8,
                "d ln gamma[{i}]/dn[{p}]: analytic {} vs finite {approx}",
                d[i * 2 + p]
            );
        }
    }
}

/// The same numbers, with every matrix resolved from the databank by name.
///
/// The test above types the water/ethanol parameters here, which is how they were
/// checked while `parse_kij` could not read the columns they come from. This one resolves
/// them - `HVTYPE` into the pair flags, `HVALPHA`/`HVGIJ`/`HVGJI` into the fitted
/// parameters, `HVGIJT`/`HVGJIT` into the temperature coefficient - and asserts the same
/// `ln gamma`, so a resolver that reads the wrong column or reverses a pair fails here
/// rather than quietly changing an answer somewhere downstream.
#[test]
fn the_databank_resolves_the_same_water_ethanol_parameters() {
    let params =
        azoth_eos::databank::huron_vidal_parameters(&["water", "ethanol"], Cubic::Pr, None)
            .expect("the pair resolves");
    let (a, b) = water_ethanol();
    let x = [0.5, 0.5];

    assert_eq!(params.hv_pairs, vec![false, true, true, false]);
    assert_eq!(params.hv_gij, vec![0.0, -2612.51, 2207.03, 0.0]);
    assert_eq!(params.hv_gij_t, vec![0.0, 7.3, -4.6, 0.0]);
    assert_eq!(params.hv_alpha, vec![0.0, 0.2245, 0.2245, 0.0]);

    let ln_gamma = hv_ln_gamma(
        &x,
        300.0,
        &a,
        &b,
        &params.kij,
        &params.hv_gij,
        &params.hv_gij_t,
        &params.hv_alpha,
        &params.hv_pairs,
        std::f64::consts::LN_2,
    );
    assert!(
        (ln_gamma[0] - 0.703_759_864_912_056_9).abs() < 1e-12,
        "water: {}",
        ln_gamma[0]
    );
    assert!(
        (ln_gamma[1] - 0.509_936_195_107_294_5).abs() < 1e-12,
        "ethanol: {}",
        ln_gamma[1]
    );
}

/// A pair the interaction table marks `Classic` is not a fitted pair, and the flag is
/// not a default that can be left off.
///
/// `HVTYPE` is what decides it - `methane`/`nitrogen` is `Classic` and `water`/`ethanol`
/// is `HV` - so a resolver that defaulted every pair to fitted, or every pair to
/// classic, fails one of the two tests. The name order is deliberately the one the
/// compiled table does not store (`nitrogen,methane` is the row), because the resolver
/// is what has to undo that.
#[test]
fn a_pair_the_table_marks_classic_is_not_fitted() {
    let params =
        azoth_eos::databank::huron_vidal_parameters(&["methane", "nitrogen"], Cubic::Pr, None)
            .expect("the pair resolves");
    assert_eq!(
        params.hv_pairs,
        vec![false, false, false, false],
        "methane/nitrogen is not an `HV` pair"
    );
    assert!(
        params.hv_gij.iter().all(|&v| v == 0.0),
        "and carries no fitted energy: {:?}",
        params.hv_gij
    );
}

/// `CO2`/`water`: the pair where the Huron-Vidal and Wong-Sandler rules differ in *which*
/// column supplies `DijT`, and the pair whose two columns are both non-zero.
///
/// `WSGIJT` is 0.96 and `HVGIJT` is -0.842 for it, so unlike `water`/`ethanol` - where
/// both are zero and the choice is unobservable - this state tells the two apart.
/// NeqSim's own GE evaluation, read through `validation/neqsim/GeGamma.java`, gives:
///
/// | rule | NeqSim | which is |
/// |---|---|---|
/// | 3, CLASSIC_HV - no `DijT` at all | `[1.6711808766525502, 1.1266767804843818]` | `hv_gij_t = 0` |
/// | 4, Huron-Vidal with `HVGIJT` | `[1.3117413880259408, 0.8943368257405854]` | `hv_gij_t = HVGIJT` |
/// | 5, Wong-Sandler with `WSGIJT` | `[2.9589427251536327, 2.2321445185372717]` | `hv_gij_t = WSGIJT` |
///
/// The three are asserted together because the point is not that one number matches - it
/// is that each number matches the *right* rule, so a `DijT` wired to the wrong column
/// lands on one of the other two rows and fails rather than looking plausible.
#[test]
fn the_two_rules_read_their_own_temperature_coefficient() {
    let names = ["CO2", "water"];
    let (base, _) =
        azoth_eos::databank::mixture_of(&names, Cubic::Srk, None).expect("the pair resolves");
    let hv = azoth_eos::databank::huron_vidal_parameters(&names, Cubic::Srk, None)
        .expect("the pair resolves");
    let ws = azoth_eos::databank::wong_sandler_parameters(&names, None).expect("the pair resolves");

    assert_eq!(hv.hv_gij_t, vec![0.0, -0.842_025_353, -0.512_331_604, 0.0]);
    assert_eq!(ws.hv_gij_t, vec![0.0, 0.96, 3.96, 0.0]);

    let mixture = base.with_cubic(azoth_eos::Cubic::Srk);
    let reduced = mixture
        .reduced_parameters(kelvins(350.0), pascals(500_000.0))
        .expect("a state");
    let x = [0.5, 0.5];
    let lambda = std::f64::consts::LN_2;
    let kij = [0.0, 0.0, 0.0, 0.0];

    let cases: [(&str, &[f64], &[f64]); 3] = [
        (
            "CLASSIC_HV, no DijT",
            &[1.671_180_876_652_550_2, 1.126_676_780_484_381_8],
            &[0.0, 0.0, 0.0, 0.0],
        ),
        (
            "Huron-Vidal, HVGIJT",
            &[1.311_741_388_025_940_8, 0.894_336_825_740_585_4],
            &hv.hv_gij_t,
        ),
        (
            "Wong-Sandler, WSGIJT",
            &[2.958_942_725_153_632_7, 2.232_144_518_537_271_7],
            &ws.hv_gij_t,
        ),
    ];
    for (label, expected, dij_t) in cases {
        let got = hv_ln_gamma(
            &x,
            350.0,
            &reduced.a,
            &reduced.b,
            &kij,
            &hv.hv_gij,
            dij_t,
            &hv.hv_alpha,
            &hv.hv_pairs,
            lambda,
        );
        for i in 0..2 {
            assert!(
                (got[i] - expected[i]).abs() < 1e-12,
                "{label}: ln gamma[{i}] is {} but NeqSim's is {}",
                got[i],
                expected[i]
            );
        }
    }
}
