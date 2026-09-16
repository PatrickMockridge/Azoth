//! The Huron-Vidal GE model, checked against NeqSim 3.20.0's CLASSIC_HV.

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
