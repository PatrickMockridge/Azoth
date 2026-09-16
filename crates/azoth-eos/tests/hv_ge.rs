//! The Huron-Vidal GE model, checked against NeqSim 3.20.0's CLASSIC_HV.

use azoth_eos::hv_ge::hv_ln_gamma;

/// NeqSim's `SRKHuronVidal2` (CLASSIC_HV) for water/ethanol at T = 300 K, x = 0.5/0.5.
///
/// The reduced attraction and repulsion are NeqSim's SRK values for the two, the fitted
/// NRTL parameters its database carries (`HVDij` in Kelvin, `HValpha` symmetric), and the
/// pair marked HV. `lambda` is SRK's `ln(2)`.
#[test]
fn reproduces_neqsims_water_ethanol_gamma() {
    let x = [0.5, 0.5];
    let a = [0.015_684_626_494_084_11, 0.036_425_601_555_317_06];
    let b = [0.000_846_308_042_272_862_2, 0.002_417_275_020_583_037_4];
    let kij = [0.0, 0.0, 0.0, 0.0];
    let hv_gij = [0.0, -2612.51, 2207.03, 0.0];
    let hv_alpha = [0.0, 0.2245, 0.2245, 0.0];
    let hv_pairs = [false, true, true, false];

    let ln_gamma = hv_ln_gamma(
        &x,
        300.0,
        &a,
        &b,
        &kij,
        &hv_gij,
        &hv_alpha,
        &hv_pairs,
        std::f64::consts::LN_2,
    );

    assert!(
        (ln_gamma[0] - (-0.864_169_992_654_414_6)).abs() < 1e-12,
        "water: {}",
        ln_gamma[0]
    );
    assert!(
        (ln_gamma[1] - (-2.733_576_199_101_879)).abs() < 1e-12,
        "ethanol: {}",
        ln_gamma[1]
    );
}
