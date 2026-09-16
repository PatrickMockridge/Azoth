//! The BWRS mixture: linear coefficient mixing and the fugacity, checked against NeqSim.

use azoth_eos::bwrs::{BwrsCoefficients, be, bp, d_helmholtz_drho, helmholtz};
use azoth_eos::bwrs_mixture::{d_f_dn, ln_fugacity};

fn methane() -> BwrsCoefficients {
    BwrsCoefficients {
        a: [
            0.000009898,
            0.021996082,
            -0.532278800,
            20.216579620,
            -2234.398926000,
            0.000010679,
            0.000145792,
            -0.926581666,
            291.536473200,
            0.000000231,
            0.000138721,
            0.004780467,
            0.000011761,
            -0.000198209,
            -0.025128877,
            0.000009748,
            -0.000000120,
            0.000041283,
            -0.000000721,
            508.173825500,
            -91989.031920000,
            -2.732264677,
            74990.243510000,
            0.001114060,
            1.083955159,
            -0.000044909,
            -1.380337847,
            -0.000000023,
            0.000037616,
            -2.375166954e-10,
            -0.000000012,
            0.000000676,
        ],
        rhoc: 10.150000000,
    }
}

fn ethane() -> BwrsCoefficients {
    BwrsCoefficients {
        a: [
            -0.018439486,
            1.051016206,
            -16.057820303,
            848.440275620,
            -42738.409106000,
            0.000765652,
            -0.483607241,
            85.195473835,
            -16607.434721000,
            -0.000037521,
            0.028616309,
            -2.868528597,
            0.000119069,
            -0.008531571,
            3.836506384,
            0.000024986,
            0.000005797,
            -0.007164832,
            0.000125778,
            22240.102466000,
            -1480051.23280000,
            50.498054887,
            1642883.75992000,
            0.213253871,
            37.791273422,
            -0.000011857,
            -31.630780767,
            -0.000004100,
            0.001487004,
            0.000000003,
            -0.000002167,
            0.000024000,
        ],
        rhoc: 6.860000000,
    }
}

fn coeffs(t: f64) -> (Vec<[f64; 9]>, Vec<[f64; 6]>, Vec<f64>) {
    let comps = [methane(), ethane()];
    (
        comps.iter().map(|c| bp(t, &c.a)).collect(),
        comps.iter().map(|c| be(t, &c.a)).collect(),
        comps.iter().map(|c| c.rhoc).collect(),
    )
}

fn assert_close(actual: f64, expected: f64) {
    let tol = (expected.abs() * 1e-9).max(1e-14);
    assert!(
        (actual - expected).abs() <= tol,
        "actual {actual} != expected {expected} (diff {})",
        (actual - expected).abs()
    );
}

/// `d(nF)/dn_i` reproduces NeqSim's `getdFdN` at the states its (buggy) density solver
/// happens to reach.
#[test]
fn d_f_dn_reproduces_neqsim_getd_fd_n() {
    // (x_methane, T, rho, dFdN[0], dFdN[1]) from NeqSim's SystemBWRSEos.
    let cases: [(f64, f64, f64, f64, f64); 4] = [
        (
            0.5,
            300.0,
            4.453182016992e-01,
            -1.183767265094e-01,
            -2.824323899085e-01,
        ),
        (
            0.5,
            400.0,
            3.092199181508e-01,
            -3.235056523386e-02,
            -7.891283041695e-02,
        ),
        (
            0.7,
            300.0,
            4.286580639483e-01,
            -8.269902124568e-02,
            -2.406651807189e-01,
        ),
        (
            0.5,
            300.0,
            2.691966518347e+00,
            -6.400962669140e-01,
            -1.587540549470e+00,
        ),
    ];
    for (xm, t, rho, d0, d1) in cases {
        let x = [xm, 1.0 - xm];
        let (bps, bes, rhocs) = coeffs(t);
        let d = d_f_dn(t, rho, &x, &bps, &bes, &rhocs);
        assert_close(d[0], d0);
        assert_close(d[1], d1);
    }
}

/// For a pure component `ln phi = F + rho dF/drho - ln Z`, so the numerical fugacity
/// reduces to the closed-form residual fugacity.
#[test]
fn pure_component_fugacity_reduces_to_closed_form() {
    for c in [methane(), ethane()] {
        let t = 300.0;
        let rho = 0.8;
        let b = bp(t, &c.a);
        let e = be(t, &c.a);
        let gamma = c.gamma();
        let f = helmholtz(t, rho, &b, &e, gamma);
        let drho = d_helmholtz_drho(t, rho, &b, &e, gamma);
        let ln_z = (1.0 + rho * drho).ln();
        let closed = f + rho * drho - ln_z;

        let phi = ln_fugacity(t, rho, &[1.0], &[c]);
        assert_close(phi[0], closed);
    }
}
