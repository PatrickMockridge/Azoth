//! The BWRS (MBWR-32) kernel, checked against NeqSim 3.20.0's `ComponentBWRS`.

// The Taylor-series reference construction below is index arithmetic by nature, so the
// range loops it is written with are the clear form, not a needless index.
#![allow(clippy::needless_range_loop)]

use azoth_eos::bwrs::{
    BwrsCoefficients, R_MPA, R_SI, be, be_dt, be_dt_dt, bp, bp_dt, bp_dt_dt, d_helmholtz_drho,
    d_helmholtz_dt, d2_helmholtz_drho2, d2_helmholtz_dt2, d2_helmholtz_dtdrho, d3_helmholtz_drho3,
    departure, helmholtz, pressure, solve_density,
};

/// The methane coefficients, verbatim from NeqSim's `MBWR32param` table.
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

/// The ethane coefficients, verbatim from NeqSim's `MBWR32param` table.
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

/// The model assembles the kernel, the mixture and the departure into the phase state
/// the spec's worked example names.
#[test]
fn bwrs_phase_model_matches_the_worked_example() {
    let coeffs = vec![methane(), ethane()];
    let r = azoth_eos::bwrs_phase(
        &coeffs,
        azoth_core::units::kelvins(300.0),
        azoth_core::units::pascals(1_000_000.0),
        &[0.5, 0.5],
    )
    .unwrap();
    assert_close(r.z_factor, 0.900_272_824_234_523_2);
    assert_close(r.ln_phi[0], -0.013_319_302_953_322_648);
    assert_close(r.ln_phi[1], -0.177_374_966_351_585_8);
    assert_close(r.h_res.value, -914.887_389_316_650_1);
    assert_close(r.s_res.value, -2.256_864_493_554_733_7);
    assert_close(r.cp_res.value, 7.100_809_443_311_604);
}

/// The databank lookup returns the same coefficients the oracle tests hardcode.
#[test]
fn databank_lookup_matches_the_verbatim_coefficients() {
    assert_eq!(
        azoth_eos::databank::bwrs_coefficients("methane"),
        Some(methane())
    );
    assert_eq!(
        azoth_eos::databank::bwrs_coefficients("ethane"),
        Some(ethane())
    );
    assert_eq!(azoth_eos::databank::bwrs_coefficients("water"), None);
}

fn assert_close(actual: f64, expected: f64) {
    let tol = (expected.abs() * 1e-11).max(1e-18);
    assert!(
        (actual - expected).abs() <= tol,
        "actual {actual} != expected {expected} (diff {})",
        (actual - expected).abs()
    );
}

/// `bp` and `be` reproduce NeqSim's `getBP`/`getBE` arrays.
///
/// NeqSim's `BP[0]` uses the SI gas constant, so it is a thousand times this module's
/// `BP[0] = R_MPA * T`; the other polynomial coefficients and all of `BE` are the same.
#[test]
fn bp_and_be_reproduce_neqsim() {
    let c = methane();
    let a = &c.a;
    // NeqSim's arrays at T = 300 K.
    let neqsim_bp = [
        2494.338630000,
        -0.1057641400855,
        0.003500180593333,
        0.000223955890000,
        0.000011761000000,
        -0.000000939906411111,
        0.000000032493333333,
        0.000000000058700,
        -0.000000000008011111,
    ];
    let neqsim_be = [
        0.002239374656667,
        -0.00002110044165679,
        0.00000005252493181481,
        -0.0000000006694009687654,
        0.00000000000113762962963,
        -0.000000000000003000062047654,
    ];

    let b = bp(300.0, a);
    let e = be(300.0, a);
    assert_close(b[0], neqsim_bp[0] / 1000.0);
    for i in 1..9 {
        assert_close(b[i], neqsim_bp[i]);
    }
    for i in 0..6 {
        assert_close(e[i], neqsim_be[i]);
    }
}

/// The dimensionless residual Helmholtz reproduces NeqSim's `getF`, which is correct
/// even though its `calcPressure2` is not.
#[test]
fn helmholtz_reproduces_neqsim_get_f() {
    // (name, T, rho, F) from NeqSim's SystemBWRSEos at the given temperature and pressure.
    let cases: [(&str, f64, f64, f64); 5] = [
        ("methane", 300.0, 4.078009498897e-01, -1.709816845295e-02),
        ("methane", 400.0, 3.020550316125e-01, -4.628221868127e-03),
        ("methane", 300.0, 2.180940796827e+00, -8.676404931204e-02),
        ("ethane", 300.0, 5.038801034558e-01, -2.064893946421e-01),
        ("ethane", 400.0, 3.171176529670e-01, -5.258985174915e-02),
    ];
    for (name, t, rho, expected) in cases {
        let c = match name {
            "methane" => methane(),
            _ => ethane(),
        };
        let b = bp(t, &c.a);
        let e = be(t, &c.a);
        let f = helmholtz(t, rho, &b, &e, c.gamma());
        assert_close(f, expected);
    }
}

/// The pressure is the physical value - the ideal-gas limit at low density - and is tied
/// to the Helmholtz by `P = rho R T (1 + rho dF/drho)`.
///
/// NeqSim's `calcPressure2` divides by 100 instead of 1000, so it reports ten times this;
/// its `getF` compensates and is correct, which is what the previous test checks.
#[test]
fn pressure_is_physical_and_consistent_with_helmholtz() {
    let c = methane();
    let b = bp(300.0, &c.a);
    let e = be(300.0, &c.a);
    let gamma = c.gamma();

    // At rho = 0.4078 mol/L the pressure must be about the ideal-gas value (~1.0 MPa for
    // 10 bar), not NeqSim's 10.17 MPa.
    let rho = 0.4078009498897;
    let p = pressure(rho, &b, &e, gamma);
    let ideal = R_MPA * 300.0 * rho;
    assert!(
        (p - ideal).abs() < 0.05 * ideal,
        "pressure {p} MPa far from ideal-gas {ideal} MPa"
    );

    // The thermodynamic identity ties the pressure to the Helmholtz derivative exactly.
    let drho = d_helmholtz_drho(300.0, rho, &b, &e, gamma);
    let from_helmholtz = rho * R_MPA * 300.0 * (1.0 + rho * drho);
    assert_close(p, from_helmholtz);
}

// ---------------------------------------------------------------------------
// Automatic differentiation: a truncated Taylor series, used as the reference
// construction for every derivative below.
// ---------------------------------------------------------------------------

/// A value and its first three derivatives, `c[k] = f^(k)/k!`.
#[derive(Clone, Copy)]
struct Taylor {
    c: [f64; 4],
}

impl Taylor {
    fn k(v: f64) -> Self {
        Taylor {
            c: [v, 0.0, 0.0, 0.0],
        }
    }
    fn var(v: f64) -> Self {
        Taylor {
            c: [v, 1.0, 0.0, 0.0],
        }
    }
    fn d1(self) -> f64 {
        self.c[1]
    }
    fn d2(self) -> f64 {
        2.0 * self.c[2]
    }
    fn d3(self) -> f64 {
        6.0 * self.c[3]
    }
}

impl std::ops::Add for Taylor {
    type Output = Taylor;
    fn add(self, o: Taylor) -> Taylor {
        let mut c = [0.0; 4];
        for i in 0..4 {
            c[i] = self.c[i] + o.c[i];
        }
        Taylor { c }
    }
}

impl std::ops::Sub for Taylor {
    type Output = Taylor;
    fn sub(self, o: Taylor) -> Taylor {
        let mut c = [0.0; 4];
        for i in 0..4 {
            c[i] = self.c[i] - o.c[i];
        }
        Taylor { c }
    }
}

impl std::ops::Mul for Taylor {
    type Output = Taylor;
    fn mul(self, o: Taylor) -> Taylor {
        let mut c = [0.0; 4];
        for i in 0..4 {
            for j in 0..=i {
                c[i] += self.c[j] * o.c[i - j];
            }
        }
        Taylor { c }
    }
}

impl std::ops::Div for Taylor {
    type Output = Taylor;
    fn div(self, o: Taylor) -> Taylor {
        // 1/o, then multiply: b_k = -sum_{j=1..k} o_j b_{k-j} / o_0.
        let mut b = [0.0; 4];
        b[0] = 1.0 / o.c[0];
        for k in 1..4 {
            let mut s = 0.0;
            for j in 1..=k {
                s += o.c[j] * b[k - j];
            }
            b[k] = -s / o.c[0];
        }
        self * Taylor { c: b }
    }
}

/// `exp(x)` for a Taylor series: `y' = x' y`, so `y_k = (1/k) sum_{j=1..k} j x_j y_{k-j}`.
fn texp(x: Taylor) -> Taylor {
    let mut b = [0.0; 4];
    b[0] = x.c[0].exp();
    for k in 1..4 {
        let mut s = 0.0;
        for j in 1..=k {
            s += (j as f64) * x.c[j] * b[k - j];
        }
        b[k] = s / (k as f64);
    }
    Taylor { c: b }
}

/// `sqrt(x)` for a Taylor series, from `y^2 = x` coefficient by coefficient.
fn tsqrt(x: Taylor) -> Taylor {
    let mut b = [0.0; 4];
    b[0] = x.c[0].sqrt();
    for k in 1..4 {
        let mut s = 0.0;
        for j in 1..k {
            s += b[j] * b[k - j];
        }
        b[k] = (x.c[k] - s) / (2.0 * b[0]);
    }
    Taylor { c: b }
}

/// `bp(t)` evaluated in Taylor arithmetic, for the temperature-derivative checks.
fn bp_taylor(t: Taylor, a: &[f64; 32]) -> [Taylor; 9] {
    let sr = tsqrt(t);
    let t2 = t * t;
    let mut b = [Taylor::k(0.0); 9];
    b[0] = t * Taylor::k(R_MPA);
    b[1] = t * Taylor::k(a[0])
        + sr * Taylor::k(a[1])
        + Taylor::k(a[2])
        + Taylor::k(a[3]) / t
        + Taylor::k(a[4]) / t2;
    b[2] = t * Taylor::k(a[5]) + Taylor::k(a[6]) + Taylor::k(a[7]) / t + Taylor::k(a[8]) / t2;
    b[3] = t * Taylor::k(a[9]) + Taylor::k(a[10]) + Taylor::k(a[11]) / t;
    b[4] = Taylor::k(a[12]);
    b[5] = Taylor::k(a[13]) / t + Taylor::k(a[14]) / t2;
    b[6] = Taylor::k(a[15]) / t;
    b[7] = Taylor::k(a[16]) / t + Taylor::k(a[17]) / t2;
    b[8] = Taylor::k(a[18]) / t2;
    b
}

/// `be(t)` evaluated in Taylor arithmetic.
fn be_taylor(t: Taylor, a: &[f64; 32]) -> [Taylor; 6] {
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t2 * t2;
    [
        Taylor::k(a[19]) / t2 + Taylor::k(a[20]) / t3,
        Taylor::k(a[21]) / t2 + Taylor::k(a[22]) / t4,
        Taylor::k(a[23]) / t2 + Taylor::k(a[24]) / t3,
        Taylor::k(a[25]) / t2 + Taylor::k(a[26]) / t4,
        Taylor::k(a[27]) / t2 + Taylor::k(a[28]) / t3,
        Taylor::k(a[29]) / t2 + Taylor::k(a[30]) / t3 + Taylor::k(a[31]) / t4,
    ]
}

/// The exponential integral `T_i(rho)` as a plain number, for the temperature-derivative
/// checks where `rho` is held constant.
fn integral(rho: f64, gamma: f64, i: usize) -> f64 {
    const FACT: [f64; 6] = [1.0, 1.0, 2.0, 6.0, 24.0, 120.0];
    let g2 = gamma * rho * rho;
    let el = (-g2).exp();
    let mut s = 0.0;
    let mut pk = 1.0;
    for k in 0..=i {
        s += pk / FACT[k];
        pk *= g2;
    }
    let mut gpow = 1.0;
    for _ in 0..=i {
        gpow *= gamma;
    }
    FACT[i] / (2.0 * gpow) * (1.0 - el * s)
}

/// `helmholtz` as a function of `rho`, in Taylor arithmetic.
fn helmholtz_rho(rho: Taylor, t: f64, b: &[f64; 9], e: &[f64; 6], gamma: f64) -> Taylor {
    let mut pol = Taylor::k(0.0);
    let mut rp = rho;
    for i in 1..9 {
        pol = pol + rp * Taylor::k(b[i] / (i as f64));
        rp = rp * rho;
    }
    let g2 = Taylor::k(gamma) * rho * rho;
    let el = texp(Taylor::k(0.0) - g2);
    let mut exp = Taylor::k(0.0);
    let mut gpow = 1.0;
    for i in 0..6 {
        gpow *= gamma;
        const FACT: [f64; 6] = [1.0, 1.0, 2.0, 6.0, 24.0, 120.0];
        let mut s = Taylor::k(0.0);
        let mut pk = Taylor::k(1.0);
        for k in 0..=i {
            s = s + pk * Taylor::k(1.0 / FACT[k]);
            pk = pk * g2;
        }
        exp = exp + Taylor::k(e[i] * FACT[i] / (2.0 * gpow)) * (Taylor::k(1.0) - el * s);
    }
    (pol + exp) * Taylor::k(1.0 / (R_MPA * t))
}

/// The rho-derivatives match the Taylor series of the Helmholtz.
#[test]
fn rho_derivatives_match_ad() {
    let c = ethane();
    let t = 320.0;
    let rho = 1.7;
    let b = bp(t, &c.a);
    let e = be(t, &c.a);
    let gamma = c.gamma();

    let series = helmholtz_rho(Taylor::var(rho), t, &b, &e, gamma);
    assert_close(d_helmholtz_drho(t, rho, &b, &e, gamma), series.d1());
    assert_close(d2_helmholtz_drho2(t, rho, &b, &e, gamma), series.d2());
    assert_close(d3_helmholtz_drho3(t, rho, &b, &e, gamma), series.d3());
}

/// `helmholtz` as a function of `t`, in Taylor arithmetic, so the temperature derivatives
/// include the temperature dependence of `bp` and `be`.
fn helmholtz_t(t: Taylor, rho: f64, a: &[f64; 32], gamma: f64) -> Taylor {
    let b = bp_taylor(t, a);
    let e = be_taylor(t, a);
    let mut pol = Taylor::k(0.0);
    let mut rp = rho;
    for i in 1..9 {
        pol = pol + b[i] * Taylor::k(rp / (i as f64));
        rp *= rho;
    }
    let mut exp = Taylor::k(0.0);
    for i in 0..6 {
        exp = exp + e[i] * Taylor::k(integral(rho, gamma, i));
    }
    (pol + exp) / (t * Taylor::k(R_MPA))
}

/// The temperature derivatives match the Taylor series of the Helmholtz.
#[test]
fn t_derivatives_match_ad() {
    let c = ethane();
    let t = 320.0;
    let rho = 1.7;
    let b = bp(t, &c.a);
    let e = be(t, &c.a);
    let bt = bp_dt(t, &c.a);
    let et = be_dt(t, &c.a);
    let btt = bp_dt_dt(t, &c.a);
    let ett = be_dt_dt(t, &c.a);
    let gamma = c.gamma();

    let series = helmholtz_t(Taylor::var(t), rho, &c.a, gamma);
    assert_close(d_helmholtz_dt(t, rho, &b, &e, &bt, &et, gamma), series.d1());
    assert_close(
        d2_helmholtz_dt2(t, rho, &b, &e, &bt, &et, &btt, &ett, gamma),
        series.d2(),
    );
}

/// The mixed derivative `d^2 F/dT drho` matches the temperature derivative of `dF/drho`.
#[test]
fn mixed_derivative_matches_ad() {
    let c = methane();
    let t = 350.0;
    let rho = 2.0;
    let b = bp(t, &c.a);
    let e = be(t, &c.a);
    let bt = bp_dt(t, &c.a);
    let et = be_dt(t, &c.a);
    let gamma = c.gamma();

    // dF/drho as a function of t, in Taylor arithmetic, then differentiated at t.
    let series = {
        let tv = Taylor::var(t);
        let b = bp_taylor(tv, &c.a);
        let e = be_taylor(tv, &c.a);
        let mut pol = Taylor::k(0.0);
        let mut rp = 1.0;
        for i in 1..9 {
            pol = pol + b[i] * Taylor::k(rp);
            rp *= rho;
        }
        let mut g = Taylor::k(0.0);
        let mut rp = rho;
        for i in 0..6 {
            g = g + e[i] * Taylor::k(rp);
            rp *= rho * rho;
        }
        let el = Taylor::k((-gamma * rho * rho).exp());
        (pol + el * g) / (tv * Taylor::k(R_MPA))
    };

    assert_close(
        d2_helmholtz_dtdrho(t, rho, &b, &e, &bt, &et, gamma),
        series.d1(),
    );
}

fn departure_at(t: f64, rho: f64, c: &BwrsCoefficients) -> azoth_eos::bwrs::BwrsDeparture {
    let b = bp(t, &c.a);
    let bt = bp_dt(t, &c.a);
    let btt = bp_dt_dt(t, &c.a);
    let e = be(t, &c.a);
    let et = be_dt(t, &c.a);
    let ett = be_dt_dt(t, &c.a);
    departure(t, rho, &b, &bt, &btt, &e, &et, &ett, c.gamma())
}

/// The energy departures reproduce NeqSim's `getAresTV`, `getHresTP`, `getSresTP` and
/// `getGresTP`.
///
/// The heat capacities `getCvres` and `getCpres` are *not* reproduced: they read NeqSim's
/// hand-written `getFexpdTdT` (and its cross-derivative kin), which is ~0.5-3% off
/// depending on the state, whereas the closed-form `F_TT` here is checked against forward
/// AD and a finite difference in the tests below.
#[test]
fn departure_reproduces_neqsim() {
    // (name, T, rho, Ares, Hres, Sres, Gres).
    let cases: [(&str, f64, f64, f64, f64, f64, f64); 3] = [
        (
            "methane",
            300.0,
            4.078009498897e-01,
            -4.264862207444e+01,
            -1.563684951278e+02,
            -3.802675737903e-01,
            -4.228822299074e+01,
        ),
        (
            "methane",
            400.0,
            3.020550316125e-01,
            -1.539247012517e+01,
            -9.137062615646e+01,
            -1.900316887898e-01,
            -1.535795064055e+01,
        ),
        (
            "ethane",
            300.0,
            5.038801034558e-01,
            -5.150544737412e+02,
            -1.875376932375e+03,
            -4.736018062536e+00,
            -4.545715136145e+02,
        ),
    ];
    for (name, t, rho, a, h, s, g) in cases {
        let c = match name {
            "methane" => methane(),
            _ => ethane(),
        };
        let d = departure_at(t, rho, &c);
        assert_close(d.a_res, a);
        assert_close(d.h_res, h);
        assert_close(d.s_res, s);
        assert_close(d.g_res, g);
        // The Gibbs-Helmholtz identity ties the four energy departures together.
        assert_close(d.h_res - t * d.s_res, d.g_res);
    }
}

/// `Cv` and `Cp` reproduce the second temperature derivative by finite difference, the
/// independent check that NeqSim's hand-written `getFexpdTdT` recurrence fails.
#[test]
fn cv_and_cp_match_finite_difference() {
    let c = ethane();
    let t = 320.0;
    let rho = 1.7;
    let gamma = c.gamma();
    let d = departure_at(t, rho, &c);

    let h = 1e-4;
    let ft = |temp: f64| {
        d_helmholtz_dt(
            temp,
            rho,
            &bp(temp, &c.a),
            &be(temp, &c.a),
            &bp_dt(temp, &c.a),
            &be_dt(temp, &c.a),
            gamma,
        )
    };
    let f_tt_fd = (ft(t + h) - ft(t - h)) / (2.0 * h);
    let f_t = ft(t);
    // Cv = -R (2 T F_T + T^2 F_TT), the residual isochoric heat capacity. The finite
    // difference carries a truncation error the `t^2` factor amplifies, so the tolerance
    // here is looser than the closed-form checks.
    assert!(
        (d.cv_res - (-R_SI * (2.0 * t * f_t + t * t * f_tt_fd))).abs() < 1e-5 * d.cv_res.abs(),
        "cv {} vs finite difference {}",
        d.cv_res,
        -R_SI * (2.0 * t * f_t + t * t * f_tt_fd)
    );

    // Cp - Cv = R (Z + rho T F_rhoT)^2 / (1 + 2 rho F_rho + rho^2 F_rhorho), the isobaric
    // relation; here it is checked against the closed departure result.
    let phi_rho = rho * d_helmholtz_drho(t, rho, &bp(t, &c.a), &be(t, &c.a), gamma);
    let phi_rho_rho = rho * rho * d2_helmholtz_drho2(t, rho, &bp(t, &c.a), &be(t, &c.a), gamma);
    let phi_rho_t = rho
        * t
        * d2_helmholtz_dtdrho(
            t,
            rho,
            &bp(t, &c.a),
            &be(t, &c.a),
            &bp_dt(t, &c.a),
            &be_dt(t, &c.a),
            gamma,
        );
    let z = 1.0 + phi_rho;
    let cp_minus_cv = R_SI * (z + phi_rho_t).powi(2) / (1.0 + 2.0 * phi_rho + phi_rho_rho);
    assert!(
        ((d.cp_res - d.cv_res) - (cp_minus_cv - R_SI)).abs() < 1e-6 * (cp_minus_cv - R_SI).abs(),
        "cp - cv {} vs closed {}",
        d.cp_res - d.cv_res,
        cp_minus_cv - R_SI
    );
}

/// The density solver recovers the molar density NeqSim's (correct) volume solver reaches
/// for the same target pressure.
#[test]
fn solve_density_reproduces_neqsim_molar_density() {
    // (name, T, P in MPa, rho) — the target pressures are 10 and 50 bara.
    let cases: [(&str, f64, f64, f64); 3] = [
        ("methane", 300.0, 1.0, 4.078009498897e-01),
        ("methane", 300.0, 5.0, 2.180940796827e+00),
        ("ethane", 300.0, 1.0, 5.038801034558e-01),
    ];
    for (name, t, p, rho_expected) in cases {
        let c = match name {
            "methane" => methane(),
            _ => ethane(),
        };
        let b = bp(t, &c.a);
        let e = be(t, &c.a);
        let rho = solve_density(t, p, &b, &e, c.gamma());
        assert_close(rho, rho_expected);
        // The solved density is a fixed point of the pressure equation.
        assert_close(pressure(rho, &b, &e, c.gamma()), p);
    }
}
