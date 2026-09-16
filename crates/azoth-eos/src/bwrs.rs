//! The Benedict-Webb-Rubin-Starling (MBWR-32) equation of state.
//!
//! The 32-term modified Benedict-Webb-Rubin equation of Younglove and Ely
//! (J. Phys. Chem. Ref. Data 16, 577 (1987)), as NeqSim's `ComponentBWRS` carries it.
//! It is density-based rather than cubic: the pressure is a polynomial in the molar
//! density plus an exponential tail,
//!
//! ```text
//! P = sum_{i=0..8} BP[i] rho^(i+1) + exp(-gamma rho^2) sum_{i=0..5} BE[i] rho^(3+2i)
//! ```
//!
//! with `BP`, `BE` temperature polynomials built from 32 coefficients `a[0..31]`, and
//! `gamma = 1/rhoc^2` from the critical density. The 32 coefficients are calibrated in
//! native MBWR units - density in mol/L, pressure in MPa, `R = 8.3144621e-3` L.MPa/(mol.K)
//! - so this module works in those units throughout and the mixture layer converts to SI.
//!
//! The residual Helmholtz free energy is integrated from the pressure in closed form
//! rather than carried as NeqSim's derivative recurrence:
//!
//! ```text
//! F_pol = 1/(R T) sum_{i=1..8} BP[i]/i rho^i
//! F_exp = 1/(R T) sum_{i=0..5} BE[i] i!/(2 gamma^(i+1)) [1 - EL sum_{k=0..i} (gamma rho^2)^k/k!]
//! ```
//!
//! with `EL = exp(-gamma rho^2)`. The derivative set - three rho derivatives and the
//! temperature derivatives - follows from this one closed form by the chain rule, so each
//! is verified against automatic differentiation in the test rather than transcribed from
//! NeqSim's `getFexpdTdT` and its kin.
//!
//! NeqSim's `ComponentBWRS` uses the SI gas constant `R = 8.3144621` in `BP[0]` and divides
//! `calcPressure2` by 100, which makes its pressure ten times the physical value; its
//! Helmholtz `getF` compensates with a `1e3` factor and is correct. This module uses the
//! native `R = 8.3144621e-3` throughout and therefore reproduces NeqSim's Helmholtz while
//! producing the physical pressure.

/// The gas constant in the units the MBWR-32 coefficients are calibrated in, `L.MPa/(mol.K)`.
pub const R_MPA: f64 = 8.3144621e-3;

/// The 32 MBWR-32 coefficients and the critical density of one substance.
///
/// `a[0..31]` and `rhoc` come verbatim from NeqSim's `MBWR32param` table; `rhoc` is in
/// mol/L. Only methane and ethane have entries.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BwrsCoefficients {
    /// The fitted coefficients `a0..a31`.
    pub a: [f64; 32],
    /// The critical density, mol/L.
    pub rhoc: f64,
}

impl BwrsCoefficients {
    /// `gamma = 1/rhoc^2`.
    #[must_use]
    pub fn gamma(self) -> f64 {
        1.0 / (self.rhoc * self.rhoc)
    }
}

/// The polynomial coefficients `BP(T)[0..8]`, native units.
#[must_use]
pub fn bp(t: f64, a: &[f64; 32]) -> [f64; 9] {
    let sr = t.sqrt();
    let t2 = t * t;
    let mut b = [0.0; 9];
    b[0] = R_MPA * t;
    b[1] = a[0] * t + a[1] * sr + a[2] + a[3] / t + a[4] / t2;
    b[2] = a[5] * t + a[6] + a[7] / t + a[8] / t2;
    b[3] = a[9] * t + a[10] + a[11] / t;
    b[4] = a[12];
    b[5] = a[13] / t + a[14] / t2;
    b[6] = a[15] / t;
    b[7] = a[16] / t + a[17] / t2;
    b[8] = a[18] / t2;
    b
}

/// The exponential coefficients `BE(T)[0..5]`, native units.
#[must_use]
pub fn be(t: f64, a: &[f64; 32]) -> [f64; 6] {
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t2 * t2;
    [
        a[19] / t2 + a[20] / t3,
        a[21] / t2 + a[22] / t4,
        a[23] / t2 + a[24] / t3,
        a[25] / t2 + a[26] / t4,
        a[27] / t2 + a[28] / t3,
        a[29] / t2 + a[30] / t3 + a[31] / t4,
    ]
}

/// `d BP[i] / dT`, the first temperature derivative.
#[must_use]
pub fn bp_dt(t: f64, a: &[f64; 32]) -> [f64; 9] {
    let sr = t.sqrt();
    let t2 = t * t;
    let t3 = t2 * t;
    let mut b = [0.0; 9];
    b[0] = R_MPA;
    b[1] = a[0] + a[1] / (2.0 * sr) - a[3] / t2 - 2.0 * a[4] / t3;
    b[2] = a[5] - a[7] / t2 - 2.0 * a[8] / t3;
    b[3] = a[9] - a[11] / t2;
    b[4] = 0.0;
    b[5] = -a[13] / t2 - 2.0 * a[14] / t3;
    b[6] = -a[15] / t2;
    b[7] = -a[16] / t2 - 2.0 * a[17] / t3;
    b[8] = -2.0 * a[18] / t3;
    b
}

/// `d BE[i] / dT`, the first temperature derivative.
#[must_use]
pub fn be_dt(t: f64, a: &[f64; 32]) -> [f64; 6] {
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t2 * t2;
    let t5 = t4 * t;
    [
        -2.0 * a[19] / t3 - 3.0 * a[20] / t4,
        -2.0 * a[21] / t3 - 4.0 * a[22] / t5,
        -2.0 * a[23] / t3 - 3.0 * a[24] / t4,
        -2.0 * a[25] / t3 - 4.0 * a[26] / t5,
        -2.0 * a[27] / t3 - 3.0 * a[28] / t4,
        -2.0 * a[29] / t3 - 3.0 * a[30] / t4 - 4.0 * a[31] / t5,
    ]
}

/// `d^2 BP[i] / dT^2`, the second temperature derivative.
#[must_use]
pub fn bp_dt_dt(t: f64, a: &[f64; 32]) -> [f64; 9] {
    let t15 = t * t.sqrt(); // t^(3/2)
    let t3 = t * t * t;
    let t4 = t3 * t;
    let mut b = [0.0; 9];
    b[0] = 0.0;
    b[1] = -a[1] / (4.0 * t15) + 2.0 * a[3] / t3 + 6.0 * a[4] / t4;
    b[2] = 2.0 * a[7] / t3 + 6.0 * a[8] / t4;
    b[3] = 2.0 * a[11] / t3;
    b[4] = 0.0;
    b[5] = 2.0 * a[13] / t3 + 6.0 * a[14] / t4;
    b[6] = 2.0 * a[15] / t3;
    b[7] = 2.0 * a[16] / t3 + 6.0 * a[17] / t4;
    b[8] = 6.0 * a[18] / t4;
    b
}

/// `d^2 BE[i] / dT^2`, the second temperature derivative.
#[must_use]
pub fn be_dt_dt(t: f64, a: &[f64; 32]) -> [f64; 6] {
    let t2 = t * t;
    let t4 = t2 * t2;
    let t5 = t4 * t;
    let t6 = t5 * t;
    [
        6.0 * a[19] / t4 + 12.0 * a[20] / t5,
        6.0 * a[21] / t4 + 20.0 * a[22] / t6,
        6.0 * a[23] / t4 + 12.0 * a[24] / t5,
        6.0 * a[25] / t4 + 20.0 * a[26] / t6,
        6.0 * a[27] / t4 + 12.0 * a[28] / t5,
        6.0 * a[29] / t4 + 12.0 * a[30] / t5 + 20.0 * a[31] / t6,
    ]
}

/// The pressure in MPa, `sum BP[i] rho^(i+1) + EL sum BE[i] rho^(3+2i)`.
#[must_use]
pub fn pressure(rho: f64, b: &[f64; 9], e: &[f64; 6], gamma: f64) -> f64 {
    let mut p = 0.0;
    let mut rp = rho;
    for &bi in b {
        p += bi * rp;
        rp *= rho;
    }
    let el = (-gamma * rho * rho).exp();
    let mut tail = 0.0;
    let mut rp = rho * rho * rho;
    for &ei in e {
        tail += ei * rp;
        rp *= rho * rho;
    }
    p + el * tail
}

/// The factorial pre-computed for the six exponential terms.
const FACT: [f64; 6] = [1.0, 1.0, 2.0, 6.0, 24.0, 120.0];

/// The exponential integrals `T_i = i!/(2 gamma^(i+1)) [1 - EL sum_{k=0..i} (gamma rho^2)^k/k!]`.
///
/// The `rho`-dependent part of the exponential Helmholtz term, factored out so the value
/// and the temperature derivatives share one evaluation.
fn exp_integrals(rho: f64, gamma: f64) -> [f64; 6] {
    let g2 = gamma * rho * rho;
    let el = (-g2).exp();
    let mut out = [0.0; 6];
    let mut gpow = 1.0;
    for i in 0..6 {
        gpow *= gamma;
        let mut s = 0.0;
        let mut pk = 1.0;
        for &fact in FACT.iter().take(i + 1) {
            s += pk / fact;
            pk *= g2;
        }
        out[i] = FACT[i] / (2.0 * gpow) * (1.0 - el * s);
    }
    out
}

/// `G(rho) = sum BE[i] rho^(2i+1)` and its first two rho-derivatives.
///
/// The exponential Helmholtz term's rho-derivatives collapse onto this one polynomial,
/// since `dT_i/drho = EL rho^(2i+1)`.
fn g_sum(rho: f64, e: &[f64; 6]) -> (f64, f64, f64) {
    let mut g = 0.0;
    let mut gp = 0.0;
    let mut gpp = 0.0;
    for (i, &ei) in e.iter().enumerate() {
        let n = (2 * i + 1) as f64;
        g += ei * rho.powi(2 * i as i32 + 1);
        gp += ei * n * rho.powi(2 * i as i32);
        // The i = 0 term is `n (n - 1) = 0`, so it never reads rho^(-1).
        if i >= 1 {
            gpp += ei * n * (n - 1.0) * rho.powi(2 * i as i32 - 1);
        }
    }
    (g, gp, gpp)
}

/// The polynomial part of the residual Helmholtz, `F_pol`.
#[must_use]
pub fn helmholtz_pol(t: f64, rho: f64, b: &[f64; 9]) -> f64 {
    let mut sum = 0.0;
    let mut rp = rho;
    for (i, &bi) in b.iter().enumerate().skip(1) {
        sum += bi / (i as f64) * rp;
        rp *= rho;
    }
    sum / (R_MPA * t)
}

/// The exponential part of the residual Helmholtz, `F_exp`.
#[must_use]
pub fn helmholtz_exp(t: f64, rho: f64, e: &[f64; 6], gamma: f64) -> f64 {
    let integrals = exp_integrals(rho, gamma);
    integrals.iter().zip(e).map(|(ti, ei)| ti * ei).sum::<f64>() / (R_MPA * t)
}

/// The dimensionless residual Helmholtz `A^res/(R T)`, `F_pol + F_exp`.
#[must_use]
pub fn helmholtz(t: f64, rho: f64, b: &[f64; 9], e: &[f64; 6], gamma: f64) -> f64 {
    helmholtz_pol(t, rho, b) + helmholtz_exp(t, rho, e, gamma)
}

/// `dF/drho`, the first density derivative at constant temperature.
#[must_use]
pub fn d_helmholtz_drho(t: f64, rho: f64, b: &[f64; 9], e: &[f64; 6], gamma: f64) -> f64 {
    let mut pol = 0.0;
    let mut rp = 1.0;
    for &bi in b.iter().skip(1) {
        pol += bi * rp;
        rp *= rho;
    }
    let (g, _, _) = g_sum(rho, e);
    let el = (-gamma * rho * rho).exp();
    (pol + el * g) / (R_MPA * t)
}

/// `d^2 F/drho^2`, the second density derivative at constant temperature.
#[must_use]
pub fn d2_helmholtz_drho2(t: f64, rho: f64, b: &[f64; 9], e: &[f64; 6], gamma: f64) -> f64 {
    let mut pol = 0.0;
    let mut rp = 1.0;
    for (i, &bi) in b.iter().enumerate().skip(2) {
        pol += bi * (i as f64 - 1.0) * rp;
        rp *= rho;
    }
    let (g, gp, _) = g_sum(rho, e);
    let el = (-gamma * rho * rho).exp();
    (pol + el * (gp - 2.0 * gamma * rho * g)) / (R_MPA * t)
}

/// `d^3 F/drho^3`, the third density derivative at constant temperature.
#[must_use]
pub fn d3_helmholtz_drho3(t: f64, rho: f64, b: &[f64; 9], e: &[f64; 6], gamma: f64) -> f64 {
    let mut pol = 0.0;
    let mut rp = 1.0;
    for (i, &bi) in b.iter().enumerate().skip(3) {
        pol += bi * (i as f64 - 1.0) * (i as f64 - 2.0) * rp;
        rp *= rho;
    }
    let (g, gp, gpp) = g_sum(rho, e);
    let el = (-gamma * rho * rho).exp();
    (pol + el
        * ((-2.0 * gamma + 4.0 * gamma * gamma * rho * rho) * g - 4.0 * gamma * rho * gp + gpp))
        / (R_MPA * t)
}

/// `dF/dT`, the first temperature derivative at constant density.
#[must_use]
pub fn d_helmholtz_dt(
    t: f64,
    rho: f64,
    b: &[f64; 9],
    e: &[f64; 6],
    bt: &[f64; 9],
    et: &[f64; 6],
    gamma: f64,
) -> f64 {
    let mut pol = 0.0;
    let mut rp = rho;
    for i in 1..9 {
        pol += (bt[i] - b[i] / t) / (i as f64) * rp;
        rp *= rho;
    }
    let integrals = exp_integrals(rho, gamma);
    let mut exp = 0.0;
    for i in 0..6 {
        exp += (et[i] - e[i] / t) * integrals[i];
    }
    (pol + exp) / (R_MPA * t)
}

/// `d^2 F/dT^2`, the second temperature derivative at constant density.
#[must_use]
#[allow(clippy::too_many_arguments)] // the signature carries the full coefficient set
pub fn d2_helmholtz_dt2(
    t: f64,
    rho: f64,
    b: &[f64; 9],
    e: &[f64; 6],
    bt: &[f64; 9],
    et: &[f64; 6],
    btt: &[f64; 9],
    ett: &[f64; 6],
    gamma: f64,
) -> f64 {
    let t2 = t * t;
    let mut pol = 0.0;
    let mut rp = rho;
    for i in 1..9 {
        pol += (btt[i] - 2.0 * bt[i] / t + 2.0 * b[i] / t2) / (i as f64) * rp;
        rp *= rho;
    }
    let integrals = exp_integrals(rho, gamma);
    let mut exp = 0.0;
    for i in 0..6 {
        exp += (ett[i] - 2.0 * et[i] / t + 2.0 * e[i] / t2) * integrals[i];
    }
    (pol + exp) / (R_MPA * t)
}

/// `d^2 F/dT drho`, the mixed temperature-density derivative.
#[must_use]
pub fn d2_helmholtz_dtdrho(
    t: f64,
    rho: f64,
    b: &[f64; 9],
    e: &[f64; 6],
    bt: &[f64; 9],
    et: &[f64; 6],
    gamma: f64,
) -> f64 {
    let mut pol = 0.0;
    let mut rp = 1.0;
    for i in 1..9 {
        pol += (bt[i] - b[i] / t) * rp;
        rp *= rho;
    }
    let mut g_t = 0.0;
    let mut rp = rho;
    for i in 0..6 {
        g_t += (et[i] - e[i] / t) * rp;
        rp *= rho * rho;
    }
    let el = (-gamma * rho * rho).exp();
    (pol + el * g_t) / (R_MPA * t)
}
