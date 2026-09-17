//! The ammonia reference equation of state (Ammonia2023).
//!
//! A multiparameter Helmholtz formulation: the reduced Helmholtz energy
//! `alpha = alpha^0 + alpha^r` in the reduced variables `tau = Tc/T` and
//! `delta = rho/rho_c`, with an ideal-gas part of Planck-Einstein terms and a residual
//! part of polynomial, exponential, Gaussian and Gao-B terms, as NeqSim's
//! `thermo.util.referenceequations.Ammonia2023` carries it. Every property follows from
//! the Helmholtz derivatives, so this module derives them in one place and the property
//! set is a rearrangement of them.
//!
//! Unlike [`crate::bwrs`], this works entirely in SI - density in mol/m³, pressure in Pa,
//! `R = 8.31446261815324` J/(mol.K) - because that is the convention the reference
//! equation is published in.

/// The gas constant, J/(mol.K) - the CODATA value the reference equation uses.
pub const R: f64 = 8.314_462_618_153_24;
/// The molar mass of ammonia, kg/mol.
pub const MOLAR_MASS: f64 = 0.017_030_52;
/// The critical temperature, K.
pub const T_CRIT: f64 = 405.56;
/// The critical molar density, mol/m³.
pub const RHO_CRIT: f64 = 13_696.0;

// Ideal-part parameters.
const A1: f64 = -6.594_060_939_438_86;
const A2: f64 = 5.601_011_519_879;
const C0: f64 = 3.0;
const IDEAL_N: [f64; 3] = [2.224, 3.148, 0.9579];
const IDEAL_T: [f64; 3] = [
    4.058_585_659_335_240_5,
    9.776_605_187_888_352,
    17.829_667_620_080_876,
];

// Residual-part parameters. `N`, `T` and `D` hold all 18 terms; indices 0..4 are the
// polynomial terms, 5..7 the exponential, 8..17 the Gaussian.
const N: [f64; 18] = [
    0.006_132_232,
    1.739_586_6,
    -2.226_179_2,
    -0.301_275_53,
    0.089_670_23,
    -0.076_387_037,
    -0.840_639_63,
    -0.270_263_27,
    6.212_578,
    -5.784_435_7,
    2.481_754_2,
    -2.373_916_8,
    0.014_936_97,
    -3.774_926_4,
    0.000_625_434_8,
    -0.000_017_359,
    -0.134_620_33,
    0.077_490_728_39,
];
const T: [f64; 18] = [
    1.0, 0.382, 1.0, 1.0, 0.677, 2.915, 3.51, 1.063, 0.655, 1.3, 3.1, 1.4395, 1.623, 0.643, 1.13,
    4.5, 1.0, 4.0,
];
const D: [f64; 18] = [
    4.0, 1.0, 1.0, 2.0, 3.0, 3.0, 2.0, 3.0, 1.0, 1.0, 1.0, 2.0, 2.0, 1.0, 3.0, 3.0, 1.0, 1.0,
];
// Exponential terms (indices 5..7 of the three arrays above).
const L: [f64; 3] = [2.0, 2.0, 1.0];
const G: [f64; 3] = [1.0, 1.0, 1.0];
// Gaussian terms (indices 8..17).
const ETA: [f64; 10] = [
    0.42776, 0.6424, 0.8175, 0.7995, 0.91, 0.3574, 1.21, 4.14, 22.56, 22.68,
];
const BETA: [f64; 10] = [
    1.708, 1.4865, 2.0915, 2.43, 0.488, 1.1, 0.85, 1.14, 945.64, 993.85,
];
const GAMMA: [f64; 10] = [
    1.036, 1.2777, 1.083, 1.2906, 0.928, 0.934, 0.919, 1.852, 1.05897, 1.05277,
];
const EPS: [f64; 10] = [
    -0.0726, -0.1274, 0.7527, 0.57, 2.2, -0.243, 2.96, 3.02, 0.9574, 0.9576,
];
// Gao-B terms.
const GAOB_N: [f64; 2] = [-1.690_985_8, 0.937_390_74];
const GAOB_T: [f64; 2] = [4.3315, 4.015];
const GAOB_D: [f64; 2] = [1.0, 1.0];
const GAOB_ETA: [f64; 2] = [-2.8452, -2.8342];
const GAOB_BETA: [f64; 2] = [0.3696, 0.2962];
const GAOB_GAMMA: [f64; 2] = [1.108, 1.313];
const GAOB_EPS: [f64; 2] = [0.4478, 0.44689];
const GAOB_B: [f64; 2] = [1.244, 0.6826];

// Saturated-liquid density ancillary, used only to seed the dense root.
const SAT_LIQ_N: [f64; 6] = [2.447, 5.8341, -25.944, 53.383, -54.411, 22.771];
const SAT_LIQ_T: [f64; 6] = [0.384, 1.65, 2.2, 2.75, 3.35, 4.0];

/// The ideal-gas Helmholtz energy and its `tau` derivatives.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IdealHelmholtz {
    /// `alpha^0`.
    pub alpha: f64,
    /// `d alpha^0 / d tau`.
    pub alpha_tau: f64,
    /// `d^2 alpha^0 / d tau^2`.
    pub alpha_tau_tau: f64,
}

/// The residual Helmholtz energy and its derivatives.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResidualHelmholtz {
    /// `alpha^r`.
    pub alpha: f64,
    /// `d alpha^r / d delta`.
    pub alpha_delta: f64,
    /// `d^2 alpha^r / d delta^2`.
    pub alpha_delta_delta: f64,
    /// `d alpha^r / d tau`.
    pub alpha_tau: f64,
    /// `d^2 alpha^r / d tau^2`.
    pub alpha_tau_tau: f64,
    /// `d^2 alpha^r / d delta d tau`.
    pub alpha_delta_tau: f64,
}

/// The ideal-gas Helmholtz energy and its derivatives,
/// `alpha^0 = ln(delta) + A1 + A2 tau + C0 ln(tau) + sum N_i ln(1 - exp(-T_i tau))`.
#[must_use]
pub fn ideal(delta: f64, tau: f64) -> IdealHelmholtz {
    let mut alpha = delta.ln() + A1 + A2 * tau + C0 * tau.ln();
    let mut dalpha_dtau = A2 + C0 / tau;
    let mut d2alpha_dtau2 = -C0 / (tau * tau);
    for i in 0..3 {
        let exp_t = (-IDEAL_T[i] * tau).exp();
        let denom = 1.0 - exp_t;
        alpha += IDEAL_N[i] * denom.ln();
        dalpha_dtau += IDEAL_N[i] * IDEAL_T[i] * exp_t / denom;
        d2alpha_dtau2 -= IDEAL_N[i] * IDEAL_T[i] * IDEAL_T[i] * exp_t / (denom * denom);
    }
    IdealHelmholtz {
        alpha,
        alpha_tau: dalpha_dtau,
        alpha_tau_tau: d2alpha_dtau2,
    }
}

/// The residual Helmholtz energy and its derivatives.
#[must_use]
pub fn residual(delta: f64, tau: f64) -> ResidualHelmholtz {
    let mut alpha = 0.0;
    let mut da_dd = 0.0;
    let mut d2a_dd2 = 0.0;
    let mut da_dt = 0.0;
    let mut d2a_dt2 = 0.0;
    let mut d2a_dd_dt = 0.0;

    // Polynomial terms.
    for i in 0..5 {
        let del_pow = delta.powf(D[i]);
        let tau_pow = tau.powf(T[i]);
        let term = N[i] * del_pow * tau_pow;
        alpha += term;
        da_dd += term * D[i] / delta;
        d2a_dd2 += term * D[i] * (D[i] - 1.0) / (delta * delta);
        da_dt += term * T[i] / tau;
        d2a_dt2 += term * T[i] * (T[i] - 1.0) / (tau * tau);
        d2a_dd_dt += term * D[i] * T[i] / (delta * tau);
    }

    // Exponential terms.
    for i in 5..8 {
        let li = i - 5;
        let del_pow = delta.powf(D[i]);
        let tau_pow = tau.powf(T[i]);
        let term = N[i] * del_pow * tau_pow * (-G[li] * delta.powf(L[li])).exp();
        let b = D[i] / delta - G[li] * L[li] * delta.powf(L[li] - 1.0);
        alpha += term;
        da_dd += term * b;
        d2a_dd2 += term
            * (b * b
                - D[i] / (delta * delta)
                - G[li] * L[li] * (L[li] - 1.0) * delta.powf(L[li] - 2.0));
        da_dt += term * T[i] / tau;
        d2a_dt2 += term * T[i] * (T[i] - 1.0) / (tau * tau);
        d2a_dd_dt += term * b * T[i] / tau;
    }

    // Gaussian terms.
    for i in 8..18 {
        let gi = i - 8;
        let del_pow = delta.powf(D[i]);
        let tau_pow = tau.powf(T[i]);
        let a = delta - EPS[gi];
        let b = tau - GAMMA[gi];
        let term = N[i] * del_pow * tau_pow * (-ETA[gi] * a * a - BETA[gi] * b * b).exp();
        let b_delta = D[i] / delta - 2.0 * ETA[gi] * a;
        let b_tau = T[i] / tau - 2.0 * BETA[gi] * b;
        alpha += term;
        da_dd += term * b_delta;
        d2a_dd2 += term * (b_delta * b_delta - D[i] / (delta * delta) - 2.0 * ETA[gi]);
        da_dt += term * b_tau;
        d2a_dt2 += term * (b_tau * b_tau - T[i] / (tau * tau) - 2.0 * BETA[gi]);
        d2a_dd_dt += term * b_delta * b_tau;
    }

    // Gao-B terms.
    for i in 0..2 {
        let del_pow = delta.powf(GAOB_D[i]);
        let tau_pow = tau.powf(GAOB_T[i]);
        let a = delta - GAOB_EPS[i];
        let y = GAOB_BETA[i] * (tau - GAOB_GAMMA[i]).powi(2) + GAOB_B[i];
        let term = GAOB_N[i] * del_pow * tau_pow * (GAOB_ETA[i] * a * a + 1.0 / y).exp();
        let b_delta = GAOB_D[i] / delta + 2.0 * GAOB_ETA[i] * a;
        let b_t = GAOB_T[i] / tau - (2.0 * GAOB_BETA[i] * (tau - GAOB_GAMMA[i])) / (y * y);
        alpha += term;
        da_dd += term * b_delta;
        d2a_dd2 += term * (b_delta * b_delta - GAOB_D[i] / (delta * delta) + 2.0 * GAOB_ETA[i]);
        da_dt += term * b_t;
        let tau_offset = tau - GAOB_GAMMA[i];
        d2a_dt2 += term
            * (b_t * b_t - GAOB_T[i] / (tau * tau) - 2.0 * GAOB_BETA[i] / (y * y)
                + 8.0 * GAOB_BETA[i] * GAOB_BETA[i] * tau_offset * tau_offset / (y * y * y));
        d2a_dd_dt += term * b_delta * b_t;
    }

    ResidualHelmholtz {
        alpha,
        alpha_delta: da_dd,
        alpha_delta_delta: d2a_dd2,
        alpha_tau: da_dt,
        alpha_tau_tau: d2a_dt2,
        alpha_delta_tau: d2a_dd_dt,
    }
}

/// The molar thermodynamic properties of one ammonia state, in SI.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Properties {
    /// Pressure, Pa.
    pub pressure: f64,
    /// Compressibility factor.
    pub z: f64,
    /// Internal energy, J/mol.
    pub u: f64,
    /// Enthalpy, J/mol.
    pub h: f64,
    /// Entropy, J/(mol.K).
    pub s: f64,
    /// Isochoric heat capacity, J/(mol.K).
    pub cv: f64,
    /// Isobaric heat capacity, J/(mol.K).
    pub cp: f64,
    /// Sound speed, m/s.
    pub sound: f64,
    /// Gibbs energy, J/mol.
    pub g: f64,
    /// Joule-Thomson coefficient, K/Pa.
    pub jt: f64,
    /// Isothermal compressibility, 1/Pa.
    pub kappa: f64,
}

/// The molar properties at a temperature and molar density, from the Helmholtz derivatives.
#[must_use]
pub fn properties(t: f64, rho: f64) -> Properties {
    let delta = rho / RHO_CRIT;
    let tau = T_CRIT / t;
    let id = ideal(delta, tau);
    let res = residual(delta, tau);

    let cv = -R * tau * tau * (id.alpha_tau_tau + res.alpha_tau_tau);
    let numer = 1.0 + delta * res.alpha_delta - delta * tau * res.alpha_delta_tau;
    let denom = 1.0 + 2.0 * delta * res.alpha_delta + delta * delta * res.alpha_delta_delta;
    let cp = cv + R * numer * numer / denom;

    let u = R * t * tau * (id.alpha_tau + res.alpha_tau);
    let h = R * t * (1.0 + tau * (id.alpha_tau + res.alpha_tau) + delta * res.alpha_delta);
    let s = R * (tau * (id.alpha_tau + res.alpha_tau) - id.alpha - res.alpha);
    // `g = h - T s`, so `g/(R T) = 1 + alpha + delta alpha^r_delta`. NeqSim's Ammonia2023
    // subtracts an extra `tau alpha_tau`, which is `u/(R T)`, so its `g` is the correct
    // value less the internal energy - reproduced here as the physical Gibbs energy.
    let g = R * t * (1.0 + id.alpha + res.alpha + delta * res.alpha_delta);

    let dpdrho = R * t * denom;
    let dpdt = rho * R * numer;
    let kappa = 1.0 / (rho * dpdrho);
    let dv_dt_p = (dpdt / dpdrho) / (rho * rho);
    let jt = (t * dv_dt_p - 1.0 / rho) / cp;

    let sound_sq = R * t / MOLAR_MASS
        * (denom - numer * numer / (tau * tau * (id.alpha_tau_tau + res.alpha_tau_tau)));
    let sound = sound_sq.max(0.0).sqrt();

    let pressure = rho * R * t * (1.0 + delta * res.alpha_delta);
    let z = 1.0 + delta * res.alpha_delta;

    Properties {
        pressure,
        z,
        u,
        h,
        s,
        cv,
        cp,
        sound,
        g,
        jt,
        kappa,
    }
}

/// The saturated-liquid density guess, mol/m³, used only to seed the dense root.
fn saturated_liquid_density_guess(t: f64) -> f64 {
    if t >= T_CRIT {
        return RHO_CRIT;
    }
    let theta = 1.0 - t / T_CRIT;
    let mut reduced = 1.0;
    for i in 0..6 {
        reduced += SAT_LIQ_N[i] * theta.powf(SAT_LIQ_T[i]);
    }
    RHO_CRIT * reduced
}

/// Solve `P(rho) = p` for the molar density by Newton's method from the ideal-gas guess
/// (or the saturated-liquid guess for the dense root). `p` in Pa, `rho` in mol/m³.
#[must_use]
pub fn solve_density(t: f64, p: f64, liquid: bool) -> f64 {
    let mut rho = if liquid {
        saturated_liquid_density_guess(t)
    } else {
        p / (R * t)
    };
    for _ in 0..100 {
        let delta = rho / RHO_CRIT;
        let tau = T_CRIT / t;
        let res = residual(delta, tau);
        let p_calc = rho * R * t * (1.0 + delta * res.alpha_delta);
        let dpdrho =
            R * t * (1.0 + 2.0 * delta * res.alpha_delta + delta * delta * res.alpha_delta_delta);
        let step = (p_calc - p) / dpdrho;
        rho -= step;
        if step.abs() < 1e-8 * rho.abs() {
            break;
        }
    }
    rho
}
