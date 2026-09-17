//! The Vega reference equation of state for helium.
//!
//! A multiparameter Helmholtz formulation in the reduced variables `tau = Tc/T` and
//! `delta = rho/rho_c`, with a monatomic ideal-gas part and a residual of polynomial,
//! exponential and Gaussian terms, as NeqSim's `thermo.util.Vega.Vega` (the NIST helium
//! EOS, IR 8474) carries it. Every property follows from the Helmholtz derivatives, so
//! this module derives them in one place.
//!
//! NeqSim works this equation in mol/L and kPa; this module works in SI - density in
//! mol/m³, pressure in Pa - with the equation's own `R = 8.314472` J/(mol.K).
//!
//! NeqSim's `AlpharVega` carries two derivative defects this module does not reproduce:
//! its Gaussian tau-derivative drops the `(de/dtau)^2` term (a ~5e-9 error in `Cv`), and
//! its exponential delta-derivative flips the sign of `d^2 e/ddelta^2` (a ~2e-5 error in
//! `Cp`). The closed-form log-derivatives here are checked against a finite difference in
//! the test.

/// The gas constant, J/(mol.K) - the value the NIST helium EOS publishes.
pub const R: f64 = 8.314472;
/// The molar mass of helium, kg/mol.
pub const MOLAR_MASS: f64 = 0.004_002_602;
/// The critical temperature, K.
pub const T_CRIT: f64 = 5.1953;
/// The critical molar density, mol/m³ (17.3837 mol/L).
pub const RHO_CRIT: f64 = 17_383.7;

// Ideal-part parameters: `alpha^0 = a1 + a2 tau + ln(delta) + 1.5 ln(tau)`, the monatomic
// ideal gas with `Cv = 1.5 R`.
const A1: f64 = 0.173_348_793_283_576_4;
const A2: f64 = 0.467_452_220_155_081_5;

// Residual term counts.
const N_POLY: usize = 6;
const N_EXP: usize = 6;
const N_GAUSS: usize = 11;
const N_TOTAL: usize = N_POLY + N_EXP + N_GAUSS; // 23

const N: [f64; N_TOTAL] = [
    0.015_559_018,
    3.063_893_2,
    -4.242_084_4,
    0.054_418_088,
    -0.189_719_04,
    0.087_856_262,
    2.283_356_6,
    -0.533_315_95,
    -0.532_965_02,
    0.994_449_15,
    -0.300_788_96,
    -1.643_256_3,
    0.802_910_2,
    0.026_838_669,
    0.046_876_78,
    -0.148_327_66,
    0.030_162_11,
    -0.019_986_041,
    0.142_835_14,
    0.007_418_269,
    -0.229_897_93,
    0.792_248_29,
    -0.049_386_338,
];
const T: [f64; N_TOTAL] = [
    1.0, 0.425, 0.63, 0.69, 1.83, 0.575, 0.925, 1.585, 1.69, 1.51, 2.9, 0.8, 1.26, 3.51, 2.785,
    1.0, 4.22, 0.83, 1.575, 3.447, 0.73, 1.634, 6.13,
];
const D: [f64; N_TOTAL] = [
    4.0, 1.0, 1.0, 2.0, 2.0, 3.0, 1.0, 1.0, 3.0, 2.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 1.0, 3.0, 2.0,
    2.0, 3.0, 2.0, 2.0,
];
// Exponential terms (indices 6..12): `exp(-delta^l)`.
const L: [f64; N_TOTAL] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 1.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
];
// Gaussian terms (indices 12..23).
const ETA: [f64; N_TOTAL] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5497, 9.245, 4.76323, 6.3826,
    8.7023, 0.255, 0.3523, 0.1492, 0.05, 0.1668, 42.2358,
];
const BETA: [f64; N_TOTAL] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2471, 0.0983, 0.1556, 2.6782,
    2.7077, 0.6621, 0.1775, 0.4821, 0.3069, 0.1758, 1357.6577,
];
const GAMMA: [f64; N_TOTAL] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.15, 2.54505, 1.2513, 1.9416,
    0.5984, 2.2282, 1.606, 3.815, 1.61958, 0.6407, 1.076,
];
const EPS: [f64; N_TOTAL] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.596, 0.3423, 0.761, 0.9747,
    0.5868, 0.5627, 2.5346, 3.6763, 4.5245, 5.039, 0.959,
];

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

/// The ideal-gas Helmholtz energy and its derivatives: `a1 + a2 tau + ln(delta) + 1.5 ln(tau)`.
#[must_use]
pub fn ideal(delta: f64, tau: f64) -> IdealHelmholtz {
    IdealHelmholtz {
        alpha: A1 + A2 * tau + delta.ln() + 1.5 * tau.ln(),
        alpha_tau: A2 + 1.5 / tau,
        alpha_tau_tau: -1.5 / (tau * tau),
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
    for i in 0..N_POLY {
        let term = N[i] * delta.powf(D[i]) * tau.powf(T[i]);
        alpha += term;
        da_dd += term * D[i] / delta;
        d2a_dd2 += term * D[i] * (D[i] - 1.0) / (delta * delta);
        da_dt += term * T[i] / tau;
        d2a_dt2 += term * T[i] * (T[i] - 1.0) / (tau * tau);
        d2a_dd_dt += term * D[i] * T[i] / (delta * tau);
    }

    // Exponential terms `n delta^d tau^t exp(-delta^l)`.
    for i in N_POLY..N_POLY + N_EXP {
        let term = N[i] * delta.powf(D[i]) * tau.powf(T[i]) * (-delta.powf(L[i])).exp();
        let u = D[i] / delta - L[i] * delta.powf(L[i] - 1.0);
        let u_prime = -D[i] / (delta * delta) - L[i] * (L[i] - 1.0) * delta.powf(L[i] - 2.0);
        alpha += term;
        da_dd += term * u;
        d2a_dd2 += term * (u * u + u_prime);
        da_dt += term * T[i] / tau;
        d2a_dt2 += term * T[i] * (T[i] - 1.0) / (tau * tau);
        d2a_dd_dt += term * u * T[i] / tau;
    }

    // Gaussian terms `n delta^d tau^t exp(-eta (delta-eps)^2 - beta (tau-gamma)^2)`.
    for i in N_POLY + N_EXP..N_TOTAL {
        let dr = delta - EPS[i];
        let tr = tau - GAMMA[i];
        let term = N[i]
            * delta.powf(D[i])
            * tau.powf(T[i])
            * (-ETA[i] * dr * dr - BETA[i] * tr * tr).exp();
        let u = D[i] / delta - 2.0 * ETA[i] * dr;
        let u_prime = -D[i] / (delta * delta) - 2.0 * ETA[i];
        let v = T[i] / tau - 2.0 * BETA[i] * tr;
        let v_prime = -T[i] / (tau * tau) - 2.0 * BETA[i];
        alpha += term;
        da_dd += term * u;
        d2a_dd2 += term * (u * u + u_prime);
        da_dt += term * v;
        d2a_dt2 += term * (v * v + v_prime);
        d2a_dd_dt += term * u * v;
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

/// The molar thermodynamic properties of one helium state, in SI.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Properties {
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
    /// Gibbs energy, J/mol.
    pub g: f64,
}

/// The molar properties at a temperature and molar density.
#[must_use]
pub fn properties(t: f64, rho: f64) -> Properties {
    let delta = rho / RHO_CRIT;
    let tau = T_CRIT / t;
    let id = ideal(delta, tau);
    let res = residual(delta, tau);

    let at = id.alpha_tau + res.alpha_tau;
    let numer = 1.0 + delta * res.alpha_delta - delta * tau * res.alpha_delta_tau;
    let denom = 1.0 + 2.0 * delta * res.alpha_delta + delta * delta * res.alpha_delta_delta;

    let z = 1.0 + delta * res.alpha_delta;
    let u = R * t * tau * at;
    let h = R * t * (1.0 + delta * res.alpha_delta + tau * at);
    let s = R * (tau * at - id.alpha - res.alpha);
    let cv = -R * tau * tau * (id.alpha_tau_tau + res.alpha_tau_tau);
    let cp = cv + R * numer * numer / denom;
    let g = R * t * (1.0 + delta * res.alpha_delta + id.alpha + res.alpha);

    Properties {
        z,
        u,
        h,
        s,
        cv,
        cp,
        g,
    }
}

/// Solve `P(rho) = p` for the molar density by Newton's method from the ideal-gas guess.
/// `p` in Pa, `rho` in mol/m³.
#[must_use]
pub fn solve_density(t: f64, p: f64) -> f64 {
    let mut rho = p / (R * t);
    for _ in 0..100 {
        let delta = rho / RHO_CRIT;
        let tau = T_CRIT / t;
        let res = residual(delta, tau);
        let p_calc = rho * R * t * (1.0 + delta * res.alpha_delta);
        let dpdrho =
            R * t * (1.0 + 2.0 * delta * res.alpha_delta + delta * delta * res.alpha_delta_delta);
        let step = (p_calc - p) / dpdrho;
        rho -= step;
        if step.abs() < 1e-12 * rho.abs() {
            break;
        }
    }
    rho
}
