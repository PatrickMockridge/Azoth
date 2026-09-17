//! The Leachman reference equation of state for hydrogen.
//!
//! A multiparameter Helmholtz formulation in the reduced variables `tau = Tc/T` and
//! `delta = rho/rho_c`, with a Planck-Einstein ideal-gas part and a residual of
//! polynomial, exponential and Gaussian terms, as NeqSim's `thermo.util.leachman.Leachman`
//! carries it for the three hydrogen spin-isomers (normal, para, ortho). Every property
//! follows from the Helmholtz derivatives, so this module derives them in one place.
//!
//! Works in SI - density in mol/m³, pressure in Pa - with the equation's own
//! `R = 8.31451` J/(mol.K), the same rounded value Span-Wagner and Vega use.

// The coefficient tables are fitted data; one of them (0.6366) happens to sit near 2/pi.
#![allow(clippy::approx_constant)]

/// The gas constant, J/(mol.K).
pub const R: f64 = 8.31451;
/// The molar mass of hydrogen, kg/mol.
pub const MOLAR_MASS: f64 = 0.002_015_88;

/// The hydrogen spin-isomer the equation is parameterised for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HydrogenType {
    /// Normal hydrogen (the 3:1 ortho:para equilibrium mixture).
    Normal,
    /// Para-hydrogen.
    Para,
    /// Ortho-hydrogen.
    Ortho,
}

impl std::str::FromStr for HydrogenType {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim().to_lowercase().as_str() {
            "normal" => Ok(HydrogenType::Normal),
            "para" => Ok(HydrogenType::Para),
            "ortho" => Ok(HydrogenType::Ortho),
            other => Err(format!(
                "unknown hydrogen type `{other}`; expected `normal`, `para` or `ortho`"
            )),
        }
    }
}

impl HydrogenType {
    /// The critical temperature, K.
    #[must_use]
    pub const fn tc(self) -> f64 {
        match self {
            HydrogenType::Normal => 33.145,
            HydrogenType::Para => 32.938,
            HydrogenType::Ortho => 33.22,
        }
    }

    /// The critical molar density, mol/m³.
    #[must_use]
    pub const fn rhoc(self) -> f64 {
        match self {
            HydrogenType::Normal => 15_508.0,
            HydrogenType::Para => 15_538.0,
            HydrogenType::Ortho => 15_445.0,
        }
    }
}

// Ideal-part parameters: `alpha^0 = ln(delta) + 1.5 ln(tau) + a0[0] + a0[1] tau +
// sum_k a0[k] ln(1 - exp(b0[k] tau))`, the Planck-Einstein ideal gas.
const NORMAL_A0: [f64; 7] = [
    -1.457_985_647_5,
    1.888_076_782,
    1.616,
    -0.4117,
    -0.792,
    0.758,
    1.217,
];
const NORMAL_B0: [f64; 7] = [
    0.0,
    0.0,
    -16.020_515_914_9,
    -22.658_017_800_6,
    -60.009_051_138_9,
    -74.943_430_381_7,
    -206.939_206_516_8,
];
const PARA_A0: [f64; 9] = [
    -1.448_589_113_4,
    1.884_521_239,
    4.302_56,
    13.0289,
    -47.7365,
    50.0013,
    -18.6261,
    0.993_973,
    0.536_078,
];
const PARA_B0: [f64; 9] = [
    0.0,
    0.0,
    -15.149_675_147_2,
    -25.092_598_214_8,
    -29.473_556_378_7,
    -35.405_914_141_7,
    -40.724_998_482,
    -163.792_579_998_8,
    -309.217_317_384_2,
];
const ORTHO_A0: [f64; 6] = [
    -1.467_544_233_6,
    1.884_506_886_2,
    2.541_51,
    -2.3661,
    1.003_65,
    1.224_47,
];
const ORTHO_B0: [f64; 6] = [
    0.0,
    0.0,
    -25.767_609_873_6,
    -43.467_790_487_7,
    -66.044_551_475_0,
    -209.753_160_746_5,
];

// Residual term counts: 7 polynomial, 2 exponential, 5 Gaussian.
const N_POLY: usize = 7;
const N_EXP: usize = 2;
const N_GAUSS: usize = 5;
const N_TOTAL: usize = N_POLY + N_EXP + N_GAUSS; // 14

const NORMAL_N: [f64; N_TOTAL] = [
    -6.93643, 0.01, 2.1101, 4.52059, 0.732564, -1.34086, 0.130985, -0.777414, 0.351944, -0.0211716,
    0.0226312, 0.032187, -0.0231752, 0.0557346,
];
const NORMAL_T: [f64; N_TOTAL] = [
    0.6844, 1.0, 0.989, 0.489, 0.803, 1.1444, 1.409, 1.754, 1.311, 4.187, 5.646, 0.791, 7.249,
    2.986,
];
const NORMAL_D: [f64; N_TOTAL] = [
    1.0, 4.0, 1.0, 1.0, 2.0, 2.0, 3.0, 1.0, 3.0, 2.0, 1.0, 3.0, 1.0, 1.0,
];
const NORMAL_PHI: [f64; N_TOTAL] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.685, -0.489, -0.103, -2.506, -1.607,
];
const NORMAL_BETA: [f64; N_TOTAL] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.171, -0.2245, -0.1304, -0.2785, -0.3967,
];
const NORMAL_GAMMA: [f64; N_TOTAL] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7164, 1.3444, 1.4517, 0.7204, 1.5445,
];
const NORMAL_EPS: [f64; N_TOTAL] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.506, 0.156, 1.736, 0.67, 1.662,
];

const PARA_N: [f64; N_TOTAL] = [
    -7.33375, 0.01, 2.60375, 4.66279, 0.68239, -1.47078, 0.135801, -1.05327, 0.328239, -0.057783,
    0.044974, 0.070346, -0.040176, 0.11951,
];
const PARA_T: [f64; N_TOTAL] = [
    0.6855, 1.0, 1.0, 0.489, 0.774, 1.133, 1.386, 1.619, 1.162, 3.96, 5.276, 0.99, 6.791, 3.19,
];
const PARA_D: [f64; N_TOTAL] = [
    1.0, 4.0, 1.0, 1.0, 2.0, 2.0, 3.0, 1.0, 3.0, 2.0, 1.0, 3.0, 1.0, 1.0,
];
const PARA_PHI: [f64; N_TOTAL] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.7437, -0.5516, -0.0634, -2.1341, -1.777,
];
const PARA_BETA: [f64; N_TOTAL] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.194, -0.2019, -0.0301, -0.2383, -0.3253,
];
const PARA_GAMMA: [f64; N_TOTAL] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8048, 1.5248, 0.6648, 0.6832, 1.493,
];
const PARA_EPS: [f64; N_TOTAL] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5487, 0.1785, 1.28, 0.6319, 1.7104,
];

const ORTHO_N: [f64; N_TOTAL] = [
    -6.83148, 0.01, 2.11505, 4.38353, 0.211292, -1.00939, 0.142086, -0.87696, 0.804927, -0.710775,
    0.0639688, 0.0710858, -0.087654, 0.647088,
];
const ORTHO_T: [f64; N_TOTAL] = [
    0.7333, 1.0, 1.1372, 0.5136, 0.5638, 1.6248, 1.829, 2.404, 2.105, 4.1, 7.658, 1.259, 7.589,
    3.946,
];
const ORTHO_D: [f64; N_TOTAL] = [
    1.0, 4.0, 1.0, 1.0, 2.0, 2.0, 3.0, 1.0, 3.0, 2.0, 1.0, 3.0, 1.0, 1.0,
];
const ORTHO_PHI: [f64; N_TOTAL] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.169, -0.894, -0.04, -2.072, -1.306,
];
const ORTHO_BETA: [f64; N_TOTAL] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.4555, -0.4046, -0.0869, -0.4415, -0.5743,
];
const ORTHO_GAMMA: [f64; N_TOTAL] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5444, 0.6627, 0.763, 0.6587, 1.4327,
];
const ORTHO_EPS: [f64; N_TOTAL] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6366, 0.3876, 0.9437, 0.3976, 0.9626,
];

/// The residual coefficient arrays for a hydrogen type, as `(n, t, d, phi, beta, gamma, eps)`.
type ResidualParams = (
    &'static [f64; N_TOTAL],
    &'static [f64; N_TOTAL],
    &'static [f64; N_TOTAL],
    &'static [f64; N_TOTAL],
    &'static [f64; N_TOTAL],
    &'static [f64; N_TOTAL],
    &'static [f64; N_TOTAL],
);

fn residual_params(ht: HydrogenType) -> ResidualParams {
    match ht {
        HydrogenType::Normal => (
            &NORMAL_N,
            &NORMAL_T,
            &NORMAL_D,
            &NORMAL_PHI,
            &NORMAL_BETA,
            &NORMAL_GAMMA,
            &NORMAL_EPS,
        ),
        HydrogenType::Para => (
            &PARA_N,
            &PARA_T,
            &PARA_D,
            &PARA_PHI,
            &PARA_BETA,
            &PARA_GAMMA,
            &PARA_EPS,
        ),
        HydrogenType::Ortho => (
            &ORTHO_N,
            &ORTHO_T,
            &ORTHO_D,
            &ORTHO_PHI,
            &ORTHO_BETA,
            &ORTHO_GAMMA,
            &ORTHO_EPS,
        ),
    }
}

fn ideal_params(ht: HydrogenType) -> (&'static [f64], &'static [f64]) {
    match ht {
        HydrogenType::Normal => (&NORMAL_A0, &NORMAL_B0),
        HydrogenType::Para => (&PARA_A0, &PARA_B0),
        HydrogenType::Ortho => (&ORTHO_A0, &ORTHO_B0),
    }
}

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

/// The ideal-gas Helmholtz energy and its derivatives.
#[must_use]
pub fn ideal(delta: f64, tau: f64, ht: HydrogenType) -> IdealHelmholtz {
    let (a0, b0) = ideal_params(ht);
    let mut alpha = delta.ln() + 1.5 * tau.ln() + a0[0] + a0[1] * tau;
    let mut da_dt = 1.5 / tau + a0[1];
    let mut d2a_dt2 = -1.5 / (tau * tau);
    for k in 2..a0.len() {
        let e = (b0[k] * tau).exp();
        let denom = 1.0 - e;
        alpha += a0[k] * denom.ln();
        da_dt += -a0[k] * b0[k] * e / denom;
        d2a_dt2 += -a0[k] * b0[k] * b0[k] * e / (denom * denom);
    }
    IdealHelmholtz {
        alpha,
        alpha_tau: da_dt,
        alpha_tau_tau: d2a_dt2,
    }
}

/// The residual Helmholtz energy and its derivatives.
#[must_use]
pub fn residual(delta: f64, tau: f64, ht: HydrogenType) -> ResidualHelmholtz {
    let (n, t, d, phi, beta, gamma, eps) = residual_params(ht);
    let mut alpha = 0.0;
    let mut da_dd = 0.0;
    let mut d2a_dd2 = 0.0;
    let mut da_dt = 0.0;
    let mut d2a_dt2 = 0.0;
    let mut d2a_dd_dt = 0.0;

    // Polynomial terms.
    for i in 0..N_POLY {
        let term = n[i] * delta.powf(d[i]) * tau.powf(t[i]);
        alpha += term;
        da_dd += term * d[i] / delta;
        d2a_dd2 += term * d[i] * (d[i] - 1.0) / (delta * delta);
        da_dt += term * t[i] / tau;
        d2a_dt2 += term * t[i] * (t[i] - 1.0) / (tau * tau);
        d2a_dd_dt += term * d[i] * t[i] / (delta * tau);
    }

    // Exponential terms `n delta^d tau^t exp(-delta)`.
    for i in N_POLY..N_POLY + N_EXP {
        let term = n[i] * delta.powf(d[i]) * tau.powf(t[i]) * (-delta).exp();
        let u = d[i] / delta - 1.0;
        let u_prime = -d[i] / (delta * delta);
        alpha += term;
        da_dd += term * u;
        d2a_dd2 += term * (u * u + u_prime);
        da_dt += term * t[i] / tau;
        d2a_dt2 += term * t[i] * (t[i] - 1.0) / (tau * tau);
        d2a_dd_dt += term * u * t[i] / tau;
    }

    // Gaussian terms `n delta^d tau^t exp(phi (delta-eps)^2 + beta (tau-gamma)^2)`.
    for i in N_POLY + N_EXP..N_TOTAL {
        let dr = delta - eps[i];
        let tr = tau - gamma[i];
        let term =
            n[i] * delta.powf(d[i]) * tau.powf(t[i]) * (phi[i] * dr * dr + beta[i] * tr * tr).exp();
        let u = d[i] / delta + 2.0 * phi[i] * dr;
        let u_prime = -d[i] / (delta * delta) + 2.0 * phi[i];
        let v = t[i] / tau + 2.0 * beta[i] * tr;
        let v_prime = -t[i] / (tau * tau) + 2.0 * beta[i];
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

/// The molar thermodynamic properties of one hydrogen state, in SI.
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
pub fn properties(t: f64, rho: f64, ht: HydrogenType) -> Properties {
    let delta = rho / ht.rhoc();
    let tau = ht.tc() / t;
    let id = ideal(delta, tau, ht);
    let res = residual(delta, tau, ht);

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
pub fn solve_density(t: f64, p: f64, ht: HydrogenType) -> f64 {
    let mut rho = p / (R * t);
    for _ in 0..100 {
        let delta = rho / ht.rhoc();
        let tau = ht.tc() / t;
        let res = residual(delta, tau, ht);
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
