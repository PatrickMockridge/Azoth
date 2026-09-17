//! The Span-Wagner reference equation of state for carbon dioxide.
//!
//! A multiparameter Helmholtz formulation in the reduced variables `tau = Tc/T` and
//! `delta = rho/rho_c`, with an ideal-gas part of Planck-Einstein terms and a residual
//! part of power-exponential and Gaussian terms, as NeqSim's
//! `thermo.util.spanwagner.NeqSimSpanWagner` carries it. Every property follows from the
//! Helmholtz derivatives, so this module derives them in one place.
//!
//! Works entirely in SI - density in mol/m³, pressure in Pa - with the equation's own
//! `R = 8.31451` J/(mol.K), which is the value Span and Wagner publish rather than the
//! full CODATA constant.

/// The gas constant, J/(mol.K) - Span-Wagner's published value.
pub const R: f64 = 8.31451;
/// The molar mass of carbon dioxide, kg/mol.
pub const MOLAR_MASS: f64 = 0.044_009_8;
/// The critical temperature, K.
pub const T_CRIT: f64 = 304.1282;
/// The critical molar density, mol/m³.
pub const RHO_CRIT: f64 = 10624.9063;

// Ideal-part parameters. The leading and offset constants are combined.
const A1: f64 = -6.124_871_062_431_9; // 8.37304456 - 14.4979156224319
const A2: f64 = 5.115_596_318_014_53; // -3.70454304 + 8.82013935801453
const C0: f64 = 2.5;
const IDEAL_N: [f64; 5] = [
    1.994_270_42,
    0.621_052_48,
    0.411_952_93,
    1.040_289_22,
    0.083_276_78,
];
const IDEAL_T: [f64; 5] = [3.15163, 6.1119, 6.77708, 11.32384, 27.08792];

// Residual power-exponential terms (34), where `L[i] = 0` means no exponential factor.
const N: [f64; 34] = [
    0.388_568_232_032,
    2.938_547_594_27,
    -5.586_718_853_5,
    -0.767_531_995_925,
    0.317_290_055_804,
    0.548_033_158_978,
    0.122_794_112_203,
    2.165_896_154_32,
    1.584_173_510_97,
    -0.231_327_054_055,
    0.058_116_916_431_4,
    -0.553_691_372_054,
    0.489_466_159_094,
    -0.024_275_739_843_5,
    0.062_494_790_501_7,
    -0.121_758_602_252,
    -0.370_556_852_701,
    -0.016_775_879_700_4,
    -0.119_607_366_38,
    -0.045_619_362_508_8,
    0.035_612_789_270_3,
    -0.007_442_772_713_21,
    -0.001_739_570_490_24,
    -0.021_810_121_289_5,
    0.024_332_166_559_2,
    -0.037_440_133_423_5,
    0.143_387_157_569,
    -0.134_919_690_833,
    -0.023_151_225_053_5,
    0.012_363_125_492_9,
    0.002_105_832_197_29,
    -0.000_339_585_190_264,
    0.005_599_365_177_16,
    -0.000_303_351_180_556,
];
const D: [f64; 34] = [
    1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 1.0, 2.0, 4.0, 5.0, 5.0, 5.0, 6.0, 6.0, 6.0, 1.0, 1.0, 4.0,
    4.0, 4.0, 7.0, 8.0, 2.0, 3.0, 3.0, 5.0, 5.0, 6.0, 7.0, 8.0, 10.0, 4.0, 8.0,
];
const T: [f64; 34] = [
    0.0, 0.75, 1.0, 2.0, 0.75, 2.0, 0.75, 1.5, 1.5, 2.5, 0.0, 1.5, 2.0, 0.0, 1.0, 2.0, 3.0, 6.0,
    3.0, 6.0, 8.0, 6.0, 0.0, 7.0, 12.0, 16.0, 22.0, 24.0, 16.0, 24.0, 8.0, 2.0, 28.0, 14.0,
];
const L: [f64; 34] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0,
    2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 5.0, 6.0,
];
// Gaussian terms.
const GN: [f64; 5] = [
    -213.654_886_883,
    26_641.569_149_3,
    -24_027.212_204_6,
    -283.416_034_24,
    212.472_844_002,
];
const GD: [f64; 5] = [2.0, 2.0, 2.0, 3.0, 3.0];
const GT: [f64; 5] = [1.0, 0.0, 1.0, 3.0, 3.0];
const GBETA: [f64; 5] = [325.0, 300.0, 300.0, 275.0, 275.0];
const GGAMMA: [f64; 5] = [1.16, 1.19, 1.19, 1.25, 1.22];
const GEPS: [f64; 5] = [1.0, 1.0, 1.0, 1.0, 1.0];
const GETA: [f64; 5] = [25.0, 25.0, 25.0, 15.0, 20.0];

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
pub fn ideal(delta: f64, tau: f64) -> IdealHelmholtz {
    let mut alpha = delta.ln() + A1 + A2 * tau + C0 * tau.ln();
    let mut dalpha_dtau = A2 + C0 / tau;
    let mut d2alpha_dtau2 = -C0 / (tau * tau);
    for i in 0..5 {
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

    // Power-exponential terms.
    for i in 0..34 {
        let expc = if L[i] == 0.0 {
            1.0
        } else {
            (-delta.powf(L[i])).exp()
        };
        let common = N[i] * delta.powf(D[i]) * tau.powf(T[i]) * expc;
        let dterm = D[i] / delta
            - if L[i] == 0.0 {
                0.0
            } else {
                L[i] * delta.powf(L[i] - 1.0)
            };
        alpha += common;
        da_dd += common * dterm;
        d2a_dd2 += common
            * (dterm * dterm
                - D[i] / (delta * delta)
                - if L[i] == 0.0 {
                    0.0
                } else {
                    L[i] * (L[i] - 1.0) * delta.powf(L[i] - 2.0)
                });
        let tterm = T[i] / tau;
        da_dt += common * tterm;
        d2a_dt2 += common * tterm * (tterm - 1.0 / tau);
        d2a_dd_dt += common * dterm * tterm;
    }

    // Gaussian terms.
    for i in 0..5 {
        let dr = delta - GEPS[i];
        let tr = tau - GGAMMA[i];
        let common = GN[i]
            * delta.powf(GD[i])
            * tau.powf(GT[i])
            * (-GETA[i] * dr * dr - GBETA[i] * tr * tr).exp();
        let dterm = GD[i] / delta - 2.0 * GETA[i] * dr;
        alpha += common;
        da_dd += common * dterm;
        d2a_dd2 += common * (dterm * dterm - GD[i] / (delta * delta) - 2.0 * GETA[i]);
        let tterm = GT[i] / tau - 2.0 * GBETA[i] * tr;
        da_dt += common * tterm;
        d2a_dt2 += common * (tterm * tterm - GT[i] / (tau * tau) - 2.0 * GBETA[i]);
        d2a_dd_dt += common * dterm * tterm;
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

/// The molar thermodynamic properties of one CO2 state, in SI.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Properties {
    /// Molar density, mol/m³.
    pub rho: f64,
    /// Compressibility factor.
    pub z: f64,
    /// Enthalpy, J/mol.
    pub h: f64,
    /// Entropy, J/(mol.K).
    pub s: f64,
    /// Isobaric heat capacity, J/(mol.K).
    pub cp: f64,
    /// Isochoric heat capacity, J/(mol.K).
    pub cv: f64,
    /// Internal energy, J/mol.
    pub u: f64,
    /// Gibbs energy, J/mol.
    pub g: f64,
    /// Sound speed, m/s.
    pub sound: f64,
    /// Fugacity coefficient.
    pub phi: f64,
    /// Joule-Thomson coefficient, K/Pa.
    pub jt: f64,
}

/// The molar properties at a temperature and molar density.
#[must_use]
pub fn properties(t: f64, rho: f64) -> Properties {
    let delta = rho / RHO_CRIT;
    let tau = T_CRIT / t;
    let id = ideal(delta, tau);
    let res = residual(delta, tau);

    let z = 1.0 + delta * res.alpha_delta;
    let at = id.alpha_tau + res.alpha_tau;
    let att = id.alpha_tau_tau + res.alpha_tau_tau;
    let h = R * t * (1.0 + tau * at + delta * res.alpha_delta);
    let s = R * (tau * at - id.alpha - res.alpha);
    let cv = -R * tau * tau * att;
    let numer = 1.0 + delta * res.alpha_delta - delta * tau * res.alpha_delta_tau;
    let denom = 1.0 + 2.0 * delta * res.alpha_delta + delta * delta * res.alpha_delta_delta;
    let cp = cv + R * numer * numer / denom;
    let u = R * t * tau * at;
    let g = R * t * (id.alpha + res.alpha + 1.0 + delta * res.alpha_delta);
    let dpdrho = R * t * denom;
    let sound = (cp / cv * dpdrho / MOLAR_MASS).max(0.0).sqrt();
    let dpdt = R * rho * numer;
    let jt = (t / (rho * rho) * (dpdt / dpdrho) - 1.0 / rho) / cp;
    let ln_phi = res.alpha + delta * res.alpha_delta - z.ln();

    Properties {
        rho,
        z,
        h,
        s,
        cp,
        cv,
        u,
        g,
        sound,
        phi: ln_phi.exp(),
        jt,
    }
}

/// Solve `P(rho) = p` for the molar density by Newton's method. The gas root starts from
/// the ideal-gas guess; the liquid root from `delta = 3`. `p` in Pa, `rho` in mol/m³.
#[must_use]
pub fn solve_density(t: f64, p: f64, liquid: bool) -> f64 {
    let mut delta = if liquid { 3.0 } else { p / (R * t) / RHO_CRIT };
    let tau = T_CRIT / t;
    for _ in 0..100 {
        let res = residual(delta, tau);
        let f = R * t * RHO_CRIT * delta * (1.0 + delta * res.alpha_delta) - p;
        let df = R
            * t
            * RHO_CRIT
            * (1.0 + 2.0 * delta * res.alpha_delta + delta * delta * res.alpha_delta_delta);
        let new_delta = delta - f / df;
        if (new_delta - delta).abs() < 1e-12 {
            return new_delta * RHO_CRIT;
        }
        delta = new_delta;
    }
    delta * RHO_CRIT
}
