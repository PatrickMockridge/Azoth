//! PC-SAFT's residual Helmholtz energy, per mole, and the layers it is built from.
//!
//! Gross and Sadowski's perturbed-chain statistical associating fluid theory as NeqSim
//! carries it (`PhasePCSAFT`): a hard-sphere term with the Carnahan-Starling-Leland
//! compressibility, a chain term from the Boublik hard-sphere pair correlation, and the
//! two dispersion terms of the perturbed-chain expansion with `C1` from the
//! hard-sphere-chain compressibility.
//!
//! ```text
//! A^R/(RT) = m_bar a_hs - (m_bar - 1) ln g_hs  -  2 pi rho S1 I1  -  pi m_bar rho S2 I2 C1
//! ```
//!
//! **Everything is per mole and dimensionless**, like the rest of this crate: the
//! number density is `N_A/v` rather than NeqSim's `N_A N/V`, so `A^R/(RT)` here is the
//! molar quantity and a caller multiplies by the mole number. The volume is the state's,
//! reached by the caller's own solve - this module has no solver in it, which is what
//! lets it be checked layer by layer against `validation/neqsim/PcsaftProbe.java`.
//!
//! # The parameters
//!
//! `m` segments, each of diameter `sigma` and energy `epsilon/k`. The segment diameter
//! is temperature-dependent - `d_i = sigma_i (1 - 0.12 exp(-3 (epsilon_i/k)/T))` - and
//! the mixture's are weighted by `x_i m_i`, not by `x_i`: a molecule with more segments
//! contributes more of the packing.
//!
//! # The interaction parameter
//!
//! `k_ij` multiplies the pair's dispersion, once in `S1` and squared in `S2`. It is the
//! `KIJPCSAFT` column, a third fit rather than a convention for one of the others, and
//! NeqSim reads it in the only configuration its PC-SAFT can run in.

use azoth_core::{AzothError, Result};

/// Avogadro's constant, **NeqSim's value**, in 1/mol.
///
/// `6.023e23`, which is the pre-2019 definition and not CODATA's `6.02214076e23` - a
/// relative `1.4e-4` apart. Using the modern one would be a silent `1.4e-4` divergence
/// from every number the oracle prints, because the packing fraction and the dispersion
/// terms are built on it, so this is NeqSim's and the difference is recorded rather than
/// corrected. `validation/neqsim/PcsaftProbe.java`'s `f1vol = N_A/v` is the check: it
/// reproduces to the printed digits with `6.023e23` and to four with the modern value.
pub const AVOGADRO: f64 = 6.023e23;

/// One component's PC-SAFT parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PcsaftComponent {
    /// The number of segments, dimensionless.
    pub m: f64,
    /// The temperature-independent segment diameter, in metres.
    pub sigma: f64,
    /// The segment energy over Boltzmann's constant, in K.
    pub epsik: f64,
}

impl PcsaftComponent {
    /// Whether the table gave this component a PC-SAFT set.
    ///
    /// **Zero means absent**, which is how the compiled table spells it - 47 of the 286
    /// rows carry zeros in all three columns and none is blank - so a model refuses a
    /// component without a set rather than solving for a fluid with no segments.
    #[must_use]
    pub fn has_parameters(&self) -> bool {
        self.m > 0.0 && self.sigma > 0.0 && self.epsik > 0.0
    }

    /// The segment diameter at a temperature, in metres.
    #[must_use]
    pub fn d(self, t: f64) -> f64 {
        self.sigma * (1.0 - 0.12 * (-3.0 * self.epsik / t).exp())
    }

    /// `d d_i/dT`, in metres per kelvin.
    ///
    /// The derivative of the line above: `-0.36 sigma_i eps_i exp(-3 eps_i/T)/T^2`.
    /// Carried here rather than differenced because the departure functions need it to
    /// full precision and it is one line.
    #[must_use]
    pub fn d_d_t(self, t: f64) -> f64 {
        -0.36 * self.sigma * self.epsik * (-3.0 * self.epsik / t).exp() / (t * t)
    }
}

/// The `C1`/`I1` polynomial constants, Gross-Sadowski 2001: `[order in m][power of eta]`.
const A_CONST: [[f64; 7]; 3] = [
    [
        0.9105631445,
        0.6361281449,
        2.6861347891,
        -26.547362491,
        97.759208784,
        -159.59154087,
        91.297774084,
    ],
    [
        -0.3084016918,
        0.1860531159,
        -2.5030047259,
        21.419793629,
        -65.255885330,
        83.318680481,
        -33.746922930,
    ],
    [
        -0.0906148351,
        0.4527842806,
        0.5962700728,
        -1.7241829131,
        -4.1302112531,
        13.776631870,
        -8.6728470368,
    ],
];

/// The `I2` polynomial constants, same layout as [`A_CONST`].
const B_CONST: [[f64; 7]; 3] = [
    [
        0.7240946941,
        2.2382791861,
        -4.0025849485,
        -21.003576815,
        26.855641363,
        206.55133841,
        -355.60235612,
    ],
    [
        -0.5755498075,
        0.6995095521,
        3.8925673390,
        -17.215471648,
        192.67226447,
        -161.826465,
        -165.20769346,
    ],
    [
        0.0976883116,
        -0.2557574982,
        -9.1558561530,
        20.642075974,
        -38.804430052,
        93.626774077,
        -29.666905585,
    ],
];

/// Every layer PC-SAFT's Helmholtz energy is built from, at one state.
#[derive(Debug, Clone, PartialEq)]
pub struct PcsaftState {
    /// Each component's segment diameter at this temperature, in metres.
    pub d: Vec<f64>,
    /// `sum_i x_i m_i`, the mixture's segment number.
    pub m_bar: f64,
    /// `sum_i x_i (m_i - 1)`, the chain term's weight.
    pub m_minus_1: f64,
    /// `sum_i x_i m_i d_i^3`, in m^3/mol - what the packing fraction is built from.
    pub md3: f64,
    /// The packing fraction.
    pub eta: f64,
    /// The hard-sphere compressibility `(4 eta - 3 eta^2)/(1 - eta)^2`.
    pub a_hs: f64,
    /// The hard-sphere pair correlation `(1 - eta/2)/(1 - eta)^3`.
    pub g_hs: f64,
    /// The first dispersion sum, `sum_i sum_j x_i x_j m_i m_j (1-k_ij) sqrt(e_i e_j)/T sigma_ij^3`.
    pub s1: f64,
    /// The second dispersion sum, the same with `(e_i e_j)/T^2` and `(1-k_ij)^2`.
    pub s2: f64,
    /// The hard-sphere-chain compressibility `C1`.
    pub c1: f64,
    /// `I1(eta, m_bar)`.
    pub i1: f64,
    /// `I2(eta, m_bar)`.
    pub i2: f64,
    /// The hard-sphere and chain contribution to `A^R/(RT)`, per mole.
    pub f_hc: f64,
    /// The first dispersion contribution, per mole.
    pub f_disp1: f64,
    /// The second dispersion contribution, per mole.
    pub f_disp2: f64,
}

impl PcsaftState {
    /// `A^R/(RT)` for the mixture, per mole.
    #[must_use]
    pub fn f(&self) -> f64 {
        self.f_hc + self.f_disp1 + self.f_disp2
    }
}

/// Every layer, at a temperature, a molar volume and a composition.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the lengths disagree, or the composition does not
///   sum to one.
/// * [`AzothError::OutOfRange`] if `t`, `v` or any component's parameters are not
///   positive, or if the packing fraction reaches one, where the hard-sphere terms
///   diverge.
pub fn state(
    components: &[PcsaftComponent],
    kij: &[f64],
    x: &[f64],
    t: f64,
    v: f64,
) -> Result<PcsaftState> {
    let n = components.len();
    if x.len() != n || kij.len() != n * n {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "{n} components need {n} mole fractions and {} interaction parameters, \
                 but got {} and {}",
                n * n,
                x.len(),
                kij.len()
            ),
        ));
    }
    if !t.is_finite() || t <= 0.0 {
        return Err(AzothError::OutOfRange {
            field: "T".to_string(),
            value: t,
            detail: "PC-SAFT's segment diameter is a function of temperature, and its \
                     dispersion terms carry epsilon/kT"
                .to_string(),
        });
    }
    if !v.is_finite() || v <= 0.0 {
        return Err(AzothError::OutOfRange {
            field: "v".to_string(),
            value: v,
            detail: "the packing fraction is a volume fraction".to_string(),
        });
    }
    for (i, component) in components.iter().enumerate() {
        if !component.has_parameters() {
            return Err(AzothError::invalid_input(
                "components",
                format!(
                    "component {i} carries no PC-SAFT set (m = {}, sigma = {}, epsilon/k = \
                     {}). The table spells an absent set as zeros rather than a blank, so \
                     this is a fluid with no segments and not one to solve for",
                    component.m, component.sigma, component.epsik
                ),
            ));
        }
    }
    let total: f64 = x.iter().sum();
    if (total - 1.0).abs() > 1.0e-9 {
        return Err(AzothError::invalid_input(
            "x",
            format!("the mole fractions sum to {total}, not to one"),
        ));
    }

    let d: Vec<f64> = components.iter().map(|c| c.d(t)).collect();
    let m_bar: f64 = (0..n).map(|i| x[i] * components[i].m).sum();
    let m_minus_1: f64 = (0..n).map(|i| x[i] * (components[i].m - 1.0)).sum();
    let md3: f64 = (0..n).map(|i| x[i] * components[i].m * d[i].powi(3)).sum();

    // `eta = (pi/6) N_A (N/V) sum_i x_i m_i d_i^3`, which per mole is `(pi/6) N_A md3 / v`.
    let eta = std::f64::consts::PI / 6.0 * AVOGADRO * md3 / v;
    if eta >= 1.0 {
        return Err(AzothError::OutOfRange {
            field: "v".to_string(),
            value: v,
            detail: format!(
                "the packing fraction at this volume is {eta}, and the hard-sphere terms \
                 diverge at one. The volume is below the segments' own excluded volume, \
                 which is not a state"
            ),
        });
    }

    let a_hs = (4.0 * eta - 3.0 * eta * eta) / ((1.0 - eta) * (1.0 - eta));
    let g_hs = (1.0 - eta / 2.0) / (1.0 - eta).powi(3);

    // The two dispersion sums. `sigma_ij` is the arithmetic mean and the energy the
    // geometric one, which is the combining rule the whole term is written in.
    let mut s1 = 0.0;
    let mut s2 = 0.0;
    for i in 0..n {
        for j in 0..n {
            let sigma_ij = 0.5 * (components[i].sigma + components[j].sigma);
            let e_ij = (components[i].epsik / t * (components[j].epsik / t)).sqrt();
            let weight = x[i] * x[j] * components[i].m * components[j].m * sigma_ij.powi(3);
            let one_minus_k = 1.0 - kij[i * n + j];
            s1 += weight * e_ij * one_minus_k;
            s2 += weight * e_ij * e_ij * one_minus_k * one_minus_k;
        }
    }

    // `C1`, the hard-sphere-chain compressibility the second dispersion term is scaled by.
    let terms = c1_terms(eta);
    let c1 = 1.0 / (1.0 + m_bar * terms.a + (1.0 - m_bar) * terms.b);

    let i1 = series(&A_CONST, m_bar, eta).value;
    let i2 = series(&B_CONST, m_bar, eta).value;

    let rho = AVOGADRO / v;
    let f_hc = m_bar * a_hs - m_minus_1 * g_hs.ln();
    let f_disp1 = -2.0 * std::f64::consts::PI * rho * s1 * i1;
    let f_disp2 = -std::f64::consts::PI * m_bar * rho * s2 * i2 * c1;

    Ok(PcsaftState {
        d,
        m_bar,
        m_minus_1,
        md3,
        eta,
        a_hs,
        g_hs,
        s1,
        s2,
        c1,
        i1,
        i2,
        f_hc,
        f_disp1,
        f_disp2,
    })
}

/// `C1`'s denominator term by term: `W = 1 + m_bar A + (1 - m_bar) B`.
struct C1Terms {
    a: f64,
    a_d_eta: f64,
    a_d2_eta: f64,
    b: f64,
    b_d_eta: f64,
    b_d2_eta: f64,
}

/// `A`, `B` and their first two derivatives in `eta`.
///
/// Written once because two callers need them - the volume derivative and the composition
/// derivative - and `C1`'s dependence on `m_bar` reaches both. Each has the shape
/// `(num)/(1 - eta)^k` or `num/P^2`; the derivatives are the same quotient rule twice.
fn c1_terms(eta: f64) -> C1Terms {
    let one = 1.0 - eta;

    // `A = (8 eta - 2 eta^2)/(1 - eta)^4`, so `u = 8 + 20 eta - 4 eta^2`.
    let u_a = 8.0 + 20.0 * eta - 4.0 * eta * eta;
    let a_d_eta = u_a / one.powi(5);
    let a_d2_eta = ((20.0 - 8.0 * eta) * one + 5.0 * u_a) / one.powi(6);

    // `B = num/P^2` with `P = (1 - eta)(2 - eta)`, so `B' = f/P^3` and `B'' =
    // (f' P - 3 f P')/P^4`, where `f = num' P - 2 num P'`.
    let p = (1.0 - eta) * (2.0 - eta);
    let p_d_eta = 2.0 * eta - 3.0;
    let num = 20.0 * eta - 27.0 * eta * eta + 12.0 * eta.powi(3) - 2.0 * eta.powi(4);
    let num_d_eta = 20.0 - 54.0 * eta + 36.0 * eta * eta - 8.0 * eta.powi(3);
    let num_d2_eta = -54.0 + 72.0 * eta - 24.0 * eta * eta;
    let f = num_d_eta * p - 2.0 * num * p_d_eta;
    let f_d_eta = num_d2_eta * p - num_d_eta * p_d_eta - 4.0 * num;

    C1Terms {
        a: (8.0 * eta - 2.0 * eta * eta) / one.powi(4),
        a_d_eta,
        a_d2_eta,
        b: num / p.powi(2),
        b_d_eta: f / p.powi(3),
        b_d2_eta: (f_d_eta * p - 3.0 * f * p_d_eta) / p.powi(4),
    }
}

/// One of the two series and every derivative of it a caller needs, at one state.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Series {
    /// `sum_i a_i(m_bar) eta^i`, Gross-Sadowski's `I1` or `I2`.
    value: f64,
    /// `dI/deta`.
    d_eta: f64,
    /// `d^2I/deta^2`.
    d2_eta: f64,
    /// `dI/d m_bar`, which carries `a_i`'s own `m_bar` dependence.
    d_m_bar: f64,
}

/// `sum_{i=0..6} a_i(m_bar) eta^i` and its derivatives. `a_i` is linear in
/// `(m_bar-1)/m_bar` and `(m_bar-1)(m_bar-2)/m_bar^2`, both of which carry `m_bar`.
///
/// The four sums are one loop because they are four weights on the same `a_i`: writing
/// them apart is how a coefficient gets differentiated one way in one and another way in
/// the next, which is what `PhasePCSAFT.dF_HC_SAFTdVdV` did to a `dn/dV`.
fn series(order: &[[f64; 7]; 3], m_bar: f64, eta: f64) -> Series {
    let one = (m_bar - 1.0) / m_bar;
    let two = one * (m_bar - 2.0) / m_bar;
    let mut out = Series {
        value: 0.0,
        d_eta: 0.0,
        d2_eta: 0.0,
        d_m_bar: 0.0,
    };
    for (i, (a0, (a1, a2))) in order[0]
        .iter()
        .zip(order[1].iter().zip(order[2].iter()))
        .enumerate()
    {
        let a_i = a0 + one * a1 + two * a2;
        let power = eta.powi(i as i32);
        out.value += a_i * power;
        out.d_eta += (i as f64) * a_i * eta.powi(i as i32 - 1);
        out.d2_eta += (i as f64) * (i as f64 - 1.0) * a_i * eta.powi(i as i32 - 2);
        // `d/d m_bar` of the two ratios, against the coefficients they carry.
        out.d_m_bar += (a1 / (m_bar * m_bar) + a2 * (3.0 * m_bar - 4.0) / m_bar.powi(3)) * power;
    }
    out
}

/// `(F_eta, F_etaeta)`, the first two derivatives of `A^R/(RT)` in the packing fraction.
///
/// The volume enters the energy only through `eta`, and the number density is
/// proportional to it, so `rho/eta = 6/(pi md3)` is a constant along the whole volume
/// path: every term below is in `eta` alone.
///
/// Three of the layers have the shape `(num)/(1 - eta)^k`, so with
/// `u = (num)' (1 - eta) + k num` the first derivative is `u/(1 - eta)^(k+1)` and the
/// second is `(u' (1 - eta) + (k+1) u)/(1 - eta)^(k+2)`. The first draft of the first
/// derivative dropped the `(1 - eta)` in `u`, so they are written to that shape for all
/// three rather than re-derived per term.
fn eta_derivatives(s: &PcsaftState, v: f64) -> (f64, f64) {
    let eta = s.eta;
    let one = 1.0 - eta;

    // `a_hs = (4 eta - 3 eta^2)/(1 - eta)^2`, so `u = 4 - 2 eta`.
    let a_hs_d_eta = (4.0 - 2.0 * eta) / one.powi(3);
    let a_hs_d2_eta = (10.0 - 4.0 * eta) / one.powi(4);

    // `g_hs = (1 - eta/2)/(1 - eta)^3`, so `u = 5/2 - eta`.
    let g_hs_d_eta = (2.5 - eta) / one.powi(4);
    let g_hs_d2_eta = (9.0 - 3.0 * eta) / one.powi(5);

    // `C1 = 1/W` with `W = 1 + m_bar A + (1 - m_bar) B`, so `C1' = -C1^2 W'` and
    // `C1'' = 2 C1^3 W'^2 - C1^2 W''`.
    let t = c1_terms(eta);
    let w_d_eta = s.m_bar * t.a_d_eta + (1.0 - s.m_bar) * t.b_d_eta;
    let w_d2_eta = s.m_bar * t.a_d2_eta + (1.0 - s.m_bar) * t.b_d2_eta;
    let c1_d_eta = -s.c1 * s.c1 * w_d_eta;
    let c1_d2_eta = 2.0 * s.c1.powi(3) * w_d_eta * w_d_eta - s.c1 * s.c1 * w_d2_eta;

    let a_series = series(&A_CONST, s.m_bar, eta);
    let b_series = series(&B_CONST, s.m_bar, eta);
    let (i1_d_eta, i1_d2_eta) = (a_series.d_eta, a_series.d2_eta);
    let (i2_d_eta, i2_d2_eta) = (b_series.d_eta, b_series.d2_eta);

    // The two dispersion brackets are `(I1 + eta I1')` and `(D + eta D')` with
    // `D = I2 C1`, and a bracket `(X + eta X')` differentiates to `(2 X' + eta X'')`.
    let d_d_eta = i2_d_eta * s.c1 + s.i2 * c1_d_eta;
    let d_d2_eta = i2_d2_eta * s.c1 + 2.0 * i2_d_eta * c1_d_eta + s.i2 * c1_d2_eta;

    let pi = std::f64::consts::PI;
    let rho_over_eta = AVOGADRO / v / eta;
    let f_eta = s.m_bar * a_hs_d_eta
        - s.m_minus_1 * g_hs_d_eta / s.g_hs
        - 2.0 * pi * rho_over_eta * s.s1 * (s.i1 + eta * i1_d_eta)
        - pi * s.m_bar * rho_over_eta * s.s2 * (s.i2 * s.c1 + eta * d_d_eta);

    let f_d2_eta = s.m_bar * a_hs_d2_eta
        - s.m_minus_1 * (g_hs_d2_eta / s.g_hs - (g_hs_d_eta / s.g_hs).powi(2))
        - 2.0 * pi * rho_over_eta * s.s1 * (2.0 * i1_d_eta + eta * i1_d2_eta)
        - pi * s.m_bar * rho_over_eta * s.s2 * (2.0 * d_d_eta + eta * d_d2_eta);

    (f_eta, f_d2_eta)
}

/// The state's pressure at a trial molar volume, over `R T`.
///
/// `P/(RT) = 1/v - F_V`, with `F_V = -eta F_eta / v` - the packing fraction is the only
/// place the volume enters, and every other quantity is a function of the temperature and
/// the composition alone.
///
/// # Errors
/// As [`state`].
pub fn pressure_over_rt(
    components: &[PcsaftComponent],
    kij: &[f64],
    x: &[f64],
    t: f64,
    v: f64,
) -> Result<f64> {
    let s = state(components, kij, x, t, v)?;
    let (f_eta, _) = eta_derivatives(&s, v);
    Ok(1.0 / v + s.eta * f_eta / v)
}

/// `d(P/(RT))/dv` at a trial molar volume, in mol/m^6.
///
/// This is what a volume solve divides by. `P/(RT) = 1/v - f'(v)`, and with
/// `f' = -eta F_eta/v` and `eta` proportional to `1/v`,
/// `f'' = (2 eta F_eta + eta^2 F_etaeta)/v^2`.
///
/// # Errors
/// As [`state`].
pub fn d_pressure_over_rt_dv(
    components: &[PcsaftComponent],
    kij: &[f64],
    x: &[f64],
    t: f64,
    v: f64,
) -> Result<f64> {
    let s = state(components, kij, x, t, v)?;
    let (f_eta, f_d2_eta) = eta_derivatives(&s, v);
    let eta = s.eta;
    Ok(-(1.0 + 2.0 * eta * f_eta + eta * eta * f_d2_eta) / (v * v))
}

/// `ln phi_i = d(nF)/dn_i - ln Z` for every component, at a state.
///
/// The derivative is of the extensive `A^R/(RT)` at fixed temperature and **volume**, so
/// it is the composition and the volume together that move: with `n` the mole numbers,
/// `k` one component and `sum n = 1`,
///
/// ```text
/// d(nF)/dn_k = f + rho f_rho + n_k^eta f_eta
///                  + f_m_bar (m_k - m_bar) + f_m1 ((m_k - 1) - m1)
///                  + f_S1 (S1_k - 2 S1) + f_S2 (S2_k - 2 S2)
/// ```
///
/// where `f_rho` and the `f_eta` beside it are partials at held packing fraction and held
/// density respectively, and `n_k^eta = (pi/6) N_A m_k d_k^3/v` is the packing fraction's
/// own derivative. **That last one is not `eta/v`**: `rho` and `eta` both scale with `1/v`,
/// which is why the volume derivative can use a single path, but at fixed volume the
/// composition moves `eta` through `sum n_i m_i d_i^3` and `rho` through `sum n_i` alone.
/// `S1` and `S2` are quadratic in the mole fractions, so their derivatives are row sums.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if the state's compressibility is not positive, where the
///   logarithm is not defined.
/// * As [`state`].
pub fn ln_fugacity_coefficients(
    components: &[PcsaftComponent],
    kij: &[f64],
    x: &[f64],
    t: f64,
    v: f64,
) -> Result<Vec<f64>> {
    let n = components.len();
    let s = state(components, kij, x, t, v)?;
    let eta = s.eta;
    let one = 1.0 - eta;
    let pi = std::f64::consts::PI;
    let rho = AVOGADRO / v;

    let terms = c1_terms(eta);
    let a_series = series(&A_CONST, s.m_bar, eta);
    let b_series = series(&B_CONST, s.m_bar, eta);
    let w_d_eta = s.m_bar * terms.a_d_eta + (1.0 - s.m_bar) * terms.b_d_eta;
    let c1_d_eta = -s.c1 * s.c1 * w_d_eta;
    // `C1`'s `m_bar` derivative is its denominator's, and `dA/dm_bar - dB/dm_bar` is
    // `A - B` because only the coefficients carry the segment number.
    let c1_d_m_bar = -s.c1 * s.c1 * (terms.a - terms.b);

    // The partials of `f` in each thing the composition moves. Every one of them is a sum
    // over the layers the state already holds.
    let f_rho = -2.0 * pi * s.s1 * s.i1 - pi * s.m_bar * s.s2 * s.i2 * s.c1;
    let f_eta = s.m_bar * (4.0 - 2.0 * eta) / one.powi(3)
        - s.m_minus_1 * (2.5 - eta) / one.powi(4) / s.g_hs
        - 2.0 * pi * rho * s.s1 * a_series.d_eta
        - pi * s.m_bar * rho * s.s2 * (b_series.d_eta * s.c1 + s.i2 * c1_d_eta);
    let f_m_bar = s.a_hs
        - 2.0 * pi * rho * s.s1 * a_series.d_m_bar
        - pi * rho * s.s2 * (s.i2 * s.c1 + s.m_bar * (b_series.d_m_bar * s.c1 + s.i2 * c1_d_m_bar));
    let f_m_minus_1 = -s.g_hs.ln();
    let f_s1 = -2.0 * pi * rho * s.i1;
    let f_s2 = -pi * s.m_bar * rho * s.i2 * s.c1;

    let z = pressure_over_rt(components, kij, x, t, v)? * v;
    if !z.is_finite() || z <= 0.0 {
        return Err(AzothError::OutOfRange {
            field: "Z".to_string(),
            value: z,
            detail: "the fugacity coefficients are `d(nF)/dn_i - ln Z`, and a \
                     compressibility at or below zero is not a state"
                .to_string(),
        });
    }

    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let component = &components[k];
        let eta_k = pi / 6.0 * AVOGADRO * component.m * s.d[k].powi(3) / v;

        let mut s1_k = -2.0 * s.s1;
        let mut s2_k = -2.0 * s.s2;
        for j in 0..n {
            let sigma_ij = 0.5 * (component.sigma + components[j].sigma);
            let e_ij = (component.epsik / t * (components[j].epsik / t)).sqrt();
            let weight = x[j] * component.m * components[j].m * sigma_ij.powi(3);
            let one_minus_k = 1.0 - kij[k * n + j];
            s1_k += 2.0 * weight * e_ij * one_minus_k;
            s2_k += 2.0 * weight * e_ij * e_ij * one_minus_k * one_minus_k;
        }

        let d = s.f()
            + rho * f_rho
            + eta_k * f_eta
            + f_m_bar * (component.m - s.m_bar)
            + f_m_minus_1 * ((component.m - 1.0) - s.m_minus_1)
            + f_s1 * s1_k
            + f_s2 * s2_k;
        out.push(d - z.ln());
    }
    Ok(out)
}

/// `T d(A^R/(R T))/dT` at constant volume, per mole.
///
/// **NeqSim publishes this and its PC-SAFT value cannot be used.** `PhasePCSAFT.getdDSAFTdT`
/// multiplies the chain rule for the segment diameter by an extra `3 d_i^2`, so every
/// `eta`-dependent part of its `dFdT` is about `1e-18` too small: the volume's own
/// derivative is right and the temperature's is not. The write-up is at
/// `~/Desktop/neqsim-pcsaft-hard-chain-temperature-derivative.md`; the check this function
/// is held to is a **finite difference of NeqSim's own `F` at fixed volume**, which uses
/// none of the derivative code the defect reaches.
///
/// # The three terms
///
/// `eta` is the only place the volume and the temperature meet, so at constant `v`
/// everything moves through `T d eta/dT = eta (T d(md3)/dT)/md3` with
/// `md3 = sum_i x_i m_i d_i^3`. Then
///
/// * the hard-sphere and chain term is `m_bar a_hs - (m_bar - 1) ln g_hs`, differentiated
///   in `eta`;
/// * the first dispersion sum carries `sqrt(eps_i eps_j)/T`, so `T dS1/dT = -S1`, and `I1`
///   moves with `eta`;
/// * the second is the same with `C1` beside it, whose own `eta` derivative is
///   `-C1^2 (m_bar A' + (1 - m_bar) B')`.
#[must_use]
pub fn t_d_helmholtz_rt_dt(
    components: &[PcsaftComponent],
    x: &[f64],
    t: f64,
    v: f64,
    state: &PcsaftState,
) -> f64 {
    let n = components.len();
    // `T d(md3)/dT`, which is the whole of the temperature's reach through the packing
    // fraction.
    let t_d_md3 = (0..n)
        .map(|i| x[i] * components[i].m * 3.0 * state.d[i].powi(2) * (t * components[i].d_d_t(t)))
        .sum::<f64>();
    let t_d_eta = state.eta * t_d_md3 / state.md3;

    let eta = state.eta;
    let one = 1.0 - eta;
    // The two hard-sphere layers' own derivatives in `eta`.
    let a_hs_d_eta = ((4.0 - 6.0 * eta) * one + 2.0 * (4.0 * eta - 3.0 * eta * eta)) / one.powi(3);
    let g_hs_d_eta = (2.5 - eta) / one.powi(4);

    let t_d_f_hc =
        state.m_bar * a_hs_d_eta * t_d_eta - state.m_minus_1 * g_hs_d_eta / state.g_hs * t_d_eta;

    let terms = c1_terms(eta);
    let t_d_c1 = -state.c1
        * state.c1
        * (state.m_bar * terms.a_d_eta + (1.0 - state.m_bar) * terms.b_d_eta)
        * t_d_eta;

    let i1_d_eta = series(&A_CONST, state.m_bar, eta).d_eta;
    let i2_d_eta = series(&B_CONST, state.m_bar, eta).d_eta;
    let rho = AVOGADRO / v;

    // `-2 pi rho S1 I1` with `T dS1/dT = -S1` and `I1` moving through `eta`.
    let t_d_f_disp1 =
        -2.0 * std::f64::consts::PI * rho * state.s1 * (i1_d_eta * t_d_eta - state.i1);
    // `-pi m_bar rho S2 I2 C1`, the same with `C1`'s own derivative beside it.
    let t_d_f_disp2 = -std::f64::consts::PI
        * state.m_bar
        * rho
        * state.s2
        * (state.i2 * (-state.c1 + t_d_c1) + state.c1 * (i2_d_eta * t_d_eta - state.i2));

    t_d_f_hc + t_d_f_disp1 + t_d_f_disp2
}
