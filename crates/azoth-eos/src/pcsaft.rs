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
    let c1 = 1.0
        / (1.0
            + m_bar * (8.0 * eta - 2.0 * eta * eta) / (1.0 - eta).powi(4)
            + (1.0 - m_bar)
                * (20.0 * eta - 27.0 * eta * eta + 12.0 * eta.powi(3) - 2.0 * eta.powi(4))
                / ((1.0 - eta) * (2.0 - eta)).powi(2));

    let i1 = series(&A_CONST, m_bar, eta);
    let i2 = series(&B_CONST, m_bar, eta);

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

/// `sum_{i=0..6} a_i(m_bar) eta^i`, Gross-Sadowski's `I1` or `I2`.
fn series(order: &[[f64; 7]; 3], m_bar: f64, eta: f64) -> f64 {
    let one = (m_bar - 1.0) / m_bar;
    let two = one * (m_bar - 2.0) / m_bar;
    let mut out = 0.0;
    for (i, (a0, (a1, a2))) in order[0]
        .iter()
        .zip(order[1].iter().zip(order[2].iter()))
        .enumerate()
    {
        let a_i = a0 + one * a1 + two * a2;
        out += a_i * eta.powi(i as i32);
    }
    out
}

/// `dI/deta` for one of the two series, `sum_i i a_i(m_bar) eta^(i-1)`.
fn series_d_eta(order: &[[f64; 7]; 3], m_bar: f64, eta: f64) -> f64 {
    let one = (m_bar - 1.0) / m_bar;
    let two = one * (m_bar - 2.0) / m_bar;
    let mut out = 0.0;
    for (i, (a0, (a1, a2))) in order[0]
        .iter()
        .zip(order[1].iter().zip(order[2].iter()))
        .enumerate()
        .skip(1)
    {
        let a_i = a0 + one * a1 + two * a2;
        out += (i as f64) * a_i * eta.powi(i as i32 - 1);
    }
    out
}

/// The state's pressure at a trial molar volume, over `R T`.
///
/// `P/(RT) = 1/v - F_V`, with `F_V = -eta F_eta / v` because the packing fraction is the
/// only place the volume enters - `eta` is proportional to `1/v` and every other quantity
/// is a function of the temperature and the composition alone.
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
    let eta = s.eta;
    let d = 1.0 - eta;

    // Each of the three has the shape `(num)/(1 - eta)^k`, so its derivative is
    // `((num)' (1 - eta) + k num)/(1 - eta)^(k+1)` and the denominator never has to be
    // squared into a separate variable.
    let a_hs_d_eta = ((4.0 - 6.0 * eta) * d + 2.0 * (4.0 * eta - 3.0 * eta * eta)) / d.powi(3);
    let g_hs_d_eta = (2.5 - eta) / d.powi(4);

    // `C1 = 1/(1 + m_bar A + (1-m_bar) B)`, so `C1' = -C1^2 (m_bar A' + (1-m_bar) B')`.
    let a_d_eta = ((8.0 - 4.0 * eta) * d + 4.0 * (8.0 * eta - 2.0 * eta * eta)) / d.powi(5);
    let b_num = 20.0 * eta - 27.0 * eta * eta + 12.0 * eta.powi(3) - 2.0 * eta.powi(4);
    let b_den = ((1.0 - eta) * (2.0 - eta)).powi(2);
    let b_num_d_eta = 20.0 - 54.0 * eta + 36.0 * eta * eta - 8.0 * eta.powi(3);
    let b_den_d_eta = 2.0 * (1.0 - eta) * (2.0 - eta) * (2.0 * eta - 3.0);
    let b_d_eta = (b_num_d_eta * b_den - b_num * b_den_d_eta) / (b_den * b_den);
    let c1_d_eta = -s.c1 * s.c1 * (s.m_bar * a_d_eta + (1.0 - s.m_bar) * b_d_eta);

    let i1_d_eta = series_d_eta(&A_CONST, s.m_bar, eta);
    let i2_d_eta = series_d_eta(&B_CONST, s.m_bar, eta);

    // Every term below is proportional to `rho`, which is proportional to `eta`: so
    // `d(rho X)/deta = (rho/eta) d(eta X)/deta`, and the `rho` never has to be written.
    let rho_over_eta = AVOGADRO / v / eta;
    let f_eta = s.m_bar * a_hs_d_eta
        - s.m_minus_1 * g_hs_d_eta / s.g_hs
        - 2.0 * std::f64::consts::PI * rho_over_eta * s.s1 * (s.i1 + eta * i1_d_eta)
        - std::f64::consts::PI
            * s.m_bar
            * rho_over_eta
            * s.s2
            * (s.i2 * s.c1 + eta * i2_d_eta * s.c1 + eta * s.i2 * c1_d_eta);

    Ok(1.0 / v + eta * f_eta / v)
}
