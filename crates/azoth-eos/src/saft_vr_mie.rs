//! SAFT-VR-Mie's parameters and the effective segment diameter.
//!
//! The Mie potential `u(r) = C epsilon ((sigma/r)^lambda_r - (sigma/r)^lambda_a)` replaces
//! the square well of the original SAFT, so a component carries five numbers rather than
//! three: the two exponents, the segment number, the diameter and the energy. Everything
//! downstream is built on the **temperature-dependent diameter** here, which is the
//! Barker-Henderson integral of the potential's Boltzmann factor.
//!
//! The layers the Helmholtz energy is built from come next; this module is the parameters
//! and the diameter, so each can be checked against `validation/neqsim/SaftVrMieProbe.java`
//! before anything is built on it.

use azoth_core::{AzothError, Result};

/// One component's SAFT-VR-Mie parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MieComponent {
    /// The number of segments, dimensionless.
    pub m: f64,
    /// The repulsive exponent, dimensionless.
    pub lambda_r: f64,
    /// The attractive exponent, dimensionless.
    pub lambda_a: f64,
    /// The temperature-independent segment diameter, in metres.
    pub sigma: f64,
    /// The segment energy over Boltzmann's constant, in K.
    pub epsik: f64,
}

impl MieComponent {
    /// Whether the table gave this component a SAFT-VR-Mie set.
    ///
    /// **`m` is the marker and the exponents are not.** The table carries the standard
    /// `12`/`6` on every one of its 286 rows, whether or not the row has a set, so
    /// `lambda_r > 0` says nothing; `m`, `sigma` and `epsik` are zero together on the 274
    /// rows without one.
    #[must_use]
    pub fn has_parameters(&self) -> bool {
        self.m > 0.0 && self.sigma > 0.0 && self.epsik > 0.0
    }

    /// The effective segment diameter at a temperature, in metres.
    ///
    /// # Errors
    /// * [`AzothError::OutOfRange`] if `t` is not positive, or the exponents do not make a
    ///   potential - `lambda_a` must be above three, or the integral diverges as the
    ///   attractive well reaches its limit, and `lambda_r` must exceed it.
    pub fn d(self, t: f64) -> Result<f64> {
        if !t.is_finite() || t <= 0.0 {
            return Err(AzothError::OutOfRange {
                field: "T".to_string(),
                value: t,
                detail: "the effective diameter is a Boltzmann average of the potential, \
                         so it is a function of temperature"
                    .to_string(),
            });
        }
        if self.lambda_r <= self.lambda_a {
            return Err(AzothError::invalid_input(
                "lambda_r",
                format!(
                    "the repulsive exponent {} is not above the attractive one {}, so there \
                     is no potential well for this integral to average",
                    self.lambda_r, self.lambda_a
                ),
            ));
        }
        let theta = mie_prefactor(self.lambda_r, self.lambda_a) * self.epsik / t;
        Ok(self.sigma * barker_henderson(theta, self.lambda_r, self.lambda_a))
    }
}

/// `C = lambda_r/(lambda_r - lambda_a) (lambda_r/lambda_a)^(lambda_a/(lambda_r-lambda_a))`,
/// the constant that puts the Mie well's minimum at `-epsilon`.
#[must_use]
pub fn mie_prefactor(lambda_r: f64, lambda_a: f64) -> f64 {
    lambda_r / (lambda_r - lambda_a) * (lambda_r / lambda_a).powf(lambda_a / (lambda_r - lambda_a))
}

/// The 10-point Gauss-Legendre nodes on `[0, 1]`, **NeqSim's own rounded literals**.
///
/// Taken as NeqSim writes them rather than as the exact nodes, for the reason every other
/// constant in this port is: where the two would disagree, NeqSim wins, and this one is
/// multiplied by the diameter and then read back to fifteen digits.
const GAUSS_NODES: [f64; 10] = [
    0.013_046_735_74,
    0.067_468_316_65,
    0.160_295_215_85,
    0.283_302_302_94,
    0.425_562_830_50,
    0.574_437_169_50,
    0.716_697_697_06,
    0.839_704_784_15,
    0.932_531_683_35,
    0.986_953_264_26,
];

/// The matching weights, same provenance as [`GAUSS_NODES`].
///
/// **Ten points carry this integrand poorly and that is the model's choice, not the
/// port's.** Measured against a converged integral the weak branch is `1e-4` out at
/// `theta = 0.9` and `4e-3` at `0.05`, because the integrand is a step-like
/// `1 - exp(-u*)` and ten points cannot resolve the knee. Reproducing NeqSim means
/// reproducing the rule; improving it would be porting a different model.
const GAUSS_WEIGHTS: [f64; 10] = [
    0.033_335_672_15,
    0.074_725_674_58,
    0.109_543_181_26,
    0.134_633_359_65,
    0.147_762_112_36,
    0.147_762_112_36,
    0.134_633_359_65,
    0.109_543_181_26,
    0.074_725_674_58,
    0.033_335_672_15,
];

/// `integral_0^1 (1 - exp(-u*(x))) dx` with `u* = theta (x^-lambda_r - x^-lambda_a)`, the
/// Barker-Henderson integral.
///
/// **Two branches, and the split is NeqSim's rather than an improvement on it.** For
/// `theta <= 1` the integrand is smooth on `[0, 1]` and ten Gauss-Legendre points carry
/// it. Above that the well is deep enough that `exp(-u*)` underflows near the origin -
/// `1 - exp(-u*)` is one to machine precision over a growing sub-interval - so NeqSim
/// finds where that stops being true, `x_cut = exp(-ln(theta + 20)/lambda_r)` clamped to
/// `[0.5, 0.999]`, treats `[0, x_cut]` as exactly `x_cut`, and puts the same ten points on
/// what is left. The approximation is the model's, not the port's.
#[must_use]
pub fn barker_henderson(theta: f64, lambda_r: f64, lambda_a: f64) -> f64 {
    let reduced = |x: f64| theta * (x.powf(-lambda_r) - x.powf(-lambda_a));

    if theta <= 1.0 {
        let mut sum = 0.0;
        for (node, weight) in GAUSS_NODES.iter().zip(GAUSS_WEIGHTS.iter()) {
            if *node < 1.0e-20 {
                sum += weight;
                continue;
            }
            sum += weight * (1.0 - (-reduced(*node)).exp());
        }
        return sum;
    }

    let x_cut = (-(theta + 20.0).ln() / lambda_r).exp().clamp(0.5, 0.999);
    let (low, high) = (x_cut, 1.0);
    let (half, middle) = ((high - low) / 2.0, (high + low) / 2.0);
    let mut sum = x_cut;
    for (node, weight) in GAUSS_NODES.iter().zip(GAUSS_WEIGHTS.iter()) {
        let x = middle + half * (2.0 * node - 1.0);
        if x < 1.0e-20 {
            sum += weight * (high - low);
            continue;
        }
        sum += weight * (1.0 - (-reduced(x)).exp()) * (high - low);
    }
    sum
}

/// Lafitte 2013's effective-packing-fraction coefficients, `[coefficient][power of 1/lambda]`.
///
/// `c_i(lambda) = A[i][0] + A[i][1]/lambda + A[i][2]/lambda^2 + A[i][3]/lambda^3`, from
/// that paper's Table 5. NeqSim's own rounded literals, like the Gauss nodes above.
const ETA_EFF_COEFFS: [[f64; 4]; 4] = [
    [0.81096, 1.7888, -37.578, 92.284],
    [1.0205, -19.341, 151.26, -463.50],
    [-1.9057, 22.845, -228.14, 973.92],
    [1.08850, -6.1962, 106.98, -677.64],
];

/// `eta_eff = c1 eta + c2 eta^2 + c3 eta^3 + c4 eta^4`, the packing fraction a Mie
/// exponent's own hard-sphere reference is evaluated at.
///
/// The whole chain term is built on this: the bare `aS1` a segment's pair correlation
/// carries is the Carnahan-Starling contact value at `eta_eff`, not at `eta`, because the
/// softness of the potential lets neighbours interpenetrate.
#[must_use]
pub fn eta_effective(eta: f64, lambda: f64) -> f64 {
    let inv = 1.0 / lambda;
    let mut out = 0.0;
    let mut power = eta;
    for coefficient in &ETA_EFF_COEFFS {
        let c = coefficient[0]
            + coefficient[1] * inv
            + coefficient[2] * inv * inv
            + coefficient[3] * inv * inv * inv;
        out += c * power;
        power *= eta;
    }
    out
}

/// The bare `aS1(eta, lambda)` a Mie exponent contributes: the contact value at its own
/// effective packing fraction, scaled by `1/(lambda - 3)`.
#[must_use]
pub fn a_s1_bare(eta: f64, lambda: f64) -> f64 {
    let effective = eta_effective(eta, lambda);
    let one = 1.0 - effective;
    let g_cs = (1.0 - effective / 2.0) / one.powi(3);
    -g_cs / (lambda - 3.0)
}

/// The bare `B(eta, lambda, x0)` a Mie exponent contributes, `x0 = sigma/d`.
#[must_use]
pub fn b_bare(eta: f64, lambda: f64, x0: f64) -> f64 {
    let x0_3l = x0.powf(3.0 - lambda);
    let cap_i = (1.0 - x0_3l) / (lambda - 3.0);
    let cap_j = (1.0 - (lambda - 3.0) * x0.powf(4.0 - lambda) + (lambda - 4.0) * x0_3l)
        / ((lambda - 3.0) * (lambda - 4.0));
    let one = 1.0 - eta;
    cap_i * (1.0 - eta / 2.0) / one.powi(3) - 9.0 * cap_j * eta * (eta + 1.0) / (2.0 * one.powi(3))
}

/// `K_HS`, the hard-sphere isothermal compressibility the second-order chain term is
/// scaled by.
#[must_use]
pub fn k_hs(eta: f64) -> f64 {
    let one = 1.0 - eta;
    one.powi(4) / (1.0 + 4.0 * eta + 4.0 * eta * eta - 4.0 * eta.powi(3) + eta.powi(4))
}

/// Lafitte's `alpha = C (1/(lambda_a - 3) - 1/(lambda_r - 3))`, the well's softness.
#[must_use]
pub fn mie_alpha(lambda_r: f64, lambda_a: f64) -> f64 {
    mie_prefactor(lambda_r, lambda_a) * (1.0 / (lambda_a - 3.0) - 1.0 / (lambda_r - 3.0))
}

/// The Carnahan-Starling compressibility `(4 eta - 3 eta^2)/(1 - eta)^2`, which SAFT-VR-Mie
/// shares with PC-SAFT: the hard-sphere term is the same one, and only the diameter it is
/// evaluated at differs.
#[must_use]
pub fn a_hs(eta: f64) -> f64 {
    let one = 1.0 - eta;
    (4.0 * eta - 3.0 * eta * eta) / (one * one)
}

/// The Carnahan-Starling contact value `(1 - eta/2)/(1 - eta)^3`.
///
/// **Not the mixture's `g` when `m_minus_1 > 0`** - see [`chain_contact_value`], which
/// replaces it.
#[must_use]
pub fn g_hs(eta: f64) -> f64 {
    let one = 1.0 - eta;
    (1.0 - eta / 2.0) / one.powi(3)
}

/// The Mie contact value at one component's own exponents, `g_HS0(eta, x0)`.
///
/// A quartic in `x0 = sigma/d` with `eta`-dependent coefficients: the hard-sphere value of
/// a soft potential is not the hard-sphere value of a hard one, and the difference grows
/// with the ratio of the potential's range to the segment's effective size.
#[must_use]
pub fn contact_value_0(eta: f64, x0: f64) -> f64 {
    let one = 1.0 - eta;
    let (eta2, eta3, eta4) = (eta * eta, eta * eta * eta, eta.powi(4));
    let om3 = one.powi(3);
    let k0 = -one.ln() + (42.0 * eta - 39.0 * eta2 + 9.0 * eta3 - 2.0 * eta4) / (6.0 * om3);
    let k1 = (-12.0 * eta + 6.0 * eta2 + eta4) / (2.0 * om3);
    let k2 = -3.0 * eta2 / (8.0 * one * one);
    let k3 = (3.0 * eta + 3.0 * eta2 - eta4) / (6.0 * om3);
    (k0 + k1 * x0 + k2 * x0 * x0 + k3 * x0 * x0 * x0).exp()
}

/// The step NeqSim's numerical `eta` derivatives use, `max(|eta| 1e-5, 1e-12)`, with its
/// floor on the lower point.
///
/// Written once because four derivatives use it and a second copy is a second chance to
/// differ; the floor matters, because at a small packing fraction `eta - step` would
/// otherwise go negative and the derivative would be taken across a state that is not one.
fn eta_step(eta: f64) -> (f64, f64, f64) {
    let step = (eta.abs() * 1.0e-5).max(1.0e-12);
    let high = eta + step;
    let low = (eta - step).max(1.0e-15);
    (high, low, (high - low) / 2.0)
}

/// `g1`, the first-order chain perturbation at one component's exponents.
///
/// Contains two `eta` derivatives NeqSim takes by central difference; the shape is
/// reproduced rather than differentiated, so the port lands on the same number for the same
/// reason.
#[must_use]
pub fn chain_g1(eta: f64, lambda_r: f64, lambda_a: f64, c_mie: f64, x0: f64) -> f64 {
    let bare = |e: f64, lambda: f64| a_s1_bare(e, lambda) + b_bare(e, lambda, x0);
    let (high, low, step) = eta_step(eta);

    let d_full = |lambda: f64| (high * bare(high, lambda) - low * bare(low, lambda)) / (2.0 * step);
    let da1_drho =
        c_mie * (x0.powf(lambda_a) * d_full(lambda_a) - x0.powf(lambda_r) * d_full(lambda_r));

    3.0 * da1_drho
        - c_mie
            * (lambda_a * x0.powf(lambda_a) * bare(eta, lambda_a)
                - lambda_r * x0.powf(lambda_r) * bare(eta, lambda_r))
}

/// `g2`, the second-order chain perturbation at one component's exponents.
///
/// The `gamma_c` factor is Lafitte 2013's Eq. 39 correction, which switches off the
/// second-order term as the potential softens past `alpha = 0.57`.
#[must_use]
pub fn chain_g2(
    eta: f64,
    zeta_st: f64,
    lambda_r: f64,
    lambda_a: f64,
    eps_over_kt: f64,
    c_mie: f64,
    x0: f64,
) -> f64 {
    let bare = |e: f64, lambda: f64| a_s1_bare(e, lambda) + b_bare(e, lambda, x0);
    let inner = |e: f64| {
        x0.powf(2.0 * lambda_a) * bare(e, 2.0 * lambda_a)
            - 2.0 * x0.powf(lambda_a + lambda_r) * bare(e, lambda_a + lambda_r)
            + x0.powf(2.0 * lambda_r) * bare(e, 2.0 * lambda_r)
    };

    let (high, low, step) = eta_step(eta);
    let da2_product =
        (high * k_hs(high) * inner(high) - low * k_hs(low) * inner(low)) / (2.0 * step);
    let da2_drho = 0.5 * c_mie * c_mie * da2_product;

    let g_mca2 = 3.0 * da2_drho
        - k_hs(eta)
            * c_mie
            * c_mie
            * (lambda_r * x0.powf(2.0 * lambda_r) * bare(eta, 2.0 * lambda_r)
                - (lambda_a + lambda_r)
                    * x0.powf(lambda_a + lambda_r)
                    * bare(eta, lambda_a + lambda_r)
                + lambda_a * x0.powf(2.0 * lambda_a) * bare(eta, 2.0 * lambda_a));

    let alpha = mie_alpha(lambda_r, lambda_a);
    let theta = eps_over_kt.exp() - 1.0;
    let gamma_c = 10.0
        * (-(10.0 * (0.57 - alpha)).tanh() + 1.0)
        * zeta_st
        * theta
        * (-6.7 * zeta_st - 8.0 * zeta_st * zeta_st).exp();

    (1.0 + gamma_c) * g_mca2
}

/// The mixture's chain contact value: `exp(sum_i w_i ln g_i / sum_i w_i)`.
///
/// **A pure fluid with `m = 1` keeps the Carnahan-Starling value**, and that is the model's
/// behaviour rather than a special case here: `w_i = x_i (m_i - 1)` is zero for a
/// one-segment molecule, so the weighted mean has no weight to take and NeqSim's block is
/// skipped outright. A mixture always has weight, and then the value is *not* the contact
/// value - which is why `validation/neqsim/SaftVrMieProbe.java` prints methane's `gHS` as
/// the CS number and methane/n-butane's as the blended one.
#[must_use]
pub fn chain_contact_value(
    components: &[MieComponent],
    x: &[f64],
    t: f64,
    eta: f64,
    diameters: &[f64],
) -> f64 {
    let weight: f64 = x
        .iter()
        .zip(components)
        .map(|(xi, c)| xi * (c.m - 1.0))
        .sum();
    if weight <= 1.0e-10 {
        return g_hs(eta);
    }
    let mut weighted = 0.0;
    for ((xi, component), d) in x.iter().zip(components).zip(diameters) {
        let w = xi * (component.m - 1.0);
        if w < 1.0e-30 {
            continue;
        }
        let x0 = if *d > 0.0 { component.sigma / d } else { 1.0 };
        let eps_over_kt = component.epsik / t;
        let c_mie = mie_prefactor(component.lambda_r, component.lambda_a);
        let zeta_st = eta * x0 * x0 * x0;
        let g0 = contact_value_0(eta, x0);
        let g1 = chain_g1(eta, component.lambda_r, component.lambda_a, c_mie, x0);
        let g2 = chain_g2(
            eta,
            zeta_st,
            component.lambda_r,
            component.lambda_a,
            eps_over_kt,
            c_mie,
            x0,
        );
        // The full Mie correction: the class's default blend fraction is `1.0`, and its
        // adaptive selection only ever lowers it.
        let g_mie = g0 * (eps_over_kt * (g1 + eps_over_kt * g2) / g0).exp();
        weighted += w * g_mie.ln();
    }
    (weighted / weight).exp()
}

/// Lafitte 2013's Padé coefficients, `[coefficient][function]`.
///
/// `f_i(alpha) = (T[0][i] + T[1][i] a + T[2][i] a^2 + T[3][i] a^3)
/// / (1 + T[4][i] a + T[5][i] a^2 + T[6][i] a^3)` with `a` the Mie softness
/// [`mie_alpha`]. Six functions, used by `chi` and by the third-order term.
const PHI_PADE: [[f64; 6]; 7] = [
    [7.5365557, -359.440, 1550.9, -1.199320, -1911.2800, 9236.9],
    [-37.604630, 1825.60, -5070.1, 9.063632, 21390.175, -129430.0],
    [71.745953, -3168.00, 6534.6, -17.94820, -51320.700, 357230.0],
    [-46.835520, 1884.20, -3288.7, 11.34027, 37064.540, -315530.0],
    [-2.4679820, -0.82376, -2.7171, 20.52142, 1103.7420, 1390.2],
    [-0.5027200, -3.19350, 2.0883, -56.63770, -3264.6100, -4518.2],
    [8.0956883, 3.70900, 0.0000, 40.53683, 2556.1810, 4241.6],
];

/// `f_i(alpha)`, Lafitte 2013's Padé approximant of the `i`th function.
///
/// # Panics
/// Never: the index is a compile-time constant at every call site, and the table is
/// checked to have seven rows. A caller passing an index outside the six functions would
/// be a defect rather than a runtime condition.
#[must_use]
pub fn pade_f(index: usize, alpha: f64) -> f64 {
    let a2 = alpha * alpha;
    let numerator = PHI_PADE[0][index]
        + PHI_PADE[1][index] * alpha
        + PHI_PADE[2][index] * a2
        + PHI_PADE[3][index] * a2 * alpha;
    let denominator = 1.0
        + PHI_PADE[4][index] * alpha
        + PHI_PADE[5][index] * a2
        + PHI_PADE[6][index] * a2 * alpha;
    numerator / denominator
}

/// `a1^S(eta, lambda)`, Sutherland's first-order attractive term at one exponent.
#[must_use]
pub fn a1_sutherland(eta: f64, lambda: f64, eps_over_kt: f64) -> f64 {
    let effective = eta_effective(eta, lambda);
    let one = 1.0 - effective;
    -12.0 * eps_over_kt * eta / (lambda - 3.0) * (1.0 - effective / 2.0) / one.powi(3)
}

/// `B(eta, lambda, x0)`, the correction the softness of the potential adds to Sutherland's
/// term.
///
/// The same `capI`/`capJ` pair the chain's bare `B` uses, but against the *real* contact
/// value rather than the effective packing fraction's.
#[must_use]
pub fn b_correction(eta: f64, lambda: f64, eps_over_kt: f64, x0: f64) -> f64 {
    let x0_3ml = x0.powf(3.0 - lambda);
    let x0_4ml = x0.powf(4.0 - lambda);
    let cap_i = (1.0 - x0_3ml) / (lambda - 3.0);
    let cap_j = (1.0 - (lambda - 3.0) * x0_4ml + (lambda - 4.0) * x0_3ml)
        / ((lambda - 3.0) * (lambda - 4.0));
    let one = 1.0 - eta;
    let om3 = one.powi(3);
    let contact = (1.0 - eta / 2.0) / om3;
    let b_bar = contact * cap_i - 9.0 * eta * (1.0 + eta) / (2.0 * om3) * cap_j;
    12.0 * eta * eps_over_kt * b_bar
}

/// `a_1`, the mean-attractive term at one component's exponents.
#[must_use]
pub fn a1_mie(
    eta: f64,
    lambda_r: f64,
    lambda_a: f64,
    eps_over_kt: f64,
    c_mie: f64,
    x0: f64,
) -> f64 {
    let bare = |lambda: f64| {
        a1_sutherland(eta, lambda, eps_over_kt) + b_correction(eta, lambda, eps_over_kt, x0)
    };
    c_mie * (x0.powf(lambda_a) * bare(lambda_a) - x0.powf(lambda_r) * bare(lambda_r))
}

/// `chi`, the correction the second-order term's density dependence carries.
#[must_use]
pub fn mie_chi(zeta_st: f64, lambda_r: f64, lambda_a: f64) -> f64 {
    let alpha = mie_alpha(lambda_r, lambda_a);
    pade_f(0, alpha) * zeta_st
        + pade_f(1, alpha) * zeta_st.powi(5)
        + pade_f(2, alpha) * zeta_st.powi(8)
}

/// `a_2`, the second-order perturbation term at one component's exponents.
#[must_use]
pub fn a2_mie(
    eta: f64,
    zeta_st: f64,
    lambda_r: f64,
    lambda_a: f64,
    eps_over_kt: f64,
    c_mie: f64,
    x0: f64,
) -> f64 {
    let bare = |lambda: f64| {
        a1_sutherland(eta, lambda, eps_over_kt) + b_correction(eta, lambda, eps_over_kt, x0)
    };
    let inner = x0.powf(2.0 * lambda_a) * bare(2.0 * lambda_a)
        - 2.0 * x0.powf(lambda_a + lambda_r) * bare(lambda_a + lambda_r)
        + x0.powf(2.0 * lambda_r) * bare(2.0 * lambda_r);
    let chi = mie_chi(zeta_st, lambda_r, lambda_a);
    0.5 * k_hs(eta) * (1.0 + chi) * eps_over_kt * c_mie * c_mie * inner
}

/// `a_3`, the third-order term, which is a function of the reduced density alone.
#[must_use]
pub fn a3_mie(zeta_st: f64, lambda_r: f64, lambda_a: f64, eps_over_kt: f64) -> f64 {
    let alpha = mie_alpha(lambda_r, lambda_a);
    -(eps_over_kt.powi(3))
        * pade_f(3, alpha)
        * zeta_st
        * (pade_f(4, alpha) * zeta_st + pade_f(5, alpha) * zeta_st * zeta_st).exp()
}

/// The three dispersion terms of a mixture, by pair summation.
///
/// Lafitte 2013 Eqs. 37-40 with the cross parameters of its Eq. 36, and **the cross
/// parameters are not averages of the pure ones**: `sigma_ij` is arithmetic, `epsilon_ij`
/// carries a `sigma^3` correction that keeps the well's volume, and the two exponents are
/// `3 + sqrt((lambda_i - 3)(lambda_j - 3))`. The diameter is then the Barker-Henderson
/// integral of the *cross* potential rather than the mean of the two pure diameters - so
/// the pair's `x0` is a genuinely new number, not an interpolation.
///
/// The weights are segment fractions, `xi_i m_i / m_bar`, because the dispersion is a sum
/// over segments rather than over molecules.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the lengths disagree.
pub fn dispersion_pair_sum(
    components: &[MieComponent],
    x: &[f64],
    t: f64,
    eta: f64,
) -> Result<(f64, f64, f64)> {
    let n = components.len();
    if x.len() != n {
        return Err(AzothError::invalid_input(
            "x",
            format!(
                "{n} components need {n} mole fractions, but got {}",
                x.len()
            ),
        ));
    }
    let m_bar: f64 = x.iter().zip(components).map(|(xi, c)| xi * c.m).sum();
    if m_bar <= 0.0 {
        return Err(AzothError::invalid_input(
            "components",
            "the mixture's segment number is zero, so no segment fraction is defined".to_string(),
        ));
    }
    let segment: Vec<f64> = x
        .iter()
        .zip(components)
        .map(|(xi, c)| xi * c.m / m_bar)
        .collect();

    let mut sum = (0.0, 0.0, 0.0);
    for (i, ci) in components.iter().enumerate() {
        for (j, cj) in components.iter().enumerate() {
            let weight = segment[i] * segment[j];
            if weight < 1.0e-30 {
                continue;
            }
            let sigma_ij = 0.5 * (ci.sigma + cj.sigma);
            let sigma3 = ci.sigma.powi(3) * cj.sigma.powi(3);
            let eps_ij = (ci.epsik * cj.epsik).sqrt() * sigma3.sqrt() / sigma_ij.powi(3);
            let lambda_r = 3.0 + ((ci.lambda_r - 3.0) * (cj.lambda_r - 3.0)).sqrt();
            let lambda_a = 3.0 + ((ci.lambda_a - 3.0) * (cj.lambda_a - 3.0)).sqrt();

            let cross = MieComponent {
                m: 1.0,
                lambda_r,
                lambda_a,
                sigma: sigma_ij,
                epsik: eps_ij,
            };
            let d_ij = cross.d(t)?;
            let c_mie = mie_prefactor(lambda_r, lambda_a);
            let x0 = if d_ij > 0.0 { sigma_ij / d_ij } else { 1.0 };
            let beta = eps_ij / t;
            let zeta = eta * x0 * x0 * x0;

            sum.0 += weight * a1_mie(eta, lambda_r, lambda_a, beta, c_mie, x0);
            sum.1 += weight * a2_mie(eta, zeta, lambda_r, lambda_a, beta, c_mie, x0);
            sum.2 += weight * a3_mie(zeta, lambda_r, lambda_a, beta);
        }
    }
    Ok(sum)
}

/// Every layer SAFT-VR-Mie's Helmholtz energy is built from, at one state.
#[derive(Debug, Clone, PartialEq)]
pub struct MieState {
    /// Each component's effective segment diameter at this temperature, in metres.
    pub d: Vec<f64>,
    /// `sum_i x_i m_i`, the mixture's segment number.
    pub m_bar: f64,
    /// `sum_i x_i (m_i - 1)`, the chain term's weight.
    pub m_minus_1: f64,
    /// `sum_i x_i m_i d_i^3`, in m^3/mol.
    pub md3: f64,
    /// The packing fraction.
    pub eta: f64,
    /// The hard-sphere compressibility at the mixture's packing fraction.
    pub a_hs: f64,
    /// The chain term's contact value - the Carnahan-Starling one for a fluid with no
    /// chain, the Mie-weighted one otherwise.
    pub g_hs: f64,
    /// The first dispersion term, per mole.
    pub a1: f64,
    /// The second.
    pub a2: f64,
    /// The third.
    pub a3: f64,
}

impl MieState {
    /// `A^R/(RT)` per mole: `m_bar a_hs - m_minus_1 ln g_hs + m_bar (a_1 + a_2 + a_3)`.
    #[must_use]
    pub fn f(&self) -> f64 {
        self.f_hc() + self.f_disp()
    }

    /// The hard-sphere and chain part.
    #[must_use]
    pub fn f_hc(&self) -> f64 {
        self.m_bar * self.a_hs - self.m_minus_1 * self.g_hs.ln()
    }

    /// The three dispersion terms together, which is what NeqSim's `F_DISP_SAFT` is.
    ///
    /// **Scaled by the segment number.** NeqSim's `F_DISP_SAFT` is
    /// `n m_bar (a_1 + a_2 + a_3)`, and the `m_bar` is easy to lose because a one-segment
    /// fluid has `m_bar = 1` and agrees either way: methane's `F_disp` is the plain sum to
    /// every printed digit, and methane/n-butane's is 1.34056 times it. The same trap as
    /// the chain RDF's, one layer down.
    #[must_use]
    pub fn f_disp(&self) -> f64 {
        self.m_bar * (self.a1 + self.a2 + self.a3)
    }
}

/// Every layer at a temperature, a molar volume and a composition.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the lengths disagree, the composition is not one, or a
///   component has no set.
/// * [`AzothError::OutOfRange`] if `t`, `v` or a parameter is not positive, or the packing
///   fraction reaches one.
pub fn state(components: &[MieComponent], x: &[f64], t: f64, v: f64) -> Result<MieState> {
    let n = components.len();
    if x.len() != n {
        return Err(AzothError::invalid_input(
            "x",
            format!(
                "{n} components need {n} mole fractions, but got {}",
                x.len()
            ),
        ));
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
                    "component {i} carries no SAFT-VR-Mie set (m = {}, sigma = {}, \
                     epsilon/k = {}). The table spells an absent set as zeros rather than a \
                     blank, so this is a fluid with no segments",
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

    let d = components
        .iter()
        .map(|c| c.d(t))
        .collect::<Result<Vec<_>>>()?;
    let m_bar: f64 = x.iter().zip(components).map(|(xi, c)| xi * c.m).sum();
    let m_minus_1: f64 = x
        .iter()
        .zip(components)
        .map(|(xi, c)| xi * (c.m - 1.0))
        .sum();
    let md3: f64 = x
        .iter()
        .zip(components)
        .zip(&d)
        .map(|((xi, c), di)| xi * c.m * di.powi(3))
        .sum();

    // `eta = (pi/6) N_A (N/V) sum_i x_i m_i d_i^3`, which per mole is `(pi/6) N_A md3 / v`.
    let eta = std::f64::consts::PI / 6.0 * super::pcsaft::AVOGADRO * md3 / v;
    if eta >= 1.0 {
        return Err(AzothError::OutOfRange {
            field: "v".to_string(),
            value: v,
            detail: format!(
                "the packing fraction at this volume is {eta}, and the hard-sphere terms \
                 diverge at one"
            ),
        });
    }

    let (a1, a2, a3) = dispersion_at(components, x, t, eta, &d)?;

    let g_hs = chain_contact_value(components, x, t, eta, &d);
    Ok(MieState {
        d,
        m_bar,
        m_minus_1,
        md3,
        eta,
        a_hs: a_hs(eta),
        g_hs,
        a1,
        a2,
        a3,
    })
}

/// The three dispersion terms at a packing fraction, by NeqSim's branch.
///
/// **One component evaluates directly and a mixture sums over pairs**, and the direct
/// expression is not the pair sum at `n = 1` - it is a different expression whose cross
/// parameters happen to be the pure ones. Taking NeqSim's branch keeps the two agreeing by
/// construction rather than by algebra.
///
/// `diameters` are the *pure* components' effective diameters, which the direct branch
/// needs; the pair sum computes its own cross diameters from the temperature.
///
/// # Errors
/// As [`dispersion_pair_sum`].
pub fn dispersion_at(
    components: &[MieComponent],
    x: &[f64],
    t: f64,
    eta: f64,
    diameters: &[f64],
) -> Result<(f64, f64, f64)> {
    if components.len() != 1 {
        return dispersion_pair_sum(components, x, t, eta);
    }
    let c = components[0];
    let d = diameters.first().copied().unwrap_or(0.0);
    let x0 = if d > 0.0 { c.sigma / d } else { 1.0 };
    let beta = c.epsik / t;
    let c_mie = mie_prefactor(c.lambda_r, c.lambda_a);
    let zeta = eta * x0 * x0 * x0;
    Ok((
        a1_mie(eta, c.lambda_r, c.lambda_a, beta, c_mie, x0),
        a2_mie(eta, zeta, c.lambda_r, c.lambda_a, beta, c_mie, x0),
        a3_mie(zeta, c.lambda_r, c.lambda_a, beta),
    ))
}

/// The state's pressure at a trial molar volume, over `R T`.
///
/// `P/(RT) = 1/v - F_V` with `F_V = -eta f_eta/v`, because the packing fraction is the only
/// place the volume enters.
///
/// **The `eta` derivatives are taken the way NeqSim takes them**, which is not one way: the
/// hard-sphere term is differentiated in closed form there and the chain contact value and
/// all three dispersion terms are central differences at a relative step of `1e-5`, with a
/// floor on the lower point. The port does the same rather than deriving them.
///
/// **At that step the difference is noise-limited on both sides**, and it is worth knowing
/// before reading a test's tolerance as a model error. `P v/(RT)` at NeqSim's own converged
/// volume comes out `0.851583259006834` against the probe's `0.851583764555351` - a
/// relative `5.9e-7` - and *shrinking* the step makes it worse rather than better
/// (`1.4e-6` at `1e-8`), which is what a difference of two nearly equal energies does.
/// Ten times larger, where the cancellation is not in charge, the same port gives
/// `0.851583766815854`: a relative `2.7e-9`. So the model is right and the `5.9e-7` is the
/// arithmetic both implementations are doing.
///
/// # Errors
/// As [`state`].
pub fn pressure_over_rt(components: &[MieComponent], x: &[f64], t: f64, v: f64) -> Result<f64> {
    let s = state(components, x, t, v)?;
    let eta = s.eta;

    // NeqSim's step, in the packing fraction: `max(|eta| 1e-5, 1e-12)` with the lower point
    // floored, and the halved span that its `etaM` correction makes.
    let step = (eta.abs() * 1.0e-5).max(1.0e-12);
    let high = eta + step;
    let low = (eta - step).max(1.0e-15);
    let step = (high - low) / 2.0;

    let g_eta = (chain_contact_value(components, x, t, high, &s.d)
        - chain_contact_value(components, x, t, low, &s.d))
        / (2.0 * step);
    let dispersion = |e: f64| -> Result<f64> {
        let (a1, a2, a3) = dispersion_at(components, x, t, e, &s.d)?;
        Ok(a1 + a2 + a3)
    };
    let dispersion_eta = (dispersion(high)? - dispersion(low)?) / (2.0 * step);

    // The hard-sphere term is the one NeqSim differentiates in closed form.
    let one = 1.0 - eta;
    let a_hs_eta = (4.0 - 2.0 * eta) / one.powi(3);

    let f_eta = s.m_bar * a_hs_eta - s.m_minus_1 * g_eta / s.g_hs + s.m_bar * dispersion_eta;
    Ok(1.0 / v + eta * f_eta / v)
}

/// The dispersion's **row sums**: `sum_l xi_l a^{il}` for each component `i`.
///
/// Not a derivative and not the total: the total is `sum_i xi_i` times this. It is what the
/// *composition* derivative of the quadratic form needs - differentiating
/// `a = sum_i sum_j xi_i xi_j a^{ij}` in one mole number brings half of this down with it -
/// and NeqSim computes exactly these, in `calcPairDispSumPerComp`, for that reason.
///
/// # Errors
/// As [`dispersion_pair_sum`].
pub fn dispersion_row_sums(
    components: &[MieComponent],
    x: &[f64],
    t: f64,
    eta: f64,
) -> Result<Vec<f64>> {
    let n = components.len();
    if x.len() != n {
        return Err(AzothError::invalid_input(
            "x",
            format!(
                "{n} components need {n} mole fractions, but got {}",
                x.len()
            ),
        ));
    }
    let m_bar: f64 = x.iter().zip(components).map(|(xi, c)| xi * c.m).sum();
    if m_bar <= 0.0 {
        return Err(AzothError::invalid_input(
            "components",
            "the mixture's segment number is zero, so no segment fraction is defined".to_string(),
        ));
    }
    let segment: Vec<f64> = x
        .iter()
        .zip(components)
        .map(|(xi, c)| xi * c.m / m_bar)
        .collect();

    let mut out = Vec::with_capacity(n);
    for ci in components {
        let mut sum = 0.0;
        for (segment_l, cj) in segment.iter().zip(components) {
            if *segment_l < 1.0e-30 {
                continue;
            }
            let sigma_ij = 0.5 * (ci.sigma + cj.sigma);
            let sigma3 = ci.sigma.powi(3) * cj.sigma.powi(3);
            let eps_ij = (ci.epsik * cj.epsik).sqrt() * sigma3.sqrt() / sigma_ij.powi(3);
            let lambda_r = 3.0 + ((ci.lambda_r - 3.0) * (cj.lambda_r - 3.0)).sqrt();
            let lambda_a = 3.0 + ((ci.lambda_a - 3.0) * (cj.lambda_a - 3.0)).sqrt();
            let cross = MieComponent {
                m: 1.0,
                lambda_r,
                lambda_a,
                sigma: sigma_ij,
                epsik: eps_ij,
            };
            let d_ij = cross.d(t)?;
            let c_mie = mie_prefactor(lambda_r, lambda_a);
            let x0 = if d_ij > 0.0 { sigma_ij / d_ij } else { 1.0 };
            let beta = eps_ij / t;
            let zeta = eta * x0 * x0 * x0;
            sum += segment_l
                * (a1_mie(eta, lambda_r, lambda_a, beta, c_mie, x0)
                    + a2_mie(eta, zeta, lambda_r, lambda_a, beta, c_mie, x0)
                    + a3_mie(zeta, lambda_r, lambda_a, beta));
        }
        out.push(sum);
    }
    Ok(out)
}
