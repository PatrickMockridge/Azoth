//! The temperature dependence of the attraction term, as a strategy object.
//!
//! NeqSim's `AttractiveTerm*` classes own `alpha(T)` and its derivatives; this is the
//! port of that seam. A cubic's shape is in [`crate::cubic::Cubic`] and its temperature
//! dependence is here, so a new alpha correlation is a new type rather than a new branch
//! in the model layer.

/// The temperature dependence of a cubic's attraction parameter.
///
/// The three functions are the ones the model layer needs:
///
/// - `alpha(Tr)` - the attraction scale itself;
/// - `psi(Tr)` - `d ln alpha / d ln T`, the logarithmic derivative the departure
///   enthalpy and entropy are written in;
/// - `psi_t(Tr)` - `T * d(psi)/dT`, already multiplied by `T` so the departure heat
///   capacity needs no absolute temperature.
pub trait AlphaTerm {
    /// The attraction scale at reduced temperature `Tr`.
    fn alpha(&self, tr: f64) -> f64;
    /// `d ln alpha / d ln T`, the logarithmic derivative.
    fn psi(&self, tr: f64) -> f64;
    /// `T * d(psi)/dT`.
    fn psi_t(&self, tr: f64) -> f64;
}

/// Which alpha correlation a mixture's attraction term uses.
///
/// The alpha *form* is shared by the Soave variants - they differ only in the `m`
/// correlation, so each is a distinct kappa calc feeding the same [`Soave`] term. RK's
/// `1/sqrt(Tr)` is the odd one out and stays coupled to [`crate::cubic::Cubic::Rk`].
/// This is NeqSim's `attractiveTermNumber`, the axis that is independent of the cubic's
/// shape.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Alpha {
    /// Peng-Robinson's `m`, from `eos.pr_kappa`.
    #[default]
    Pr,
    /// Soave's `m`, from `eos.srk_kappa`.
    Srk,
    /// The 1978 Peng-Robinson `m`, from `eos.pr78_kappa`.
    Pr78,
    /// Twu's `m`, from `eos.twu_kappa`.
    Twu,
    /// Twu-Coon's non-Soave correlation, with `m = omega`.
    TwuCoon,
    /// Gassem et al. (2001)'s non-Soave exponential correlation.
    Gassem2001,
    /// Danesh's Soave form with `m` scaled by 1.21 above the critical temperature.
    Danesh,
    /// Schwartzentruber's correlation, three fitted parameters.
    Schwartzentruber,
    /// Mollerup's correlation, three fitted parameters.
    Mollerup,
    /// Mathias-Copeman's correlation, three fitted parameters.
    MatCop,
    /// Mathias-Copeman on Peng-Robinson, three fitted parameters.
    MatCopPr,
    /// Mathias-Copeman for UMR-PRU, three fitted parameters.
    MatCopPrUmr,
    /// Five-parameter Mathias-Copeman for UMR-CPA.
    MatCop5PrUmr,
    /// Delft (1998), methane-specific, Soave otherwise.
    Delft1998,
}

impl Alpha {
    /// The short name that crosses the Python boundary.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Alpha::Pr => "pr",
            Alpha::Srk => "srk",
            Alpha::Pr78 => "pr78",
            Alpha::Twu => "twu",
            Alpha::TwuCoon => "twucoon",
            Alpha::Gassem2001 => "gassem2001",
            Alpha::Danesh => "danesh",
            Alpha::Schwartzentruber => "schwartzentruber",
            Alpha::Mollerup => "mollerup",
            Alpha::MatCop => "matcop",
            Alpha::MatCopPr => "matcop_pr",
            Alpha::MatCopPrUmr => "matcop_prumr",
            Alpha::MatCop5PrUmr => "matcop_5prumr",
            Alpha::Delft1998 => "delft1998",
        }
    }
}

impl std::str::FromStr for Alpha {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "pr" => Ok(Alpha::Pr),
            "srk" => Ok(Alpha::Srk),
            "pr78" => Ok(Alpha::Pr78),
            "twu" => Ok(Alpha::Twu),
            "twucoon" => Ok(Alpha::TwuCoon),
            "gassem2001" => Ok(Alpha::Gassem2001),
            "danesh" => Ok(Alpha::Danesh),
            "schwartzentruber" => Ok(Alpha::Schwartzentruber),
            "mollerup" => Ok(Alpha::Mollerup),
            "matcop" => Ok(Alpha::MatCop),
            "matcop_pr" => Ok(Alpha::MatCopPr),
            "matcop_prumr" => Ok(Alpha::MatCopPrUmr),
            "matcop_5prumr" => Ok(Alpha::MatCop5PrUmr),
            "delft1998" => Ok(Alpha::Delft1998),
            other => Err(format!(
                "unknown alpha `{other}`; expected `pr`, `srk`, `pr78`, `twu`, `twucoon`, \
                 `gassem2001`, `danesh`, `schwartzentruber`, `mollerup`, `matcop`, \
                 `matcop_pr`, `matcop_prumr`, `matcop_5prumr` or `delft1998`"
            )),
        }
    }
}

/// Soave's correlation: `alpha = (1 + m(1 - sqrt(Tr)))**2`.
///
/// Peng-Robinson and Soave-Redlich-Kwong share this form and differ only in the
/// coefficient `m` - PR's `0.37464 + 1.54226 w - 0.26992 w**2`, SRK's
/// `0.48 + 1.574 w - 0.176 w**2` - so the term carries `m` as a plain number and the
/// correlation that produced it lives in the `kappa` calc.
#[derive(Debug, Clone, Copy)]
pub struct Soave {
    /// The coefficient `m`, from `eos.pr_kappa` (or `eos.srk_kappa`).
    pub kappa: f64,
}

impl AlphaTerm for Soave {
    fn alpha(&self, tr: f64) -> f64 {
        let attraction = 1.0 + self.kappa * (1.0 - tr.sqrt());
        attraction * attraction
    }

    fn psi(&self, tr: f64) -> f64 {
        let sqrt_tr = tr.sqrt();
        -self.kappa * sqrt_tr / (1.0 + self.kappa * (1.0 - sqrt_tr))
    }

    fn psi_t(&self, tr: f64) -> f64 {
        let sqrt_tr = tr.sqrt();
        -self.kappa * (1.0 + self.kappa) * tr
            / (2.0 * sqrt_tr * (1.0 + self.kappa * (1.0 - sqrt_tr)).powi(2))
    }
}

/// Danesh's correlation: Soave's form with `m` scaled by 1.21 above the critical
/// temperature. The base `m` is `eos.pr78_kappa`, the same value NeqSim's
/// `AttractiveTermPrDanesh` carries after its constructor runs.
#[derive(Debug, Clone, Copy)]
pub struct Danesh {
    /// The base coefficient `m`, from `eos.pr78_kappa`.
    pub kappa: f64,
}

impl Danesh {
    /// `1.21 m` above the critical temperature, `m` otherwise.
    fn m_mod(&self, tr: f64) -> f64 {
        if tr > 1.0 {
            1.21 * self.kappa
        } else {
            self.kappa
        }
    }
}

impl AlphaTerm for Danesh {
    fn alpha(&self, tr: f64) -> f64 {
        let m_mod = self.m_mod(tr);
        let attraction = 1.0 + m_mod * (1.0 - tr.sqrt());
        attraction * attraction
    }

    fn psi(&self, tr: f64) -> f64 {
        let m_mod = self.m_mod(tr);
        let sqrt_tr = tr.sqrt();
        -m_mod * sqrt_tr / (1.0 + m_mod * (1.0 - sqrt_tr))
    }

    fn psi_t(&self, tr: f64) -> f64 {
        let m_mod = self.m_mod(tr);
        let sqrt_tr = tr.sqrt();
        -m_mod * (1.0 + m_mod) * tr / (2.0 * sqrt_tr * (1.0 + m_mod * (1.0 - sqrt_tr)).powi(2))
    }
}

/// Redlich-Kwong's original correlation: `alpha = 1/sqrt(Tr)`.
///
/// Kappa-free, which is the whole difference from Soave: the logarithmic derivative is
/// the constant `-1/2` and its temperature derivative is zero, so neither needs a
/// coefficient.
#[derive(Debug, Clone, Copy)]
pub struct RkAlpha;

impl AlphaTerm for RkAlpha {
    fn alpha(&self, tr: f64) -> f64 {
        1.0 / tr.sqrt()
    }

    fn psi(&self, _tr: f64) -> f64 {
        -0.5
    }

    fn psi_t(&self, _tr: f64) -> f64 {
        0.0
    }
}

/// Twu-Coon's correlation, a non-Soave `alpha` with the acentric factor as `m`.
///
/// `alpha = Tr^a exp(b(1 - Tr^c)) + m (Tr^d exp(e(1 - Tr^f)) - Tr^a exp(b(1 - Tr^c)))`,
/// with six constants and `m = omega`. The logarithmic derivative and its temperature
/// derivative follow from `A = Tr^a exp(b(1 - Tr^c))` and `D = Tr^d exp(e(1 - Tr^f))`.
#[derive(Debug, Clone, Copy)]
pub struct TwuCoon {
    /// The acentric factor, Twu-Coon's `m`.
    pub omega: f64,
}

impl TwuCoon {
    /// `(alpha, d alpha/d Tr, d2 alpha/d Tr^2)` at `Tr`, the three quantities the trait
    /// methods reduce to. `A` and `D` are the two power-exponential terms; `la`/`ld`
    /// are their logarithmic derivatives and `la2`/`ld2` the derivative of those.
    fn values(&self, tr: f64) -> (f64, f64, f64) {
        let tr_c = tr.powf(2.29528);
        let tr_f = tr.powf(2.63165);
        let a_term = tr.powf(-0.201158) * (0.141599 * (1.0 - tr_c)).exp();
        let d_term = tr.powf(-0.660145) * (0.500315 * (1.0 - tr_f)).exp();

        let la = -0.201158 / tr - 0.141599 * 2.29528 * tr_c / tr;
        let ld = -0.660145 / tr - 0.500315 * 2.63165 * tr_f / tr;

        let alpha = a_term + self.omega * (d_term - a_term);
        let d_alpha = a_term * la + self.omega * (d_term * ld - a_term * la);

        let la2 = 0.201158 / (tr * tr) - 0.141599 * 2.29528 * 1.29528 * tr_c / (tr * tr);
        let ld2 = 0.660145 / (tr * tr) - 0.500315 * 2.63165 * 1.63165 * tr_f / (tr * tr);
        let d2_alpha = a_term * (la * la + la2)
            + self.omega * (d_term * (ld * ld + ld2) - a_term * (la * la + la2));

        (alpha, d_alpha, d2_alpha)
    }
}

impl AlphaTerm for TwuCoon {
    fn alpha(&self, tr: f64) -> f64 {
        self.values(tr).0
    }

    fn psi(&self, tr: f64) -> f64 {
        let (alpha, d_alpha, _) = self.values(tr);
        tr * d_alpha / alpha
    }

    fn psi_t(&self, tr: f64) -> f64 {
        let (alpha, d_alpha, d2_alpha) = self.values(tr);
        tr * d_alpha / alpha + tr * tr * (d2_alpha * alpha - d_alpha * d_alpha) / (alpha * alpha)
    }
}

/// Gassem et al. (2001)'s correlation, a non-Soave `alpha` with the acentric factor in
/// the exponent.
///
/// `alpha = exp((A + B Tr)(1 - Tr^g))` with `g = C + D w + E w**2`, five constants and
/// `w` the acentric factor. The logarithmic derivative is
/// `Tr (B (1 - Tr^g) - g Tr^(g-1) (A + B Tr))`, and its temperature derivative is
/// `B Tr - g**2 A Tr^g - B (1 + g)**2 Tr^(g+1)`.
#[derive(Debug, Clone, Copy)]
pub struct Gassem2001 {
    /// The acentric factor, Gassem's `w`.
    pub omega: f64,
}

impl Gassem2001 {
    const A: f64 = 2.0;
    const B: f64 = 0.836;
    const C: f64 = 0.134;
    const D: f64 = 0.508;
    const E: f64 = -0.0467;

    /// The exponent `g = C + D w + E w**2`, a constant per component.
    fn g(&self) -> f64 {
        Self::C + Self::D * self.omega + Self::E * self.omega * self.omega
    }
}

impl AlphaTerm for Gassem2001 {
    fn alpha(&self, tr: f64) -> f64 {
        ((Self::A + Self::B * tr) * (1.0 - tr.powf(self.g()))).exp()
    }

    fn psi(&self, tr: f64) -> f64 {
        let g = self.g();
        tr * (Self::B * (1.0 - tr.powf(g)) - g * tr.powf(g - 1.0) * (Self::A + Self::B * tr))
    }

    fn psi_t(&self, tr: f64) -> f64 {
        let g = self.g();
        let tr_g = tr.powf(g);
        Self::B * tr - g * g * Self::A * tr_g - Self::B * (1.0 + g).powi(2) * tr_g * tr
    }
}

/// The `T * d(psi)/dT` of a squared form `alpha = S**2`, from `S` and its two `Tr`
/// derivatives. `psi = 2 Tr S'/S` and `psi_t = psi + 2 Tr**2 (S S'' - S'**2)/S**2`.
fn squared_psi_t(s: f64, ds: f64, dds: f64, tr: f64) -> f64 {
    2.0 * tr * ds / s + 2.0 * tr * tr * (s * dds - ds * ds) / (s * s)
}

/// The Mathias-Copeman polynomial `S = 1 + sum_k c_k (1 - sqrt(Tr))**k` and its first
/// two `Tr` derivatives.
fn matcop_polynomial(params: &[f64], tr: f64) -> (f64, f64, f64) {
    let sqrt_tr = tr.sqrt();
    let u = 1.0 - sqrt_tr;
    let du = -0.5 / sqrt_tr;
    let ddu = 0.25 / (tr * sqrt_tr);
    let mut s = 1.0;
    let mut ds = 0.0;
    let mut dds = 0.0;
    let mut u_pow = 1.0; // u**(k-1) at the start of iteration k
    let mut u_pow_minus = 0.0; // u**(k-2); unused for k = 1
    for (i, &c) in params.iter().enumerate() {
        let k = (i + 1) as f64;
        let u_k = u_pow * u; // u**k
        s += c * u_k;
        ds += c * k * u_pow * du;
        dds += c * (k * (k - 1.0) * u_pow_minus * du * du + k * u_pow * ddu);
        u_pow_minus = u_pow;
        u_pow = u_k;
    }
    (s, ds, dds)
}

/// Mathias-Copeman's base `m` for the SRK variant: `0.48 + 1.574 w - 0.175 w**2`.
pub fn matcop_kappa(omega: f64) -> f64 {
    0.48 + 1.574 * omega - 0.175 * omega * omega
}

/// The UMR-PRU quartic `m`: `0.384401 + 1.52276 w - 0.213808 w**2 + 0.034616 w**3
/// - 0.001976 w**4`.
pub fn umr_kappa(omega: f64) -> f64 {
    0.384401 + 1.52276 * omega - 0.213808 * omega * omega + 0.034616 * omega.powi(3)
        - 0.001976 * omega.powi(4)
}

/// Schwartzentruber's correlation, a squared form with three fitted parameters.
///
/// `alpha = (1 + m(1 - sqrt(Tr)) - p0 (1 - Tr)(1 + p1 Tr + p2 Tr**2))**2`, with
/// `m = 0.48508 + 1.55191 w - 0.15613 w**2`. The `Tr > 100` high-temperature form
/// NeqSim also carries is not ported; it is an extrapolation no ordinary state reaches.
#[derive(Debug, Clone)]
pub struct Schwartzentruber {
    /// The acentric factor, for `m`.
    pub omega: f64,
    /// The three fitted parameters, from the `schwartzentruber1..3` columns.
    pub params: Vec<f64>,
}

impl Schwartzentruber {
    fn kappa(&self) -> f64 {
        0.48508 + 1.55191 * self.omega - 0.15613 * self.omega * self.omega
    }

    /// `(S, dS/dTr, d2S/dTr2)` at `Tr`, where `alpha = S**2`.
    fn polynomial(&self, tr: f64) -> (f64, f64, f64) {
        let m = self.kappa();
        let p0 = self.params.first().copied().unwrap_or(0.0);
        let p1 = self.params.get(1).copied().unwrap_or(0.0);
        let p2 = self.params.get(2).copied().unwrap_or(0.0);
        let sqrt_tr = tr.sqrt();
        let q = 1.0 + (p1 - 1.0) * tr + (p2 - p1) * tr * tr - p2 * tr * tr * tr;
        let dq = p1 - 1.0 + 2.0 * (p2 - p1) * tr - 3.0 * p2 * tr * tr;
        let ddq = 2.0 * (p2 - p1) - 6.0 * p2 * tr;
        let s = 1.0 + m * (1.0 - sqrt_tr) - p0 * q;
        let ds = -m / (2.0 * sqrt_tr) - p0 * dq;
        let dds = m / (4.0 * tr * sqrt_tr) - p0 * ddq;
        (s, ds, dds)
    }
}

impl AlphaTerm for Schwartzentruber {
    fn alpha(&self, tr: f64) -> f64 {
        let (s, _, _) = self.polynomial(tr);
        s * s
    }

    fn psi(&self, tr: f64) -> f64 {
        let (s, ds, _) = self.polynomial(tr);
        2.0 * tr * ds / s
    }

    fn psi_t(&self, tr: f64) -> f64 {
        let (s, ds, dds) = self.polynomial(tr);
        squared_psi_t(s, ds, dds, tr)
    }
}

/// Mollerup's correlation, a three-parameter non-Soave form.
///
/// `alpha = 1 + p0 (1/Tr - 1) + p1 Tr ln(Tr) + p2 (Tr - 1)`. The parameters are the
/// same `schwartzentruber1..3` columns Schwartzentruber reads.
#[derive(Debug, Clone)]
pub struct Mollerup {
    /// The three fitted parameters.
    pub params: Vec<f64>,
}

impl Mollerup {
    /// `(alpha, d alpha/dTr, d2 alpha/dTr2)` at `Tr`.
    fn values(&self, tr: f64) -> (f64, f64, f64) {
        let p0 = self.params.first().copied().unwrap_or(0.0);
        let p1 = self.params.get(1).copied().unwrap_or(0.0);
        let p2 = self.params.get(2).copied().unwrap_or(0.0);
        let ln_tr = tr.ln();
        let alpha = 1.0 + p0 * (1.0 / tr - 1.0) + p1 * tr * ln_tr + p2 * (tr - 1.0);
        let d_alpha = -p0 / (tr * tr) + p1 * (ln_tr + 1.0) + p2;
        let d2_alpha = 2.0 * p0 / (tr * tr * tr) + p1 / tr;
        (alpha, d_alpha, d2_alpha)
    }
}

impl AlphaTerm for Mollerup {
    fn alpha(&self, tr: f64) -> f64 {
        self.values(tr).0
    }

    fn psi(&self, tr: f64) -> f64 {
        let (alpha, d_alpha, _) = self.values(tr);
        tr * d_alpha / alpha
    }

    fn psi_t(&self, tr: f64) -> f64 {
        let (alpha, d_alpha, d2_alpha) = self.values(tr);
        tr * d_alpha / alpha + tr * tr * (d2_alpha * alpha - d_alpha * d_alpha) / (alpha * alpha)
    }
}

/// When Mathias-Copeman's polynomial is replaced by the base Soave alpha.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatCopFallback {
    /// Never: the polynomial is used at every temperature.
    None,
    /// Above the critical temperature, or when the first coefficient is unset.
    Supercritical,
    /// Only when the first coefficient is unset.
    Unset,
    /// Only when every coefficient is unset.
    AllUnset,
}

/// Mathias-Copeman's correlation, a squared polynomial with three or five fitted
/// parameters.
///
/// `alpha = (1 + sum_k c_k (1 - sqrt(Tr))**k)**2`. The four NeqSim variants share the
/// form and differ in the base `m`, the coefficient count, and the fallback to the base
/// Soave alpha.
#[derive(Debug, Clone)]
pub struct MatCop {
    /// The base `m`, used when the first coefficient is unset and for the Soave fallback.
    pub kappa: f64,
    /// The Mathias-Copeman coefficients `c1..cn`.
    pub params: Vec<f64>,
    /// When the base Soave alpha replaces the polynomial.
    pub fallback: MatCopFallback,
}

impl MatCop {
    /// The coefficients with `c1` replaced by `kappa` when it is effectively zero.
    fn coefficients(&self) -> Vec<f64> {
        let mut c = self.params.clone();
        if let Some(c0) = c.first_mut()
            && c0.abs() < 1e-12
        {
            *c0 = self.kappa;
        }
        c
    }

    fn use_soave(&self, c: &[f64], tr: f64) -> bool {
        match self.fallback {
            MatCopFallback::None => false,
            MatCopFallback::Supercritical => tr > 1.0 || c.first().is_none_or(|p| *p < 1e-20),
            MatCopFallback::Unset => c.first().is_none_or(|p| *p < 1e-20),
            MatCopFallback::AllUnset => c.iter().all(|p| p.abs() < 1e-20),
        }
    }
}

impl AlphaTerm for MatCop {
    fn alpha(&self, tr: f64) -> f64 {
        let c = self.coefficients();
        if self.use_soave(&c, tr) {
            return Soave { kappa: self.kappa }.alpha(tr);
        }
        let (s, _, _) = matcop_polynomial(&c, tr);
        s * s
    }

    fn psi(&self, tr: f64) -> f64 {
        let c = self.coefficients();
        if self.use_soave(&c, tr) {
            return Soave { kappa: self.kappa }.psi(tr);
        }
        let (s, ds, _) = matcop_polynomial(&c, tr);
        2.0 * tr * ds / s
    }

    fn psi_t(&self, tr: f64) -> f64 {
        let c = self.coefficients();
        if self.use_soave(&c, tr) {
            return Soave { kappa: self.kappa }.psi_t(tr);
        }
        let (s, ds, dds) = matcop_polynomial(&c, tr);
        squared_psi_t(s, ds, dds, tr)
    }
}

/// Delft (1998)'s correlation: a fitted cubic for methane, the 1978 Peng-Robinson
/// Soave form for everything else.
#[derive(Debug, Clone, Copy)]
pub struct Delft1998 {
    /// The 1978 Peng-Robinson `m`, for the non-methane branch.
    pub kappa: f64,
    /// Whether this component is methane, whose fitted cubic replaces the Soave form.
    pub is_methane: bool,
}

impl AlphaTerm for Delft1998 {
    fn alpha(&self, tr: f64) -> f64 {
        if self.is_methane {
            0.969617 + 0.20089 * tr - 0.3256987 * tr * tr + 0.06653 * tr * tr * tr
        } else {
            Soave { kappa: self.kappa }.alpha(tr)
        }
    }

    fn psi(&self, tr: f64) -> f64 {
        if self.is_methane {
            let alpha = 0.969617 + 0.20089 * tr - 0.3256987 * tr * tr + 0.06653 * tr * tr * tr;
            let d_alpha = 0.20089 - 2.0 * 0.3256987 * tr + 3.0 * 0.06653 * tr * tr;
            tr * d_alpha / alpha
        } else {
            Soave { kappa: self.kappa }.psi(tr)
        }
    }

    fn psi_t(&self, tr: f64) -> f64 {
        if self.is_methane {
            let alpha = 0.969617 + 0.20089 * tr - 0.3256987 * tr * tr + 0.06653 * tr * tr * tr;
            let d_alpha = 0.20089 - 2.0 * 0.3256987 * tr + 3.0 * 0.06653 * tr * tr;
            let d2_alpha = -2.0 * 0.3256987 + 6.0 * 0.06653 * tr;
            tr * d_alpha / alpha
                + tr * tr * (d2_alpha * alpha - d_alpha * d_alpha) / (alpha * alpha)
        } else {
            Soave { kappa: self.kappa }.psi_t(tr)
        }
    }
}
