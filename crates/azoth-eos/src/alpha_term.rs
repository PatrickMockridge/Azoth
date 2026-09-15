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
            other => Err(format!(
                "unknown alpha `{other}`; expected `pr`, `srk`, `pr78`, `twu` or `twucoon`"
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
