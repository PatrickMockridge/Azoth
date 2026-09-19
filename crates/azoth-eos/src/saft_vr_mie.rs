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
