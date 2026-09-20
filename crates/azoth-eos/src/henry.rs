//! The Henry reference state, as NeqSim's `ComponentGE` carries it.
//!
//! **There is no single branch to port.** Three of the tranche's models override
//! `ComponentGE.fugcoef` and disagree about the reference state, so where a model is a
//! *solute* is that model's own arithmetic:
//!
//! | NeqSim component | a neutral solute's `phi` | an ion's `phi` |
//! |---|---|---|
//! | `ComponentGE` | `(gamma / gamma_inf) H / P` | `1e12 / P` through the cap |
//! | `ComponentGePitzer` | `gamma * H * (m / x) / P` | `1e12 / P` through the cap |
//! | `ComponentKentEisenberg` | `H / P`, with `gamma = 1` | `1e8`, a constant of its own |
//! | `ComponentDesmukhMather` | `(gamma / gamma_inf) H / P` | `1e-15`, a constant of its own |
//!
//! What they share is `H(T)` and the cap, and that is what this module is. The four
//! constants an insoluble ion gets are three different numbers - `1e12`, `1e8` and `1e-15`
//! - and each belongs to the model that chose it rather than to a rule here.
//!
//! # The data behind it is thin, and that is a finding
//!
//! Of the 348 compiled rows, 296 carry no correlation at all and the other 52 carry one
//! of **four** distinct sets:
//!
//! | rows | `h0`..`h3` | what it is |
//! |---|---|---|
//! | 40 | `900, 0, 0, 0` | `exp(900)` overflows to infinity, which the cap turns into `1e12` - a sentinel, not a constant |
//! | 5 | `2.92, 0, 0, 0` | a substance-independent constant, `33.5 bar` at 298.15 K |
//! | 4 | `0, 0, 0, 0.059565` | `9.3e7 bar` at 298.15 K, shared by `h2s`, `sf6`, `r12` and `r134a` |
//! | 3 | `94.4914, -6789.04, -11.4519, -0.010454` | the CO2 set - shared by `co2`, `co` and `cos` |
//!
//! So the table fits **one** gas, and it gives that gas's correlation to carbon monoxide
//! and carbonyl sulphide as well. The port reproduces what NeqSim computes - hydrogen
//! sulphide's Henry constant really is `9.3e7 bar`, and the model that reads it should
//! return that rather than a literature value this library invented - and the thinness is
//! recorded here rather than smoothed over.

use azoth_core::{AzothError, Result};

use crate::databank::{Entry, ION};

/// The Henry coefficient a model uses for a substance with no usable correlation.
///
/// `ComponentGE.INSOLUBLE_HENRY_COEFFICIENT`, in bar. Effectively insoluble: the
/// coefficient enters as `H / P`, so `1e12 bar` at a process pressure gives a fugacity
/// coefficient of order `1e7` and a mole fraction of order `1e-7`.
pub const INSOLUBLE_HENRY_COEFFICIENT: f64 = 1.0e12;

/// The four coefficients of NeqSim's Henry correlation, dimensionless, in bar.
///
/// `HenryCoef1`..`HenryCoef4` of `COMP.csv`, which `Component.java:2210` reads as
/// `henryCoefParameter[0..3]`. **All four are zero on 296 of the 348 rows**, which is the
/// table's marker for a substance it fits no correlation for - and a zero set evaluated as
/// a polynomial gives `1.802 bar` for every one of them, so a caller must refuse rather
/// than evaluate. [`coefficient`] is that refusal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HenryRecord {
    /// The constant term.
    pub h0: f64,
    /// The `1/T` coefficient, in K.
    pub h1: f64,
    /// The `ln T` coefficient, dimensionless.
    pub h2: f64,
    /// The `T` coefficient, in 1/K.
    pub h3: f64,
}

impl HenryRecord {
    /// Whether the table fits this substance at all.
    ///
    /// `false` for the four zeros, and the caller refuses rather than evaluating: a
    /// polynomial of zeros is `1.802 bar`, which is a number rather than an absence.
    #[must_use]
    pub fn is_fitted(&self) -> bool {
        self.h0 != 0.0 || self.h1 != 0.0 || self.h2 != 0.0 || self.h3 != 0.0
    }
}

/// `H(T)`, the correlation alone, in bar.
///
/// NeqSim's own expression (`Component.java:2214`):
///
/// ```text
/// H(T) = exp(h0 + h1/T + h2 ln T + h3 T) * 0.01802 * 100
/// ```
///
/// The `0.01802` is water's molar mass in kg/mol and the `100` is NeqSim's own factor; the
/// two together are `1.802`, and this is that product rather than a re-derivation. The
/// value is **not** capped - see [`effective_coefficient`], which is what a model takes.
///
/// A `900` constant overflows `exp` to infinity, which is the table's sentinel for a
/// substance it lists and does not fit. That is not an error here: the cap turns it into
/// [`INSOLUBLE_HENRY_COEFFICIENT`], and refusing it would refuse 40 rows NeqSim evaluates.
#[must_use]
pub fn coefficient(entry: &Entry, t: f64) -> f64 {
    let HenryRecord { h0, h1, h2, h3 } = entry.henry;
    (h0 + h1 / t + h2 * t.ln() + h3 * t).exp() * 0.01802 * 100.0
}

/// `dH/dT`, in bar per kelvin.
///
/// NeqSim's `getHenryCoefdT`, which is `H(T)` times the correlation's own logarithmic
/// derivative - so a caller wanting the derivative of `ln H` divides by [`coefficient`],
/// the way `getLnHenryCoefficientTemperatureDerivative` does.
#[must_use]
pub fn coefficient_dt(entry: &Entry, t: f64) -> f64 {
    let HenryRecord { h1, h2, h3, .. } = entry.henry;
    coefficient(entry, t) * (-h1 / (t * t) + h2 / t + h3)
}

/// Whether a correlation must fail closed to the insoluble limit.
///
/// `ComponentGE.isHenryCoefficientCapped`: a value that is not finite, not positive, or
/// above the limit, **or a substance the table classes [`ION`]**. The ion clause is what
/// makes an ion insoluble whatever its row says, and it is the reason the cap is not
/// merely a numerical guard - 40 of the 52 fitted rows reach it through the `900`
/// sentinel, and 62 more through their class.
#[must_use]
pub fn is_capped(entry: &Entry, value: f64) -> bool {
    !value.is_finite() || value <= 0.0 || value > INSOLUBLE_HENRY_COEFFICIENT || entry.class == ION
}

/// The Henry coefficient a model takes: `ComponentGE.getEffectiveHenryCoefficient`.
///
/// The database arm of NeqSim's two, and the only one this library ports. The other
/// selects the IAPWS pure-water table for a supported neutral solute in a water-bearing
/// phase, which is a second correlation over a second table and is refused by name in
/// [`effective_coefficient`]'s caller rather than approximated here.
///
/// # Errors
/// * [`AzothError::PropertyUnavailable`] if the table fits the substance no correlation -
///   the four zeros on 296 rows.
pub fn effective_coefficient(entry: &Entry, t: f64) -> Result<f64> {
    if !entry.henry.is_fitted() {
        return Err(AzothError::property_unavailable(
            entry.name.clone(),
            "a Henry coefficient".to_string(),
            "the compiled table carries `HenryCoef1..4` as zeros for it, which is the \
             table's marker for a substance it fits no correlation for. Evaluating the \
             zeros as a polynomial gives 1.802 bar for every such substance, which is a \
             number rather than an absence"
                .to_string(),
        ));
    }
    let value = coefficient(entry, t);
    Ok(if is_capped(entry, value) {
        INSOLUBLE_HENRY_COEFFICIENT
    } else {
        value
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::databank;

    /// NeqSim's own arithmetic, evaluated by hand at one state.
    ///
    /// CO2's set - the only real correlation in the table - with the `* 0.01802 * 100`
    /// kept as written rather than folded to `1.802`, so the test fails if either factor
    /// is dropped. `H(298.15) = 51.6545 bar` is what NeqSim computes; the literature value
    /// for CO2 in the mole-fraction convention is about 0.5 bar, and the difference is
    /// NeqSim's, not this port's - the task is to reproduce it, and a divergence between
    /// the two is a finding about NeqSim rather than a licence to substitute another
    /// number.
    #[test]
    fn the_co2_correlation_is_neqsims_expression() {
        let co2 = databank::entry("co2", None).expect("co2");
        let value = coefficient(&co2, 298.15);
        assert!(
            (value - 51.6545).abs() < 1.0e-3,
            "H(298.15) = {value}, and NeqSim's expression gives 51.6545"
        );
        // Monotone over the liquid-water range, which is what a Henry constant does.
        assert!(coefficient(&co2, 323.15) > value);
        assert!(coefficient(&co2, 273.15) < value);
    }

    /// **The cap is load-bearing for most of the table, not an edge case.**
    ///
    /// Forty rows carry `900` as their constant. `exp(900)` overflows to infinity, and the
    /// cap is what turns that into a finite insoluble coefficient - so a caller that
    /// skipped the cap would propagate an infinity into a fugacity coefficient rather
    /// than refusing. Methane is the case, and it is a substance every model here carries.
    #[test]
    fn the_sentinel_rows_cap_rather_than_overflow() {
        let methane = databank::entry("methane", None).expect("methane");
        assert!(methane.henry.is_fitted(), "the row lists a correlation");
        assert!(
            !coefficient(&methane, 298.15).is_finite(),
            "and it overflows, which is how the sentinel works"
        );
        assert_eq!(
            effective_coefficient(&methane, 298.15).expect("capped"),
            INSOLUBLE_HENRY_COEFFICIENT
        );

        // An ion is capped whatever its row says, which is the clause that made the cap a
        // model statement rather than a numerical guard.
        let sodium = databank::entry("na+", None).expect("na+");
        assert_eq!(sodium.class, ION);
        assert!(is_capped(&sodium, 1.0), "an ion is insoluble by class");
    }

    /// A substance the table fits nothing for is refused rather than evaluated.
    #[test]
    fn an_unfitted_substance_is_refused() {
        let ethane = databank::entry("ethane", None).expect("ethane");
        assert!(!ethane.henry.is_fitted());
        let error = effective_coefficient(&ethane, 298.15).expect_err("no correlation");
        assert!(error.to_string().contains("ethane"), "{error}");
        // And the zeros really would have given a number, which is why this is a refusal
        // rather than a default.
        assert!((coefficient(&ethane, 298.15) - 1.802).abs() < 1.0e-12);
    }

    /// **The derivative is the correlation's, and the two agree by construction.**
    ///
    /// `getHenryCoefdT` is `H` times the logarithmic derivative, so a caller wanting
    /// `d(ln H)/dT` divides - and a factor dropped between the two shows up here rather
    /// than in a phase's temperature derivative three tranches later.
    #[test]
    fn the_temperature_derivative_matches_a_finite_difference() {
        let co2 = databank::entry("co2", None).expect("co2");
        for t in [273.15, 298.15, 323.15, 373.15] {
            let step = 1.0e-4;
            let numeric =
                (coefficient(&co2, t + step) - coefficient(&co2, t - step)) / (2.0 * step);
            let analytic = coefficient_dt(&co2, t);
            assert!(
                (analytic - numeric).abs() <= 1.0e-6 * numeric.abs(),
                "at {t} K: analytic {analytic} against finite difference {numeric}"
            );
        }
    }
}
