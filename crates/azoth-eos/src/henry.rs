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
//! What they share is `H(T)` and the cap, and that is what this module is. The three
//! constants an insoluble ion gets are [`Insoluble`] - three numbers, each belonging to the
//! model that chose it rather than to a rule here - and keeping them under one name is what
//! stops a `use` choosing between `1e8` and `1e-15` by accident.
//!
//! # The data behind it is sourced, and its history is worth keeping
//!
//! Of the 348 compiled rows, 269 carry no correlation and the other 79 carry 75 distinct
//! two-term sets, `exp(h0 + h1/T) * 1.802` bar. **Upstream PR #4045 replaced the four
//! placeholder sets this module used to describe** - forty rows of `900`, the constant
//! `2.92`, the `0.059565` quartet and the one CO2 set copied onto `co` and `cos` - with
//! sourced coefficients. Methane is now `723.75 bar` at 298.15 K rather than `exp(900)`
//! through the cap, and CO2 is `28.95 bar` where the copied set gave `51.65`. The cap
//! remains the model's rule for an ion, and it still catches a value that would overflow.
//!
//! **Filed upstream as [`equinor/neqsim#4044`](https://github.com/equinor/neqsim/issues/4044)**
//! and closed by PR #4045: the four sets and the `exp(900)` sentinel were the table's, and
//! the issue said what each row count was rather than what this library did about it.

use azoth_core::{AzothError, Result};

use crate::databank::{Entry, ION};

/// What a model gives an ion it has no data for, by model.
///
/// **Three numbers, and none of them a rule.** `ComponentGE` and `ComponentGePitzer` give
/// an insoluble ion the capped Henry coefficient - `1e12` bar, which is also the number
/// `exp(900)` overflows to, so the cap and the sentinel are one constant in two roles - and
/// the other two models each carry a fugacity coefficient of their own. `ComponentKentEisenberg`'s
/// `1e8` says "sparingly present" and `ComponentDesmukhMather`'s `1e-15` says "absent",
/// which are different claims about the same ion.
///
/// **They were three loose `f64`s under two names, both spelled `INSOLUBLE_ION`**, in a
/// crate where `use` decides which one a caller gets: 1e8 against 1e-15 is a factor of
/// 1e23, it changes a vapour composition rather than crashing, and nothing compared the two
/// modules that each defined the name. A model now names the one it takes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Insoluble {
    /// The Henry coefficient an ion is capped to, **in bar** - and the cap itself, which
    /// [`is_capped`] compares against. `ComponentGE.INSOLUBLE_HENRY_COEFFICIENT`.
    HenryCoefficient,
    /// `ComponentKentEisenberg.fugcoef`'s constant, **dimensionless**: the fugacity
    /// coefficient that model gives an ion.
    KentEisenberg,
    /// `ComponentDesmukhMather.fugcoef`'s constant, **dimensionless**, and the smallest of
    /// the three.
    DesmukhMather,
}

impl Insoluble {
    /// The value the model uses.
    #[must_use]
    pub fn value(self) -> f64 {
        match self {
            Self::HenryCoefficient => 1.0e12,
            Self::KentEisenberg => 1.0e8,
            Self::DesmukhMather => 1.0e-15,
        }
    }
}

/// The four coefficients of NeqSim's Henry correlation, dimensionless, in bar.
///
/// `HenryCoef1`..`HenryCoef4` of `COMP.csv`, which `Component.java:2210` reads as
/// `henryCoefParameter[0..3]`. **All four are zero on 269 of the 348 rows**, which is the
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
/// [`Insoluble::HenryCoefficient`], and refusing it would refuse 40 rows NeqSim evaluates.
#[must_use]
pub fn coefficient(entry: &Entry, t: f64) -> f64 {
    coefficient_from_record(&entry.henry, t)
}

/// [`coefficient`] from the four numbers alone.
///
/// The GE phases carry a `HenryRecord` in their resolved parameters and no [`Entry`], so
/// the correlation is expressed over the record and the entry form delegates to it.
#[must_use]
pub fn coefficient_from_record(record: &HenryRecord, t: f64) -> f64 {
    let HenryRecord { h0, h1, h2, h3 } = *record;
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
    is_capped_value(value, entry.class == ION)
}

/// [`is_capped`] from the value and the one component fact it reads.
///
/// A GE phase's resolved parameters carry the ion flag rather than the entry.
#[must_use]
pub fn is_capped_value(value: f64, ion: bool) -> bool {
    !value.is_finite() || value <= 0.0 || value > Insoluble::HenryCoefficient.value() || ion
}

/// The Henry coefficient a model takes: `ComponentGE.getEffectiveHenryCoefficient`.
///
/// The database arm of NeqSim's two. The other is the IAPWS pure-water table, which
/// `eos.iapws_henry_law` carries; the GE phases select it through
/// [`crate::ge_phase::HenryState`], and this function is the database arm alone.
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
        Insoluble::HenryCoefficient.value()
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
    /// CO2's sourced set, with the `* 0.01802 * 100` kept as written rather than folded to
    /// `1.802`, so the test fails if either factor is dropped. `H(298.15) = 28.95 bar` is
    /// what the table's coefficients give after upstream PR #4045 replaced the placeholder
    /// correlations with sourced data.
    #[test]
    fn the_co2_correlation_is_neqsims_expression() {
        let co2 = databank::entry("co2", None).expect("co2");
        let value = coefficient(&co2, 298.15);
        assert!(
            (value - 28.95).abs() < 1.0e-3,
            "H(298.15) = {value}, and the sourced coefficients give 28.95"
        );
        // Monotone over the liquid-water range, which is what a Henry constant does.
        assert!(coefficient(&co2, 323.15) > value);
        assert!(coefficient(&co2, 273.15) < value);
    }

    /// **The `exp(900)` sentinel is gone and the cap is not.**
    ///
    /// Upstream PR #4045 replaced the four placeholder sets - forty rows of `900`, the
    /// constant `2.92`, the `0.059565` quartet and the one CO2 set copied onto `co`,
    /// `cos` - with sourced two-term coefficients, so methane's row now evaluates to a
    /// finite `723.75 bar` rather than overflowing to the cap. The cap remains the model's
    /// rule for an ion, and it still catches a value that would overflow.
    #[test]
    fn the_cap_survives_the_sourced_data() {
        let methane = databank::entry("methane", None).expect("methane");
        assert!(methane.henry.is_fitted(), "the row lists a correlation");
        assert!(
            (coefficient(&methane, 298.15) - 723.75).abs() < 1.0e-6,
            "H(298.15) = {}",
            coefficient(&methane, 298.15)
        );

        // An ion is capped whatever its row says, which is the clause that made the cap a
        // model statement rather than a numerical guard.
        let sodium = databank::entry("na+", None).expect("na+");
        assert_eq!(sodium.class, ION);
        assert!(is_capped(&sodium, 1.0), "an ion is insoluble by class");
        assert!(is_capped_value(f64::INFINITY, false));
    }

    /// A substance the table fits nothing for is refused rather than evaluated.
    #[test]
    fn an_unfitted_substance_is_refused() {
        let mdea = databank::entry("mdea", None).expect("mdea");
        assert!(!mdea.henry.is_fitted());
        let error = effective_coefficient(&mdea, 298.15).expect_err("no correlation");
        assert!(error.to_string().contains("mdea"), "{error}");
        // And the zeros really would have given a number, which is why this is a refusal
        // rather than a default.
        assert!((coefficient(&mdea, 298.15) - 1.802).abs() < 1.0e-12);
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
