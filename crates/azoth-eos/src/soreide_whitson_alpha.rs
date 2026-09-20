//! `eos.soreide_whitson_alpha` - the Soreide-Whitson alpha function for water.
//!
//! ```text
//! alpha = (1 + 0.453*(1 - Tr*(1 - 0.0103*salinity**1.1)) + 0.0034*((1/Tr)**3 - 1))**2
//! ```
//!
//! Spec: `specs/calcs/eos/soreide_whitson_alpha.toml`, which carries the provenance and
//! why the salinity is molality and the water-only scope.

use azoth_core::{Result, apply_checks};

use crate::results::SoreideWhitsonAlphaResult;
use crate::spec_gen;

/// The Soreide-Whitson alpha function for water.
///
/// `salinity` is the molality (mol NaCl / kg H2O); `Tr` the reduced temperature. For a
/// non-water component NeqSim's class delegates to the 1978 Peng-Robinson alpha, which
/// is ported separately as `pr78_kappa` and the standard alpha form.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tr <= 0` or `salinity < 0`.
///
/// # Example
/// ```
/// use azoth_eos::soreide_whitson_alpha;
///
/// let r = soreide_whitson_alpha(1.0, 0.7)?;
/// assert!((r.alpha - 1.3125796067429514).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn soreide_whitson_alpha(salinity: f64, Tr: f64) -> Result<SoreideWhitsonAlphaResult> {
    let spec = &spec_gen::SOREIDE_WHITSON_ALPHA_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "salinity" => Some(salinity),
            "Tr" => Some(Tr),
            _ => None,
        },
        &mut warnings,
    )?;

    // Through `bracket`, which the two derivatives beside it also read: the alpha and its
    // derivatives are one expression, and three copies of it would be three chances to
    // disagree.
    let alpha = bracket(salinity, Tr).powi(2);

    apply_checks(
        spec.derived_checks(),
        |name| (name == "alpha").then_some(alpha),
        &mut warnings,
    )?;

    Ok(SoreideWhitsonAlphaResult { alpha, warnings })
}

/// `A(Tr)`, the bracket the alpha squares: `AttractiveTermSoreideWhitson.alpha`'s inner term.
///
/// Factored out because both derivatives need it and the three must agree by construction
/// rather than by three copies of the same expression.
#[must_use]
pub fn bracket(salinity: f64, reduced_temperature: f64) -> f64 {
    1.0 + 0.453 * (1.0 - reduced_temperature * (1.0 - 0.0103 * salinity.powf(1.1)))
        + 0.0034 * ((1.0 / reduced_temperature).powi(3) - 1.0)
}

/// `d alpha / dT` in **kelvin**, from `AttractiveTermSoreideWhitson.diffalphaT`.
///
/// `alpha = A(Tr)^2` and `Tr = T / Tc`, so the chain rule carries one `1/Tc`:
///
/// ```text
/// dA/dTr    = -0.453 (1 - 0.0103 s^1.1) - 3 (0.0034) Tr^-4
/// d alpha/dT = 2 A (dA/dTr) / Tc
/// ```
///
/// **`Tc` is an argument and not a constant**, which is why this is not a registered
/// calculation: the one beside it takes the reduced temperature and has no use for `Tc`,
/// while a derivative in kelvin cannot drop it.
#[must_use]
pub fn diff_alpha_t(salinity: f64, reduced_temperature: f64, critical_temperature: f64) -> f64 {
    let slope = -0.453 * (1.0 - 0.0103 * salinity.powf(1.1))
        - 3.0 * 0.0034 * (1.0 / reduced_temperature).powi(4);
    2.0 * bracket(salinity, reduced_temperature) * slope / critical_temperature
}

/// `d^2 alpha / dT^2` in kelvin squared, from `AttractiveTermSoreideWhitson.diffdiffalphaT`.
///
/// ```text
/// d^2A/dTr^2    = 12 (0.0034) Tr^-5
/// d^2 alpha/dT^2 = 2 (dA/dTr)^2 / Tc^2 + 2 A (d^2A/dTr^2) / Tc^2
/// ```
#[must_use]
pub fn diff2_alpha_t(salinity: f64, reduced_temperature: f64, critical_temperature: f64) -> f64 {
    let slope = -0.453 * (1.0 - 0.0103 * salinity.powf(1.1))
        - 3.0 * 0.0034 * (1.0 / reduced_temperature).powi(4);
    let curvature = 12.0 * 0.0034 * (1.0 / reduced_temperature).powi(5);
    let scale = critical_temperature * critical_temperature;
    (2.0 * slope * slope + 2.0 * bracket(salinity, reduced_temperature) * curvature) / scale
}

#[cfg(test)]
mod derivative_tests {
    use super::*;

    /// **The oracle's own table**, from `AttractiveTermSoreideWhitson` on a bare
    /// `ComponentSoreideWhitson` with the salinity set by reflection.
    ///
    /// It has to be set that way: `SystemSoreideWhitson.calcSalinity` never reaches
    /// `setSalinityFromPhase`, so no phase can drive these at a non-zero brine. That defect
    /// is reported upstream; this test drives the term directly, which is the only way to
    /// see what it is for.
    #[test]
    fn the_derivatives_match_neqsim() {
        let tc = 647.3;
        for (salinity, temperature, want_alpha, want_d1, want_d2) in [
            (
                0.0,
                273.15,
                1.69960412576723,
                -0.003_120_455_641_454_86,
                2.183_925_893_651_75e-5,
            ),
            (
                0.0,
                298.15,
                1.62750902693538,
                -0.002_678_843_098_775_56,
                1.418_845_932_933_90e-5,
            ),
            (
                0.0,
                573.15,
                1.10963066944357,
                -0.001_528_396_712_055_57,
                1.429_525_719_826_01e-6,
            ),
            (
                1.0,
                298.15,
                1.63299712659856,
                -0.002_664_933_273_598_47,
                1.417_847_914_657_23e-5,
            ),
            (
                1.0,
                373.15,
                1.46142955319757,
                -0.002_019_602_025_552_07,
                5.093_580_671_394_74e-6,
            ),
            (
                4.0,
                373.15,
                1.48490097993714,
                -0.001_972_604_099_844_49,
                5.037_923_497_493_47e-6,
            ),
            (
                8.0,
                573.15,
                1.19701316026062,
                -0.001_432_087_459_751_77,
                1.248_148_150_590_14e-6,
            ),
        ] {
            let tr = temperature / tc;
            let alpha = bracket(salinity, tr).powi(2);
            assert!(
                (alpha - want_alpha).abs() < 1.0e-14,
                "alpha({salinity}, {temperature}) = {alpha}, and NeqSim gives {want_alpha}"
            );
            let d1 = diff_alpha_t(salinity, tr, tc);
            assert!(
                (d1 - want_d1).abs() < 1.0e-15,
                "dalpha/dT({salinity}, {temperature}) = {d1}, and NeqSim gives {want_d1}"
            );
            let d2 = diff2_alpha_t(salinity, tr, tc);
            assert!(
                (d2 - want_d2).abs() < 1.0e-18,
                "d2alpha/dT2({salinity}, {temperature}) = {d2}, and NeqSim gives {want_d2}"
            );
        }
    }

    /// **The term's two logarithmic derivatives against the calculus**, which is what the
    /// departure surface reads them as.
    ///
    /// `psi = d ln alpha / d ln Tr` and `psi_t = Tr d(psi)/dTr`, checked against a central
    /// difference of `alpha` rather than against a second transcription of the same algebra.
    /// The step is `1e-6` and the tolerance `1e-5` relative, which is what a central
    /// difference of a smooth function buys at that step.
    #[test]
    fn the_logarithmic_derivatives_are_the_calculus() {
        use crate::alpha_term::AlphaTerm;
        for salinity in [0.0, 1.0, 4.0] {
            for tr in [0.45, 0.55, 0.75] {
                let term = crate::alpha_term::SoreideWhitsonWater {
                    salinity,
                    critical_temperature: 647.3,
                };
                let h = 1.0e-6;
                let alpha = term.alpha(tr);
                let slope = (term.alpha(tr + h) - term.alpha(tr - h)) / (2.0 * h);
                let want_psi = tr * slope / alpha;
                assert!(
                    (term.psi(tr) - want_psi).abs() < 1.0e-5 * want_psi.abs().max(1.0e-3),
                    "psi({salinity}, {tr}) = {} and the difference quotient gives {want_psi}",
                    term.psi(tr)
                );
                let second = (term.psi(tr + h) - term.psi(tr - h)) / (2.0 * h);
                let want_psi_t = tr * second;
                assert!(
                    (term.psi_t(tr) - want_psi_t).abs() < 1.0e-3 * want_psi_t.abs().max(1.0e-2),
                    "psi_t({salinity}, {tr}) = {} and the difference quotient gives \
                     {want_psi_t}",
                    term.psi_t(tr)
                );
            }
        }
    }

    /// **The derivatives are in kelvin and carry `1/Tc`, so a wrong `Tc` shows up in them
    /// and not in `alpha`.**
    ///
    /// `alpha` is a function of `Tr` alone - which is why the registered calculation takes
    /// `Tr` and no `Tc` - while each derivative carries one factor of `1/Tc`. Passing a
    /// critical temperature of one K reduces the first derivative to the reduced-temperature
    /// slope, which is the check that the factor is there at all.
    #[test]
    fn the_derivatives_carry_a_critical_temperature() {
        let (salinity, tr) = (4.0, 373.15 / 647.3);
        let slope = diff_alpha_t(salinity, tr, 1.0);
        let expected =
            -0.453 * (1.0 - 0.0103 * salinity.powf(1.1)) - 3.0 * 0.0034 * (1.0 / tr).powi(4);
        let expected = 2.0 * bracket(salinity, tr) * expected;
        assert!((slope - expected).abs() < 1.0e-15);
        assert!(
            (diff_alpha_t(salinity, tr, 647.3) - slope / 647.3).abs() < 1.0e-15,
            "a 647.3 K critical temperature divides the per-kelvin slope by 647.3"
        );
        assert!(
            (diff2_alpha_t(salinity, tr, 647.3)
                - diff2_alpha_t(salinity, tr, 1.0) / 647.3f64.powi(2))
            .abs()
                < 1.0e-18,
            "and the second by its square"
        );
    }
}
