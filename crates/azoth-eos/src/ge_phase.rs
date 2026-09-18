//! The arithmetic every activity-coefficient phase shares.
//!
//! NeqSim's `ComponentGE.fugcoef` sets `phi_i = gamma_i P0_i / P` for a component whose
//! `REFERENCESTATETYPE` is `solvent`, and *every* GE phase inherits that method - none of
//! `PhaseGENRTL`, `PhaseGEUnifac`, `PhaseGEUniquac`, `PhaseGEWilson` or `PhaseGEVanLaarAcid`
//! overrides `fugcoef`. What differs between the phases is which activity model supplies
//! `gamma_i` and which correlation supplies `P0_i`; the composition is one line and is
//! here rather than written once per phase.
//!
//! What is *not* here is the reference-state branch. A component tagged otherwise takes a
//! Henry's-law coefficient in NeqSim, which this library does not implement, so the
//! resolvers that build the parameters refuse such a component before a phase is
//! evaluated - see [`crate::databank::SOLVENT`].

use azoth_core::units::{kelvins, pascals};
use azoth_core::{AzothError, Result, Warning};

use crate::antoine_vapor_pressure::{antoine_vapor_pressure, form_from_type};
use crate::databank::AntoineRecord;

/// The fugacity coefficients of an activity-coefficient phase, and the parts they are
/// made of.
///
/// The saturation pressures come back beside the coefficients rather than folded into
/// them, so a caller can check the correlation and the arithmetic separately: a wrong
/// `P0` and a wrong `gamma` produce the same kind of wrong `phi`, and only the parts tell
/// them apart.
#[derive(Debug, Clone, PartialEq)]
pub struct GeFugacities {
    /// `ln(gamma_i P0_i / P)` per component.
    pub ln_phi: Vec<f64>,
    /// The pure-component saturation pressure at the state's temperature, in Pa.
    pub p_sat: Vec<f64>,
    /// Caveats from the correlations.
    pub warnings: Vec<Warning>,
}

/// The pure-component saturation pressures, one evaluation per record.
///
/// Split out of [`ge_fugacities`] because not every GE phase takes its `P0` from
/// Antoine. `eos.ge_van_laar_acid_phase` takes it from
/// `eos.nitric_sulfuric_acid_vapor_pressure` for the three acids it models and from
/// Antoine only for the species it does not - so it shares the *composition* and not this,
/// and the composition is the part worth sharing.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if a correlation is evaluated outside its
///   form's range.
pub fn saturation(antoine: &[AntoineRecord], t: f64) -> Result<(Vec<f64>, Vec<Warning>)> {
    let mut warnings = Vec::new();
    let mut p_sat = Vec::with_capacity(antoine.len());
    for record in antoine {
        // No error for an unmapped *label*: `form_from_type` falls through to Wagner, which
        // is NeqSim's own dispatch and the reason `eos.antoine_vapor_pressure` carries the
        // same fall-through. The phase reproduces what the correlation does rather than
        // refusing a label the upstream evaluates. A non-zero `E` outranks the label,
        // which is also NeqSim's own rule.
        // The label alone does not name the form: twenty components carry DIPPR-101
        // coefficients under a `log` label, so `E` decides. See `form_from_type`.
        //
        // **But `none` is refused.** It is upstream's marker for an *unavailable*
        // correlation, and an activity-coefficient phase's standard state **is** `P0` - so a
        // component without one cannot be in this phase at all. Evaluating the zeros would
        // return the critical pressure, which is a number rather than an absence.
        let form =
            form_from_type(&record.antoine_type, record.coefficients[4]).ok_or_else(|| {
                AzothError::invalid_input(
                    "components",
                    "the component's databank row has no vapour-pressure correlation: its \
                     `AntoineVapPresLiqType` is `none`, upstream's marker for unavailable \
                     data. An activity-coefficient phase's standard state is the pure \
                     liquid's saturation pressure, so a component without one cannot be in \
                     this phase. Evaluating the row's zeros as Wagner would return `Pc`, \
                     which is a number rather than an absence",
                )
            })?;
        let [a, b, c, d, e] = record.coefficients;
        let saturated = antoine_vapor_pressure(
            a,
            b,
            c,
            d,
            e,
            form,
            kelvins(record.tc),
            pascals(record.pc),
            kelvins(t),
        )?;
        warnings.extend(saturated.warnings);
        p_sat.push(saturated.p_sat.value);
    }
    Ok((p_sat, warnings))
}

/// `ln phi_i = ln gamma_i + ln(P0_i / P)`, the composition itself.
///
/// One line, and the reason this module exists: every activity-coefficient phase in this
/// library is this expression, and the phases differ only in where `gamma` and `P0` come
/// from.
#[must_use]
pub fn combine(gamma: &[f64], p_sat: &[f64], p: f64) -> Vec<f64> {
    gamma
        .iter()
        .zip(p_sat)
        .map(|(&g, &p0)| g.ln() + (p0 / p).ln())
        .collect()
}

/// `phi_i = gamma_i P0_i / P`, at a state and composition.
///
/// `gamma` is the activity coefficients the phase's own model produced, in component
/// order, and `antoine` the per-component vapour-pressure records in the same order.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if a correlation is evaluated outside its
///   form's range.
/// * Propagates the correlations' own checks.
pub fn ge_fugacities(
    gamma: &[f64],
    antoine: &[AntoineRecord],
    t: f64,
    p: f64,
) -> Result<GeFugacities> {
    let (p_sat, warnings) = saturation(antoine, t)?;
    Ok(GeFugacities {
        ln_phi: combine(gamma, &p_sat, p),
        p_sat,
        warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::antoine_vapor_pressure::AntoineForm;

    /// An unavailable correlation, and the number it must not become.
    ///
    /// **The test is synthetic because the defect is not yet reachable.** azoth's vendored
    /// `COMP.csv` is the 2026-09-13 snapshot, in which no row carries the marker; upstream
    /// `83b64e5` (PR #3775) has since added it to 313 of its 389. So this constructs the row
    /// the refresh will bring and checks the refusal is in place before the data moves - which
    /// is the order the two have to happen in, because a refresh without it would replace one
    /// wrong number with another.
    #[test]
    fn an_unavailable_correlation_is_refused_rather_than_evaluated() {
        let record = AntoineRecord {
            antoine_type: "none".to_string(),
            coefficients: [0.0; 5],
            tc: 658.0,
            pc: 1.82e6,
        };

        // The label is a marker, not a form.
        assert!(form_from_type("none", 0.0).is_none());

        // And the number it would otherwise become: sent to Wagner, `nc12`'s five zeros give
        // `exp(0) * Pc`, so an unavailable correlation comes back as the **critical
        // pressure** - `1.82e6 Pa` where the real vapour pressure is about `42 Pa`.
        let as_wagner = antoine_vapor_pressure(
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            AntoineForm::Wagner,
            kelvins(658.0),
            pascals(1.82e6),
            kelvins(298.15),
        )
        .expect("the Wagner form evaluates");
        assert!(
            (as_wagner.p_sat.value - 1.82e6).abs() < 1.0,
            "the fall-through gives {} Pa, not Pc",
            as_wagner.p_sat.value
        );

        let error = saturation(&[record], 298.15)
            .expect_err("a component without a correlation cannot be in an activity phase");
        assert!(matches!(error, AzothError::InvalidInput { .. }));
    }

    /// The fall-through the refusal must *not* break: `loglog` and `log10` are NeqSim's own
    /// unmapped labels and still resolve, so the guard is on the marker and not on novelty.
    #[test]
    fn an_unmapped_label_still_falls_through_to_wagner() {
        assert_eq!(form_from_type("loglog", 0.0), Some(AntoineForm::Wagner));
        assert_eq!(form_from_type("log10", 0.0), Some(AntoineForm::Wagner));
        assert_eq!(form_from_type("log", 0.0), Some(AntoineForm::Exp));
    }
}
