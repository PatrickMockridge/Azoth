//! The arithmetic every activity-coefficient phase shares.
//!
//! NeqSim's `ComponentGE.fugcoef` sets `phi_i = gamma_i P0_i / P` for a component whose
//! `REFERENCESTATETYPE` is `solvent`, and *every* GE phase inherits that method - none of
//! `PhaseGENRTL`, `PhaseGEUnifac`, `PhaseGEUniquac`, `PhaseGEWilson` or `PhaseGEVanLaarAcid`
//! overrides `fugcoef`. What differs between the phases is which activity model supplies
//! `gamma_i` and which correlation supplies `P0_i`; the composition is one line and is
//! here rather than written once per phase.
//!
//! What is here is the reference-state branch: [`ReferenceState`] carries each component's
//! arm and [`ge_fugacities_with_references`] composes the two. A phase whose resolver has
//! not moved to it yet still refuses a solute before it is evaluated.

use azoth_core::units::{kelvins, pascals};
use azoth_core::{AzothError, Result, Warning};

use crate::antoine_vapor_pressure::{antoine_vapor_pressure, form_from_type};
use crate::databank::AntoineRecord;
use crate::henry::{HenryRecord, Insoluble};
use crate::iapws_henry_law::{Gas, iapws_henry_law};

/// One bar in pascals: the databank's Henry correlation is in bar, the phase's pressure in
/// pascals, and the ratio `H / P` is dimensionless in either.
const BAR_TO_PA: f64 = 1.0e5;

/// One component's reference state in an activity-coefficient phase.
///
/// `ComponentGE.fugcoef` branches on `REFERENCESTATETYPE`: a `solvent` component (and one the
/// IAPWS aqueous override does not take) is `gamma_i P0_i / P`; every other is
/// `(gamma_i / gamma_i^inf) H_i / P`.
#[derive(Debug, Clone, PartialEq)]
pub enum ReferenceState {
    /// `phi_i = gamma_i P0_i / P`, with the component's Antoine correlation.
    Solvent(AntoineRecord),
    /// `phi_i = (gamma_i / gamma_i^inf) H_i / P`.
    Henry(HenryState),
}

/// Where a Henry-arm component's `H` comes from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HenryState {
    /// The IAPWS G7-04 pure-water table, which overrides the databank row for its 14 gases
    /// in a phase that carries water.
    Iapws(Gas),
    /// The databank's `HenryCoef1..4` correlation, capped to [`Insoluble::HenryCoefficient`].
    Database {
        /// The four coefficients, uncapped.
        record: HenryRecord,
        /// Whether the component is an ion, which the cap takes whatever the row says.
        ion: bool,
    },
}

impl HenryState {
    /// `H(T)` in pascals, with the cap applied.
    ///
    /// The IAPWS arm fails closed to the insoluble limit outside its fitted window, which is
    /// `IapwsHenryLaw.isUsable`'s false; the database arm caps a non-finite, non-positive or
    /// over-limit value and every ion.
    ///
    /// # Errors
    /// * [`AzothError::PropertyUnavailable`] if a non-ion database row carries the four
    ///   zeros the table uses for a substance it fits no correlation for. The resolver
    ///   refuses that state before a phase is evaluated; this is the guard behind it.
    pub fn coefficient(&self, t: f64) -> Result<f64> {
        match *self {
            Self::Iapws(gas) => {
                // `IapwsHenryLaw.isUsable` is the row's own fitted window, and both of its
                // false cases fail closed to the insoluble limit rather than extrapolate.
                // The window also keeps the call inside liquid water, which is the
                // correlation's domain.
                let row = gas.row();
                if t < row.minimum_temperature || t > row.maximum_temperature {
                    return Ok(Insoluble::HenryCoefficient.value() * BAR_TO_PA);
                }
                Ok(iapws_henry_law(gas, kelvins(t))?.henry.value)
            }
            Self::Database { record, ion } => {
                let raw = crate::henry::coefficient_from_record(&record, t);
                if crate::henry::is_capped_value(raw, ion) {
                    Ok(Insoluble::HenryCoefficient.value() * BAR_TO_PA)
                } else if record.is_fitted() {
                    Ok(raw * BAR_TO_PA)
                } else {
                    Err(AzothError::property_unavailable(
                        "a component".to_string(),
                        "a Henry coefficient".to_string(),
                        "the compiled table carries `HenryCoef1..4` as zeros for it, which is \
                         the table's marker for a substance it fits no correlation for"
                            .to_string(),
                    ))
                }
            }
        }
    }
}

/// The composition a Henry arm's `gamma^inf` is evaluated at.
///
/// The first index is always the solute, and `x` runs over `indices` only.
#[derive(Debug, Clone, PartialEq)]
pub struct ReferenceComposition {
    /// The full-list indices, solute first.
    pub indices: Vec<usize>,
    /// The composition over `indices`.
    pub x: Vec<f64>,
}

/// `PhaseGE.getActivityCoefficientInfDilWater` and `getActivityCoefficientInfDil`.
///
/// With water, the reference is the two-component phase `Phase.initRefPhases` builds: the
/// solute at `1e-10` mol and water at `10.0`, so `x = [1e-11, 1 - 1e-11]`. Without water,
/// the phase's own composition is cloned with the solute at `x_k = 1e-10` and the rest
/// renormalised: `x'_k = 1e-10 / denom` and `x'_j = x_j / denom`, where
/// `denom = sum_{j != k} x_j + 1e-10`.
///
/// The no-water form drops components at zero: they contribute nothing to any of the
/// activity models, and UNIFAC divides by `x_i`.
#[must_use]
pub fn reference_composition(k: usize, x: &[f64], water: Option<usize>) -> ReferenceComposition {
    if let Some(w) = water {
        return ReferenceComposition {
            indices: vec![k, w],
            x: vec![1.0e-11, 1.0 - 1.0e-11],
        };
    }
    let others: f64 = x
        .iter()
        .enumerate()
        .filter(|(j, _)| *j != k)
        .map(|(_, value)| value)
        .sum();
    let denominator = others + 1.0e-10;
    let mut indices = vec![k];
    let mut values = vec![1.0e-10 / denominator];
    for (j, &value) in x.iter().enumerate() {
        if j != k && value > 0.0 {
            indices.push(j);
            values.push(value / denominator);
        }
    }
    ReferenceComposition { indices, x: values }
}

/// The fugacity coefficients of a GE phase, from the reference state the resolver chose.
///
/// `reference_gamma(solute, composition)` evaluates the phase's own activity model on the
/// reference composition and returns the solute's `gamma^inf`. A model that divides by
/// `x_i` (UNIFAC, UNIQUAC) must evaluate on the pair or subset, not the full list with the
/// other components zeroed.
///
/// `p_sat` and `henry` are zero where the component's arm reads none, and `gamma_inf` is one
/// where it does not divide by one - the convention `eos.pitzer_phase` already uses for the
/// two Henry quantities.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if a solvent component's Antoine correlation is evaluated
///   outside its form's range.
/// * Propagates [`HenryState::coefficient`] and the reference evaluation.
pub fn ge_fugacities_with_references(
    gamma: &[f64],
    references: &[ReferenceState],
    x: &[f64],
    water: Option<usize>,
    t: f64,
    p: f64,
    mut reference_gamma: impl FnMut(usize, &ReferenceComposition) -> Result<f64>,
) -> Result<GeFugacities> {
    let n = references.len();
    let mut ln_phi = Vec::with_capacity(n);
    let mut p_sat = Vec::with_capacity(n);
    let mut henry = Vec::with_capacity(n);
    let mut gamma_inf = Vec::with_capacity(n);
    let mut warnings = Vec::new();

    for (i, reference) in references.iter().enumerate() {
        match reference {
            ReferenceState::Solvent(record) => {
                let (value, caveats) = saturation(std::slice::from_ref(record), t)?;
                warnings.extend(caveats);
                ln_phi.push(gamma[i].ln() + (value[0] / p).ln());
                p_sat.push(value[0]);
                henry.push(0.0);
                gamma_inf.push(1.0);
            }
            ReferenceState::Henry(state) => {
                let coefficient = state.coefficient(t)?;
                let composition = reference_composition(i, x, water);
                let infinite = reference_gamma(i, &composition)?;
                ln_phi.push((gamma[i] / infinite).ln() + (coefficient / p).ln());
                p_sat.push(0.0);
                henry.push(coefficient);
                gamma_inf.push(infinite);
            }
        }
    }

    Ok(GeFugacities {
        ln_phi,
        p_sat,
        henry,
        gamma_inf,
        warnings,
    })
}

/// The fugacity coefficients of an activity-coefficient phase, and the parts they are
/// made of.
///
/// The saturation pressures come back beside the coefficients rather than folded into
/// them, so a caller can check the correlation and the arithmetic separately: a wrong
/// `P0` and a wrong `gamma` produce the same kind of wrong `phi`, and only the parts tell
/// them apart.
#[derive(Debug, Clone, PartialEq)]
pub struct GeFugacities {
    /// `ln(gamma_i P0_i / P)` per component on the solvent arm and
    /// `ln(gamma_i / gamma_i^inf) + ln(H_i / P)` on the Henry arm.
    pub ln_phi: Vec<f64>,
    /// The pure-component saturation pressure at the state's temperature, in Pa, and **zero
    /// where the component's arm reads a Henry coefficient instead**.
    pub p_sat: Vec<f64>,
    /// The Henry coefficient each component's arm read, in Pa, and **zero where the arm
    /// reads none**.
    pub henry: Vec<f64>,
    /// The infinite-dilution activity coefficient each Henry arm divided by, and **one
    /// where the arm does not divide by one**.
    pub gamma_inf: Vec<f64>,
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
    let n = gamma.len();
    Ok(GeFugacities {
        ln_phi: combine(gamma, &p_sat, p),
        p_sat,
        henry: vec![0.0; n],
        gamma_inf: vec![1.0; n],
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
