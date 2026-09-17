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
use azoth_core::{Result, Warning};

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
        // No error for an unmapped label: `form_from_type` falls through to Wagner, which
        // is NeqSim's own dispatch and the reason `eos.antoine_vapor_pressure` carries the
        // same fall-through. The phase reproduces what the correlation does rather than
        // refusing a label the upstream evaluates. A non-zero `E` outranks the label,
        // which is also NeqSim's own rule.
        // The label alone does not name the form: twenty components carry DIPPR-101
        // coefficients under a `log` label, so `E` decides. See `form_from_type`.
        let form = form_from_type(&record.antoine_type, record.coefficients[4]);
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
