//! `hydraulics.crane_k_factors` - fitting losses by the equivalent-length
//! method.
//!
//! ```text
//! K = f_t * sum(n_ld for each fitting)
//! ```
//!
//! Spec: `specs/calcs/hydraulics/crane_k_factors.yaml`
//!
//! # This calc cannot validate its own inputs, and says so
//!
//! The *method* is standard: a fitting's resistance coefficient is its
//! equivalent length ratio times the friction factor. The **coefficients** come
//! from `data/fittings/crane_k_factors.csv`, where every row is currently an
//! estimated dummy value - a placeholder of plausible magnitude, not from Crane
//! TP-410 or any other standard. See that file's header.
//!
//! The consequence is unusual and worth stating plainly: because the
//! coefficients are placeholders, there is no correct value for this calc to be
//! checked against, and **no test in this repository can detect a wrong
//! coefficient**. The tests validate the arithmetic and the registry lookup.
//! That is why every result built from estimated rows carries an
//! [`azoth_core::WarningCode::EstimatedData`] warning: a pressure drop
//! computed from this data can be wrong by a factor of two and still look
//! entirely reasonable.

use crate::fittings::{self, Fitting};
use crate::provenance::VerifyStatus;
use crate::results::{KComponent, KFactorsResult};
use crate::spec_gen;
use azoth_core::{Result, Warning, WarningCode, apply_checks};

/// The provenance warnings for a set of fitting statuses.
///
/// Extracted from [`crane_k_factors`] so both branches are reachable from a
/// test. The embedded registry currently holds only placeholders, so the
/// cited-but-unconfirmed branch would otherwise be code that no test ever
/// executes - which is how a warning quietly stops working.
///
/// Two levels, because there are two different things to say. A placeholder is
/// not engineering data at all; a cited-but-unconfirmed value is a real
/// published figure that nobody has checked against an authoritative copy.
/// Collapsing them would either overstate the first or understate the second.
fn provenance_warnings(statuses: &[(VerifyStatus, &str)], total: usize) -> Vec<Warning> {
    let collect = |want: VerifyStatus| -> Vec<&str> {
        statuses
            .iter()
            .filter(|(status, _)| *status == want)
            .map(|(_, id)| *id)
            .collect()
    };

    let mut warnings = Vec::new();

    let estimated = collect(VerifyStatus::EstimatedDummy);
    if !estimated.is_empty() {
        warnings.push(Warning::new(
            WarningCode::EstimatedData,
            format!(
                "{} of {} fitting(s) use ESTIMATED DUMMY coefficients that are not \
                 engineering data ({}). This resistance coefficient is a placeholder \
                 and must not be used to size equipment. Populate \
                 data/fittings/crane_k_factors.csv from the primary standard and set \
                 verify_status=verified.",
                estimated.len(),
                total,
                estimated.join(", ")
            ),
        ));
    }

    let unverified = collect(VerifyStatus::Unverified);
    if !unverified.is_empty() {
        warnings.push(Warning::new(
            WarningCode::UnverifiedSource,
            format!(
                "{} of {} fitting(s) use coefficients that are cited but NOT CONFIRMED \
                 by a named verifier against an authoritative copy of the source ({}). \
                 They may be correct; nobody has checked. Treat this resistance \
                 coefficient as provisional and confirm it before sizing equipment.",
                unverified.len(),
                total,
                unverified.join(", ")
            ),
        ));
    }

    warnings
}

/// Total resistance coefficient for a list of fittings.
///
/// `f_t` is the Darcy friction factor used as the basis. Crane specifies `f_T`,
/// the fully turbulent friction factor at the nominal fitting size; using the
/// actual friction factor at the flow Reynolds number instead is a common and
/// slightly more accurate variant, and is what the `azoth pipe` CLI does.
/// Both are defensible and they give different answers, which is why the value
/// is an explicit input rather than something this function guesses.
///
/// # Errors
/// * [`azoth_core::AzothError::UnknownFitting`] for an id not in the registry. An error
///   rather than a skip: treating an unknown fitting as zero loss would
///   under-report pressure drop.
/// * [`azoth_core::AzothError::OutOfRange`] if `f_t` is not positive, or if the list is
///   empty - an empty list would silently return zero loss, which reads as "no
///   fittings" when the caller may have meant to supply some.
///
/// # Example
/// ```
/// use azoth_hydraulics::crane_k_factors;
/// use azoth_core::CalcResult; // for has_warning
///
/// let r = crane_k_factors(&["90_elbow", "gate_valve_open"], 0.018)?;
/// assert!((r.k_total - 0.684).abs() < 1e-12);
/// assert_eq!(r.components.len(), 2);
/// // The registry holds placeholder data, so the result says so.
/// assert!(r.has_warning(azoth_core::WarningCode::EstimatedData));
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn crane_k_factors(fittings: &[&str], f_t: f64) -> Result<KFactorsResult> {
    let spec = &spec_gen::CRANE_K_FACTORS_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |q| (q == "f_t").then_some(f_t),
        &mut warnings,
    )?;

    let count = fittings.len() as f64;
    apply_checks(
        spec.derived_checks(),
        |q| (q == "n_fittings").then_some(count),
        &mut warnings,
    )?;

    let mut components = Vec::with_capacity(fittings.len());
    let mut statuses: Vec<(VerifyStatus, &'static str)> = Vec::with_capacity(fittings.len());
    for id in fittings {
        let fitting: &'static Fitting = fittings::find(id)?;
        statuses.push((fitting.status, fitting.id.as_str()));
        components.push(KComponent {
            fitting_id: fitting.id.clone(),
            n_ld: fitting.n_ld,
            k: f_t * fitting.n_ld,
        });
    }

    warnings.extend(provenance_warnings(&statuses, fittings.len()));

    let k_total = components.iter().map(|c| c.k).sum();

    Ok(KFactorsResult {
        k_total,
        f_t,
        components,
        warnings,
    })
}

/// Every fitting id in the registry, for CLI completion and diagnostics.
///
/// # Errors
/// Propagates a malformed embedded registry.
pub fn known_fittings() -> Result<Vec<&'static str>> {
    Ok(fittings::registry()?
        .iter()
        .map(|f| f.id.as_str())
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn codes(statuses: &[(VerifyStatus, &str)], total: usize) -> Vec<WarningCode> {
        provenance_warnings(statuses, total)
            .iter()
            .map(|w| w.code)
            .collect()
    }

    #[test]
    fn placeholders_and_unconfirmed_values_say_different_things() {
        // The whole point of two levels: a placeholder is not data at all, a
        // cited value is data nobody has checked. A reader needs to tell them
        // apart, and only one of them is "this might be right".
        assert_eq!(
            codes(&[(VerifyStatus::EstimatedDummy, "a")], 1),
            vec![WarningCode::EstimatedData]
        );
        assert_eq!(
            codes(&[(VerifyStatus::Unverified, "a")], 1),
            vec![WarningCode::UnverifiedSource]
        );
    }

    #[test]
    fn verified_data_is_silent() {
        // The only state that raises nothing. If this ever fires, `verified`
        // has stopped meaning anything.
        assert!(codes(&[(VerifyStatus::Verified, "a")], 1).is_empty());
    }

    #[test]
    fn a_mixed_list_raises_both_codes() {
        let statuses = [
            (VerifyStatus::EstimatedDummy, "placeholder"),
            (VerifyStatus::Unverified, "cited"),
            (VerifyStatus::Verified, "checked"),
            (VerifyStatus::Unverified, "also_cited"),
        ];
        assert_eq!(
            codes(&statuses, statuses.len()),
            vec![WarningCode::EstimatedData, WarningCode::UnverifiedSource]
        );
    }

    #[test]
    fn the_count_and_the_named_ids_are_accurate() {
        let statuses = [
            (VerifyStatus::Unverified, "one"),
            (VerifyStatus::Unverified, "two"),
            (VerifyStatus::Verified, "three"),
        ];
        let warnings = provenance_warnings(&statuses, 3);
        assert_eq!(warnings.len(), 1);
        let message = &warnings[0].message;
        // "2 of 3" and both ids, since a reader has to know which rows to look
        // at rather than only that some rows are suspect.
        assert!(message.contains("2 of 3"), "{message}");
        assert!(
            message.contains("one") && message.contains("two"),
            "{message}"
        );
        assert!(
            !message.contains("three"),
            "verified rows must not be named: {message}"
        );
    }
}
