//! `eos.kent_eisenberg_phase` - the fugacity coefficients of a Kent-Eisenberg phase.
//!
//! Spec: `specs/models/eos/kent_eisenberg_phase.toml`. A *direct* model: no iteration, so
//! no `algorithm` block.
//!
//! NeqSim's `PhaseKentEisenberg` is a `PhaseGENRTL` whose `getActivityCoefficient` is
//! overridden to return `1.0` for **every** component, and whose components are
//! `ComponentKentEisenberg`, a `ComponentGeNRTL` that overrides `fugcoef`. So the model's
//! whole content is the branch that method takes:
//!
//! | the component | `phi_i` |
//! |---|---|
//! | a `solvent` reference state | `P0_i(T) / P` |
//! | a neutral with any other reference state | `H_i(T) / P` |
//! | an ion | `1e8`, a constant |
//!
//! **`gamma` is identically one**, so this is an ideal-solution activity model with a
//! non-ideal *reference state* - and the constant an ion gets is the model's own. NeqSim
//! gives `1e12` in `ComponentGePitzer`, `1e8` here, and `1e-15` in
//! `ComponentDesmukhMather`; the three are not variations on one number.
//!
//! # What the model does not carry
//!
//! `ComponentKentEisenberg` extends `ComponentGeNRTL`, so the NRTL parameters are loaded -
//! and unused, because the activity coefficient is overridden. Nothing here reads them.

use azoth_core::units::{kelvins, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::KentEisenbergPhaseResult;

/// The fugacity coefficient an ion gets, dimensionless: `ComponentKentEisenberg.fugcoef`.
///
/// A constant and not a correlation - an ion is treated as insoluble rather than as a
/// dissolved species, so its `phi` says "not in the vapour" rather than "poorly in it".
pub const INSOLUBLE_ION: f64 = 1.0e8;

/// The fugacity coefficients of a Kent-Eisenberg phase.
///
/// `gamma_i = 1` for every component, so `phi_i` is the component's reference state and
/// nothing else: `P0_i(T) / P` for a `solvent`, `H_i(T) / P` for a neutral solute, and
/// [`INSOLUBLE_ION`] for an ion. `P0_i` is `eos.antoine_vapor_pressure` and `H_i` is the
/// Henry correlation under `crate::henry`, so neither is recomputed here.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if a component is not in the databank, if `x` is not a
///   composition, or if a `solvent`-reference component carries no Antoine correlation.
/// * [`AzothError::PropertyUnavailable`] if a component's Henry row is the all-zero filler
///   NeqSim's database gives 296 of its 348 substances, which is an absence rather than a
///   value.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
///
/// # Example
/// ```
/// use azoth_eos::kent_eisenberg_phase::kent_eisenberg_phase;
///
/// let r = kent_eisenberg_phase(
///     &["water", "Na+", "Cl-", "CO2"],
///     313.15,
///     500_000.0,
///     &[0.89, 0.04, 0.04, 0.03],
/// )?;
/// assert_eq!(r.gamma, vec![1.0; 4]);
/// assert!((r.ln_phi[1] - 1.0e8f64.ln()).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T`, `P` and `x` are the symbols in the chemistry
pub fn kent_eisenberg_phase(
    names: &[&str],
    T: f64,
    P: f64,
    x: &[f64],
) -> Result<KentEisenbergPhaseResult> {
    let spec = &crate::model_gen::KENT_EISENBERG_PHASE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T),
            "P" => Some(P),
            _ => None,
        },
        &mut warnings,
    )?;

    let n = names.len();
    if x.len() != n {
        return Err(AzothError::invalid_input(
            "x",
            format!("a mixture of {n} components has {} mole fractions", x.len()),
        ));
    }
    if let Some(bad) = x.iter().position(|&value| value < 0.0) {
        return Err(AzothError::invalid_input(
            "x",
            format!(
                "x[{bad}] is {} but a mole fraction cannot be negative",
                x[bad]
            ),
        ));
    }
    let sum: f64 = x.iter().sum();
    if (sum - 1.0).abs() > 1.0e-09 {
        return Err(AzothError::invalid_input(
            "x",
            format!(
                "the mole fractions sum to {sum}, not to one. Renormalising them here \
                 would make a composition error invisible in every number downstream, so \
                 it is refused instead"
            ),
        ));
    }

    let mut ln_phi = Vec::with_capacity(n);
    for name in names {
        let entry = crate::databank::entry(name, None)?;
        let coefficient = if entry.reference_state == crate::databank::SOLVENT {
            vapour_pressure(&entry, T, &mut warnings)? / P
        } else if entry.ionic_charge == 0.0 {
            // `henry` reports in bar and the phase's `P` is in pascals, so the conversion
            // is here rather than in the correlation - which is NeqSim's own scale.
            crate::henry::effective_coefficient(&entry, T)? * 1.0e5 / P
        } else {
            INSOLUBLE_ION
        };
        if coefficient <= 0.0 || !coefficient.is_finite() {
            return Err(AzothError::property_unavailable(
                entry.name.clone(),
                "fugacity coefficient".to_string(),
                format!(
                    "the branch its reference state selects gives {coefficient}, which is \
                     not a coefficient. NeqSim returns it without comment"
                ),
            ));
        }
        ln_phi.push(coefficient.ln());
    }

    Ok(KentEisenbergPhaseResult {
        // **Identically one**, which is the model's defining property rather than a
        // computed column - `PhaseKentEisenberg.getActivityCoefficient` returns it.
        gamma: vec![1.0; n],
        ln_gamma: vec![0.0; n],
        ln_phi,
        warnings,
    })
}

/// `P0_i(T)` from the component's own Antoine row, in pascals.
fn vapour_pressure(
    entry: &crate::databank::Entry,
    temperature: f64,
    warnings: &mut Vec<azoth_core::Warning>,
) -> Result<f64> {
    let Some((coefficients, label)) = entry.antoine.as_ref() else {
        return Err(AzothError::property_unavailable(
            entry.name.clone(),
            "Antoine vapour-pressure coefficients".to_string(),
            "its reference state is `solvent`, so `ComponentKentEisenberg.fugcoef` gives it \
             `P0_i(T) / P` and there is no correlation to evaluate. A keycard supplies the \
             parameters a cubic reads and not these"
                .to_string(),
        ));
    };
    let Some(form) = crate::antoine_vapor_pressure::form_from_type(label, coefficients[4]) else {
        return Err(AzothError::property_unavailable(
            entry.name.clone(),
            "Antoine vapour-pressure coefficients".to_string(),
            format!(
                "its row is marked `{label}`, which upstream retracted - the 93 rows that \
                 carried one repeated filler tuple. NeqSim evaluates the filler and returns \
                 a number; this library refuses it"
            ),
        ));
    };
    let pressure = crate::antoine_vapor_pressure::antoine_vapor_pressure(
        coefficients[0],
        coefficients[1],
        coefficients[2],
        coefficients[3],
        coefficients[4],
        form,
        kelvins(entry.tc),
        pascals(entry.pc),
        kelvins(temperature),
    )?;
    warnings.extend(pressure.warnings);
    Ok(pressure.p_sat.value)
}
