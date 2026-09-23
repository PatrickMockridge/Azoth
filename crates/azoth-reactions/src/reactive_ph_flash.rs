//! The reactive PH flash as a model: `reactions.reactive_ph_flash`.
//!
//! The class's outer loop is [`crate::ph_flash_loop::solve_enthalpy_spec`], and everything it
//! needs comes from the pieces this crate already carries: [`crate::reactive_tp_flash`] for the
//! inner flash at a trial temperature, and `azoth-eos`' `molar_enthalpy_entropy` for the state's
//! sensible enthalpy and its heat capacity. **SRK and the vapour root on every phase**, for the
//! reason the TP model gives - the oracle's systems are `SystemSrkEos` and NeqSim's reactive
//! phases are all `gas`-typed.
//!
//! # What the specification is, and where the class reads it from
//!
//! NeqSim's process-stream enthalpy is a *sensible* one and excludes the ideal-gas formation
//! enthalpies, which a reactive calculation cannot: the composition moves, so the heat that
//! moving it releases belongs in the balance. The class's constructor therefore takes the
//! caller's sensible enthalpy and adds the inventory `sum_i n_i dHf_i` **at whatever
//! composition the system happens to hold when it is constructed** - and its own test constructs
//! it after a flash, so that is the flashed composition and not the feed.
//!
//! A model has no such system to read. This takes the **thermochemical** specification directly
//! and adds nothing, which makes the caller responsible for the same convention on both sides of
//! the comparison - and makes the number the case pins unambiguous. A caller with a sensible
//! enthalpy adds the inventory at their own composition before calling, which is the class's own
//! arithmetic with the composition stated rather than assumed.

use azoth_core::units::{kelvins, pascals};
use azoth_core::warning::Warning;
use azoth_core::{AzothError, CalcResult, Result};
use azoth_eos::databank::mixture_of;
use azoth_eos::{Cubic, IdealGasModel, Mixture, RootSide, molar_enthalpy_entropy};

use crate::databank::formation_properties;
use crate::ph_flash_loop::{PhState, solve_enthalpy_spec};
use crate::reactive_tp_flash::{ReactiveTpFlashResult, reactive_tp_flash};

/// The cubic this model flashes with, from the oracle's `SystemSrkEos`.
pub const CUBIC: Cubic = Cubic::Srk;

/// What the model answers with.
///
/// The loop's own [`crate::ph_flash_loop::PhFlashOutcome`] plus the warnings every registered
/// model carries; the fields are the same four because the loop's answer *is* the model's.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactivePhFlashResult {
    /// The temperature the loop stopped at.
    pub temperature: f64,
    /// `isConverged`. Also true where the bracket closed rather than the residual.
    pub converged: bool,
    /// `getOuterIterations`: the temperature steps taken, a path quantity.
    pub outer_iterations: u32,
    /// `getTotalInnerIterations`: every inner flash's passes, summed.
    pub total_inner_iterations: u32,
    /// Caveats. This model has no range check of its own, so it is empty wherever the
    /// caller's state was accepted.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ReactivePhFlashResult {
    const CALC_ID: &'static str = "reactions.reactive_ph_flash";
    const FIELDS: &'static [&'static str] = &[
        "temperature",
        "converged",
        "outer_iterations",
        "total_inner_iterations",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// `reactions.reactive_ph_flash`: the temperature at which a reactive fluid's thermochemical
/// enthalpy matches a specification, at fixed pressure.
///
/// `initial_temperature` is where the search starts, and it matters: the loop is a secant, so a
/// different start is a different path to the same answer - and where the fluid is stable the
/// answer may be one the path never reaches.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a charged component, a shape disagreement, or a component
/// the databank does not carry; whatever the inner flash and the enthalpy raise.
#[allow(clippy::too_many_arguments)]
pub fn reactive_ph_flash(
    components: &[String],
    initial_temperature: f64,
    pressure: f64,
    moles: &[f64],
    thermochemical_enthalpy: f64,
    max_phases: usize,
) -> Result<ReactivePhFlashResult> {
    if components.len() != moles.len() {
        return Err(AzothError::InvalidInput {
            field: "moles".to_string(),
            reason: format!(
                "{} entry(ies) against {} component(s)",
                moles.len(),
                components.len()
            ),
        });
    }
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, ideal) = mixture_of(&names, CUBIC, None)?;

    // The formation column, which the thermochemical convention is built on. Read here rather
    // than beside the enthalpy because it is the caller's number that carries it.
    let mut formation = Vec::with_capacity(components.len());
    for component in components {
        let properties =
            formation_properties(component)?.ok_or_else(|| AzothError::InvalidInput {
                field: "components".to_string(),
                reason: format!(
                    "the component databank carries no formation row for `{component}`"
                ),
            })?;
        formation.push(properties.enthalpy_of_formation);
    }

    let components = components.to_vec();
    let moles = moles.to_vec();
    let mut inner = |temperature: f64| -> Result<PhState> {
        let outcome = reactive_tp_flash(&components, temperature, pressure, &moles, max_phases)?;
        let (enthalpy, cp) = state_enthalpy(
            &mixture,
            &ideal,
            &outcome,
            temperature,
            pressure,
            &formation,
        )?;
        Ok(PhState {
            iterations: outcome.total_iterations,
            thermochemical_enthalpy: enthalpy,
            cp,
        })
    };

    let outcome = solve_enthalpy_spec(initial_temperature, thermochemical_enthalpy, &mut inner)?;
    Ok(ReactivePhFlashResult {
        temperature: outcome.temperature,
        converged: outcome.converged,
        outer_iterations: outcome.outer_iterations,
        total_inner_iterations: outcome.total_inner_iterations,
        warnings: Vec::new(),
    })
}

/// One flashed state's thermochemical enthalpy and its heat capacity.
///
/// Each phase's molar enthalpy is taken at its own composition and the cubic's vapour root, and
/// weighted by its moles; the formation inventory is added over the same rows. The class's own
/// `getFormationEnthalpyInventory` and `getThermochemicalEnthalpy` are this, one phase object at
/// a time.
fn state_enthalpy(
    mixture: &Mixture,
    ideal: &IdealGasModel,
    outcome: &ReactiveTpFlashResult,
    temperature: f64,
    pressure: f64,
    formation: &[f64],
) -> Result<(f64, f64)> {
    let reduced = mixture.reduced_parameters(kelvins(temperature), pascals(pressure))?;
    let mut sensible = 0.0_f64;
    let mut heat_capacity = 0.0_f64;
    let mut inventory = 0.0_f64;

    for row in &outcome.phase_moles {
        let total: f64 = row.iter().sum();
        if total <= 0.0 {
            continue;
        }
        let composition: Vec<f64> = row.iter().map(|moles| moles / total).collect();
        let z = mixture
            .phase_state(&reduced, &composition, RootSide::Vapour)?
            .z;
        let state = molar_enthalpy_entropy(
            mixture,
            ideal,
            kelvins(temperature),
            pascals(pressure),
            &composition,
            z,
        )?;
        sensible += total * state.h.value;
        heat_capacity += total * state.cp.value;
        for (moles, d_hf) in row.iter().zip(formation) {
            inventory += moles * d_hf;
        }
    }
    Ok((sensible + inventory, heat_capacity))
}
