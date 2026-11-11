//! `eos.pitzer_hydrate_formation_pressure` - the pressure a Pitzer brine's hydrate appears at.

use azoth_core::units::{ThermodynamicTemperature, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::cubic::Cubic;
use crate::hydrate::{self, HydrateModel};
use crate::model_gen;
use crate::pitzer_hydrate::{dataset_of, henry_floor, requested_structure, saturate, structure_at};
use crate::results::PitzerHydrateFormationPressureResult;

/// The class's default pressure window, Pa.
const MINIMUM_PRESSURE: f64 = 1.0e5;
const MAXIMUM_PRESSURE: f64 = 1000.0e5;
/// The class's temperature window, K.
const MINIMUM_TEMPERATURE: f64 = 273.15;
const MAXIMUM_TEMPERATURE: f64 = 323.15;

/// The pressure at which a Pitzer brine's hydrate appears, at a temperature.
///
/// # Errors
/// As `eos.pitzer_hydrate_formation_temperature`, with the temperature window's refusal in
/// place of the pressure window's.
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn pitzer_hydrate_formation_pressure(
    names: &[&str],
    cubic: Cubic,
    moles: &[f64],
    t: ThermodynamicTemperature,
    structure: Option<&str>,
    minimum_pressure: Option<f64>,
    maximum_pressure: Option<f64>,
    pitzer_dataset: Option<&str>,
) -> Result<PitzerHydrateFormationPressureResult> {
    let spec = &model_gen::PITZER_HYDRATE_FORMATION_PRESSURE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "minimum_pressure" => minimum_pressure,
            "maximum_pressure" => maximum_pressure,
            _ => None,
        },
        &mut warnings,
    )?;
    let requested = requested_structure(structure.unwrap_or("auto"))?;
    let dataset = dataset_of(pitzer_dataset)?;
    let algorithm = crate::algorithm_of(spec)?;
    if !(MINIMUM_TEMPERATURE..=MAXIMUM_TEMPERATURE).contains(&t.value) {
        return Err(AzothError::out_of_range(
            "T",
            t.value,
            "the Pitzer hydrate calculation requires 273.15-323.15 K",
        ));
    }
    let hydration = hydrate::hydration_of(names, None, HydrateModel::Pvtsim)?;
    let floor = henry_floor(names, moles, &hydration);
    if t.value < floor {
        return Err(AzothError::out_of_range(
            "T",
            t.value,
            format!("the temperature is below the guests' Henry reference limit of {floor} K"),
        ));
    }
    let lower = minimum_pressure
        .unwrap_or(MINIMUM_PRESSURE)
        .clamp(MINIMUM_PRESSURE, MAXIMUM_PRESSURE);
    let upper = maximum_pressure
        .unwrap_or(MAXIMUM_PRESSURE)
        .clamp(MINIMUM_PRESSURE, MAXIMUM_PRESSURE);

    let (pressure, value, iterations) = saturate(
        names,
        cubic,
        moles,
        &hydration,
        t.value,
        lower,
        upper,
        true,
        requested,
        dataset,
        algorithm.tolerance,
        algorithm.max_iterations as usize,
    )?;
    Ok(PitzerHydrateFormationPressureResult {
        pressure: pascals(pressure),
        structure: structure_at(value.structure),
        water_activity: value.water_activity,
        iterations,
        residual: value.residual,
        warnings,
    })
}
