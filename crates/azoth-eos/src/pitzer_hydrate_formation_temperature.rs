//! `eos.pitzer_hydrate_formation_temperature` - the temperature a Pitzer brine's hydrate appears at.

use azoth_core::units::{Pressure, kelvins};
use azoth_core::{AzothError, Result, apply_checks};

use crate::cubic::Cubic;
use crate::hydrate::{self, HydrateModel};
use crate::model_gen;
use crate::pitzer_hydrate::{dataset_of, henry_floor, requested_structure, saturate, structure_at};
use crate::results::PitzerHydrateFormationTemperatureResult;

/// The pressure at which the class's temperature window ends, Pa.
const MAXIMUM_PRESSURE: f64 = 1000.0e5;
/// The class's default temperature window, K.
const MINIMUM_TEMPERATURE: f64 = 273.15;
const MAXIMUM_TEMPERATURE: f64 = 323.15;

/// The temperature at which a Pitzer brine's hydrate appears, at a pressure.
///
/// # Errors
/// * [`AzothError::InvalidInput`] for a structure name outside the three, or a feed with no
///   water or no hydrate former.
/// * [`AzothError::OutOfRange`] for a pressure or temperature outside the class's window, or
///   a water activity outside `(0, 1]`.
/// * [`AzothError::SolverNotConverged`] if the scan finds no sign change or the bisection
///   hits its cap.
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn pitzer_hydrate_formation_temperature(
    names: &[&str],
    cubic: Cubic,
    moles: &[f64],
    p: Pressure,
    structure: Option<&str>,
    minimum_temperature: Option<f64>,
    maximum_temperature: Option<f64>,
    pitzer_dataset: Option<&str>,
) -> Result<PitzerHydrateFormationTemperatureResult> {
    let spec = &model_gen::PITZER_HYDRATE_FORMATION_TEMPERATURE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "P" => Some(p.value),
            "minimum_temperature" => minimum_temperature,
            "maximum_temperature" => maximum_temperature,
            _ => None,
        },
        &mut warnings,
    )?;
    let requested = requested_structure(structure.unwrap_or("auto"))?;
    let dataset = dataset_of(pitzer_dataset)?;
    let algorithm = crate::algorithm_of(spec)?;
    if !(p.value > 0.0 && p.value <= MAXIMUM_PRESSURE) {
        return Err(AzothError::out_of_range(
            "P",
            p.value,
            "the Pitzer hydrate calculation requires 0 < P <= 1000 bara",
        ));
    }
    let hydration = hydrate::hydration_of(names, None, HydrateModel::Pvtsim)?;
    let floor = henry_floor(names, moles, &hydration);
    let lower = minimum_temperature
        .unwrap_or(MINIMUM_TEMPERATURE)
        .clamp(MINIMUM_TEMPERATURE, MAXIMUM_TEMPERATURE)
        .max(floor);
    let upper = maximum_temperature
        .unwrap_or(MAXIMUM_TEMPERATURE)
        .clamp(MINIMUM_TEMPERATURE, MAXIMUM_TEMPERATURE);
    if upper < floor {
        return Err(AzothError::invalid_input(
            "maximum_temperature",
            format!(
                "the temperature bounds do not overlap the guests' Henry range, whose floor is \
                 {floor} K"
            ),
        ));
    }

    let (temperature, value, iterations) = saturate(
        names,
        cubic,
        moles,
        &hydration,
        p.value,
        lower,
        upper,
        false,
        requested,
        dataset,
        algorithm.tolerance,
        algorithm.max_iterations as usize,
    )?;
    Ok(PitzerHydrateFormationTemperatureResult {
        temperature: kelvins(temperature),
        structure: structure_at(value.structure),
        water_activity: value.water_activity,
        iterations,
        residual: value.residual,
        warnings,
    })
}
