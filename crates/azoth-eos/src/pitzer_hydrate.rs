//! The Pitzer-brine hydrate saturation's shared scan and bisection.
//!
//! NeqSim's `PitzerHydrateFlash`: a scan and a bisection on `ln(f_hydrate^water / f_aq^water)`
//! where the aqueous side is a `SystemPitzer`'s hybrid EoS/GE flash and the hydrate side is
//! `ComponentHydratePitzer`'s water fugacity coefficient. The inner flash is
//! [`crate::hybrid_eos_ge_flash`], the id that seam already carries, and the coefficient is
//! [`crate::hydrate::pitzer_structure`] with the pure-water Antoine reference.
//!
//! The two public models - [`crate::pitzer_hydrate_formation_temperature`] and
//! [`crate::pitzer_hydrate_formation_pressure`] - are thin wrappers over [`saturate`].

use azoth_core::{AzothError, Result};

use crate::antoine_vapor_pressure::saturation_pressure;
use crate::cubic::Cubic;
use crate::databank;
use crate::hybrid_eos_ge_flash::hybrid_eos_ge_flash_with_dataset;
use crate::hydrate::{self, Hydration};
use crate::pitzer_phase::PitzerDatasetSelection;
use crate::results::HydrateStructure;

/// The class's temperature grid step, K.
const TEMPERATURE_STEP: f64 = 5.0;
/// The pressure grid's additive step, Pa - the class's `lower + 0.5` bara.
const PRESSURE_STEP: f64 = 0.5e5;
/// The pressure grid's multiplicative step.
const PRESSURE_RATIO: f64 = 1.3;
/// The floor the Henry reference starts from, K.
const MINIMUM_TEMPERATURE: f64 = 273.15;
/// The water-activity bound the class checks after every trial.
const WATER_ACTIVITY_TOLERANCE: f64 = 1.0e-8;

/// The trial state's residual and the two diagnostics the class exposes.
pub(crate) struct Trial {
    pub(crate) residual: f64,
    pub(crate) structure: usize,
    pub(crate) water_activity: f64,
}

/// The guests' Henry-reference floor: the largest fitted minimum over the present guests.
///
/// `PitzerHydrateFlash.minimumHenryTemperature()` starts at 273.15 K and raises it for every
/// molecular guest with a positive inventory that the IAPWS table carries.
pub(crate) fn henry_floor(names: &[&str], moles: &[f64], hydration: &Hydration) -> f64 {
    let mut floor = MINIMUM_TEMPERATURE;
    for (index, name) in names.iter().enumerate() {
        if moles.get(index).copied().unwrap_or(0.0) <= 0.0 || !hydration.guests[index].former {
            continue;
        }
        if let Some(gas) = crate::iapws_henry_law::gas_from_name(name) {
            floor = floor.max(gas.row().minimum_temperature);
        }
    }
    floor
}

/// One trial: the hybrid flash, the guests' aqueous fugacities, the hydrate coefficient.
#[allow(clippy::too_many_arguments)] // The class's own list of state and controls.
pub(crate) fn trial(
    names: &[&str],
    cubic: Cubic,
    moles: &[f64],
    hydration: &Hydration,
    temperature: f64,
    pressure: f64,
    requested: Option<usize>,
    dataset: PitzerDatasetSelection,
) -> Result<Trial> {
    let flash =
        hybrid_eos_ge_flash_with_dataset(names, cubic, temperature, pressure, moles, dataset)?;
    // The active roles keep `[gas, oil, aqueous]` order and the aqueous one is last; a
    // gas/brine state has two rows, and a feed the solver refuses has none.
    let aqueous = flash.x.len().checked_sub(1).ok_or_else(|| {
        AzothError::invalid_input(
            "components",
            "the hybrid flash returned no roles, so there is no aqueous phase to read".to_string(),
        )
    })?;
    let composition = &flash.x[aqueous];
    let ln_phi = &flash.ln_phi[aqueous];
    let fugacity = |index: usize| composition[index] * ln_phi[index].exp() * pressure;

    let water_index = hydration
        .water_index
        .expect("the resolver refuses no water");
    let f_water = fugacity(water_index);
    if f_water.is_nan() || f_water <= 0.0 {
        return Err(AzothError::out_of_range(
            "components",
            f_water,
            "the aqueous phase carries no water, so there is no solvent fugacity to match",
        ));
    }
    let guests: Vec<f64> = hydration
        .guests
        .iter()
        .enumerate()
        .map(|(index, guest)| if guest.former { fugacity(index) } else { 0.0 })
        .collect();

    let mut warnings = Vec::new();
    let p_sat = saturation_pressure(&databank::entry("water", None)?, temperature, &mut warnings)?;
    let (structure, coefficient) = hydrate::pitzer_structure(
        &hydration.guests,
        &guests,
        temperature,
        pressure,
        p_sat,
        requested,
    )?;
    let water_activity = f_water / p_sat;
    if !(water_activity > 0.0 && water_activity <= 1.0 + WATER_ACTIVITY_TOLERANCE) {
        return Err(AzothError::out_of_range(
            "water_activity",
            water_activity,
            "the Pitzer aqueous phase's water activity is outside (0, 1]; downstream NeqSim \
             refuses it as a non-physical solvent reference",
        ));
    }
    Ok(Trial {
        residual: (coefficient * pressure / f_water).ln(),
        structure,
        water_activity,
    })
}

/// The scan and the bisection, shared by the two solves.
///
/// `variable` is the quantity being solved: `true` pressure at fixed temperature, `false`
/// temperature at fixed pressure.
#[allow(clippy::too_many_arguments)] // The class's own list of state, bounds and controls.
pub(crate) fn saturate(
    names: &[&str],
    cubic: Cubic,
    moles: &[f64],
    hydration: &Hydration,
    fixed: f64,
    lower: f64,
    upper: f64,
    solve_pressure: bool,
    requested: Option<usize>,
    dataset: PitzerDatasetSelection,
    tolerance: f64,
    max_bisections: usize,
) -> Result<(f64, Trial, u32)> {
    let evaluate = |value: f64| -> Result<Trial> {
        if solve_pressure {
            trial(
                names, cubic, moles, hydration, fixed, value, requested, dataset,
            )
        } else {
            trial(
                names, cubic, moles, hydration, value, fixed, requested, dataset,
            )
        }
    };
    let step = |current: f64| -> f64 {
        if solve_pressure {
            (current * PRESSURE_RATIO)
                .max(current + PRESSURE_STEP)
                .min(upper)
        } else {
            (current + TEMPERATURE_STEP).min(upper)
        }
    };

    let mut iterations = 1u32;
    let mut last = evaluate(lower)?;
    if last.residual.abs() <= tolerance {
        return Ok((lower, last, iterations));
    }
    let mut current = lower;
    let mut current_residual = last.residual;
    while current < upper {
        let next = step(current);
        let value = evaluate(next)?;
        iterations += 1;
        if value.residual.abs() <= tolerance {
            return Ok((next, value, iterations));
        }
        if current_residual.signum() != value.residual.signum() {
            return bisect(
                &evaluate,
                current,
                next,
                current_residual,
                iterations,
                tolerance,
                max_bisections,
            );
        }
        current = next;
        current_residual = value.residual;
        last = value;
    }
    Err(AzothError::SolverNotConverged {
        iterations,
        residual: last.residual,
        tolerance,
    })
}

/// The class's plain midpoint halving, capped at seventy evaluations.
pub(crate) fn bisect(
    evaluate: &impl Fn(f64) -> Result<Trial>,
    mut lower: f64,
    mut upper: f64,
    mut lower_residual: f64,
    mut iterations: u32,
    tolerance: f64,
    max_bisections: usize,
) -> Result<(f64, Trial, u32)> {
    for _ in 0..max_bisections {
        let middle = 0.5 * (lower + upper);
        let value = evaluate(middle)?;
        iterations += 1;
        if value.residual.abs() <= tolerance {
            return Ok((middle, value, iterations));
        }
        if lower_residual.signum() != value.residual.signum() {
            upper = middle;
        } else {
            lower = middle;
            lower_residual = value.residual;
        }
    }
    Err(AzothError::SolverNotConverged {
        iterations,
        residual: 0.5 * (lower + upper),
        tolerance,
    })
}

/// The structure request the spec's enum names.
/// The Pitzer dataset the spec's enum names.
pub(crate) fn dataset_of(name: Option<&str>) -> Result<PitzerDatasetSelection> {
    match name {
        None | Some("auto") => Ok(PitzerDatasetSelection::Auto),
        Some("phreeqc_co2_chloride") => Ok(PitzerDatasetSelection::PhreeqcCo2Chloride),
        Some(other) => Err(AzothError::invalid_input(
            "pitzer_dataset",
            format!("`{other}` is not one of auto, phreeqc_co2_chloride"),
        )),
    }
}

pub(crate) fn requested_structure(structure: &str) -> Result<Option<usize>> {
    match structure {
        "auto" => Ok(None),
        "si" => Ok(Some(0)),
        "sii" => Ok(Some(1)),
        other => Err(AzothError::invalid_input(
            "structure",
            format!("`{other}` is not one of auto, si, sii"),
        )),
    }
}

pub(crate) fn structure_at(index: usize) -> HydrateStructure {
    if index == 0 {
        HydrateStructure::StructureI
    } else {
        HydrateStructure::StructureIi
    }
}
