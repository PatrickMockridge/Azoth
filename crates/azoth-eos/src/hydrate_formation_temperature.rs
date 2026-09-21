//! `eos.hydrate_formation_temperature` - the temperature at which a fluid's hydrate appears.
//!
//! ```text
//! solve  f_w^hydrate(T, P, fugacities) / f_w^fluid(T, P) - 1 = 0   for T
//! ```
//!
//! NeqSim's `HydrateFormationTemperatureFlash`. Each trial is a full flash of the fluid,
//! because the guests' fugacities come from the gas it produces - which is what NeqSim's
//! `setFug` copies in before every occupancy evaluation - and then the hydrate's water
//! fugacity is built from those by [`crate::hydrate`].
//!
//! # What settles the answer
//!
//! Water. The hydrate and the fluid meet where water's fugacity is the same in both, and the
//! two structures are evaluated at every trial because the **lower** coefficient is the stable
//! hydrate - so the structure is an output of the same comparison the temperature is.
//!
//! # The reference, and why the fluid's equation is an input
//!
//! The coefficient's third term is `ln(f_w^ref/f_w^fluid)`, and `f_w^ref` is the fugacity of
//! pure water on **the same cubic the fluid runs** - NeqSim builds it from a one-component
//! phase of the host's own class. It is computed here from a water-only mixture rather than
//! taken, so the model needs to know the equation; it reads it from the mixture itself.

use azoth_core::units::{Pressure, kelvins};
use azoth_core::{AzothError, Result, apply_checks};

use crate::databank;
use crate::hydrate;
use crate::mixture::{Mixture, RootSide};
use crate::model_gen;
use crate::pt_flash::pt_flash;
use crate::results::{HydrateFormationTemperatureResult, HydrateStructure, PtFlashResult};

/// The fugacity of pure water at a state on the cubic the fluid runs, in Pa.
///
/// NeqSim's `refPhase` is a one-component phase of the *host's* class - an SRK phase for an
/// SRK fluid, not `ComponentWater` - and this builds the same thing: the same cubic, water
/// alone, on the liquid root.
///
/// # Errors
/// * [`AzothError`] from the water-only mixture's own state.
fn reference_water_fugacity(mixture: &Mixture, t: f64, p: f64) -> Result<f64> {
    let (water, _) = databank::mixture_of(&["water"], mixture.cubic(), None)?;
    let reduced = water.reduced_parameters(kelvins(t), azoth_core::units::pascals(p))?;
    let state = water.phase_state(&reduced, &[1.0], RootSide::Liquid)?;
    Ok(state.ln_phi[0].exp() * p)
}

/// The guests' fugacities at a flash of the fluid, in Pa.
///
/// From the **vapour** phase's composition and coefficients, which is the phase NeqSim's
/// `setFug` reads. A single-phase feed has no vapour to read, so the feed's own composition is
/// its phase composition - and the fugacity coefficients come from whichever phase the flash
/// says it is.
fn guest_fugacities(flash: &PtFlashResult, z: &[f64], p: f64) -> Vec<f64> {
    let composition = match flash.phase {
        crate::results::Phase::TwoPhase => flash.y.clone(),
        _ => z.to_vec(),
    };
    let ln_phi = match flash.phase {
        crate::results::Phase::AllLiquid => &flash.ln_phi_liquid,
        _ => &flash.ln_phi_vapour,
    };
    composition
        .iter()
        .zip(ln_phi)
        .map(|(fraction, coefficient)| fraction * coefficient.exp() * p)
        .collect()
}

/// The temperature at which a fluid's hydrate appears, at a pressure.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the mixture was not resolved for a hydrate - build it
///   with [`hydrate::hydrate_mixture_of`], which is what attaches the guests' tables.
/// * [`AzothError::OutOfRange`] if `P` is not positive, or if a trial temperature's cavity sum
///   has no value.
/// * [`AzothError::SolverNotConverged`] if the scan finds no sign change or the bisection hits
///   its cap.
pub fn hydrate_formation_temperature(
    mixture: &Mixture,
    p: Pressure,
    z: &[f64],
) -> Result<HydrateFormationTemperatureResult> {
    let spec = &model_gen::HYDRATE_FORMATION_TEMPERATURE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| (quantity == "P").then_some(p.value),
        &mut warnings,
    )?;

    let hydration = mixture.hydration().ok_or_else(|| {
        AzothError::invalid_input(
            "components",
            "this mixture was not resolved for a hydrate; `hydrate_mixture_of` is what \
             attaches the guests' tables"
                .to_string(),
        )
    })?;
    let water_index = hydration
        .water_index
        .expect("the resolver refuses no water");

    let algorithm = crate::algorithm_of(spec)?;
    let bracket = algorithm.bracket.ok_or_else(|| {
        AzothError::invalid_input(
            "algorithm",
            "this procedure needs a `[algorithm.bracket]` to scan for the sign change, and \
             the spec declares none"
                .to_string(),
        )
    })?;

    let residual = |t: f64| -> Result<f64> {
        let flash = pt_flash(mixture, kelvins(t), p, z)?;
        let fugacities = guest_fugacities(&flash, z, p.value);
        let reference = reference_water_fugacity(mixture, t, p.value)?;
        let (_structure, coefficient) = hydrate::stable_structure(
            &hydration.guests,
            &fugacities,
            hydration.model,
            t,
            p.value,
            reference,
        )?;
        Ok(coefficient * p.value / fugacities[water_index] - 1.0)
    };

    // The scan, then the bisection: a hydrate's formation temperature is monotone in pressure
    // at fixed composition, and a linear scan is what the spec declares.
    let steps = bracket.steps.max(2);
    let mut iterations = 0;
    let mut lower: Option<(f64, f64)> = None;
    let mut previous: Option<(f64, f64)> = None;
    let mut bracket_found = None;
    for step in 0..steps {
        let fraction = step as f64 / (steps - 1) as f64;
        let t = bracket.lower + (bracket.upper - bracket.lower) * fraction;
        let value = residual(t)?;
        iterations += 1;
        if let Some((previous_t, previous_value)) = previous
            && previous_value * value <= 0.0
        {
            bracket_found = Some((previous_t, t));
            break;
        }
        previous = Some((t, value));
        lower = Some((t, value));
    }
    let (mut lo, mut hi) = bracket_found.ok_or(AzothError::SolverNotConverged {
        iterations,
        residual: lower.map(|(_, value)| value).unwrap_or(f64::NAN),
        tolerance: algorithm.tolerance,
    })?;
    let mut lo_value = residual(lo)?;
    iterations += 1;

    for _ in 0..algorithm.max_iterations {
        let mid = 0.5 * (lo + hi);
        let value = residual(mid)?;
        iterations += 1;
        if value.abs() < algorithm.tolerance || hi - lo < 1.0e-9 {
            return finish(mixture, p, z, hydration, mid, value, iterations, warnings);
        }
        if lo_value * value <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
            lo_value = value;
        }
    }

    Err(AzothError::SolverNotConverged {
        iterations,
        residual: 0.5 * (lo + hi),
        tolerance: algorithm.tolerance,
    })
}

/// The answer, with the structure the same comparison chose.
#[allow(clippy::too_many_arguments)]
fn finish(
    mixture: &Mixture,
    p: Pressure,
    z: &[f64],
    hydration: &hydrate::Hydration,
    temperature: f64,
    residual: f64,
    iterations: u32,
    warnings: Vec<azoth_core::Warning>,
) -> Result<HydrateFormationTemperatureResult> {
    let flash = pt_flash(mixture, kelvins(temperature), p, z)?;
    let fugacities = guest_fugacities(&flash, z, p.value);
    let reference = reference_water_fugacity(mixture, temperature, p.value)?;
    let (structure, _coefficient) = hydrate::stable_structure(
        &hydration.guests,
        &fugacities,
        hydration.model,
        temperature,
        p.value,
        reference,
    )?;

    Ok(HydrateFormationTemperatureResult {
        temperature: kelvins(temperature),
        structure: if structure == 0 {
            HydrateStructure::StructureI
        } else {
            HydrateStructure::StructureIi
        },
        iterations,
        residual,
        warnings,
    })
}
