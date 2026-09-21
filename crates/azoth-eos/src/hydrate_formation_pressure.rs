//! `eos.hydrate_formation_pressure` - the pressure at which a fluid's hydrate appears.
//!
//! ```text
//! solve  f_w^hydrate(T, P, fugacities) / f_w^fluid(T, P) - 1 = 0   for P
//! ```
//!
//! NeqSim's `HydrateFormationPressureFlash`, and the same equilibrium
//! [`crate::hydrate_formation_temperature`] reads the other way round: both solve for the
//! state where the two water fugacities are equal, one at a fixed pressure and this one at a
//! fixed temperature. Everything the kernel does is shared - the guests' fugacities come from
//! the vapour phase of a full flash, the structure is the lower of the two coefficients, and
//! the reference water phase is the same cubic the fluid runs.
//!
//! # The divergence
//!
//! **NeqSim iterates `P <- P (f_w^hydrate/f_w^fluid)` and this brackets and bisects**, which
//! is the same divergence the formation temperature states and for the same reason: the fixed
//! point converges only where its own derivative permits, and the pressure enters the
//! residual twice - once through the hydrate's chemical potential and once through every
//! fugacity the fluid contributes - so the derivative is not one a caller can bound. The
//! equality is monotone in the pressure at a fixed composition, so a bracket exists and is
//! what the spec declares.
//!
//! Measured, the two agree on the root and differ on where to stop. At 288.15 K NeqSim lands
//! on `47.5586513856629` bara and this on `4755816.196713141` Pa, `1.03e-5` relative lower -
//! and the residual at *NeqSim's* pressure is `-1.151e-6` while at this one it is `-1.650e-9`.
//! That is its own objective at the state it stopped at (`-1.24e-06`, which it logs): it
//! scales the pressure until the ratio is one to `1e-8`, and a part in `1e-6` of the ratio is
//! a part in `1e-5` of the pressure.

use azoth_core::units::{ThermodynamicTemperature, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::databank;
use crate::hydrate;
use crate::mixture::{Mixture, RootSide};
use crate::model_gen;
use crate::pt_flash::pt_flash;
use crate::results::{HydrateFormationPressureResult, HydrateStructure, PtFlashResult};

/// The fugacity of pure water at a state on the cubic the fluid runs, in Pa.
///
/// # Errors
/// * [`AzothError`] from the water-only mixture's own state.
fn reference_water_fugacity(mixture: &Mixture, t: f64, p: f64) -> Result<f64> {
    let (water, _) = databank::mixture_of(&["water"], mixture.cubic(), None)?;
    let reduced =
        water.reduced_parameters(azoth_core::units::kelvins(t), azoth_core::units::pascals(p))?;
    let state = water.phase_state(&reduced, &[1.0], RootSide::Liquid)?;
    Ok(state.ln_phi[0].exp() * p)
}

/// The guests' fugacities at a flash of the fluid, in Pa.
///
/// From the **vapour** phase's, which is the phase NeqSim's `setFug` reads. A single-phase
/// feed has no vapour to read, so the feed's own composition is its phase composition.
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

/// The pressure at which a fluid's hydrate appears, at a temperature.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the mixture was not resolved for a hydrate, or the spec
///   declares no bracket.
/// * [`AzothError::OutOfRange`] if `T` is not positive, or a trial's cavity sum has no value.
/// * [`AzothError::SolverNotConverged`] if the scan finds no sign change or the bisection hits
///   its cap.
pub fn hydrate_formation_pressure(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    z: &[f64],
) -> Result<HydrateFormationPressureResult> {
    let spec = &model_gen::HYDRATE_FORMATION_PRESSURE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| (quantity == "T").then_some(t.value),
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

    let residual = |p: f64| -> Result<f64> {
        let flash = pt_flash(mixture, t, pascals(p), z)?;
        let fugacities = guest_fugacities(&flash, z, p);
        let reference = reference_water_fugacity(mixture, t.value, p)?;
        let (_structure, coefficient) = hydrate::stable_structure(
            &hydration.guests,
            &fugacities,
            hydration.model,
            t.value,
            p,
            reference,
        )?;
        Ok(coefficient * p / fugacities[water_index] - 1.0)
    };

    // **The scan steps in the logarithm of the pressure**, because that is the variable a
    // hydrate curve is drawn on and the one the bracket's decades are stated in. A linear
    // scan over `1e4` to `2e7` Pa would put every point of interest in its first two steps.
    let steps = bracket.steps.max(2);
    let mut iterations = 0;
    let mut last: Option<(f64, f64)> = None;
    let mut previous: Option<(f64, f64)> = None;
    let mut bracket_found = None;
    for step in 0..steps {
        let fraction = step as f64 / (steps - 1) as f64;
        let p = bracket.lower * (bracket.upper / bracket.lower).powf(fraction);
        let value = residual(p)?;
        iterations += 1;
        if let Some((previous_p, previous_value)) = previous
            && previous_value * value <= 0.0
        {
            bracket_found = Some((previous_p, p));
            break;
        }
        previous = Some((p, value));
        last = Some((p, value));
    }
    let (mut lo, mut hi) = bracket_found.ok_or(AzothError::SolverNotConverged {
        iterations,
        residual: last.map(|(_, value)| value).unwrap_or(f64::NAN),
        tolerance: algorithm.tolerance,
    })?;
    let mut lo_value = residual(lo)?;
    iterations += 1;

    for _ in 0..algorithm.max_iterations {
        let mid = 0.5 * (lo + hi);
        let value = residual(mid)?;
        iterations += 1;
        if value.abs() < algorithm.tolerance || hi - lo < 1.0e-3 {
            return finish(mixture, t, z, hydration, mid, value, iterations, warnings);
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
    t: ThermodynamicTemperature,
    z: &[f64],
    hydration: &hydrate::Hydration,
    pressure: f64,
    residual: f64,
    iterations: u32,
    warnings: Vec<azoth_core::Warning>,
) -> Result<HydrateFormationPressureResult> {
    let flash = pt_flash(mixture, t, pascals(pressure), z)?;
    let fugacities = guest_fugacities(&flash, z, pressure);
    let reference = reference_water_fugacity(mixture, t.value, pressure)?;
    let (structure, _coefficient) = hydrate::stable_structure(
        &hydration.guests,
        &fugacities,
        hydration.model,
        t.value,
        pressure,
        reference,
    )?;

    Ok(HydrateFormationPressureResult {
        pressure: pascals(pressure),
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
