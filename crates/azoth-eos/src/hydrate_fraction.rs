//! `eos.hydrate_fraction` - how much of a feed is hydrate at a state.
//!
//! ```text
//! f_i = z_i - beta w_i                   the fluid the hydrate leaves, as mole numbers
//! solve  beta = min_i z_i / w_i          the fraction at which a component is exhausted
//! ```
//!
//! NeqSim's `TPHydrateFlash` computes the same quantity from the same construction, and the
//! two agree: measured at 288.15 K and 100 bara its balance error is `7.3e-13` on methane and
//! `0.0` on water, and its fraction moves with temperature as this one's does. **The solve is
//! what differs** - it optimises over the water extent where this takes the bound directly.
//!
//! # The objective is flat, and the answer is where the water runs out
//!
//! Measured in the fraction, at 288.15 K and 100 bara on the probe's feed: the aqueous phase
//! pins water's fugacity, so `ln(f_w^hydrate/f_w^fluid)` sits at `-0.071993842` at a fraction
//! of zero and has moved only to `-0.069796728` by the time the fluid holds a trace - a
//! change of `2e-3` across the whole range, with the fluid's own water fugacity identical to
//! nine digits. Past that trace the aqueous phase is gone and the objective climbs to `+inf`
//! at the state where the last water molecule has left. **So it never crosses zero: the
//! hydrate is stable at every fraction up to the one where the water is exhausted, and that
//! fraction is the answer.**
//!
//! This is why the two earlier readings of the probe looked like physics and were not. Its
//! objective appears flat in the fraction - it is, and the flatness is real; what was an
//! artefact is where NeqSim's flat part *ends*, because its composition is wrong. And its
//! fraction appears temperature-independent, which it is, because its bound is the constant
//! `z_water / (46/54)`.
//!
//! # What decides zero
//!
//! The same objective measured at the **feed**, where the hydrate has taken nothing. A feed
//! above its formation temperature has no hydrate at all, and that is exactly the condition
//! [`crate::hydrate_formation_temperature`] solves for - so the two models are one statement
//! read in two directions, and the residual here is the quantity the other drives to zero.
//!
//! # The composition is the cages' own
//!
//! From both cavity types of the stable structure ([`crate::hydrate::composition`]).
//! The balance closes by construction: the fluid is `z_i - beta w_i`, which sums to `1 - beta` for any `w` that
//! sums to one, so no phase can hold matter the feed did not have.
//!
//! # The two loops
//!
//! The composition and the fluid are mutually dependent: the cages fill by the fluid's
//! fugacities, and the fluid is what the cages have already taken from. The fraction and the
//! composition are one fixed point - `beta = min_i z_i/w_i` is a function of `w`, and the
//! fluid it leaves is a function of both - so the iteration is over that pair, and each step
//! flashes a state the hydrate has already taken the water out of.

use azoth_core::units::{Pressure, ThermodynamicTemperature, kelvins};
use azoth_core::{AzothError, Result, apply_checks};

use crate::Cubic;
use crate::databank;
use crate::hydrate::{self, Hydration};
use crate::mixture::{Mixture, RootSide};
use crate::model_gen;
use crate::pt_flash::pt_flash;
use crate::results::{HydrateFractionResult, HydrateStructure, Phase, PtFlashResult};

/// The fugacity of pure water at a state on the cubic the fluid runs, in Pa.
///
/// # Errors
/// * [`AzothError`] from the water-only mixture's own state.
fn reference_water_fugacity(eos: Cubic, t: f64, p: f64) -> Result<f64> {
    let (water, _) = databank::mixture_of(&["water"], eos, None)?;
    let reduced = water.reduced_parameters(kelvins(t), azoth_core::units::pascals(p))?;
    let state = water.phase_state(&reduced, &[1.0], RootSide::Liquid)?;
    Ok(state.ln_phi[0].exp() * p)
}

/// The guests' fugacities at a flash, in Pa, from the phase NeqSim's `setFug` reads.
fn guest_fugacities(flash: &PtFlashResult, fluid: &[f64], p: f64) -> Vec<f64> {
    let composition = if flash.phase == Phase::TwoPhase {
        flash.y.clone()
    } else {
        fluid.to_vec()
    };
    let ln_phi = if flash.phase == Phase::AllLiquid {
        &flash.ln_phi_liquid
    } else {
        &flash.ln_phi_vapour
    };
    composition
        .iter()
        .zip(ln_phi)
        .map(|(fraction, coefficient)| fraction * coefficient.exp() * p)
        .collect()
}

/// The hydrate the fluid at one state builds, and the water fugacity that comes with it.
struct Trial {
    /// The hydrate's mole fractions, from the cages the fluid fills.
    composition: Vec<f64>,
    /// The stable structure, as the index [`hydrate`] uses.
    structure: usize,
    /// The guests' fugacities the cages were filled from, in Pa.
    fugacities: Vec<f64>,
    /// The hydrate's water fugacity divided by the pressure.
    coefficient: f64,
}

/// One trial: flash the fluid, fill the cages by what the flash leaves, keep the structure.
fn trial(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    p: Pressure,
    fluid: &[f64],
    hydration: &Hydration,
) -> Result<Trial> {
    let water_index = hydration
        .water_index
        .expect("the resolver refuses no water");
    let flash = pt_flash(mixture, t, p, fluid)?;
    let fugacities = guest_fugacities(&flash, fluid, p.value);
    let reference = reference_water_fugacity(mixture.cubic(), t.value, p.value)?;
    let (structure, coefficient) = hydrate::stable_structure(
        &hydration.guests,
        &fugacities,
        hydration.model,
        t.value,
        p.value,
        reference,
    )?;
    Ok(Trial {
        composition: hydrate::composition(
            &hydration.guests,
            &fugacities,
            hydration.model,
            structure,
            t.value,
            water_index,
        ),
        structure,
        fugacities,
        coefficient,
    })
}

/// `ln(f_w^hydrate / f_w^fluid)` at a trial, in Pa over Pa.
fn objective(trial: &Trial, p: f64, water_index: usize) -> f64 {
    (trial.coefficient * p / trial.fugacities[water_index]).ln()
}

/// The largest fraction a feed's own composition allows: `min_i z_i / w_i`.
///
/// A fraction past it asks some component for more than the feed holds, so this is where a
/// component is exhausted - and for a water-limited feed the water is what is.
fn maximum_fraction(z: &[f64], in_hydrate: &[f64]) -> f64 {
    z.iter()
        .zip(in_hydrate)
        .filter(|(_, taken)| **taken > 1.0e-12)
        .map(|(total, taken)| total / taken)
        .fold(f64::INFINITY, f64::min)
}

/// The fluid the hydrate leaves at one fraction, as mole fractions.
///
/// `f_i = z_i - beta w_i` is that fluid's own mole numbers, summing to `1 - beta`. A
/// component the hydrate has taken all of comes out at zero rather than negative, which is
/// the state the bound is *defined* by: it is what makes the water-free trial a state the
/// hydrate has really left rather than one it has overrun.
fn fluid_at(z: &[f64], in_hydrate: &[f64], beta: f64) -> Vec<f64> {
    let mut fluid: Vec<f64> = z
        .iter()
        .zip(in_hydrate)
        .map(|(total, taken)| (total - beta * taken).max(0.0))
        .collect();
    let sum: f64 = fluid.iter().sum();
    for fraction in &mut fluid {
        *fraction /= sum;
    }
    fluid
}

/// The fraction of a feed that is hydrate at a temperature and pressure.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the mixture was not resolved for a hydrate.
/// * [`AzothError::OutOfRange`] if a trial's cavity sum has no value.
/// * [`AzothError::SolverNotConverged`] if the fraction and the composition do not settle.
pub fn hydrate_fraction(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<HydrateFractionResult> {
    let spec = &model_gen::HYDRATE_FRACTION_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "P" => Some(p.value),
            _ => None,
        },
        &mut warnings,
    )?;
    let algorithm = crate::algorithm_of(spec)?;

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

    // **The feed's own state decides whether there is any hydrate at all.** The objective is
    // monotonically increasing in the fraction - the fluid's water fugacity falls as the
    // water leaves it - so a feed whose objective is already positive has no root anywhere
    // and no hydrate, which is the same statement as its being above its formation
    // temperature. Measured, the two models cross zero together: the residual here is
    // `-0.003202` at 293 K and `+0.010844` at 294 K, and the formation temperature is
    // `293.2277` K.
    let mut iterations = 1;
    let feed = trial(mixture, t, p, z, hydration)?;
    let residual = objective(&feed, p.value, water_index);
    if residual > 0.0 {
        return finish(mixture, t, p, z, &feed, 0.0, residual, iterations, warnings);
    }

    // **Below it the hydrate takes the water there is**, and the fraction and the cages are
    // one fixed point: the fraction is `min_i z_i/w_i` from the cages' own composition, and
    // the cages are filled from the fluid that fraction leaves. The start is the feed's own
    // fluid - a hydrate in equilibrium with nothing - and each step flashes a state with the
    // water already taken out, which is the state at the answer rather than one near it.
    let mut hydrated = feed;
    let mut beta = maximum_fraction(z, &hydrated.composition);
    for step in 2..=algorithm.max_iterations {
        if !beta.is_finite() {
            break;
        }
        let fluid = fluid_at(z, &hydrated.composition, beta);
        hydrated = trial(mixture, t, p, &fluid, hydration)?;
        iterations = step;
        let next = maximum_fraction(z, &hydrated.composition);
        if (next - beta).abs() < algorithm.tolerance {
            return finish(
                mixture, t, p, z, &hydrated, next, residual, iterations, warnings,
            );
        }
        beta = next;
    }

    Err(AzothError::SolverNotConverged {
        iterations,
        residual: beta,
        tolerance: algorithm.tolerance,
    })
}

/// The answer, with the material balance it exists to keep measured and reported.
///
/// The fluid is `z_i - beta w_i` and the hydrate's fraction is `beta w_i`, so the
/// recombination is exact by construction whatever the composition is. What is *measured* is
/// the state the flash leaves: `sum_p beta_p x_ip` over the fluid's own phases and the
/// hydrate, against the feed's `z_i`. That is the same quantity NeqSim's state misses on
/// water by `0.0874`, and it tests this library's flash as much as this model.
#[allow(clippy::too_many_arguments)]
fn finish(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    hydrated: &Trial,
    beta: f64,
    residual: f64,
    iterations: u32,
    warnings: Vec<azoth_core::Warning>,
) -> Result<HydrateFractionResult> {
    let fluid = fluid_at(z, &hydrated.composition, beta);
    let flash = pt_flash(mixture, t, p, &fluid)?;
    // `beta` is the **vapour** fraction, so the fluid's own composition is
    // `(1 - split) x + split y`. At the answer the fluid is the water-free one the hydrate
    // leaves, which flashes single phase, so this split is the general form rather than the
    // one the bound is read at.
    let split = flash.beta.unwrap_or(0.0);
    let (x, y) = match flash.phase {
        Phase::TwoPhase => (flash.x.clone(), flash.y.clone()),
        _ => (fluid.clone(), fluid.clone()),
    };
    let mut balance_error: f64 = 0.0;
    for i in 0..z.len() {
        let left = ((1.0 - split) * x[i] + split * y[i]) * (1.0 - beta);
        let total = beta * hydrated.composition[i] + left;
        balance_error = balance_error.max((total - z[i]).abs());
    }

    Ok(HydrateFractionResult {
        beta,
        structure: if hydrated.structure == 1 {
            HydrateStructure::StructureIi
        } else {
            HydrateStructure::StructureI
        },
        balance_error,
        iterations: iterations + 1,
        residual,
        warnings,
    })
}
