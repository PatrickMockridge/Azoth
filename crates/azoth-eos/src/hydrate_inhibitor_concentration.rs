//! `eos.hydrate_inhibitor_concentration` - how much inhibitor holds a hydrate temperature down.
//!
//! NeqSim's `HydrateInhibitorConcentrationFlash`. A **secant on the inhibitor's moles** whose
//! residual is `T_hydrate - T_target`: each trial adds MEG or methanol to the feed and asks
//! what the hydrate temperature is now, and the answer is the inventory that lands on the
//! target.
//!
//! # Why the feed is stated in moles and not as a composition
//!
//! The secant walks an *absolute* amount. NeqSim's `addComponent` grows the system's total, and
//! both its first three steps (`error * 0.01` moles) and its secant step
//! (`-error / dError/dC * 0.5` moles) are in moles - so a port that carried a normalised
//! composition would reproduce the equation and not the path, and its iteration count and
//! stopping point would be different numbers. The answer is a mole number, so this takes one.
//!
//! # Nothing here reads a phase type
//!
//! Unlike its `...wtFlash` sibling, which reads the aqueous phase's own composition, this one
//! only asks [`crate::hydrate_formation_temperature`] for a temperature. So it needs no
//! aqueous-phase rule, which is the other half of the family's work and is not this model's.
//!
//! # The divergence it carries
//!
//! The inner equilibrium is on **water's fugacity**, so an associating mixture and a classical
//! one give different answers: measured on NeqSim's own fluid at 100 bara, its
//! `SystemSrkCPAstatoil` needs `0.326183190896948` mol of MEG to reach 270.9 K and the same
//! composition on a plain `SystemSrkEos` needs `1.66321547941976` - the association moves
//! water's fugacity and the hydrate equilibrium follows it. This library can flash neither
//! associating water nor MEG, so the plain-cubic column is the one it reproduces and the CPA
//! one is recorded in the case's `source`.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result};

use crate::hydrate_formation_temperature::hydrate_formation_temperature;
use crate::mixture::Mixture;
use crate::results::HydrateInhibitorConcentrationResult;

/// The residual the secant stops on, in kelvin. NeqSim's own `1e-3`.
const TOLERANCE: f64 = 1.0e-3;

/// The step cap and the floor, NeqSim's own `iter < 100` and `|| iter < 3`.
const MAXIMUM_STEPS: u32 = 100;
const MINIMUM_STEPS: u32 = 3;

/// The moles the moles vector is indexed by, as a total.
fn total(moles: &[f64]) -> f64 {
    moles.iter().sum()
}

/// The mass fraction of `inhibitor` in the inhibitor-and-water pair, which is what a dosing
/// figure means and what NeqSim's own entry point reports.
fn weight_fraction(
    mixture: &Mixture,
    moles: &[f64],
    inhibitor: usize,
    water: usize,
) -> Result<f64> {
    let molar_mass = |index: usize| -> Result<f64> {
        mixture.components()[index].molar_mass.ok_or_else(|| {
            AzothError::invalid_input(
                "components",
                "a dosing fraction is a mass fraction of the inhibitor and the water, and \
                     one of them carries no molar mass"
                    .to_string(),
            )
        })
    };
    let inhibitor_mass = moles[inhibitor] * molar_mass(inhibitor)?;
    let water_mass = moles[water] * molar_mass(water)?;
    let total_mass = inhibitor_mass + water_mass;
    if total_mass <= 0.0 {
        return Err(AzothError::invalid_input(
            "moles",
            "the feed has no inhibitor and no water, so there is no fraction to report".to_string(),
        ));
    }
    Ok(inhibitor_mass / total_mass)
}

/// The moles of inhibitor that hold a hydrate temperature down to a target.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `moles` does not match the mixture, if `inhibitor` is not
///   one of its components, if the feed has no water, or if `P` is not positive.
/// * [`AzothError`] from any trial's hydrate equilibrium - which is the model's own refusal,
///   not this loop's.
/// * [`AzothError::SolverNotConverged`] if the secant reaches its step cap without landing.
pub fn hydrate_inhibitor_concentration(
    mixture: &Mixture,
    inhibitor: &str,
    moles: &[f64],
    t_target: ThermodynamicTemperature,
    p: Pressure,
) -> Result<HydrateInhibitorConcentrationResult> {
    let spec = &crate::model_gen::HYDRATE_INHIBITOR_CONCENTRATION_SPEC;
    let mut warnings = Vec::new();
    azoth_core::apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "P" => Some(p.value),
            "T_target" => Some(t_target.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let components = mixture.components();
    if moles.len() != components.len() {
        return Err(AzothError::invalid_input(
            "moles",
            format!(
                "{} component(s) and {} mole number(s), and the secant walks them side by side",
                components.len(),
                moles.len()
            ),
        ));
    }
    let index = mixture.index_of(inhibitor).ok_or_else(|| {
        AzothError::invalid_input(
            "inhibitor",
            format!(
                "`{inhibitor}` is not one of the feed's components, so there is nothing to add \
                 to it"
            ),
        )
    })?;
    let water = mixture.index_of("water").ok_or_else(|| {
        AzothError::invalid_input(
            "components",
            "no water: an inhibitor holds a hydrate temperature down by moving water's \
             fugacity, so a feed without it has nothing to inhibit"
                .to_string(),
        )
    })?;

    let mut moles = moles.to_vec();
    let mut error = 1.0;
    let mut old_error = 1.0;
    let mut old_c = moles[index];
    let mut iterations: u32 = 0;

    // The loop breaks with the temperature it stopped at, which is what the residual is read
    // from - so there is no binding outside it to seed with a value nothing reads.
    let hydrate_temperature = loop {
        iterations += 1;
        let c = moles[index];
        // NeqSim's own secant. The first step's denominator is zero - both amounts are the
        // feed's - so its derivative is `NaN`, which the first three steps never read.
        let derivative = (error - old_error) / (c - old_c);
        old_error = error;
        old_c = c;

        let step = if iterations < MINIMUM_STEPS + 1 {
            error * 0.01
        } else {
            -error / derivative * 0.5
        };
        moles[index] += step;
        if moles[index] < 0.0 {
            // A secant may overshoot through zero, and NeqSim's `addComponent` would take the
            // amount negative and carry on with a feed that cannot exist.
            return Err(AzothError::out_of_range(
                "moles",
                moles[index],
                "the secant stepped the inhibitor's moles below zero, which is not a feed",
            ));
        }

        let sum = total(&moles);
        if sum <= 0.0 {
            return Err(AzothError::out_of_range(
                "moles",
                sum,
                "the feed has no moles left to flash",
            ));
        }
        let fractions: Vec<f64> = moles.iter().map(|value| value / sum).collect();
        let solved = hydrate_formation_temperature(mixture, p, &fractions)?;
        let reached = solved.temperature.value;
        error = reached - t_target.value;

        if !((error.abs() > TOLERANCE && iterations < MAXIMUM_STEPS) || iterations < MINIMUM_STEPS)
        {
            break reached;
        }
    };

    if error.abs() > TOLERANCE {
        return Err(AzothError::SolverNotConverged {
            iterations,
            residual: error,
            tolerance: TOLERANCE,
        });
    }

    Ok(HydrateInhibitorConcentrationResult {
        inhibitor_moles: moles[index],
        weight_fraction: weight_fraction(mixture, &moles, index, water)?,
        hydrate_temperature: azoth_core::units::kelvins(hydrate_temperature),
        iterations,
        residual: error,
        warnings,
    })
}
