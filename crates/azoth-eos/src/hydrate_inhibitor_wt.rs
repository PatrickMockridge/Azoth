//! `eos.hydrate_inhibitor_wt` - the inhibitor dose that reaches a target weight fraction.
//!
//! NeqSim's `HydrateInhibitorwtFlash`. A **secant on the inhibitor's moles**, as its
//! `...ConcentrationFlash` sibling, but its residual is the aqueous phase's mass fraction
//! rather than a temperature:
//!
//! ```text
//! wtp   = x_inh MW_inh / (x_inh MW_inh + x_water MW_water)   over the AQUEOUS phase
//! error = -(wtp - wt_target)
//! ```
//!
//! # It is the one that reads a phase's label
//!
//! `system.getPhase(PhaseType.AQUEOUS)` is where the composition comes from, and a cubic phase
//! gets that label from [`label`] - NeqSim's `PhaseEos.init`, whose rule is three branches and
//! one number. The concentration flash needs none of this because its inner step is a hydrate
//! equilibrium, which reads the fluid as a whole.
//!
//! # The association barely matters here, and it is the opposite of its sibling
//!
//! Measured on NeqSim's own fluid at 100 bara, the CPA and plain-SRK columns agree to `3e-5`
//! relative - `0.124372418372090` against `0.124375874492093` mol of MEG at a 0.30 mass
//! fraction. Its `...ConcentrationFlash` sibling is five times apart on the same feed, because
//! *its* equilibrium is on water's fugacity and this one's inner step is an ordinary flash.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result};

use crate::mixture::{Mixture, ReducedParameters};
use crate::pt_flash::pt_flash;
use crate::results::{HydrateInhibitorWtResult, Phase};

/// The residual the secant stops on, in mass fraction. NeqSim's own `1e-5`.
const TOLERANCE: f64 = 1.0e-5;

/// The step cap and the floor, NeqSim's own `iter < 100` and `|| iter < 3`.
const MAXIMUM_STEPS: u32 = 100;
const MINIMUM_STEPS: u32 = 3;

/// The bar NeqSim compares `getVolume() / getB()` against.
///
/// **It is the ratio the name says, with no unit slip in it.** Both halves are in NeqSim's own
/// `1e-5 m3/mol`: methane at 300 K and 50 bar prints `46.0917080848140` and
/// `2.98485123264092`, whose quotient is `15.4418778332323`. That is `Z / b_r`, where `b_r` is
/// the reduced covolume `b P/(R T)` the reduced parameters carry - which is the same number as
/// `V_m/b` in any consistent unit, and is how this computes it.
const GAS_VOLUME_OVER_B: f64 = 1.75;

/// Which of NeqSim's three labels a phase carries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseLabel {
    /// The volume ratio is past the bar.
    Gas,
    /// Not a gas, and its hydrocarbons outweigh its aqueous components.
    Oil,
    /// Not a gas, and they do not.
    Aqueous,
}

/// NeqSim's `PhaseEos.init` rule, on one phase.
///
/// ```text
/// if (getVolume() / getB() > 1.75)        GAS
/// else if (sumHydrocarbons > sumAqueous)  OIL
/// else                                    AQUEOUS
/// ```
///
/// **"Hydrocarbon" is wider than the name: `COMPTYPE == "HC"` *or* `"inert"`.** NeqSim's
/// `isInert()` is a separate predicate and it is counted on the same side, which is what puts
/// a pure CO2 liquid on the OIL branch rather than the aqueous one - measured, CO2's row is
/// `inert`, and a probe state at 260 K and 60 bar comes back `OIL`.
///
/// NeqSim also asks `isTBPfraction()` and `isPlusFraction()`, flags a pseudo-component carries
/// and a databank row does not; this library has no pseudo-component mixtures, so the two
/// types are the whole of it here. Water is excluded from the hydrocarbon sum by name
/// upstream, and its class is `water` rather than either, which is the same exclusion.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if a component carries no `COMPTYPE`, which a card's
///   substance does not - and a label that silently assumed one would decide between an oil
///   and an aqueous phase on nothing.
pub fn label(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    composition: &[f64],
    z_factor: f64,
) -> Result<PhaseLabel> {
    let (_, b_mix) = mixture.mixture_parameters(reduced, composition);
    if b_mix <= 0.0 {
        return Err(AzothError::out_of_range(
            "composition",
            b_mix,
            "a phase with no covolume has no volume ratio to compare",
        ));
    }
    // `Z / b_r` is `V_m/b`: the reduced covolume already carries `P/(R T)`, so no unit enters
    // here at all. The caller's `b_mix` is the *reduced* one `mixture_parameters` returns.
    if z_factor / b_mix > GAS_VOLUME_OVER_B {
        return Ok(PhaseLabel::Gas);
    }

    let mut hydrocarbons = 0.0;
    let mut aqueous = 0.0;
    for (component, fraction) in mixture.components().iter().zip(composition) {
        if component.class.is_empty() {
            return Err(AzothError::invalid_input(
                "components",
                "a component carries no `COMPTYPE`, and a phase's label is decided by whether \
                 its hydrocarbons outweigh its aqueous components"
                    .to_string(),
            ));
        }
        if matches!(component.class.as_str(), "hc" | "inert") && *fraction > 0.0 {
            hydrocarbons += fraction;
        } else {
            aqueous += fraction;
        }
    }
    Ok(if hydrocarbons > aqueous {
        PhaseLabel::Oil
    } else {
        PhaseLabel::Aqueous
    })
}

/// The phases a two-phase flash reported, as `(amount, composition, z factor)`.
fn phases_of(flash: &crate::results::PtFlashResult) -> Vec<(f64, Vec<f64>, f64)> {
    match flash.phase {
        Phase::TwoPhase => {
            // **Vapour first, which is NeqSim's own order** and the order the probe prints.
            // This model does not care - it looks for the aqueous one - but a capture is read
            // side by side with it, and an unstated order is one more thing to check.
            let beta = flash.beta.unwrap_or(0.0);
            vec![
                (beta, flash.y.clone(), flash.z_vapour),
                (1.0 - beta, flash.x.clone(), flash.z_liquid),
            ]
        }
        Phase::AllLiquid => vec![(1.0, flash.x.clone(), flash.z_liquid)],
        _ => vec![(1.0, flash.y.clone(), flash.z_vapour)],
    }
}

/// The mass fraction of the inhibitor in the inhibitor-and-water pair of one phase.
fn weight_fraction(
    mixture: &Mixture,
    composition: &[f64],
    inhibitor: usize,
    water: usize,
) -> Result<f64> {
    let molar_mass = |index: usize| -> Result<f64> {
        mixture.components()[index].molar_mass.ok_or_else(|| {
            AzothError::invalid_input(
                "components",
                "a dosing fraction is a mass fraction of the inhibitor and the water, and one \
                 of them carries no molar mass"
                    .to_string(),
            )
        })
    };
    let inhibitor_mass = composition[inhibitor] * molar_mass(inhibitor)?;
    let water_mass = composition[water] * molar_mass(water)?;
    let total = inhibitor_mass + water_mass;
    if total <= 0.0 {
        return Err(AzothError::invalid_input(
            "components",
            "this phase holds no inhibitor and no water, so it has no dosing fraction".to_string(),
        ));
    }
    Ok(inhibitor_mass / total)
}

/// The moles of inhibitor that put the aqueous phase at a target mass fraction.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `moles` does not match the mixture, if `inhibitor` is not
///   one of its components, if the feed has no water, or if `P` is not positive.
/// * [`AzothError::SolverNotConverged`] if no trial produces an aqueous phase, or the secant
///   reaches its step cap without landing.
pub fn hydrate_inhibitor_wt(
    mixture: &Mixture,
    inhibitor: &str,
    moles: &[f64],
    wt_target: f64,
    t: ThermodynamicTemperature,
    p: Pressure,
) -> Result<HydrateInhibitorWtResult> {
    let spec = &crate::model_gen::HYDRATE_INHIBITOR_WT_SPEC;
    let mut warnings = Vec::new();
    azoth_core::apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "P" => Some(p.value),
            "T" => Some(t.value),
            "wt_target" => Some(wt_target),
            _ => None,
        },
        &mut warnings,
    )?;

    if moles.len() != mixture.components().len() {
        return Err(AzothError::invalid_input(
            "moles",
            format!(
                "{} component(s) and {} mole number(s), and the secant walks them side by side",
                mixture.components().len(),
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
            "no water: this reads an aqueous phase's own inhibitor fraction, so a feed without \
             one has nothing to hold"
                .to_string(),
        )
    })?;

    let mut moles = moles.to_vec();
    let mut error = 1.0;
    let mut old_error = 1.0;
    let mut old_c = moles[index];
    let mut iterations: u32 = 0;

    let (reached, reached_fraction) = loop {
        iterations += 1;
        let c = moles[index];
        // NeqSim's own secant, and its first step is degenerate for the same reason its
        // sibling's is: both amounts are the feed's.
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
            return Err(AzothError::out_of_range(
                "moles",
                moles[index],
                "the secant stepped the inhibitor's moles below zero, which is not a feed",
            ));
        }

        let sum: f64 = moles.iter().sum();
        if sum <= 0.0 {
            return Err(AzothError::out_of_range(
                "moles",
                sum,
                "the feed has no moles left to flash",
            ));
        }
        let fractions: Vec<f64> = moles.iter().map(|value| value / sum).collect();
        let flash = pt_flash(mixture, t, p, &fractions)?;
        let reduced = mixture.reduced_parameters(t, p)?;

        // The aqueous phase's own fraction, which is what the residual is on.
        let mut found = None;
        for (_, composition, z_factor) in phases_of(&flash) {
            if label(mixture, &reduced, &composition, z_factor)? == PhaseLabel::Aqueous {
                found = Some(weight_fraction(mixture, &composition, index, water)?);
                break;
            }
        }
        let Some(wtp) = found else {
            return Err(AzothError::SolverNotConverged {
                iterations,
                residual: error,
                tolerance: TOLERANCE,
            });
        };
        error = -(wtp - wt_target);

        if !((error.abs() > TOLERANCE && iterations < MAXIMUM_STEPS) || iterations < MINIMUM_STEPS)
        {
            break (flash, wtp);
        }
    };

    if error.abs() > TOLERANCE {
        return Err(AzothError::SolverNotConverged {
            iterations,
            residual: error,
            tolerance: TOLERANCE,
        });
    }

    Ok(HydrateInhibitorWtResult {
        inhibitor_moles: moles[index],
        weight_fraction: reached_fraction,
        phases: phases_of(&reached).len() as u32,
        iterations,
        residual: error,
        warnings,
    })
}
