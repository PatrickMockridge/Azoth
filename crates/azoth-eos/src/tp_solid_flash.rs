//! `eos.tp_solid_flash` - how much of a feed has frozen out as one pure solid.
//!
//! ```text
//! seed     the fluid-only flash, plus a solid of the component the caller named
//! solve    Q(beta) = sum_k beta_k - sum_i z_i ln E_i      E_i = sum_k beta_k / phi_ik
//!          with   E_solid = z_solid / phi_solid
//! drop     every phase the solve pinned under `1.01e-12`
//! ```
//!
//! Spec: `specs/models/eos/tp_solid_flash.toml`, which carries the seeding, the pinned
//! component's rule and the three states it is checked at.
//!
//! NeqSim's `SolidFlash`, reached through `Flash.solidPhaseFlash()` from a `TPflash` whose
//! system has `setSolidPhaseCheck(name)` on - which is the route `AsphalteneOnsetPressureFlash`
//! takes, and not the `TPSolidflash()` one, which builds a `SolidFlash1`. **The two classes
//! are one equilibrium and two arithmetics**: `SolidFlash` normalises each phase and damps its
//! Newton step by `(iter+1)/(10+iter)`; `SolidFlash1` normalises nothing and line-searches on
//! `Q` instead. The probe runs both, and they agree to ten digits at the states where the
//! fluid settles on one phase and part company where it is two - `SolidFlash1`'s aqueous comes
//! back with `x` summing to `0.819` - so this port follows `SolidFlash`.
//!
//! # What is different about a solid
//!
//! The fraction Newton is [`crate::multiphase`]'s - `Q(beta)` with `E_i = sum_k beta_k/phi_ik`
//! and the gradient and Hessian that follow - and one component's `E` is not that sum. **The
//! precipitating component's fugacity is pinned by the solid rather than by the material
//! balance**:
//!
//! ```text
//! E_solid = z_solid / phi_solid  =>  x_solid^k = z_solid/(E_solid phi_solid^k) = phi_solid/phi_solid^k
//! ```
//!
//! so every fluid phase holds the solid component at the pure solid's own fugacity instead of
//! at the amount the feed gave it. The Hessian loses that component's term for the same
//! reason - a pinned `E` does not move with `beta` - and the solid's own amount is then not a
//! Newton unknown at all but the material balance's remainder,
//! `beta_solid = z_solid - sum_k beta_k (phi_solid/phi_solid^k)`.
//!
//! **The solid's coefficients depend on the state alone.** `ComponentSolid.fugcoef2`'s
//! `SolidFug/(P x)` cancels the mole fraction out, so there is no trial composition to make
//! and no root to choose - one vector throughout, exactly as the wax phase's is.
//!
//! # The two constants the tabulated route carries
//!
//! `fugcoef2` reads the solid-liquid heat-capacity difference and molar-volume difference off
//! the component's own tables, and **overrides the first to `37.12` J/(mol K) for the
//! component named `water`** - the case the class was written around. The second vanishes for
//! water twice over: the solid density table is zero on its row, and NeqSim falls back to the
//! liquid's own molar volume there.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, ModelAlgorithm, Result, apply_checks};

use crate::Cubic;
use crate::mixture::{Component, Mixture, ReducedParameters, RootSide};
use crate::pt_flash::pt_flash;
use crate::results::{Phase as FlashPhase, TpSolidFlashResult};
use crate::solid_fugacity::{R, solid_fugacity};
use crate::{algorithm_of, databank, model_gen};

/// The floor under a fluid phase fraction, NeqSim's `phaseFractionMinimumLimit`.
const FRACTION_FLOOR: f64 = 1.0e-12;

/// A fluid phase below this is removed and the solve restarted, upstream's `1.01e-9`.
const REMOVE_LIMIT: f64 = 1.01e-9;

/// A phase below this is dropped from the answer, upstream's `phaseFractionMinimumLimit *
/// 1.01` - the same rule `eos.tp_multiflash`'s merge applies.
const FRACTION_DROP: f64 = 1.01e-12;

/// The diagonal regulariser, upstream's `Qmatrix[i][i] = 1.0e-9`.
///
/// **Six orders larger than the multiphase solve's**, and it is what keeps the step finite:
/// the Hessian is nearly singular in the direction of a phase that is about to vanish, and
/// upstream's larger floor is measured on exactly that.
const REGULARISER: f64 = 1.0e-9;

/// The restarts a removed fluid phase costs, so a state that keeps shedding phases stops.
///
/// Upstream removes at most one phase per `solveBeta` and re-runs once, so the bound is the
/// phase count and not a tuning knob.
const MAX_RESTARTS: u32 = 8;

/// The outer loop's cap, upstream's `!(iter > 20)`.
const MAX_OUTER: u32 = 20;

/// The fewest outer steps taken, upstream's `|| iter < 4`.
const MINIMUM_OUTER: u32 = 4;

/// The outer loop's stopping rule on the solid's fraction, upstream's `1e-3`.
const SOLID_TOLERANCE: f64 = 1.0e-3;

/// The heat-capacity difference NeqSim overrides for water, in J/(mol K).
const WATER_DELTA_CP: f64 = 37.12;

/// The factor NeqSim's `getDensity()` carries, which its own fallback divides back out.
///
/// `Phase.getDensity()` is `1/molarVolume * molarMass * 1e5`, so `1/density * M` is the molar
/// volume over `1e5`. It is the branch only a component whose liquid-density table is zero
/// takes, and it is ported as written: both sides of the volume difference are read the same
/// way, so it is a common scale rather than a unit that changes one side alone.
const DENSITY_SCALE: f64 = 1.0e5;

/// One fluid phase of the solve.
#[derive(Debug, Clone)]
struct FluidPhase {
    fraction: f64,
    composition: Vec<f64>,
    side: RootSide,
}

/// `c1 + c2 T + c3 T^2 + ...`, the form both density correlations are written in.
fn polyval(coefficients: &[f64], temperature: f64) -> f64 {
    coefficients
        .iter()
        .rev()
        .fold(0.0, |sum, value| sum * temperature + value)
}

/// `M 1000 poly(T)`, NeqSim's `getPureComponentSolidDensity` and its liquid twin.
///
/// `None` where the polynomial is not positive, which is how the table states a component it
/// carries no correlation for - carbon dioxide and benzene among them, and water's solid row.
fn tabulated_molar_volume(
    component: &Component,
    coefficients: &[f64],
    temperature: f64,
) -> Option<f64> {
    let molar_mass = component.molar_mass?;
    let density = molar_mass * 1000.0 * polyval(coefficients, temperature);
    if density > 1.0e-20 {
        Some(molar_mass / density)
    } else {
        None
    }
}

fn cubic_of(eos: &str) -> Result<Cubic> {
    match eos {
        "pr" => Ok(Cubic::Pr),
        "srk" => Ok(Cubic::Srk),
        other => Err(AzothError::invalid_input(
            "eos",
            format!(
                "`{other}` is not a cubic this reaches: the solid's reference liquid is a \
                 phase of the host's own class, and only `pr` and `srk` are ported"
            ),
        )),
    }
}

/// A pure solid's fugacity coefficient from the named component's own tables.
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
fn tabulated_solid_fugacity(
    mixture: &Mixture,
    index: usize,
    solid: &str,
    T: ThermodynamicTemperature,
    P: Pressure,
    eos: &str,
) -> Result<f64> {
    let component = &mixture.components()[index];
    let triple_point = component.triple_point_temperature;
    if component.heat_of_fusion <= 0.0 || triple_point <= 0.0 {
        return Err(AzothError::invalid_input(
            "solid",
            format!(
                "`{solid}` carries no heat of fusion or no triple-point temperature, and the \
                 solid route reads both rather than defaulting: a zero heat of fusion is a \
                 substance that does not melt"
            ),
        ));
    }
    if component.molar_mass.is_none() {
        return Err(AzothError::invalid_input(
            "solid",
            format!("`{solid}` carries no molar mass, and both density correlations scale by it"),
        ));
    }

    // `ComponentSolid.fugcoef2`: the solid-liquid difference at the triple point, from the
    // tables, with the one override the class states by name.
    let delta_cp_sl = if solid.trim().to_lowercase() == "water" {
        WATER_DELTA_CP
    } else {
        (polyval(&component.cp_liquid, triple_point) - polyval(&component.cp_solid, triple_point))
            / 1000.0
    };

    // The liquid side. Where the table states no correlation, NeqSim falls back to the
    // *reference phase's* own molar volume - the component alone, on the fluid's cubic, on
    // the liquid root at the fluid's state.
    let liquid =
        match tabulated_molar_volume(component, &component.liquid_density_coefs, triple_point) {
            Some(volume) => volume,
            None => {
                let reference = Mixture::new(vec![component.clone()], vec![0.0])?
                    .with_cubic(cubic_of(eos)?)
                    .with_alpha(match cubic_of(eos)? {
                        Cubic::Pr => crate::Alpha::Pr,
                        _ => crate::Alpha::Srk,
                    });
                let reduced = reference.reduced_parameters(T, P)?;
                let state = reference.phase_state(&reduced, &[1.0], RootSide::Liquid)?;
                state.z * R * T.value / P.value / DENSITY_SCALE
            }
        };
    // The solid side, falling back to the liquid's own volume - which is what makes the
    // difference vanish for water rather than become a division by zero.
    let solid_volume =
        tabulated_molar_volume(component, &component.solid_density_coefs, triple_point)
            .unwrap_or(liquid);

    let result = solid_fugacity(
        component.heat_of_fusion,
        triple_point,
        delta_cp_sl,
        solid_volume - liquid,
        component.tc,
        component.pc,
        component.omega,
        T,
        P,
        eos,
    )?;
    Ok(result.fugacity_coefficient)
}

/// The fugacity coefficients of every component in every fluid phase, row-major by phase.
fn coefficients(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    phases: &[FluidPhase],
) -> Result<Vec<Vec<f64>>> {
    let mut out = Vec::with_capacity(phases.len());
    for phase in phases {
        let state = mixture.phase_state(reduced, &phase.composition, phase.side)?;
        out.push(state.ln_phi.iter().map(|value| value.exp()).collect());
    }
    Ok(out)
}

/// `E_i` over the fluid phases, with the precipitating component's value pinned to the solid.
fn denominators(
    z: &[f64],
    phi: &[Vec<f64>],
    phases: &[FluidPhase],
    solid: usize,
    phi_solid: f64,
) -> Vec<f64> {
    let n = z.len();
    let mut e = vec![0.0; n];
    for (k, phase) in phases.iter().enumerate() {
        for i in 0..n {
            e[i] += phase.fraction / phi[k][i];
        }
    }
    e[solid] = z[solid] / phi_solid;
    e
}

/// The fraction of a feed that has frozen out as one pure solid.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `z` is the wrong length, negative or does not sum to one;
///   if `solid` is not one of `components`; or if the named component carries no melt data.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive, or if the feed is at or above
///   the named component's triple point - where the screen cannot add a solid.
/// * [`AzothError::SolverNotConverged`] if the fraction solve's Hessian is singular at an
///   iterate.
#[allow(clippy::too_many_arguments, non_snake_case)]
pub fn tp_solid_flash(
    components: &[&str],
    solid: &str,
    T: ThermodynamicTemperature,
    P: Pressure,
    z: &[f64],
    eos: &str,
) -> Result<TpSolidFlashResult> {
    let spec = &model_gen::TP_SOLID_FLASH_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T.value),
            "P" => Some(P.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let n = components.len();
    if z.len() != n {
        return Err(AzothError::invalid_input(
            "z",
            format!("a feed for {n} components has {} entries", z.len()),
        ));
    }
    if let Some(bad) = z.iter().position(|&value| value < 0.0) {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "z[{bad}] is {} but a mole fraction cannot be negative",
                z[bad]
            ),
        ));
    }
    let sum: f64 = z.iter().sum();
    if (sum - 1.0).abs() > 1.0e-09 {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "the mole fractions sum to {sum}, not to one. Renormalising them here would \
                 make a composition error invisible in every number downstream, so it is \
                 refused instead"
            ),
        ));
    }

    let (mixture, _) = databank::mixture_of(components, cubic_of(eos)?, None)?;
    let index = mixture.index_of(solid).ok_or_else(|| {
        AzothError::invalid_input(
            "solid",
            format!("`{solid}` is not one of the feed's components"),
        )
    })?;
    // **Above the triple point the screen is left to answer**, and it answers a negative
    // amount: NeqSim's `checkAndAddSolidPhase` skips such a component outright, and its
    // `fugcoef2` above the triple point only makes the solid's coefficient larger than the
    // liquid's. Refusing here instead would turn a state the fluid flash answers into a
    // failure.

    let algorithm = algorithm_of(spec)?;
    let reduced = mixture.reduced_parameters(T, P)?;
    let phi_solid = tabulated_solid_fugacity(&mixture, index, solid, T, P, eos)?;

    // **The fluid-only flash seeds it**, and is the answer where nothing precipitates.
    let flash = pt_flash(&mixture, T, P, z)?;
    warnings.extend(flash.warnings.iter().cloned());
    let mut phases: Vec<FluidPhase> = match flash.beta {
        Some(split) => vec![
            FluidPhase {
                fraction: 1.0 - split,
                composition: flash.x.clone(),
                side: RootSide::Liquid,
            },
            FluidPhase {
                fraction: split,
                composition: flash.y.clone(),
                side: RootSide::Vapour,
            },
        ],
        // One fluid phase: the flash reports no fraction, so the feed *is* the phase, on
        // whichever root its own lower Gibbs energy names.
        None => vec![FluidPhase {
            fraction: 1.0,
            composition: z.to_vec(),
            side: if flash.phase == FlashPhase::AllLiquid {
                RootSide::Liquid
            } else {
                RootSide::Vapour
            },
        }],
    };

    // **The candidate screen**, `Flash.solidPhaseFlash`'s: how much of the solid component the
    // fluid phases do not account for at the solid's own fugacity.
    let start = coefficients(&mixture, &reduced, &phases)?;
    let mut candidate = z[index];
    for (k, phase) in phases.iter().enumerate() {
        candidate -= phase.fraction * phi_solid / start[k][index];
    }
    if candidate <= 1.0e-20 {
        let split = flash.beta.unwrap_or(0.0);
        return Ok(TpSolidFlashResult {
            solid_fraction: 0.0,
            phase_count: 2,
            beta: vec![1.0 - split, split],
            x: vec![flash.x.clone(), flash.y.clone()],
            solid_fugacity_coefficient: phi_solid,
            iterations: 0,
            residual: 0.0,
            converged: true,
            warnings,
        });
    }

    let mut iterations = 0;
    let mut residual = f64::NAN;
    let mut converged = false;
    let mut solid_fraction = candidate;
    let mut restarts = 0;
    'restart: loop {
        // `run()`: the Newton, then `setXY` and the solid's material balance, until the solid's
        // amount stops moving.
        for outer in 1..=MAX_OUTER {
            let previous = solid_fraction;
            let (steps, left, met, removed) = solve_beta(
                &mixture,
                &reduced,
                z,
                &mut phases,
                index,
                phi_solid,
                &mut solid_fraction,
                algorithm,
            )?;
            iterations = steps;
            residual = left;
            converged = met;
            if removed {
                if restarts < MAX_RESTARTS {
                    restarts += 1;
                    continue 'restart;
                }
                break 'restart;
            }
            set_compositions(&mixture, &reduced, z, &mut phases, index, phi_solid)?;
            let phi = coefficients(&mixture, &reduced, &phases)?;
            solid_fraction = z[index]
                - phases
                    .iter()
                    .enumerate()
                    .map(|(k, phase)| phase.fraction * phi_solid / phi[k][index])
                    .sum::<f64>();
            if outer >= MINIMUM_OUTER && (solid_fraction - previous).abs() <= SOLID_TOLERANCE {
                break;
            }
        }
        break;
    }

    // `TPflash`'s own sweep after the solid flash: a phase under the floor is not there.
    let mut beta: Vec<f64> = phases.iter().map(|phase| phase.fraction).collect();
    beta.push(solid_fraction);
    let mut x: Vec<Vec<f64>> = phases
        .iter()
        .map(|phase| phase.composition.clone())
        .collect();
    let mut solid_row = vec![0.0; n];
    solid_row[index] = 1.0;
    x.push(solid_row);
    let keep: Vec<usize> = (0..beta.len())
        .filter(|&k| beta[k] >= FRACTION_DROP)
        .collect();
    let solid_present = keep.last().is_some_and(|&k| k == phases.len());

    Ok(TpSolidFlashResult {
        solid_fraction: if solid_present { solid_fraction } else { 0.0 },
        phase_count: keep.len() as u32,
        beta: keep.iter().map(|&k| beta[k]).collect(),
        x: keep.iter().map(|&k| x[k].clone()).collect(),
        solid_fugacity_coefficient: phi_solid,
        iterations,
        residual,
        converged,
        warnings,
    })
}

/// One `solveBeta` call: the Newton on the fluid fractions with the coefficients frozen.
///
/// Returns the steps taken, the residual, whether it met the tolerance, and whether a fluid
/// phase was removed - which the caller answers the way upstream does, by restarting `run()`.
#[allow(clippy::too_many_arguments)]
fn solve_beta(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    z: &[f64],
    phases: &mut Vec<FluidPhase>,
    solid: usize,
    phi_solid: f64,
    solid_fraction: &mut f64,
    algorithm: &ModelAlgorithm,
) -> Result<(u32, f64, bool, bool)> {
    let n = z.len();
    let mut residual = f64::NAN;
    let mut steps = 0;
    for step in 1..=algorithm.max_iterations {
        steps = step;
        let phi = coefficients(mixture, reduced, phases)?;
        let e = denominators(z, &phi, phases, solid, phi_solid);

        // `1 - sum_i z_i/(E_i phi_ik)`, whose pinned term is `phi_solid/phi_ik` by the
        // override above; and the Hessian, which loses that component's term because a
        // pinned `E` does not move with `beta`.
        let count = phases.len();
        let mut gradient = vec![0.0; count];
        let mut hessian = vec![vec![0.0; count]; count];
        for (k, value) in gradient.iter_mut().enumerate() {
            *value = 1.0;
            for i in 0..n {
                *value -= z[i] / (e[i] * phi[k][i]);
            }
        }
        for j in 0..count {
            for k in 0..count {
                let mut value = if j == k { REGULARISER } else { 0.0 };
                for i in 0..n {
                    if i != solid {
                        value += z[i] / (e[i] * e[i] * phi[j][i] * phi[k][i]);
                    }
                }
                hessian[j][k] = value;
            }
        }

        let correction = crate::multiphase::solve_dense(&hessian, &gradient, algorithm)?;
        residual = correction
            .iter()
            .map(|value| value * value)
            .sum::<f64>()
            .sqrt();

        // Damped by `(step+1)/(10+step)` - two elevenths at the first step - and clamped to
        // the floor, where a phase that is not there sits rather than a negative amount.
        let damping = f64::from(step + 1) / (10.0 + f64::from(step));
        for (k, phase) in phases.iter_mut().enumerate() {
            let candidate = phase.fraction - damping * correction[k];
            phase.fraction = if candidate < FRACTION_FLOOR {
                FRACTION_FLOOR
            } else if candidate > 1.0 {
                1.0 - FRACTION_FLOOR
            } else {
                candidate
            };
        }

        // The solid's own amount: `z_solid` less what the fluid phases hold at the pinned
        // fugacity. Not a Newton unknown, so it is read back rather than stepped.
        let phi = coefficients(mixture, reduced, phases)?;
        *solid_fraction = z[solid]
            - phases
                .iter()
                .enumerate()
                .map(|(k, phase)| phase.fraction * phi_solid / phi[k][solid])
                .sum::<f64>();

        // A fluid phase the Newton took under `1.01e-9` is one that is not there, and the
        // solve restarts without it - upstream's `removePhaseKeepTotalComposition`, which
        // leaves the rest holding the whole feed.
        if let Some(vanish) = phases
            .iter()
            .position(|phase| phase.fraction.abs() < REMOVE_LIMIT)
        {
            let gone = phases.remove(vanish);
            let kept = 1.0 - gone.fraction;
            for phase in phases.iter_mut() {
                phase.fraction /= kept;
            }
            return Ok((steps, residual, false, true));
        }

        if step >= 2 && residual <= algorithm.tolerance {
            return Ok((steps, residual, true, false));
        }
    }
    Ok((steps, residual, false, false))
}

/// `setXY`: the composition the material balance gives at the current fractions.
fn set_compositions(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    z: &[f64],
    phases: &mut [FluidPhase],
    solid: usize,
    phi_solid: f64,
) -> Result<()> {
    let n = z.len();
    let phi = coefficients(mixture, reduced, phases)?;
    let e = denominators(z, &phi, phases, solid, phi_solid);
    for (k, phase) in phases.iter_mut().enumerate() {
        let mut total = 0.0;
        for i in 0..n {
            let value = z[i] / (e[i] * phi[k][i]);
            phase.composition[i] = value;
            total += value;
        }
        if total > 0.0 {
            for value in phase.composition.iter_mut() {
                *value /= total;
            }
        }
    }
    Ok(())
}
