//! `eos.hybrid_eos_ge_flash` - NeqSim's fixed-role EoS-gas / EoS-oil / GE-aqueous flash.
//!
//! Spec: `specs/models/eos/hybrid_eos_ge_flash.toml`, which records the fixed topology, the
//! two rules that make an ion aqueous, and the acceptance contract the answer is held to.
//!
//! Ported from `TPHybridEosGeFlash`, which extends `TPmultiflash` and overrides exactly four
//! things: `calcE`, `calcQ`, `setXY` and `limitBetaStepScale`. The fraction Newton is the
//! parent's and this module's is that Newton written once more, because the four overrides
//! change both what a coefficient is and how a composition is built - see the spec's
//! assumptions for why sharing [`crate::multiphase::solve_phase_fractions`] would have meant
//! parameterising it past what it is.

use azoth_core::{AzothError, ModelAlgorithm, Result};

use crate::cubic::Cubic;
use crate::databank::HYDROCARBON;
use crate::mixture::{Mixture, ReducedParameters, RootSide};
use crate::multiphase::{MultiphasePhase, PhaseKind, coefficients};
use crate::results::HybridEosGeFlashResult;

/// The roles, in the solver's own phase order.
///
/// `SystemEosGE.restoreHybridEosGePhaseRoles` maps slot 0 to the gas, slot 2 to the oil and
/// slot 1 to the GE liquid, and then *indexes them* `[gas, oil, aqueous]` - so the creation
/// order is not the array order, and the array order is what every vector here is in.
pub const ROLE_ORDER: [&str; 3] = ["gas", "oil", "aqueous"];

/// A phase fraction is kept strictly inside `(0, 1)`, upstream's `phaseFractionMinimumLimit`.
const FRACTION_FLOOR: f64 = 1.0e-12;

/// `E_i`'s floor, upstream's.
const MINIMUM_E: f64 = 1.0e-100;

/// The ion concentration a non-aqueous phase is given instead of a computed one.
///
/// Not the same number as the `1e-40` the acceptance contract asserts: this is the *floor a
/// composition is built on*, and that is the tolerance the answer is read at.
const ION_CONFINEMENT: f64 = 1.0e-50;

/// The whole-number cap on the beta step where the inventory carries an ion.
const IONIC_BETA_STEP_SCALE: f64 = 0.1;

/// The ratio a phase fraction may move by in one projection, upstream's
/// `HYBRID_AQUEOUS_FRACTION_STEP_LIMIT`.
const AQUEOUS_STEP_LIMIT: f64 = 2.0;

/// The aqueous fraction's margin over the ionic inventory, in units of the floor.
const ION_CAPACITY_MARGIN: f64 = 100.0 * FRACTION_FLOOR;

/// The diagonal regulariser, `TPHybridEosGeFlash.calcQ`'s own.
const REGULARISER: f64 = 1.0e-10;

/// The gradient norm the iteration must also reach, upstream's second convergence test.
const GRADIENT_TOLERANCE: f64 = 1.0e-10;

/// A fraction at or below this is a phase that is not there, in the acceptance contract.
const MATERIAL_PHASE_FRACTION: f64 = 1.0e-10;

/// The fixed-topology solve's own cap: `MAXIMUM_HYBRID_ITERATIONS`, one pass per call of the
/// parent's `solveBeta` that the outer loop makes.
const MAXIMUM_HYBRID_ITERATIONS: u32 = 50;

/// The fewest outer passes, upstream's `iteration >= 4`.
///
/// The outer loop exists to re-seed the composition from the *settled* fractions, so its
/// first passes are not converged work: `iteration >= 4` is what stops a state that happens
/// to look converged on its first pass from being reported.
const MINIMUM_OUTER_ITERATIONS: u32 = 4;

/// A gas-oil-aqueous split whose roles are fixed before the fractions are solved.
///
/// The feed is in **moles**, one per component, and `ion` says which components may not
/// leave the aqueous phase. Returns the fractions and compositions of the three roles in
/// [`ROLE_ORDER`]'s order, plus the two residuals the acceptance contract is stated in.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `moles` is not a positive vector.
/// * [`AzothError::SolverNotConverged`] if a Newton correction cannot be solved at an
///   iterate.
/// * [`AzothError::PropertyUnavailable`] if a phase model cannot answer at a trial
///   composition the iteration reaches.
pub fn hybrid_eos_ge_flash(
    names: &[&str],
    cubic: Cubic,
    t: f64,
    p: f64,
    moles: &[f64],
) -> Result<HybridEosGeFlashResult> {
    let (mixture, _) = crate::databank::mixture_of_with_ions(names, cubic, None)?;
    let spec = &crate::model_gen::HYBRID_EOS_GE_FLASH_SPEC;
    let mut warnings = Vec::new();
    azoth_core::apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t),
            "P" => Some(p),
            _ => None,
        },
        &mut warnings,
    )?;
    let algorithm = crate::algorithm_of(spec)?;
    let reduced =
        mixture.reduced_parameters(azoth_core::units::kelvins(t), azoth_core::units::pascals(p))?;
    // The ions are the databank's own classification, which is where NeqSim reads them from
    // too: `class == ION` is the compiled `COMPTYPE`, and a charge is not needed beside it.
    let ion: Vec<bool> = names
        .iter()
        .map(|name| {
            crate::databank::entry(name, None).map(|entry| entry.class == crate::databank::ION)
        })
        .collect::<Result<_>>()?;
    let mut result = solve_fixed_topology(&mixture, &reduced, moles, &ion, algorithm)?;
    result.warnings = warnings;
    Ok(result)
}

/// The fixed-topology solve over a mixture that already carries its ions.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `moles` is not a positive vector or `ion` does not match
///   it.
/// * [`AzothError::SolverNotConverged`] if a Newton correction cannot be solved at an
///   iterate, or if the aqueous role cannot hold the ionic inventory.
/// * [`AzothError::PropertyUnavailable`] if a phase model cannot answer at a trial
///   composition the iteration reaches.
pub(crate) fn solve_fixed_topology(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    feed: &[f64],
    ion: &[bool],
    algorithm: &ModelAlgorithm,
) -> Result<HybridEosGeFlashResult> {
    let n = feed.len();
    if ion.len() != n {
        return Err(AzothError::invalid_input(
            "components",
            format!("a feed of {n} components has {} ion flags", ion.len()),
        ));
    }
    let total: f64 = feed.iter().sum();
    // `!(x > 0.0)` would say this in one clause, and clippy is right that it reads as a
    // negated comparison; the two reject a `NaN` as well as a zero or a negative, which is
    // the point.
    if !total.is_finite() || total <= 0.0 {
        return Err(AzothError::invalid_input(
            "moles",
            format!("the feed sums to {total}, which is not a mole count"),
        ));
    }
    let z: Vec<f64> = feed.iter().map(|value| value / total).collect();

    let mut phases = seed_phases(mixture, &z, ion)?;
    solve_fixed_topology_from(mixture, reduced, feed, ion, algorithm, &mut phases)
}

/// The fixed-topology Newton **from a split the caller already has**.
///
/// Split out of [`solve_fixed_topology`] because the reactive coupling re-enters it: that loop
/// re-equilibrates the brine's chemistry, writes the adjusted species inventory back into the
/// roles and solves the fractions again, and each pass starts from the fractions the last one
/// converged on rather than from the seed. Upstream's `run` is arranged the same way - its
/// `system` carries the split between passes - and re-seeding each pass would throw away the
/// only thing the passes share.
///
/// `phases` is the split, in and out, in `[gas, oil, aqueous]` order and with the kinds
/// [`seed_phases`] gives it.
///
/// # Errors
/// The same as [`solve_fixed_topology`], less the feed checks - the caller has made those.
pub fn solve_fixed_topology_from(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    feed: &[f64],
    ion: &[bool],
    algorithm: &ModelAlgorithm,
    phases: &mut [MultiphasePhase],
) -> Result<HybridEosGeFlashResult> {
    let n = feed.len();
    let total: f64 = feed.iter().sum();
    if !total.is_finite() || total <= 0.0 {
        return Err(AzothError::invalid_input(
            "moles",
            format!("the feed sums to {total}, which is not a mole count"),
        ));
    }
    let z: Vec<f64> = feed.iter().map(|value| value / total).collect();
    let aqueous = 2;

    let mut iterations = 0;
    let mut residual = f64::NAN;
    let mut gradient_norm = f64::INFINITY;
    // **Two nested caps, and both are load-bearing.** `solveFixedTopologyPhaseEquilibrium`
    // runs up to fifty of `solveBeta`, and each `solveBeta` runs up to fifty steps - and the
    // step scale is capped at a tenth while the inventory carries an ion, so a single pass
    // of fifty moves the fractions by about five Newton steps' worth. One cap alone would
    // report the fifth outer pass's mid-flight state as the answer.
    for outer in 0..MAXIMUM_HYBRID_ITERATIONS {
        for step in 1..=algorithm.max_iterations {
            iterations += 1;
            let phi = coefficients(mixture, reduced, phases)?;
            let inverted: Vec<Vec<f64>> = phi
                .iter()
                .enumerate()
                .map(|(k, row)| {
                    row.iter()
                        .enumerate()
                        .map(|(i, value)| {
                            if (ion[i] && k != aqueous) || !value.is_finite() || *value <= 0.0 {
                                // An ion is excluded from the EoS roles **exactly**, rather than
                                // by a large finite penalty, and a coefficient the phase model
                                // could not answer with is excluded the same way.
                                return 0.0;
                            }
                            1.0 / value
                        })
                        .collect()
                })
                .collect();

            let mut e = vec![0.0; n];
            for (k, phase) in phases.iter().enumerate() {
                for i in 0..n {
                    e[i] += phase.fraction * inverted[k][i];
                }
            }
            for value in e.iter_mut() {
                if !value.is_finite() || *value < MINIMUM_E {
                    *value = MINIMUM_E;
                }
            }

            let mut gradient = vec![0.0; 3];
            let mut hessian = vec![vec![0.0; 3]; 3];
            for (k, value) in gradient.iter_mut().enumerate() {
                *value = 1.0;
                for i in 0..n {
                    *value -= z[i] * inverted[k][i] / e[i];
                }
            }
            for (j, row) in hessian.iter_mut().enumerate() {
                for (k, value) in row.iter_mut().enumerate() {
                    for i in 0..n {
                        *value += z[i] * inverted[j][i] * inverted[k][i] / (e[i] * e[i]);
                    }
                    if j == k {
                        *value += REGULARISER;
                    }
                }
            }
            gradient_norm = gradient
                .iter()
                .map(|value| value * value)
                .sum::<f64>()
                .sqrt();
            let correction = crate::multiphase::solve_dense(&hessian, &gradient, algorithm)?;
            residual = correction
                .iter()
                .map(|value| value * value)
                .sum::<f64>()
                .sqrt();

            // The step is damped by `n/(n+3)`, except that a phase carrying the whole of an ionic
            // inventory may not move by more than a tenth of its correction: an unconstrained
            // proposal can push the aqueous fraction below the ions it has to hold.
            // `calcQ` captures the settled fraction *before* the proposal is written, and the
            // projection is stated against it: the band a trial may move in is a band around
            // where the iteration already was, not around where this step just put it.
            let previous_aqueous = phases[aqueous].fraction;
            let scale = if ion.iter().any(|flag| *flag) {
                (f64::from(step) / (f64::from(step) + 3.0)).min(IONIC_BETA_STEP_SCALE)
            } else {
                f64::from(step) / (f64::from(step) + 3.0)
            };
            let mut total = 0.0;
            for (k, phase) in phases.iter_mut().enumerate() {
                let candidate = phase.fraction - scale * correction[k];
                phase.fraction = candidate.clamp(FRACTION_FLOOR, 1.0 - FRACTION_FLOOR);
                total += phase.fraction;
            }
            for phase in phases.iter_mut() {
                phase.fraction /= total;
            }
            enforce_aqueous_fraction_bounds(&mut *phases, ion, &z, aqueous, previous_aqueous)?;

            set_compositions(&mut *phases, &z, &e, &inverted, ion, aqueous)?;

            if step >= 2 && residual <= algorithm.tolerance && gradient_norm <= GRADIENT_TOLERANCE {
                break;
            }
        }
        if outer >= MINIMUM_OUTER_ITERATIONS
            && residual <= algorithm.tolerance
            && gradient_norm <= GRADIENT_TOLERANCE
        {
            break;
        }
    }

    let beta: Vec<f64> = phases.iter().map(|phase| phase.fraction).collect();
    let composition: Vec<Vec<f64>> = phases
        .iter()
        .map(|phase| phase.composition.clone())
        .collect();
    let ln_phi: Vec<Vec<f64>> = coefficients(mixture, reduced, phases)?
        .iter()
        .map(|row| row.iter().map(|value| value.ln()).collect())
        .collect();

    let (balance, fugacity) = acceptance(&beta, &composition, &ln_phi, &z, ion, reduced)?;
    let min_t_over_tc = reduced
        .reduced_temperatures
        .iter()
        .fold(f64::INFINITY, |smallest, value| smallest.min(*value));

    Ok(HybridEosGeFlashResult {
        beta,
        x: composition,
        ln_phi,
        iterations,
        residual: residual.max(gradient_norm),
        max_material_balance_residual: balance,
        max_log_fugacity_residual: fugacity,
        min_t_over_tc,
        warnings: Vec::new(),
    })
}

/// The seed `SystemEosGE.prepareHybridEosGeFlash` gives the three roles.
///
/// **Per component class and not per state**, which is what makes the fixed topology work: a
/// heavy hydrocarbon starts in the oil, water and the polar solvents start in the aqueous
/// phase, and an ion starts in the aqueous phase with the `1e-50` floor in the other two.
/// The fractions are the feed's own shares of those three groups, floored at `1e-5`.
pub fn seed_phases(mixture: &Mixture, z: &[f64], ion: &[bool]) -> Result<Vec<MultiphasePhase>> {
    let n = z.len();
    let names = mixture.names().ok_or_else(|| {
        AzothError::invalid_input(
            "components",
            "the roles are seeded from the components' own names and classes, and this \
             mixture carries neither"
                .to_string(),
        )
    })?;
    let mut gas = vec![0.0; n];
    let mut oil = vec![0.0; n];
    let mut aqueous = vec![0.0; n];
    let (mut gas_feed, mut oil_feed, mut aqueous_feed) = (0.0, 0.0, 0.0);

    for i in 0..n {
        let component = &mixture.components()[i];
        let fraction = z[i].max(ION_CONFINEMENT);
        let is_water_like = is_aqueous_feed_component(&names[i]);
        let heavy = component.molar_mass.unwrap_or_default() > 0.045;
        if ion[i] {
            gas[i] = ION_CONFINEMENT;
            oil[i] = ION_CONFINEMENT;
            aqueous[i] = fraction;
            aqueous_feed += z[i].max(0.0);
        } else if is_water_like {
            gas[i] = (fraction * 1.0e-3).min(1.0e-8);
            oil[i] = (fraction * 1.0e-4).min(1.0e-10);
            aqueous[i] = fraction;
            aqueous_feed += z[i].max(0.0);
        } else if component.class == HYDROCARBON {
            if heavy {
                gas[i] = (fraction * 1.0e-2).max(1.0e-16);
                oil[i] = fraction;
                oil_feed += z[i].max(0.0);
            } else {
                gas[i] = fraction;
                oil[i] = (fraction * 5.0e-2).max(1.0e-16);
                gas_feed += z[i].max(0.0);
            }
            aqueous[i] = (fraction * 1.0e-8).max(1.0e-30);
        } else {
            gas[i] = fraction;
            oil[i] = (fraction * 1.0e-1).max(1.0e-16);
            aqueous[i] = (fraction * 1.0e-2).max(1.0e-20);
            gas_feed += z[i].max(0.0);
        }
    }
    for values in [&mut gas, &mut oil, &mut aqueous] {
        let sum: f64 = values.iter().sum();
        if sum > 0.0 {
            for value in values.iter_mut() {
                *value /= sum;
            }
        }
    }

    let floor = 1.0e-5;
    let gas_feed = gas_feed.max(floor);
    let oil_feed = oil_feed.max(floor);
    let aqueous_feed = aqueous_feed.max(floor);
    let total = gas_feed + oil_feed + aqueous_feed;

    Ok(vec![
        MultiphasePhase {
            fraction: gas_feed / total,
            composition: gas,
            kind: PhaseKind::Cubic(RootSide::Vapour),
        },
        MultiphasePhase {
            fraction: oil_feed / total,
            composition: oil,
            kind: PhaseKind::Cubic(RootSide::Liquid),
        },
        MultiphasePhase {
            fraction: aqueous_feed / total,
            composition: aqueous,
            kind: PhaseKind::Ge,
        },
    ])
}

/// `SystemEosGE.isAqueousComponent`: an ion, water, or one of the polar production solvents
/// that partition to a brine rather than to a hydrocarbon.
///
/// **A name test and not a class test**, because that is what upstream is: `meg`, `teg`,
/// `deg`, `methanol` and `ethanol` are the five it lists beside water, and a glycol that
/// arrives under another spelling seeds to the wrong role.
fn is_aqueous_feed_component(name: &str) -> bool {
    matches!(
        name.trim().to_lowercase().as_str(),
        "water" | "meg" | "teg" | "deg" | "methanol" | "ethanol"
    )
}

/// `enforceAqueousFractionBounds`: the aqueous fraction projected onto its admissible band.
///
/// The band's floor is the ionic inventory plus a margin - an aqueous phase smaller than the
/// ions it must hold is not a state - and its ceiling is the phase floor's complement. Both
/// are also limited to a factor of [`AQUEOUS_STEP_LIMIT`] from the fraction the last
/// correction started at, which bounds how fast a trial may move.
///
/// # Errors
/// * [`AzothError::SolverNotConverged`] if the other phases cannot give up enough fraction to
///   reach the floor. Upstream throws the same state.
fn enforce_aqueous_fraction_bounds(
    phases: &mut [MultiphasePhase],
    ion: &[bool],
    z: &[f64],
    aqueous: usize,
    previous: f64,
) -> Result<()> {
    let inventory: f64 = z
        .iter()
        .zip(ion)
        .filter(|(_, is_ion)| **is_ion)
        .map(|(value, _)| value.max(0.0))
        .sum();
    if !inventory.is_finite() || inventory <= 0.0 {
        return Ok(());
    }

    let maximum = 1.0 - 2.0 * FRACTION_FLOOR;
    let stepped_minimum = if previous > 0.0 {
        previous / AQUEOUS_STEP_LIMIT
    } else {
        inventory + ION_CAPACITY_MARGIN
    };
    let minimum = maximum.min((inventory + ION_CAPACITY_MARGIN).max(stepped_minimum));
    let stepped_maximum = (previous * AQUEOUS_STEP_LIMIT).min(maximum);
    let current = phases[aqueous].fraction;
    if current >= minimum && current <= stepped_maximum {
        return Ok(());
    }

    if current > stepped_maximum {
        let excess = current - stepped_maximum;
        let others: f64 = (0..3)
            .filter(|&k| k != aqueous)
            .map(|k| phases[k].fraction)
            .sum();
        for k in (0..3).filter(|&k| k != aqueous) {
            let share = if others > 0.0 {
                phases[k].fraction / others
            } else {
                0.5
            };
            phases[k].fraction += excess * share;
        }
        phases[aqueous].fraction = stepped_maximum;
        return Ok(());
    }

    let required = minimum - current;
    let adjustable: f64 = (0..3)
        .filter(|&k| k != aqueous)
        .map(|k| (phases[k].fraction - FRACTION_FLOOR).max(0.0))
        .sum();
    if adjustable + FRACTION_FLOOR < required {
        return Err(AzothError::SolverNotConverged {
            iterations: 0,
            residual: inventory,
            tolerance: minimum,
        });
    }
    let mut remaining = required;
    let mut last = None;
    for k in (0..3).filter(|&k| k != aqueous) {
        let available = (phases[k].fraction - FRACTION_FLOOR).max(0.0);
        if available <= 0.0 {
            continue;
        }
        last = Some(k);
        let transfer = (required * available / adjustable).min(remaining);
        phases[k].fraction -= transfer;
        remaining -= transfer;
    }
    if remaining > 0.0
        && let Some(k) = last
    {
        phases[k].fraction -= remaining;
    }
    let others: f64 = (0..3)
        .filter(|&k| k != aqueous)
        .map(|k| phases[k].fraction)
        .sum();
    phases[aqueous].fraction = 1.0 - others;
    Ok(())
}

/// `TPHybridEosGeFlash.setXY`: the compositions the material balance gives at the fractions.
///
/// Three rules, not one. An **ion** is `z_i / beta` in the aqueous role and
/// [`ION_CONFINEMENT`] in the other two. A **neutral** is the plain `z_i / (E_i phi_ik)`,
/// floored. And the aqueous phase's neutrals are then scaled so that they and the ions
/// together sum to one - which is what makes the ions' share of the aqueous phase exact
/// rather than diluted by the normalisation the other two roles get.
///
/// # Errors
/// * [`AzothError::SolverNotConverged`] if the aqueous role cannot hold the ionic inventory
///   at this fraction. Upstream throws the same state.
fn set_compositions(
    phases: &mut [MultiphasePhase],
    z: &[f64],
    e: &[f64],
    inverted: &[Vec<f64>],
    ion: &[bool],
    aqueous: usize,
) -> Result<()> {
    let n = z.len();
    let count = phases.len();
    for k in 0..count {
        let fraction = phases[k].fraction.max(FRACTION_FLOOR);
        let mut ion_sum = 0.0;
        let mut neutral_sum = 0.0;
        for i in 0..n {
            let mut value = if ion[i] {
                if k == aqueous {
                    z[i] / fraction
                } else {
                    ION_CONFINEMENT
                }
            } else {
                z[i] / e[i] * inverted[k][i]
            };
            if !value.is_finite() || value <= 0.0 {
                value = ION_CONFINEMENT;
            }
            phases[k].composition[i] = value;
            if ion[i] {
                ion_sum += value;
            } else {
                neutral_sum += value;
            }
        }
        if k == aqueous {
            if ion_sum.is_nan() || ion_sum >= 1.0 || neutral_sum.is_nan() || neutral_sum <= 0.0 {
                return Err(AzothError::SolverNotConverged {
                    iterations: 0,
                    residual: ion_sum,
                    tolerance: 1.0,
                });
            }
            let neutral_total = 1.0 - ion_sum;
            for (value, is_ion) in phases[k].composition.iter_mut().zip(ion) {
                if !*is_ion {
                    *value = neutral_total * *value / neutral_sum;
                }
            }
        } else {
            let sum: f64 = phases[k].composition.iter().sum();
            if sum > 0.0 {
                for value in phases[k].composition.iter_mut() {
                    *value /= sum;
                }
            }
        }
    }
    Ok(())
}

/// The acceptance contract's two residuals, which is what the answer is *for*.
///
/// The material balance is `z_i - sum_k beta_k x_ik`; the equilibrium is the spread of
/// `ln(x_i phi_i P)` over the phases that are there, taken over the components the contract
/// admits - a zero feed fraction or an ion takes no part, because an ion's coefficient is a
/// model constant rather than an equilibrium quantity.
fn acceptance(
    beta: &[f64],
    composition: &[Vec<f64>],
    ln_phi: &[Vec<f64>],
    z: &[f64],
    ion: &[bool],
    reduced: &ReducedParameters,
) -> Result<(f64, f64)> {
    let n = z.len();
    let mut balance: f64 = 0.0;
    for i in 0..n {
        let split: f64 = (0..beta.len()).map(|k| beta[k] * composition[k][i]).sum();
        balance = balance.max((z[i] - split).abs());
    }
    let mut fugacity: f64 = 0.0;
    for i in 0..n {
        if z[i] <= 1.0e-30 || ion[i] {
            continue;
        }
        let mut first: Option<f64> = None;
        for k in 0..beta.len() {
            if beta[k] <= MATERIAL_PHASE_FRACTION {
                continue;
            }
            let value = composition[k][i] * ln_phi[k][i].exp() * reduced.pressure;
            if !value.is_finite() || value <= 0.0 {
                continue;
            }
            let logged = value.ln();
            match first {
                None => first = Some(logged),
                Some(reference) => fugacity = fugacity.max((reference - logged).abs()),
            }
        }
    }
    if !balance.is_finite() || !fugacity.is_finite() {
        return Err(AzothError::SolverNotConverged {
            iterations: 0,
            residual: balance.max(fugacity),
            tolerance: 0.0,
        });
    }
    Ok((balance, fugacity))
}
