//! `eos.co2_brine_phase_equilibrium` - NeqSim's constrained CO2/brine phase selection.
//!
//! Spec: `specs/models/eos/co2_brine_phase_equilibrium.toml`, which records the supported
//! feed, the normalised tangent-plane trial, the ion-confined split and the four checks the
//! answer is accepted on.
//!
//! Ported from `CO2BrinePhaseEquilibrium`, a class of 383 lines that owns its whole solve.
//! **It never calls `TPflash`**: the reference is one aqueous phase at the feed's own
//! composition, each trial is a single phase on a requested root, and the split is its own
//! Rachford-Rice bisection in which an ion's `K` is zero. That is why this is a model of its
//! own rather than a call into [`crate::tp_multiflash`] - the coefficients are fugacity
//! coefficients at compositions the class has already fixed, not a flash's.
//!
//! # The Gibbs comparison, and why it is a mixing term
//!
//! Upstream orders the candidate splits by `SystemInterface.getGibbsEnergy()`. This library
//! has no absolute Gibbs energy for this phase - [`crate::electrolyte_cpa_statoil_phase`]'s
//! spec says so, because the electrolyte terms' departure is not assembled - but **the
//! ordering does not need one**. At fixed `T`, `P` and inventory,
//!
//! ```text
//! G = sum_i n_i mu_i^0(T) + R T [ n_total ln P + sum_phases sum_i n_iphi ln(x_iphi phi_iphi) ]
//! ```
//!
//! and the first two terms are the same for every candidate, so ordering by the bracket
//! orders the candidates exactly as `G` does. What is *not* ported is the class's
//! `gibbs > aqueousGibbs + max(1e-7, |G| 1e-12)` guard, which is written in the absolute
//! energy this library cannot form; the spec names that and what qualifies a split instead.

use azoth_core::Result;
use azoth_core::units::{Pressure, ThermodynamicTemperature};

use crate::mixture::{Mixture, ReducedParameters, RootSide};
use crate::model_gen;
use crate::results::Co2BrinePhaseEquilibriumResult;

/// The class's `FLOOR`: what an ion's mole fraction is written as outside the aqueous phase,
/// and the floor every composition update takes.
const FLOOR: f64 = 1.0e-50;

/// The class's `TOLERANCE`, which decides whether a trial is unstable and how far apart the
/// two phases' molecular fugacities may be.
const TOLERANCE: f64 = 1.0e-8;

/// The composition change that ends the trial and the split, the class's `1e-12`.
const CHANGE_TOLERANCE: f64 = 1.0e-12;

/// The class's cap on the normalised trial.
const MAX_TRIAL_STEPS: usize = 100;

/// The class's cap on the split.
const MAX_SPLIT_STEPS: usize = 150;

/// The class's bisection depth for `beta`.
const BISECTION_STEPS: usize = 80;

/// The class's seed for the rich phase's fraction, which is the value the first split pass
/// measures its change against.
const SEED_BETA: f64 = 0.1;

/// The class's upper end of the `beta` bracket. One is the single-aqueous answer and is
/// excluded, because the question is whether a *second* phase is there.
const BETA_CEILING: f64 = 1.0 - 1.0e-12;

/// The worst `|z_i - sum_k beta_k x_ik|` an accepted answer may have, upstream's `1e-10`.
const CONSERVATION_TOLERANCE: f64 = 1.0e-10;

/// One phase's state, as the solve reads it: a composition and its log fugacity coefficients.
struct Layer {
    x: Vec<f64>,
    ln_phi: Vec<f64>,
    /// `ln(x_i phi_i)`, the term the mixing Gibbs is built from.
    ln_x_phi: Vec<f64>,
}

impl Layer {
    fn of(
        mixture: &Mixture,
        reduced: &ReducedParameters,
        x: &[f64],
        side: RootSide,
    ) -> Result<Self> {
        let state = mixture.phase_state(reduced, x, side)?;
        let ln_x_phi = x
            .iter()
            .zip(&state.ln_phi)
            .map(|(fraction, ln_phi)| fraction.ln() + ln_phi)
            .collect();
        Ok(Self {
            x: x.to_vec(),
            ln_phi: state.ln_phi,
            ln_x_phi,
        })
    }

    /// `sum_i x_i ln(x_i phi_i)`, one phase's part of the mixing Gibbs.
    ///
    /// A component the phase does not carry contributes nothing rather than `0 * -inf`.
    fn mixing_term(&self) -> f64 {
        self.x
            .iter()
            .zip(&self.ln_x_phi)
            .filter(|(x, _)| **x > 0.0)
            .map(|(x, l)| x * l)
            .sum()
    }
}

/// The feed's mixing Gibbs over `R T`, per mole of feed.
fn mixing_gibbs(layers: &[&Layer], fractions: &[f64]) -> f64 {
    layers
        .iter()
        .zip(fractions)
        .map(|(layer, beta)| beta * layer.mixing_term())
        .sum()
}

/// The constrained Rachford-Rice residual, `sum_i z_i (K_i - 1) / (1 + beta (K_i - 1))`.
fn rachford_rice(z: &[f64], ratios: &[f64], beta: f64) -> f64 {
    z.iter()
        .zip(ratios)
        .map(|(overall, ratio)| {
            let difference = ratio - 1.0;
            overall * difference / (1.0 + beta * difference)
        })
        .sum()
}

/// The constrained equilibrium of a CO2/water/ion feed with the electrolyte-CPA equation of
/// state, at a temperature and pressure.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if `moles` is not one entry per component or is
///   not a feed, a name is not in the databank, or the feed is outside the class's scope -
///   which is CO2, water and ions alone, with more water than CO2.
/// * [`azoth_core::AzothError::PropertyUnavailable`] if a component carries no heat-capacity
///   coefficients.
/// * [`azoth_core::AzothError::OutOfRange`] if `T` or `P` is not positive, an unstable feed's
///   split does not converge, or the answer fails its acceptance checks.
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn co2_brine_phase_equilibrium(
    components: &[String],
    T: ThermodynamicTemperature,
    P: Pressure,
    moles: &[f64],
) -> Result<Co2BrinePhaseEquilibriumResult> {
    let spec = &model_gen::CO2_BRINE_PHASE_EQUILIBRIUM_SPEC;
    let mut warnings = Vec::new();
    azoth_core::apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T.value),
            "P" => Some(P.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let n = components.len();
    if moles.len() != n {
        return Err(azoth_core::AzothError::invalid_input(
            "moles",
            format!("a feed of {n} components has {} mole numbers", moles.len()),
        ));
    }
    if moles.iter().any(|m| !m.is_finite() || *m < 0.0) {
        return Err(azoth_core::AzothError::invalid_input(
            "moles",
            "every mole number must be finite and not negative",
        ));
    }
    let total: f64 = moles.iter().sum();
    if total <= 0.0 || !total.is_finite() {
        return Err(azoth_core::AzothError::invalid_input(
            "moles",
            format!("the feed's mole numbers total {total}, which is not a feed"),
        ));
    }
    let z: Vec<f64> = moles.iter().map(|m| m / total).collect();

    // **The supported composition, and it is a refusal and not a flag.** Upstream's
    // `hasSupportedComposition` asks for CO2, water and ions alone with more water than CO2.
    // Its other three tests are about a mutable system this library does not have - a
    // solid-phase check, a forced phase type, and a chemical-reaction state that the
    // reactive sibling is the model for.
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let charges: Vec<f64> = names
        .iter()
        .map(|name| crate::databank::entry(name, None).map(|entry| entry.ionic_charge))
        .collect::<Result<Vec<_>>>()?;
    let is_ion: Vec<bool> = charges.iter().map(|charge| *charge != 0.0).collect();
    let co2 = names
        .iter()
        .position(|name| name.eq_ignore_ascii_case("CO2"))
        .ok_or_else(|| unsupported("it has no CO2"))?;
    let water = names
        .iter()
        .position(|name| name.eq_ignore_ascii_case("water"))
        .ok_or_else(|| unsupported("it has no water"))?;
    if names.iter().enumerate().any(|(i, name)| {
        !is_ion[i]
            && !name.eq_ignore_ascii_case("CO2")
            && !name.eq_ignore_ascii_case("water")
            && moles[i] > 0.0
    }) {
        return Err(unsupported(
            "it carries a substance that is neither CO2, water nor an ion",
        ));
    }
    if moles[co2] <= 0.0 || moles[water] <= moles[co2] {
        return Err(unsupported("its water is not in excess of its CO2"));
    }

    let (mixture, _) = crate::furst_electrolyte::electrolyte_cpa_statoil_mixture_of(&names, None)?;
    let reduced = mixture.reduced_parameters(T, P)?;
    warnings.extend(reduced.warnings.iter().cloned());

    // The reference: one aqueous phase at the feed's own composition. A non-finite
    // coefficient there surfaces in the trial's own weight check rather than here.
    let reference = Layer::of(&mixture, &reduced, &z, RootSide::Liquid)?;

    let mut iterations = 0usize;
    let mut distance = f64::INFINITY;
    let mut best: Option<(RootSide, &'static str, f64, Layer, Layer, f64)> = None;
    for (side, label) in [(RootSide::Vapour, "gas"), (RootSide::Liquid, "oil")] {
        let (trial_distance, seed, trial_steps) =
            normalised_trial(&mixture, &reduced, &z, &is_ion, &reference, side)?;
        iterations += trial_steps;
        distance = distance.min(trial_distance);
        if trial_distance >= -TOLERANCE {
            continue;
        }
        let Some((beta, rich, aqueous, split_steps)) =
            solve_split(&mixture, &reduced, &z, &is_ion, &seed, side)?
        else {
            return Err(unqualified(MAX_SPLIT_STEPS, f64::NAN, CHANGE_TOLERANCE));
        };
        iterations += split_steps;
        let gibbs = mixing_gibbs(&[&rich, &aqueous], &[beta, 1.0 - beta]);
        let better = match &best {
            Some((_, _, incumbent, ..)) => gibbs < *incumbent,
            None => true,
        };
        if better {
            best = Some((side, label, gibbs, rich, aqueous, beta));
        }
    }

    let (accepted_root, phase_fractions, layers): (&str, Vec<f64>, Vec<&Layer>) = match &best {
        Some((_, label, _, rich, aqueous, beta)) => {
            (label, vec![*beta, 1.0 - *beta], vec![rich, aqueous])
        }
        // Neither root was unstable: the feed is the single aqueous phase the reference is,
        // at a fraction of one.
        None => ("aqueous", vec![1.0], vec![&reference]),
    };

    let max_material_balance_residual = (0..n)
        .map(|i| {
            let recovered: f64 = layers
                .iter()
                .zip(&phase_fractions)
                .map(|(layer, beta)| beta * layer.x[i])
                .sum();
            (z[i] - recovered).abs()
        })
        .fold(0.0_f64, f64::max);
    if max_material_balance_residual > CONSERVATION_TOLERANCE {
        return Err(unqualified(
            iterations,
            max_material_balance_residual,
            CONSERVATION_TOLERANCE,
        ));
    }

    let mut max_log_fugacity_residual: f64 = 0.0;
    if layers.len() == 2 {
        if layers[0].x[co2] <= 0.5 {
            return Err(unqualified(iterations, 0.5 - layers[0].x[co2], 0.5));
        }
        for i in 0..n {
            if is_ion[i] || z[i] <= FLOOR {
                continue;
            }
            let ratio = layers[0].ln_x_phi[i] - layers[1].ln_x_phi[i];
            max_log_fugacity_residual = max_log_fugacity_residual.max(ratio.abs());
            if ratio.abs() > TOLERANCE {
                return Err(unqualified(
                    iterations,
                    max_log_fugacity_residual,
                    TOLERANCE,
                ));
            }
        }
    }

    Ok(Co2BrinePhaseEquilibriumResult {
        phase_count: phase_fractions.len() as u32,
        accepted_root: accepted_root.to_string(),
        phase_fractions,
        x: layers.iter().map(|layer| layer.x.clone()).collect(),
        ln_phi: layers.iter().map(|layer| layer.ln_phi.clone()).collect(),
        minimum_trial_distance: distance,
        iterations: iterations as u32,
        max_material_balance_residual,
        max_log_fugacity_residual,
        warnings,
    })
}

/// The normalised tangent-plane trial on one root.
///
/// `w_i = z_i phi_i^reference / phi_i^trial` for a molecular component and zero for an ion,
/// which is the class's `f_i^aqueous / (P phi_i^trial)` with the `P` cancelled - the
/// reference's fugacity carries one and the trial's denominator carries one. `x_i` is
/// `max(FLOOR, w_i / sum w)` and the distance is `-ln(sum w)`.
///
/// The seed is the class's own: CO2-rich, water-poor, ions at the floor. The first pass
/// overwrites every entry from the weights, so it fixes the phase the first coefficients are
/// read at and nothing else.
fn normalised_trial(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    z: &[f64],
    is_ion: &[bool],
    reference: &Layer,
    side: RootSide,
) -> Result<(f64, Vec<f64>, usize)> {
    let mut x: Vec<f64> = (0..z.len())
        .map(|i| {
            if is_ion[i] {
                FLOOR
            } else if i == 0 {
                0.9999
            } else {
                0.0001
            }
        })
        .collect();

    let mut change: f64 = 0.0;
    let mut weights = vec![0.0; z.len()];
    for step in 1..=MAX_TRIAL_STEPS {
        let trial = Layer::of(mixture, reduced, &x, side)?;
        let mut sum = 0.0;
        for i in 0..z.len() {
            if is_ion[i] || z[i] <= FLOOR {
                weights[i] = 0.0;
            } else {
                weights[i] = z[i] * reference.ln_phi[i].exp() / trial.ln_phi[i].exp();
                if weights[i] <= 0.0 || !weights[i].is_finite() {
                    return Err(unqualified(step, weights[i], CHANGE_TOLERANCE));
                }
            }
            sum += weights[i];
        }
        if sum <= 0.0 || !sum.is_finite() {
            return Err(unqualified(step, sum, CHANGE_TOLERANCE));
        }
        change = 0.0;
        for i in 0..z.len() {
            let updated = (weights[i] / sum).max(FLOOR);
            change = change.max((updated - x[i]).abs());
            // `Component.setx(0)` keeps the old value upstream, so an ion-free trial needs
            // the floor stated rather than a zero.
            x[i] = updated;
        }
        if change < CHANGE_TOLERANCE {
            return Ok((-sum.ln(), x, step));
        }
    }
    Err(unqualified(MAX_TRIAL_STEPS, change, CHANGE_TOLERANCE))
}

/// The ion-confined split, or `None` where no bracket exists.
///
/// `K_i` is `phi_i^aqueous / phi_i^rich` and an ion's is zero. `beta` is a bisection on the
/// constrained Rachford-Rice - not a Newton, because the zero `K` puts the function's
/// curvature on the wrong side for one.
fn solve_split(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    z: &[f64],
    is_ion: &[bool],
    seed: &[f64],
    side: RootSide,
) -> Result<Option<(f64, Layer, Layer, usize)>> {
    let n = z.len();
    let mut rich = Layer::of(mixture, reduced, seed, side)?;
    let mut aqueous = Layer::of(mixture, reduced, z, RootSide::Liquid)?;
    let mut beta = SEED_BETA;
    let mut ratios = vec![0.0; n];

    for step in 1..=MAX_SPLIT_STEPS {
        for i in 0..n {
            if is_ion[i] {
                ratios[i] = 0.0;
                continue;
            }
            let ratio = aqueous.ln_phi[i].exp() / rich.ln_phi[i].exp();
            if ratio <= 0.0 || !ratio.is_finite() {
                return Ok(None);
            }
            ratios[i] = ratio;
        }
        // A bracket exists only where the residual is positive at zero and negative at
        // the ceiling. Both ends are checked finite first, so a NaN cannot pass as a
        // bracket the way a plain ordered comparison would let it.
        let (at_zero, at_ceiling) = (
            rachford_rice(z, &ratios, 0.0),
            rachford_rice(z, &ratios, BETA_CEILING),
        );
        if !at_zero.is_finite() || !at_ceiling.is_finite() || at_zero <= 0.0 || at_ceiling >= 0.0 {
            return Ok(None);
        }
        let (mut lower, mut upper) = (0.0, BETA_CEILING);
        for _ in 0..BISECTION_STEPS {
            let middle = 0.5 * (lower + upper);
            if rachford_rice(z, &ratios, middle) > 0.0 {
                lower = middle;
            } else {
                upper = middle;
            }
        }
        let updated_beta = 0.5 * (lower + upper);
        let mut change = (updated_beta - beta).abs();
        beta = updated_beta;

        let mut aqueous_x = Vec::with_capacity(n);
        let mut rich_x = Vec::with_capacity(n);
        for i in 0..n {
            let liquid = z[i] / (1.0 + beta * (ratios[i] - 1.0));
            let vapour = (ratios[i] * liquid).max(FLOOR);
            change = change.max((liquid - aqueous.x[i]).abs());
            change = change.max((vapour - rich.x[i]).abs());
            aqueous_x.push(liquid.max(FLOOR));
            rich_x.push(vapour);
        }
        aqueous = Layer::of(mixture, reduced, &aqueous_x, RootSide::Liquid)?;
        rich = Layer::of(mixture, reduced, &rich_x, side)?;
        if change < CHANGE_TOLERANCE {
            return Ok(Some((beta, rich, aqueous, step)));
        }
    }
    Ok(None)
}

/// The refusal a feed outside the class's scope gets.
fn unsupported(why: &str) -> azoth_core::AzothError {
    azoth_core::AzothError::invalid_input(
        "components",
        format!(
            "the feed is outside `CO2BrinePhaseEquilibrium`'s scope, which is CO2, water and \
             ions alone with water in excess: {why}"
        ),
    )
}

/// The error a state that cannot be qualified gets.
///
/// The four acceptance checks are the spec's contract, and a state that fails one is a solve
/// that did not reach an equilibrium - upstream's `IllegalStateException`. The residual
/// carries the number that failed, so a caller sees which check it was.
fn unqualified(iterations: usize, residual: f64, tolerance: f64) -> azoth_core::AzothError {
    azoth_core::AzothError::SolverNotConverged {
        iterations: iterations as u32,
        residual,
        tolerance,
    }
}
