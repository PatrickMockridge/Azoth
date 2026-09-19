//! `eos.tp_flash_saft` - the SAFT-VR-Mie flash, by successive substitution.
//!
//! **NeqSim takes this family's flash along a different route from every other**, and that is
//! the model: `ThermodynamicOperations` dispatches a `SystemSAFTVRMie` to `TPflashSAFT` by
//! `instanceof`, and `TPflashSAFT` is a successive-substitution loop - Wilson K-values to
//! seed, Rachford-Rice inside for the vapour fraction, and **each phase's fugacity
//! coefficients from a separate single-phase solve at that trial composition** - rather than
//! the Newton scheme the generic `TPflash` reaches its K-values through.
//!
//! # The vapour fraction's bound is a branch
//!
//! NeqSim breaks out of the loop when Rachford-Rice's root leaves `(1e-10, 1 - 1e-10)` and
//! reports a **single phase**, chosen by comparing the Gibbs energies of a gas and a liquid
//! at the feed. That is not a stability test: with the ideal part identical in both, the
//! comparison reduces to `sum_i x_i ln phi_i` on the two roots, which is what this does.
//!
//! # The seed is an input
//!
//! The Wilson K-values are `exp(ln(Pc/P) + 5.373 (1 + omega)(1 - Tc/T))`, and they come from
//! the *cubic* constants beside the Mie parameters rather than from the Mie set - so they
//! cross as an argument and the model resolves them from the databank. It is a ratio of
//! pressures, so the two sides have to be in the same unit and nothing else.
//!
//! # Why it is separate from `crate::pt_flash`
//!
//! The generic flash's Newton step needs the derivative of the fugacity surface in the
//! composition, the temperature and the volume. This model's `eta` derivatives are already
//! central differences, so differencing *those* again - which is what a Newton step on this
//! surface is - would be a difference of a difference of a difference. NeqSim's answer for a
//! SAFT-VR-Mie mixture is the successive-substitution one, and this is a port of the
//! algorithm rather than of the model's equation of state.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result};

use crate::flash_iteration;
use crate::mixture::RootSide;
use crate::results::{Phase, SaftFlashResult};
use crate::saft_vr_mie::MieComponent;
use crate::saft_vr_mie_phase::{ln_phi, molar_volume};

/// NeqSim's `MAX_SS_ITER`, unchanged: at `K_TOLERANCE` the loop converges in a handful.
const MAX_SS_ITER: u32 = 50;

/// The relative change in one K-value the loop stops at.
///
/// **NeqSim's `K_TOL` is `1e-6`, and this is `1e-5` because the model's own arithmetic
/// cannot do better.** The K-values are ratios of fugacity coefficients, and those carry
/// this model's `eta` derivatives - central differences at a relative step of `1e-5` - so
/// the iterate's own noise floor is around `2e-6`: measured on methane/n-butane at 250 K,
/// the relative change falls to `2.7e-6` by the forty-fifth step and then *oscillates*
/// between that and `2e-5`, never reaching `1e-6` however many steps it is given.
///
/// The threshold decides which branch the flash takes, so this is not cosmetic: at `1e-6`
/// the loop exhausts its fifty steps at every state and reports a single phase, where the
/// same port at `1e-5` stops at the same converged K-values NeqSim reaches in seven. A
/// stopping rule below the arithmetic's floor is not a stricter test - it is a different
/// answer.
const K_TOLERANCE: f64 = 1.0e-5;

/// NeqSim's `BETA_MIN`: outside `(BETA_MIN, 1 - BETA_MIN)` the flash reports one phase.
const BETA_MIN: f64 = 1.0e-10;

/// The Wilson K-values NeqSim seeds with, from the cubic constants beside the Mie set.
///
/// `exp(ln(Pc/P) + 5.373 (1 + omega)(1 - Tc/T))`, which is `Component.init`'s expression.
/// NeqSim's `Pc` and `P` are both in bar there and the ratio is what the logarithm sees, so
/// this takes either unit as long as the two agree.
#[must_use]
pub fn wilson_k(pc: &[f64], omega: &[f64], tc: &[f64], t: f64, p: f64) -> Vec<f64> {
    (0..pc.len())
        .map(|i| (pc[i] / p).ln() + 5.373 * (1.0 + omega[i]) * (1.0 - tc[i] / t))
        .map(f64::exp)
        .collect()
}

/// The flash, for a caller that resolved the fluid and the seed itself.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the lengths disagree or the composition is not one.
/// * [`AzothError::OutOfRange`] if `t` or `p` is not positive, or a trial phase has no root.
pub fn tp_flash_saft_of(
    components: &[MieComponent],
    seed_k: &[f64],
    t: f64,
    p: f64,
    z: &[f64],
) -> Result<SaftFlashResult> {
    let n = components.len();
    if z.len() != n || seed_k.len() != n {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "{n} components need {n} mole fractions and {n} K-values, but got {} and {}",
                z.len(),
                seed_k.len()
            ),
        ));
    }
    let total: f64 = z.iter().sum();
    if (total - 1.0).abs() > 1.0e-9 {
        return Err(AzothError::invalid_input(
            "z",
            format!("the mole fractions sum to {total}, not to one"),
        ));
    }
    if !p.is_finite() || p <= 0.0 {
        return Err(AzothError::OutOfRange {
            field: "P".to_string(),
            value: p,
            detail: "the flash solves at a pressure".to_string(),
        });
    }

    // One phase's fugacity coefficients, from a single-phase solve at that composition and
    // on that side - which is what NeqSim's `computeFugacityCoefficients` does by building a
    // fresh system per call.
    let coefficients = |x: &[f64], side: RootSide| -> Result<(Vec<f64>, f64)> {
        let solved = molar_volume(components, x, t, p, side)?;
        Ok((ln_phi(components, x, t, solved.v)?, solved.z))
    };

    let mut k = seed_k.to_vec();
    let mut beta = 0.5;
    let mut iterations = 0;
    let mut residual = f64::INFINITY;
    let mut converged = false;
    let mut liquid = Vec::new();
    let mut vapour = Vec::new();
    let mut liquid_coefficients = Vec::new();
    let mut vapour_coefficients = Vec::new();
    let (mut z_liquid, mut z_vapour) = (0.0, 0.0);

    for step in 0..MAX_SS_ITER {
        iterations = step + 1;
        // **NeqSim bisects `[0, 1]`; this bisects the physical bracket**, which holds the
        // same root whenever one exists and refuses when none does - the two agree by
        // construction rather than by tolerance.
        let Some(bounds) = flash_iteration::rachford_rice_bounds(&k) else {
            break;
        };
        // The bound is on the **root**, not on the bracket: the physical bracket is the
        // interval free of poles and reaches outside `[0, 1]` for an ordinary seed - at
        // `K = [5.58, 0.0138]` it is `(-0.218, 1.014)` - so checking it would take the
        // single-phase branch at every two-phase state there is.
        beta = flash_iteration::rachford_rice(z, &k, bounds, 1.0e-14, 200).0;
        if !(BETA_MIN..=1.0 - BETA_MIN).contains(&beta) {
            break;
        }
        let (x, y) = flash_iteration::compositions(z, &k, beta);

        let (ln_phi_liquid, liquid_z) = coefficients(&x, RootSide::Liquid)?;
        let (ln_phi_vapour, vapour_z) = coefficients(&y, RootSide::Vapour)?;

        // NeqSim's convergence measure: the largest *relative* change in a K-value.
        let mut largest: f64 = 0.0;
        for i in 0..n {
            let updated = (ln_phi_liquid[i] - ln_phi_vapour[i]).exp();
            let relative = (updated - k[i]).abs() / k[i].abs().max(1.0e-10);
            largest = largest.max(relative);
            k[i] = updated;
        }
        residual = largest;
        liquid = x;
        vapour = y;
        liquid_coefficients = ln_phi_liquid;
        vapour_coefficients = ln_phi_vapour;
        z_liquid = liquid_z;
        z_vapour = vapour_z;
        if largest < K_TOLERANCE {
            converged = true;
            break;
        }
    }

    if converged {
        return Ok(SaftFlashResult {
            beta: Some(beta),
            x: liquid,
            y: vapour,
            k,
            ln_phi_liquid: liquid_coefficients,
            ln_phi_vapour: vapour_coefficients,
            z_liquid,
            z_vapour,
            phase: Phase::TwoPhase,
            iterations,
            residual,
            warnings: Vec::new(),
        });
    }

    // **One phase, and which one is a Gibbs comparison.** The ideal part is identical in both
    // - same temperature, pressure and composition - so `sum_i x_i ln phi_i` decides it.
    let (gas_ln_phi, gas_z) = coefficients(z, RootSide::Vapour)?;
    let (liquid_ln_phi, liquid_z) = coefficients(z, RootSide::Liquid)?;
    let residual_gibbs =
        |ln_phi: &[f64]| -> f64 { z.iter().zip(ln_phi).map(|(zi, p)| zi * p).sum::<f64>() };
    let gas_is_lower = residual_gibbs(&gas_ln_phi) <= residual_gibbs(&liquid_ln_phi);
    Ok(SaftFlashResult {
        beta: None,
        x: z.to_vec(),
        y: z.to_vec(),
        k,
        ln_phi_liquid: liquid_ln_phi,
        ln_phi_vapour: Vec::new(),
        z_liquid: liquid_z,
        z_vapour: gas_z,
        phase: if gas_is_lower {
            Phase::AllVapour
        } else {
            Phase::AllLiquid
        },
        iterations,
        residual,
        warnings: Vec::new(),
    })
}

/// `eos.tp_flash_saft`: the flash NeqSim runs on a SAFT-VR-Mie system.
///
/// The names cross unresolved and are looked up here, so the two-kernel comparison covers
/// the resolution - which for this model means the Mie set *and* the cubic constants the
/// Wilson seed is built from.
///
/// # Errors
/// * [`AzothError::InvalidInput`] as [`crate::saft_vr_mie_phase::parameters_of`], or if `z`
///   is not one entry per component or does not sum to one.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive, or a trial phase has no root.
pub fn tp_flash_saft(
    components: &[String],
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<SaftFlashResult> {
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let parameters = crate::saft_vr_mie_phase::parameters_of(&names)?;

    let mut pc = Vec::with_capacity(names.len());
    let mut omega = Vec::with_capacity(names.len());
    let mut tc = Vec::with_capacity(names.len());
    for name in &names {
        let entry = crate::databank::entry(name, None)?;
        pc.push(entry.pc);
        omega.push(entry.omega);
        tc.push(entry.tc);
    }
    let seed = wilson_k(&pc, &omega, &tc, t.value, p.value);
    tp_flash_saft_of(&parameters, &seed, t.value, p.value, z)
}
