//! `eos.ge_flash` - the isothermal flash of a cubic vapour over a named activity-coefficient
//! liquid.
//!
//! Spec: `specs/models/eos/ge_flash.toml`
//!
//! **This is `eos.ge_nrtl_flash` with the liquid's name taken out.** NeqSim's direct
//! gamma-phi route iterates `K_i = phi_i^L / phi_i^V` (`TPflash.sucsSubsGammaPhi`), and
//! which `PhaseGE` subclass supplies `phi_i^L` is the system's decision - `SystemNRTL` pairs
//! SRK with `PhaseGENRTL` and `SystemVanLaarActivitySRK` pairs it with
//! `PhaseGEVanLaarAcid`, and `requiresDirectGammaPhiFlash()` is what picks the route. So the
//! liquid is a parameter here and the loop is unchanged.
//!
//! The four liquids are the ones NeqSim's own systems configure and this crate carries a
//! phase id for: `SystemNRTL`, `SystemUNIFAC`, `SystemGEWilson` and - the only one whose
//! `requiresDirectGammaPhiFlash()` is true - `SystemVanLaarActivitySRK`. **`uniquac` is not
//! among them**: eight systems extend `SystemEosGE` and none of them is a UNIQUAC one. The kernel takes the liquid as a **closure over composition**, the way
//! `azoth-reactions` carries a phase model it does not own, so the loop below cannot tell
//! which model it is running - which is the point: a second liquid is a match arm and not a
//! second loop.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result, apply_checks};

use crate::algorithm_of;
use crate::databank::{
    ge_nrtl_phase_parameters, ge_unifac_phase_parameters, ge_van_laar_acid_phase_parameters,
    ge_wilson_phase_parameters,
};
use crate::flash_iteration::{
    Outcome, compositions, is_trivial, negative_flash_warning, rachford_rice, rachford_rice_bounds,
    rms_delta, single_phase_warning, trivial_warning,
};
use crate::ge_nrtl_phase::ge_nrtl_phase;
use crate::ge_unifac_phase::ge_unifac_phase;
use crate::ge_van_laar_acid_phase::ge_van_laar_acid_phase;
use crate::ge_wilson_phase::ge_wilson_phase;
use crate::mixture::{Mixture, RootSide, wilson_k};
use crate::model_gen;

use crate::results::{GeFlashResult, Phase};

/// The liquid half of the flash, as a closure over composition.
///
/// `ln phi_i^L` for a composition, from whichever phase model the caller resolved. The loop
/// cannot see which - the same seam `azoth-reactions` uses to carry a phase model it does
/// not own.
pub type LiquidLogPhi<'a> = &'a mut dyn FnMut(&[f64]) -> Result<Vec<f64>>;

/// Which activity-coefficient phase is the liquid side.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LiquidModel {
    /// `eos.ge_nrtl_phase`.
    Nrtl,
    /// `eos.ge_unifac_phase`.
    Unifac,
    /// `eos.ge_wilson_phase`, whose Wilson correlation reads the cubic's critical
    /// properties as well as the databank's `lambda`.
    Wilson,
    /// `eos.ge_van_laar_acid_phase` - the one NeqSim's direct gamma-phi route is reachable
    /// through.
    VanLaarAcid,
}

impl LiquidModel {
    /// The model a spec's `liquid_model` names.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] for a name this enum does not carry - a liquid this
    /// library has no phase model for is not a member, rather than a silent default.
    pub fn parse(name: &str) -> Result<Self> {
        match name.trim().to_lowercase().as_str() {
            "nrtl" => Ok(Self::Nrtl),
            "unifac" => Ok(Self::Unifac),
            "wilson" => Ok(Self::Wilson),
            "van_laar_acid" => Ok(Self::VanLaarAcid),
            other => Err(AzothError::invalid_input(
                "liquid_model",
                format!(
                    "`{other}` is not an activity-coefficient phase this library carries;                      the four are nrtl, unifac, wilson and van_laar_acid"
                ),
            )),
        }
    }
}

/// The isothermal flash of a mixture whose liquid is a named activity-coefficient phase
/// and whose vapour is a cubic.
///
/// `z` is the overall composition, in mole fractions, and is **checked rather than
/// renormalised** - silently rescaling a caller's composition would make their error
/// invisible in a way that changes every number downstream.
///
/// `names` resolves both halves - the liquid's own parameters through `liquid_model` and the
/// vapour's `Tc`, `Pc`, `omega` and `kij` through `cubic` - so the two cannot be resolved
/// from different component lists, which is the error `eos.ge_nrtl_flash`'s signature can
/// only document.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
/// * [`AzothError::InvalidInput`] if `cubic` or `liquid_model` is not one this library
///   carries, if a component is unknown, or if `z` is the wrong length, has a negative entry
///   or does not sum to one.
/// * [`AzothError::SolverNotConverged`] if the iteration hits its cap.
/// * Propagates the phase models' range checks.
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn ge_flash(
    names: &[&str],
    cubic: &str,
    liquid_model: &str,
    T: ThermodynamicTemperature,
    P: Pressure,
    z: &[f64],
) -> Result<GeFlashResult> {
    let cubic: crate::Cubic = cubic
        .parse()
        .map_err(|reason: String| AzothError::invalid_input("cubic", reason))?;
    let (mixture, _) = crate::databank::mixture_of(names, cubic, None)?;
    let mixture = &mixture;
    let model = LiquidModel::parse(liquid_model)?;
    match model {
        LiquidModel::Nrtl => {
            let params = ge_nrtl_phase_parameters(names, None)?;
            let mut liquid = |x: &[f64]| Ok(ge_nrtl_phase(&params, T.value, P.value, x)?.ln_phi);
            ge_flash_with(&mut liquid, mixture, T, P, z)
        }
        LiquidModel::Unifac => {
            let params = ge_unifac_phase_parameters(names, None)?;
            let mut liquid = |x: &[f64]| Ok(ge_unifac_phase(&params, T.value, P.value, x)?.ln_phi);
            ge_flash_with(&mut liquid, mixture, T, P, z)
        }
        LiquidModel::Wilson => {
            let params = ge_wilson_phase_parameters(names, None)?;
            let mut liquid =
                |x: &[f64]| Ok(ge_wilson_phase(&params, mixture, T.value, P.value, x)?.ln_phi);
            ge_flash_with(&mut liquid, mixture, T, P, z)
        }
        LiquidModel::VanLaarAcid => {
            let params = ge_van_laar_acid_phase_parameters(names, None)?;
            let mut liquid =
                |x: &[f64]| Ok(ge_van_laar_acid_phase(&params, T.value, P.value, x)?.ln_phi);
            ge_flash_with(&mut liquid, mixture, T, P, z)
        }
    }
}

/// The loop, with the liquid's coefficients taken as a closure.
///
/// # Errors
/// As [`ge_flash`], whose parameter resolution happens before this is called.
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn ge_flash_with(
    liquid: LiquidLogPhi<'_>,
    mixture: &Mixture,
    T: ThermodynamicTemperature,
    P: Pressure,
    z: &[f64],
) -> Result<GeFlashResult> {
    let spec = &model_gen::GE_FLASH_SPEC;
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

    let n = mixture.len();
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
                "the feed's mole fractions sum to {sum}, not to one. Renormalising it \
                 here would make a composition error invisible in every number \
                 downstream, so it is refused instead"
            ),
        ));
    }

    let min_t_over_tc = mixture
        .components()
        .iter()
        .map(|c| T.value / c.tc.value)
        .fold(f64::INFINITY, f64::min);
    apply_checks(
        spec.derived_checks(),
        |quantity| (quantity == "min_t_over_tc").then_some(min_t_over_tc),
        &mut warnings,
    )?;

    let reduced = mixture.reduced_parameters(T, P)?;
    warnings.extend(reduced.warnings.iter().cloned());

    let algorithm = algorithm_of(spec)?;
    let inner = algorithm.inner.ok_or_else(|| AzothError::InvalidInput {
        field: "algorithm.inner".to_string(),
        reason: format!(
            "scheme `{}` solves Rachford-Rice every iteration but the spec declares \
             no inner scheme",
            algorithm.scheme
        ),
    })?;

    let mut k = wilson_k(mixture, T, P);
    let mut iterations = 0;
    let mut residual = f64::NAN;

    // How the loop finished, when it finished somewhere other than convergence.
    let mut settled: Option<Outcome> = None;

    for step in 1..=algorithm.max_iterations {
        iterations = step;

        // Both guards run *before* the Rachford-Rice solve, because both are cases
        // where the bracket is a division by zero: `K_i = 1` puts a pole at
        // infinity, and K-values all on one side of one put a bracket end there.
        if is_trivial(&k) {
            settled = Some(Outcome::Trivial);
            break;
        }
        // A K-value that is not finite and positive is the iteration diverging
        // rather than a state to diagnose. Guarded because every path below treats
        // the K-values as a composition ratio, and `NaN > 1.0` being false would
        // otherwise report a diverged iteration as a single-phase feed.
        if k.iter().any(|value| !value.is_finite() || *value <= 0.0) {
            return Err(AzothError::SolverNotConverged {
                iterations,
                residual,
                tolerance: algorithm.tolerance,
            });
        }
        let Some(bounds) = rachford_rice_bounds(&k) else {
            settled = Some(Outcome::SinglePhase(if k.iter().all(|&v| v > 1.0) {
                Phase::AllVapour
            } else {
                Phase::AllLiquid
            }));
            break;
        };

        let (beta, _) = rachford_rice(z, &k, bounds, inner.tolerance, inner.max_iterations);
        let (x, y) = compositions(z, &k, beta);

        // The liquid's coefficients are whichever phase model was resolved, and the
        // vapour's are the cubic's. This pair of lines is the whole of what makes this
        // flash a gamma-phi flash rather than `eos.pt_flash`.
        let ln_phi_liquid = liquid(&x)?;
        let vapour = mixture.phase_state(&reduced, &y, RootSide::Vapour)?;

        let ln_k_new: Vec<f64> = ln_phi_liquid
            .iter()
            .zip(&vapour.ln_phi)
            .map(|(&l, &v)| l - v)
            .collect();
        residual = rms_delta(&ln_k_new, &k);
        k = ln_k_new.iter().map(|value| value.exp()).collect();

        if residual <= algorithm.tolerance {
            break;
        }
    }

    // Classify what the loop settled on. Re-evaluated once more at the converged
    // `K`, so the returned `beta`, `x`, `y` and `K` are one consistent state rather
    // than one state and its predecessor's vapour fraction.
    let (beta, x, y, phase) = match settled {
        // The trivial solution is stated *exactly* - `x = y = z` and `K = 1` -
        // rather than as the last iterate that approached it. That iterate is a
        // function of where the bisection stopped on an identically-zero function
        // and is not reproducible; this is.
        Some(Outcome::Trivial) => {
            k = vec![1.0; n];
            warnings.push(trivial_warning());
            (None, z.to_vec(), z.to_vec(), Phase::Trivial)
        }
        // No root at all: the feed is single phase by proof, and `beta` is absent
        // because there is no vapour fraction to report.
        Some(Outcome::SinglePhase(phase)) => {
            warnings.push(single_phase_warning(phase));
            (None, z.to_vec(), z.to_vec(), phase)
        }
        None if residual > algorithm.tolerance => {
            return Err(AzothError::SolverNotConverged {
                iterations,
                residual,
                tolerance: algorithm.tolerance,
            });
        }
        None => {
            if is_trivial(&k) {
                k = vec![1.0; n];
                warnings.push(trivial_warning());
                (None, z.to_vec(), z.to_vec(), Phase::Trivial)
            } else {
                let bounds = rachford_rice_bounds(&k).ok_or(AzothError::SolverNotConverged {
                    iterations,
                    residual,
                    tolerance: algorithm.tolerance,
                })?;
                let (beta, _) = rachford_rice(z, &k, bounds, inner.tolerance, inner.max_iterations);
                let (x, y) = compositions(z, &k, beta);
                let phase = match beta {
                    value if value < 0.0 => Phase::AllLiquid,
                    value if value > 1.0 => Phase::AllVapour,
                    _ => Phase::TwoPhase,
                };
                if phase != Phase::TwoPhase {
                    warnings.push(negative_flash_warning(phase, beta));
                }
                (Some(beta), x, y, phase)
            }
        }
    };

    // The one state the result reports, evaluated at the compositions actually
    // returned rather than at the last iterate's. A feed refused by the bracket
    // before its first evaluation reaches this with `x = y = z`, which is the
    // single-phase answer it was classified as.
    let ln_phi_liquid = liquid(&x)?;
    let vapour = mixture.phase_state(&reduced, &y, RootSide::Vapour)?;

    Ok(GeFlashResult {
        beta,
        x,
        y,
        k,
        ln_phi_liquid,
        ln_phi_vapour: vapour.ln_phi,
        z_vapour: vapour.z,
        min_t_over_tc,
        phase,
        iterations,
        residual,
        warnings,
    })
}
