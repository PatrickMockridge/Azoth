//! `eos.ge_nrtl_flash` - the isothermal flash of an SRK vapour over an NRTL liquid.
//!
//! Spec: `specs/models/eos/ge_nrtl_flash.toml`
//!
//! NeqSim's `SystemNRTL` pairs `PhaseSrkEos` with `PhaseGENRTL`, and the K-value its
//! gamma-phi iteration updates is `K_i = phi_i^L / phi_i^V` - see `TPflash`'s
//! `sucsSubsGammaPhi`. For an activity-coefficient liquid `phi_i^L` is
//! `gamma_i P0_i / P`, which is [`crate::ge_nrtl_phase`] exactly, so the liquid half of
//! this model is that model and the vapour half is the cubic.
//!
//! The scaffold is [`crate::flash_iteration`]'s, shared with [`crate::pt_flash`]; what
//! is here is the one line that differs, and the phase evaluations either side of it.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result, apply_checks};

use crate::algorithm_of;
use crate::databank::GeNrtlPhaseParameters;
use crate::flash_iteration::{
    Outcome, compositions, is_trivial, negative_flash_warning, rachford_rice, rachford_rice_bounds,
    rms_delta, single_phase_warning, trivial_warning,
};
use crate::ge_nrtl_phase::ge_nrtl_phase;
use crate::mixture::{Mixture, RootSide, wilson_k};
use crate::model_gen;

use crate::results::{GeNrtlFlashResult, Phase};

/// The isothermal flash of a mixture whose liquid is an NRTL activity-coefficient
/// phase and whose vapour is a cubic.
///
/// `z` is the overall composition, in mole fractions, and is **checked rather than
/// renormalised** - silently rescaling a caller's composition would make their error
/// invisible in a way that changes every number downstream.
///
/// `params` is the resolved NRTL phase parameters, by component, from
/// [`crate::databank::ge_nrtl_phase_parameters`]. `mixture` is the same components
/// resolved for the cubic; the two are separate arguments because the vapour reads
/// `Tc`, `Pc`, `omega` and `kij` and the liquid reads the NRTL matrices and the
/// Antoine correlations, and a caller that resolved them from different component
/// lists has made an error this signature cannot see.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
/// * [`AzothError::InvalidInput`] if `z` is the wrong length, has a negative entry,
///   does not sum to one, or if `params` and `mixture` disagree in length.
/// * [`AzothError::SolverNotConverged`] if the iteration hits its cap.
/// * Propagates the phase models' range checks.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, pascals};
/// use azoth_eos::databank::{ge_nrtl_phase_parameters, mixture_of};
/// use azoth_eos::ge_nrtl_flash;
///
/// let names = ["methanol", "water"];
/// let params = ge_nrtl_phase_parameters(&names, None)?;
/// let (mixture, _) = mixture_of(&names, None)?;
/// let r = ge_nrtl_flash(&params, &mixture, kelvins(350.0), pascals(1.0e5), &[0.5, 0.5])?;
/// assert_eq!(r.phase, azoth_eos::Phase::TwoPhase);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn ge_nrtl_flash(
    params: &GeNrtlPhaseParameters,
    mixture: &Mixture,
    T: ThermodynamicTemperature,
    P: Pressure,
    z: &[f64],
) -> Result<GeNrtlFlashResult> {
    let spec = &model_gen::GE_NRTL_FLASH_SPEC;
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
    if params.antoine.len() != n {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "the phase has {} component(s) but the cubic has {n}",
                params.antoine.len()
            ),
        ));
    }
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

        // The liquid's coefficients are the NRTL phase's, and the vapour's are the
        // cubic's. This pair of lines is the whole of what makes this flash a
        // gamma-phi flash rather than `eos.pt_flash`.
        let liquid = ge_nrtl_phase(params, T.value, P.value, &x)?;
        let vapour = mixture.phase_state(&reduced, &y, RootSide::Vapour)?;

        let ln_k_new: Vec<f64> = liquid
            .ln_phi
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
    let liquid = ge_nrtl_phase(params, T.value, P.value, &x)?;
    let vapour = mixture.phase_state(&reduced, &y, RootSide::Vapour)?;

    Ok(GeNrtlFlashResult {
        beta,
        x,
        y,
        k,
        ln_phi_liquid: liquid.ln_phi,
        ln_phi_vapour: vapour.ln_phi,
        z_vapour: vapour.z,
        min_t_over_tc,
        phase,
        iterations,
        residual,
        warnings,
    })
}
