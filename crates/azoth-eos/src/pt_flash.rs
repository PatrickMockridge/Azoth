//! `eos.pt_flash` - the isothermal two-phase flash.
//!
//! Spec: `specs/models/eos/pt_flash.toml`
//!
//! The two-phase split of a mixture at a fixed temperature and pressure: Wilson
//! K-value estimates, then successive substitution, with Rachford-Rice bisected each
//! iteration. A *model* rather than a calculation - what the spec pins down is the
//! procedure, and this module reads the procedure from the generated table rather than
//! choosing it.
//!
//! The scaffold - the Rachford-Rice bracket and bisection, the convergence measure,
//! and what a trivial solution is - lives in [`crate::flash_iteration`], which
//! [`crate::ge_nrtl_flash`] shares. What is here is the part that is Peng-Robinson's:
//! both fugacity coefficients come from the cubic.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result, apply_checks};

use crate::algorithm_of;
use crate::flash_iteration::{
    Outcome, compositions, is_trivial, negative_flash_warning, rachford_rice, rachford_rice_bounds,
    rms_delta, single_phase_warning, trivial_warning,
};
use crate::flash_newton;
use crate::mixture::{Mixture, RootSide, wilson_k};
use crate::model_gen;

use crate::results::{Phase, PtFlashResult};

/// The tolerance the loop settled against, which depends on which scheme finished it.
///
/// The two schemes measure different things - the outer one the change in `ln K` and
/// the fallback the step it took - so they carry their own stopping rules, and the
/// comparison at the end has to be against the one that actually ran.
fn settled_tolerance(algorithm: &azoth_core::ModelAlgorithm, from_newton: bool) -> f64 {
    match algorithm.fallback {
        Some(fallback) if from_newton => fallback.algorithm.tolerance,
        _ => algorithm.tolerance,
    }
}

/// The isothermal two-phase flash of a mixture at a temperature and pressure.
///
/// `z` is the overall composition, in mole fractions, and is **checked rather than
/// renormalised** - silently rescaling a caller's composition would make their error
/// invisible in a way that changes every number downstream.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
/// * [`AzothError::InvalidInput`] if `z` is the wrong length, has a negative entry,
///   or does not sum to one.
/// * [`AzothError::SolverNotConverged`] if the iteration hits its cap.
/// * Propagates the kernels' range checks.
///
/// # Example
/// ```
/// use azoth_eos::Cubic;
/// use azoth_core::units::{kelvins, pascals};
/// use azoth_eos::{databank, pt_flash};
///
/// let mixture = databank::mixture_of(&["methane", "n-butane"], Cubic::Pr, None)
///     .expect("the pair resolves")
///     .0;
/// let r = pt_flash(&mixture, kelvins(330.0), pascals(2_500_000.0), &[0.6, 0.4])?;
/// assert_eq!(r.phase, azoth_eos::Phase::TwoPhase);
/// assert!((r.beta.expect("a split has a vapour fraction") - 0.8422055475803881).abs() < 1e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn pt_flash(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<PtFlashResult> {
    let spec = &model_gen::PT_FLASH_SPEC;
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
        .map(|c| t.value / c.tc.value)
        .fold(f64::INFINITY, f64::min);
    apply_checks(
        spec.derived_checks(),
        |quantity| (quantity == "min_t_over_tc").then_some(min_t_over_tc),
        &mut warnings,
    )?;

    let reduced = mixture.reduced_parameters(t, p)?;
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

    let mut k = wilson_k(mixture, t, p);
    let mut iterations = 0;
    let mut residual = f64::NAN;
    // The vapour mole numbers the fallback Newton works in, seeded from the last
    // successive-substitution iterate so the handover starts where the outer scheme
    // had got to rather than from the feed.
    let mut u = vec![0.0; n];
    let mut newton_steps = 0u32;
    let mut from_newton = false;
    // Whether the last outer iterate was a two-phase split. The fallback works in the
    // vapour mole numbers `u = beta y`, which describe two phases only while `beta` is
    // inside `(0, 1)` - a negative flash has no such `u`, and every trial the line
    // search makes at one is infeasible and halved through.
    let mut splits = false;

    // How the loop finished. `None` means it is still iterating, so a value that
    // is still `None` after the loop hit its cap is a genuine failure to converge.
    let mut settled: Option<Outcome> = None;

    for step in 1..=algorithm.max_iterations {
        iterations = step;

        // The handover. NeqSim's `TPflash` runs its own scheme until `activeNewtonLimit`
        // steps have passed and replaces it with the second-order one from then on. Two
        // things end the handover without ending the solve: a step the line search
        // refuses, and the second-order scheme's own attempt budget. Both hand the
        // iteration back to the outer scheme, which is what NeqSim's `catch` does -
        // the cap that ends a solve is the outer one, never the fallback's.
        if let Some(fallback) = algorithm
            .fallback
            .filter(|f| splits && step >= f.after && newton_steps < f.algorithm.max_iterations)
        {
            newton_steps += 1;
            if let Ok(stepped) = flash_newton::step(mixture, &reduced, z, &u, fallback.algorithm) {
                u = stepped.u;
                residual = stepped.residual;
                from_newton = true;
                // The K-values follow the second-order iterate, so that a handover back
                // at `max_iterations` resumes from the state the solve is actually in
                // rather than re-running the outer scheme from where it was when it
                // handed over.
                let (x, y) = flash_newton::split_at(z, &u);
                k = y.iter().zip(&x).map(|(&yi, &xi)| yi / xi).collect();
                if residual <= fallback.algorithm.tolerance {
                    break;
                }
                continue;
            }
            from_newton = false;
        }

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
        let liquid = mixture.phase_state(&reduced, &x, RootSide::Liquid)?;
        let vapour = mixture.phase_state(&reduced, &y, RootSide::Vapour)?;

        let ln_k_new: Vec<f64> = liquid
            .ln_phi
            .iter()
            .zip(&vapour.ln_phi)
            .map(|(&l, &v)| l - v)
            .collect();
        residual = rms_delta(&ln_k_new, &k);
        k = ln_k_new.iter().map(|value| value.exp()).collect();
        u = flash_newton::vapour_moles(beta, &y);
        splits = beta > 0.0 && beta < 1.0;

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
        // and is not reproducible; this is. See the spec's correction 2.
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
        None if residual > settled_tolerance(algorithm, from_newton) => {
            return Err(AzothError::SolverNotConverged {
                iterations,
                residual,
                tolerance: settled_tolerance(algorithm, from_newton),
            });
        }
        // The second-order scheme converged, and it converges on its own measure -
        // the relative step - rather than on the outer scheme's. Its answer is a
        // state in `u`, so the phases and the K-values are read back off it.
        None if from_newton => {
            let beta = u.iter().sum::<f64>();
            let (x, y) = flash_newton::split_at(z, &u);
            k = y.iter().zip(&x).map(|(&yi, &xi)| yi / xi).collect();
            // The trivial solution can wear a small vapour fraction: `u` near `beta z`
            // puts the two phases on top of each other with `beta` just above zero, and
            // the residual is satisfied there as it is at every other trivial point.
            // Read back as a split it would report a two-phase feed whose two phases
            // are the feed - measured at methane/n-butane 355 K, 120 bar, where the
            // second-order scheme lands on `beta = 2.8e-4` with `K` within 1e-8 of one.
            if is_trivial(&k) {
                k = vec![1.0; n];
                warnings.push(trivial_warning());
                (None, z.to_vec(), z.to_vec(), Phase::Trivial)
            } else {
                let phase = flash_newton::phase_of(beta);
                if phase != Phase::TwoPhase {
                    warnings.push(negative_flash_warning(phase, beta));
                }
                (Some(beta), x, y, phase)
            }
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

    let liquid = mixture.phase_state(&reduced, &x, RootSide::Liquid)?;
    let vapour = mixture.phase_state(&reduced, &y, RootSide::Vapour)?;

    Ok(PtFlashResult {
        beta,
        x,
        y,
        k,
        ln_phi_liquid: liquid.ln_phi,
        ln_phi_vapour: vapour.ln_phi,
        z_liquid: liquid.z,
        z_vapour: vapour.z,
        min_t_over_tc,
        phase,
        iterations,
        residual,
        warnings,
    })
}
