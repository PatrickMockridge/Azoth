//! The temperature iteration `eos.bubble_temperature` and `eos.dew_temperature` share.
//!
//! Specs: `specs/models/eos/bubble_temperature.toml` and
//! `specs/models/eos/dew_temperature.toml`
//!
//! The mirror image of [`crate::phase_boundary`]: *given one phase's composition, at what
//! temperature does the other phase appear at a fixed pressure?* The same successive
//! substitution on the K-values, with the outer loop moving temperature by a Newton step
//! on `S = sum x_i K_i` (or `sum x_i / K_i` for a dew point) whose slope is the
//! central-difference of the fugacity coefficients with respect to temperature.

use azoth_core::units::{Pressure, kelvins};
use azoth_core::{AzothError, Result, Warning};

use crate::algorithm_of;
use crate::mixture::{Mixture, RootSide, normalise};
use crate::model_gen;
use crate::phase_boundary::{Incipient, TRIVIAL_TOLERANCE, check_composition, wilson_psat};

/// The outcome of finding a saturation temperature.
pub struct PhaseBoundaryTemperature {
    /// The boundary temperature.
    pub temperature: f64,
    /// The incipient phase's composition.
    pub incipient: Vec<f64>,
    /// K-values at the converged temperature.
    pub k: Vec<f64>,
    /// The held phase's root.
    pub z_held: f64,
    /// The incipient phase's root.
    pub z_incipient: f64,
    /// Temperature updates taken.
    pub iterations: u32,
    /// `|S - 1|` at the last completed step.
    pub residual: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

/// `d ln phi_i / dT` for a phase at a composition and root, by central difference.
///
/// The temperature analogue of the K-value substitution's ingredient: the Newton step on
/// `S` needs `dK_i/dT = K_i (d ln phi_i^liquid/dT - d ln phi_i^vapour/dT)`. Taken from the
/// fugacity coefficients themselves rather than an analytic derivative, for the same
/// reason the flash-property solver takes its slopes by central difference.
fn dfugdt(mixture: &Mixture, p: Pressure, x: &[f64], side: RootSide, t: f64) -> Result<Vec<f64>> {
    let delta = (1.0e-4 * t).max(1.0e-6);
    let above =
        mixture.phase_state(&mixture.reduced_parameters(kelvins(t + delta), p)?, x, side)?;
    let below = mixture.phase_state(
        &mixture.reduced_parameters(kelvins((t - delta).max(1.0)), p)?,
        x,
        side,
    )?;
    Ok(above
        .ln_phi
        .iter()
        .zip(&below.ln_phi)
        .map(|(&a, &b)| (a - b) / (2.0 * delta))
        .collect())
}

/// The temperature at which the incipient phase appears, at a fixed pressure.
///
/// # Errors
/// * [`AzothError::OutOfRange`] on `min_t_over_tc` if the iteration converges to the
///   trivial solution.
/// * [`AzothError::SolverNotConverged`] if the iteration hits its cap.
pub fn phase_boundary_temperature(
    mixture: &Mixture,
    p: Pressure,
    held: &[f64],
    incipient: Incipient,
) -> Result<PhaseBoundaryTemperature> {
    let n = mixture.len();
    if n < 2 {
        return Err(AzothError::invalid_input(
            "components",
            "a saturation temperature needs two phases with different compositions, and \
             one component cannot have them. A pure component's is its saturation \
             pressure, which `eos.pure_saturation` computes.",
        ));
    }
    check_composition(held, n, "held")?;

    let spec = match incipient {
        Incipient::Vapour => &model_gen::BUBBLE_TEMPERATURE_SPEC,
        Incipient::Liquid => &model_gen::DEW_TEMPERATURE_SPEC,
    };
    let algorithm = algorithm_of(spec)?;
    let mut warnings = Vec::new();

    // Wilson's saturation-temperature estimate, iterated to `S = 1` the same way the
    // pressure search seeds itself, then a Raoult's-law starting temperature.
    let mut temperature = algorithm.initial_temperature.unwrap_or(300.0).max(50.0);

    // Seed the incipient phase *away* from the held one, as `phase_boundary_pressure`
    // does: the two equal would make the first step trivial.
    let psat = wilson_psat(mixture, kelvins(temperature));
    let wilson_k: Vec<f64> = (0..n).map(|i| psat[i] / p.value).collect();
    let mut other: Vec<f64> = (0..n)
        .map(|i| match incipient {
            Incipient::Vapour => held[i] * wilson_k[i],
            Incipient::Liquid => held[i] / wilson_k[i],
        })
        .collect();
    normalise(&mut other);

    let mut iterations = 0;
    let mut residual = f64::NAN;

    for step in 1..=algorithm.max_iterations {
        iterations = step;

        let reduced = mixture.reduced_parameters(kelvins(temperature), p)?;
        warnings.extend(reduced.warnings.iter().cloned());

        let (liquid, vapour) = match incipient {
            Incipient::Vapour => (held, other.as_slice()),
            Incipient::Liquid => (other.as_slice(), held),
        };
        let liquid_state = mixture.phase_state(&reduced, liquid, RootSide::Liquid)?;
        let vapour_state = mixture.phase_state(&reduced, vapour, RootSide::Vapour)?;

        let k: Vec<f64> = liquid_state
            .ln_phi
            .iter()
            .zip(&vapour_state.ln_phi)
            .map(|(&l, &v)| (l - v).exp())
            .collect();

        if k.iter().all(|value| value.ln().abs() < TRIVIAL_TOLERANCE) {
            return Err(AzothError::OutOfRange {
                field: "min_t_over_tc".to_string(),
                value: mixture
                    .components()
                    .iter()
                    .map(|c| temperature / c.tc.value)
                    .fold(f64::INFINITY, f64::min),
                detail: format!(
                    "the two phases converged onto the same composition at every \
                     temperature (max |ln K| fell below {TRIVIAL_TOLERANCE:e}), so this \
                     mixture has no {} point at this pressure.",
                    incipient.as_str()
                ),
            });
        }

        let s: f64 = match incipient {
            Incipient::Vapour => (0..n).map(|i| held[i] * k[i]).sum(),
            Incipient::Liquid => (0..n).map(|i| held[i] / k[i]).sum(),
        };
        for i in 0..n {
            other[i] = match incipient {
                Incipient::Vapour => held[i] * k[i],
                Incipient::Liquid => held[i] / k[i],
            } / s;
        }
        residual = (s - 1.0).abs();

        if residual <= algorithm.tolerance {
            return Ok(PhaseBoundaryTemperature {
                temperature,
                incipient: other.clone(),
                k,
                z_held: match incipient {
                    Incipient::Vapour => liquid_state.z,
                    Incipient::Liquid => vapour_state.z,
                },
                z_incipient: match incipient {
                    Incipient::Vapour => vapour_state.z,
                    Incipient::Liquid => liquid_state.z,
                },
                iterations,
                residual,
                warnings,
            });
        }

        // The Newton step on temperature, with the derivative of `S` assembled from the
        // per-component fugacity temperature derivatives. Re-derived from the *updated*
        // `other`, whose previous borrow ended with the update above.
        let (liquid, vapour) = match incipient {
            Incipient::Vapour => (held, other.as_slice()),
            Incipient::Liquid => (other.as_slice(), held),
        };
        let d_liquid = dfugdt(mixture, p, liquid, RootSide::Liquid, temperature)?;
        let d_vapour = dfugdt(mixture, p, vapour, RootSide::Vapour, temperature)?;
        let mut dsdt = 0.0;
        for i in 0..n {
            match incipient {
                Incipient::Vapour => dsdt += held[i] * k[i] * (d_liquid[i] - d_vapour[i]),
                Incipient::Liquid => dsdt -= held[i] / k[i] * (d_liquid[i] - d_vapour[i]),
            }
        }
        if !dsdt.is_finite() || dsdt == 0.0 {
            // Fall back to a fixed step when the slope cannot be evaluated.
            temperature += if s > 1.0 { 5.0 } else { -5.0 };
        } else {
            temperature -= (s - 1.0) / dsdt;
            if !temperature.is_finite() || temperature <= 0.0 {
                return Err(AzothError::SolverNotConverged {
                    iterations,
                    residual,
                    tolerance: algorithm.tolerance,
                });
            }
        }
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(AzothError::SolverNotConverged {
                iterations,
                residual,
                tolerance: algorithm.tolerance,
            });
        }
    }

    Err(AzothError::SolverNotConverged {
        iterations,
        residual,
        tolerance: algorithm.tolerance,
    })
}
