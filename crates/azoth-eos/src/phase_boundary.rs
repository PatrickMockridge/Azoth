//! The pressure iteration `eos.bubble_pressure` and `eos.dew_pressure` share.
//!
//! Both models ask the same question with the phases exchanged: *given one phase's
//! composition, at what pressure does the other phase appear?* Their equations differ,
//! `sum_i x_i K_i = 1` against `sum_i y_i / K_i = 1`, but the loop, the initialisation,
//! the update, the tolerance and the guard against the trivial solution are one piece
//! of code run in two directions - and writing them twice would be writing the guard
//! twice.
//!
//! # The update, and why the two directions are not symmetric
//!
//! The rule in both is one sentence: **if there is too much of the incipient phase,
//! move away from it.**
//!
//! - A bubble point's incipient phase is vapour, `S = sum_i x_i K_i` is how much
//!   vapour the liquid would give off, and `S > 1` means too much - so the pressure
//!   rises, `P <- P * S`.
//! - A dew point's is liquid, `S = sum_i y_i / K_i` is how much liquid the vapour
//!   would condense, and `S < 1` means too *little* - so the pressure rises here
//!   too, and `P / S` is what raises it.
//!
//! Measured with the dew update's sign inverted, before this was written down: the
//! dew point of methane/n-butane at 300 K ran down to `P = 3.6e-07` Pa over nine
//! iterations and then failed inside [`crate::pr_z_factor`], because the pressure had
//! gone to zero rather than to the answer.
//!
//! # The trivial solution
//!
//! `S = 1` has two kinds of solution. One is the boundary this model is looking for.
//! The other is `K_i = 1` for every `i`, which satisfies it at *every* pressure
//! because `sum_i x_i = 1` and `sum_i y_i = 1` identically - and which the iteration
//! converges to whenever the composition has no boundary at this temperature,
//! because that state is a fixed point of the update rule and the physical one is
//! not there to attract it.
//!
//! **The residual cannot detect it, and that is the point worth recording.**
//! `S - 1 = sum_i x_i (K_i - 1)` is a weighted sum whose terms cancel: measured, a
//! state with no bubble point reached `K = (1.000000375, 0.999999724, 0.999999476)`,
//! deviations of order `1e-7`, and `S - 1` came out at `1.5e-13` - inside a tolerance
//! of `1e-12`. A convergence test on `S` alone accepts the trivial solution silently
//! and returns a pressure that looks entirely plausible.
//!
//! So the guard is on the K-values. [`TRIVIAL_TOLERANCE`] carries the measurement
//! behind the constant.
//!
//! # Why this threshold is `1e-2` and the flash's is `1e-8`
//!
//! They are not inconsistent, and the difference is a consequence of the two
//! convergence tests. A flash stops on the change in `ln K` itself, so when it lands
//! on the trivial solution the K-values have stopped moving *because they are* 1 -
//! a tight guard catches it. This iteration stops on `S - 1`, which cancels, so the
//! K-values can still be `1e-6` from 1 when it fires. A `1e-8` guard here would sit
//! below the noise and never trigger.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result, Warning};

use crate::mixture::{Mixture, RootSide};
use crate::model_gen;

/// The `max_i |ln K_i|` below which the two phases have merged.
///
/// Not tight, and deliberately. The measurement is not of where the two kinds of
/// state *end up* - a degenerate one drives `max_i |ln K_i|` smoothly through every
/// value on its way down - but of the lowest value each reaches at any point in its
/// iteration, which is what a threshold has to clear:
///
/// * every state with a genuine boundary kept `max_i |ln K_i| >= 1.30` at **every**
///   step, across the binaries and the ternary in the test files;
/// * every state with no boundary drove it below `1e-06` before its residual test
///   could fire.
///
/// So the gap spans six orders of magnitude and this sits in the middle of it. A
/// tighter threshold would be nearer the degenerate cluster for no gain; a looser one
/// would start refusing genuine boundaries near the critical region, where the
/// K-values legitimately approach 1.
pub const TRIVIAL_TOLERANCE: f64 = 1.0e-02;

/// Wilson's constant, shared with [`crate::pt_flash`]'s initialisation.
///
/// The estimate here is used for a saturation-pressure guess rather than a
/// K-value, but it is the same correlation and the same constant.
const WILSON_CONSTANT: f64 = 5.373;

/// Which phase appears at the boundary being sought.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Incipient {
    /// A bubble point: a liquid is held, and vapour appears.
    Vapour,
    /// A dew point: a vapour is held, and liquid appears.
    Liquid,
}

impl Incipient {
    /// How this is spelled in a message.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Vapour => "bubble",
            Self::Liquid => "dew",
        }
    }
}

/// What a boundary solve found.
#[derive(Debug, Clone, PartialEq)]
pub struct PhaseBoundary {
    /// The boundary pressure.
    pub pressure: f64,
    /// The incipient phase's composition.
    pub incipient: Vec<f64>,
    /// K-values at the converged pressure.
    pub k: Vec<f64>,
    /// The held phase's root.
    pub z_held: f64,
    /// The incipient phase's root.
    pub z_incipient: f64,
    /// Pressure updates taken.
    pub iterations: u32,
    /// `|S - 1|` at the last completed step.
    pub residual: f64,
    /// Warnings raised on the way, from the kernels.
    pub warnings: Vec<Warning>,
}

/// Wilson's saturation-pressure estimate for one component.
///
/// Extrapolated without complaint for a component above its critical temperature,
/// where it is a large number with no physical meaning - which is what a starting
/// guess needs and why it is not reported to the caller.
fn wilson_psat(mixture: &Mixture, t: ThermodynamicTemperature) -> Vec<f64> {
    mixture
        .components()
        .iter()
        .map(|c| {
            c.pc.value * (WILSON_CONSTANT * (1.0 + c.omega) * (1.0 - c.tc.value / t.value)).exp()
        })
        .collect()
}

/// The pressure at which the incipient phase appears, at a fixed temperature.
///
/// `held` is the composition of the phase that is present; the incipient phase's is
/// solved for along with the pressure.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `held` is not one entry per component, does not
///   lie in `[0, 1]`, or does not sum to one; or if the mixture has a single
///   component, whose boundary is a saturation pressure and belongs to
///   [`crate::pure_saturation`].
/// * [`AzothError::OutOfRange`] on `min_t_over_tc` if the iteration converges to the
///   trivial solution, which is how a composition with no boundary at this
///   temperature is reported.
/// * [`AzothError::SolverNotConverged`] if the iteration hits its cap.
pub fn phase_boundary_pressure(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    held: &[f64],
    incipient: Incipient,
) -> Result<PhaseBoundary> {
    let n = mixture.len();
    if n < 2 {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "a {} point needs two phases with different compositions, and one \
                 component cannot have them. A pure component's is its saturation \
                 pressure, which `eos.pure_saturation` computes with \
                 `eos.pr_kappa`, `eos.pr_z_factor` and `eos.pr_departure` behind it.",
                incipient.as_str()
            ),
        ));
    }
    check_composition(held, n, "held")?;

    // The spec the guard's threshold and the update come from. Both schemes read
    // the same settings, so either spec's `algorithm` is the one to follow; the
    // caller has already selected the model whose spec is in force.
    let spec = match incipient {
        Incipient::Vapour => &model_gen::BUBBLE_PRESSURE_SPEC,
        Incipient::Liquid => &model_gen::DEW_PRESSURE_SPEC,
    };
    let algorithm = spec.algorithm;

    let psat = wilson_psat(mixture, t);
    let mut warnings = Vec::new();

    // Raoult's law, which is exact in the ideal limit and a fine starting guess
    // everywhere else.
    let mut pressure = match incipient {
        Incipient::Vapour => (0..n).map(|i| held[i] * psat[i]).sum::<f64>(),
        Incipient::Liquid => 1.0 / (0..n).map(|i| held[i] / psat[i]).sum::<f64>(),
    };
    // Written as an explicit finiteness-and-positivity test rather than
    // `!(pressure > 0.0)`, which reads as a double negative and is the shape clippy
    // flags. Both reject NaN and non-positive values; the explicit form says so.
    if !pressure.is_finite() || pressure <= 0.0 {
        return Err(AzothError::invalid_input(
            "held",
            "the saturation-pressure estimate for this composition is not a positive \
             finite number, so there is no pressure to start the search from",
        ));
    }

    // The incipient phase starts *away* from the held one. Beginning with the two
    // equal would make the first step evaluate both phases at one composition,
    // where the K-values are 1 unless the cubic has two roots to tell them apart -
    // which is the trivial solution wearing a converged face.
    let wilson_k: Vec<f64> = (0..n).map(|i| psat[i] / pressure).collect();
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

        let reduced = mixture.reduced_parameters(t, pascals_from(pressure))?;
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

        // Checked before the update, and on the K-values rather than on the
        // residual: `S - 1` cancels, so it is small here long before the K-values
        // are near 1. See the module documentation.
        if k.iter().all(|value| value.ln().abs() < TRIVIAL_TOLERANCE) {
            return Err(AzothError::OutOfRange {
                field: "min_t_over_tc".to_string(),
                value: mixture
                    .components()
                    .iter()
                    .map(|c| t.value / c.tc.value)
                    .fold(f64::INFINITY, f64::min),
                detail: format!(
                    "the two phases converged onto the same composition at every \
                     pressure (max |ln K| fell below {TRIVIAL_TOLERANCE:e}), so this \
                     mixture has no {} point at this temperature - it is at or above \
                     its critical condition. Unlike the flash, there is no vapour \
                     fraction to report as absent: the pressure itself is what was \
                     being solved for and it is not determined.",
                    incipient.as_str()
                ),
            });
        }

        // The amount of incipient phase the held phase would produce.
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
            return Ok(PhaseBoundary {
                pressure,
                incipient: other,
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

        // Too much incipient phase means move away from it; which direction that is
        // differs between the two, and the derivation is in the module docs.
        pressure = match incipient {
            Incipient::Vapour => pressure * s,
            Incipient::Liquid => pressure / s,
        };
        if !pressure.is_finite() || pressure <= 0.0 {
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

/// A composition must be one entry per component, in `[0, 1]`, summing to one.
fn check_composition(values: &[f64], n: usize, field: &str) -> Result<()> {
    if values.len() != n {
        return Err(AzothError::invalid_input(
            field,
            format!(
                "a composition for {n} components has {} entries",
                values.len()
            ),
        ));
    }
    if let Some(bad) = values.iter().position(|&value| value < 0.0) {
        return Err(AzothError::invalid_input(
            field,
            format!(
                "{field}[{bad}] is {} but a mole fraction cannot be negative",
                values[bad]
            ),
        ));
    }
    let sum: f64 = values.iter().sum();
    if (sum - 1.0).abs() > 1.0e-09 {
        return Err(AzothError::invalid_input(
            field,
            format!(
                "the composition sums to {sum}, not to one. Renormalising it here \
                 would make a caller's error invisible in every number downstream, \
                 so it is refused instead"
            ),
        ));
    }
    Ok(())
}

/// Rescale a vector to sum to one, in place.
fn normalise(values: &mut [f64]) {
    let sum: f64 = values.iter().sum();
    if sum > 0.0 {
        for value in values.iter_mut() {
            *value /= sum;
        }
    }
}

/// A `uom` pressure from a bare pascal magnitude.
///
/// The iteration's interior is dimensionless in everything but the pressure itself,
/// which is what the reduced parameters are reduced against. Keeping it a bare
/// `f64` and converting at the two points where a `Pressure` is needed is the same
/// division of labour the kernels use.
fn pascals_from(value: f64) -> Pressure {
    azoth_core::units::pascals(value)
}
