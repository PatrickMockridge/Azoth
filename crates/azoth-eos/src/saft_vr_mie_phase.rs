//! SAFT-VR-Mie's volume at a temperature, a pressure and a composition.
//!
//! ```text
//! solve P_calc(v) = P   for   v,   Z = P v/(R T)
//! ```
//!
//! The layers are [`crate::saft_vr_mie`]'s and the solve is [`crate::volume_solve`]'s, which
//! is what turns them into a state. This module is the two together and the floor.

use azoth_core::{AzothError, Result};

use crate::association::R;
use crate::mixture::RootSide;
use crate::pcsaft::AVOGADRO;
use crate::saft_vr_mie::{self, MieComponent};
use crate::volume_solve;

/// A solved SAFT-VR-Mie volume.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MieVolume {
    /// The molar volume, in m^3/mol.
    pub v: f64,
    /// The compressibility factor at that volume.
    pub z: f64,
    /// The packing fraction at that volume.
    pub eta: f64,
    /// The Newton steps the vapour branch took, or zero on the liquid branch, which is
    /// bracketed and bisected. Reported so a caller can see a solve that did not converge.
    pub iterations: usize,
}

/// The step the pressure's slope is differenced at, as a fraction of the volume.
///
/// **Ten times NeqSim's own `eta` step, and it is a choice between two noises.** The
/// pressure here is already a central difference - NeqSim differentiates the chain and
/// dispersion terms that way - so this differences a difference, and at `1e-5` the
/// cancellation owns the slope.
///
/// Measured on the binary state, over two decades: the converged volume lands `5.9e-6` from
/// NeqSim's at `1e-3`, `4.0e-7` at `3e-4`, `8.2e-7` at `1e-4`, `7.2e-7` at `3e-5` and
/// `3.1e-6` at `1e-5`. So there is a plateau and `1e-4` is on it; the model's answer is not
/// in question and the last digit of any root here is the arithmetic.
const SLOPE_STEP: f64 = 1.0e-4;

/// Solve `P_calc(v) = p` for the molar volume, in m^3/mol.
///
/// The packing fraction reaches one at `v = (pi/6) N_A md3`, where every hard-sphere term
/// diverges; that volume is the floor of the solve, and the Newton does walk down to it -
/// the vapour branch steps towards whatever root the isotherm has, and on a state with a
/// liquid root that means towards the floor - so the slope is taken one-sided there rather
/// than across a divergence.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `t` or `p` is not positive, or the wanted branch has no
///   root at this state.
/// * [`AzothError::InvalidInput`] as [`saft_vr_mie::state`].
pub fn molar_volume(
    components: &[MieComponent],
    x: &[f64],
    t: f64,
    p: f64,
    side: RootSide,
) -> Result<MieVolume> {
    if !p.is_finite() || p <= 0.0 {
        return Err(AzothError::OutOfRange {
            field: "P".to_string(),
            value: p,
            detail: "the volume solve divides by the pressure".to_string(),
        });
    }
    let rt = R * t;
    let ideal = rt / p;
    let md3 = saft_vr_mie::state(components, x, t, ideal)?.md3;
    let floor = std::f64::consts::PI / 6.0 * AVOGADRO * md3;

    // `P_calc(v) - p` in Pa, and its slope in v from a central difference of the pressure.
    let residual = |v: f64| -> Result<(f64, f64)> {
        let pressure = saft_vr_mie::pressure_over_rt(components, x, t, v)? * rt;
        let step = SLOPE_STEP * v;
        let (high, low) = (v + step, v - step);
        // **One-sided near the floor**, because a difference that reaches under it is a
        // difference across a divergence rather than a slope. The Newton does walk down
        // there - the vapour branch steps towards whatever root the isotherm has, and on a
        // state with a liquid root that means towards the floor - so the choice is between
        // a one-sided slope and refusing a solve that has a perfectly good answer above it.
        let slope = if low <= floor {
            let upper = saft_vr_mie::pressure_over_rt(components, x, t, high)? * rt;
            (upper - pressure) / step
        } else {
            let upper = saft_vr_mie::pressure_over_rt(components, x, t, high)? * rt;
            let lower = saft_vr_mie::pressure_over_rt(components, x, t, low)? * rt;
            (upper - lower) / (2.0 * step)
        };
        Ok((pressure - p, slope))
    };

    let (v, iterations) = volume_solve::molar_volume(residual, ideal, floor, side, |branch| {
        AzothError::OutOfRange {
            field: "P".to_string(),
            value: p,
            detail: format!(
                "the {branch} branch has no zero at {t} K and {p} Pa. The isotherm at this \
                 temperature does not have a root on that side, so there is no such phase \
                 rather than a volume to report"
            ),
        }
    })?;

    let eta = std::f64::consts::PI / 6.0 * AVOGADRO * md3 / v;
    Ok(MieVolume {
        v,
        z: p * v / rt,
        eta,
        iterations,
    })
}
