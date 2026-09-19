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

/// The fugacity coefficient of a **pure** component, as `ln phi`.
///
/// NeqSim's `ComponentSAFTVRMie.dFdN` branches on the component count, and the pure branch
/// is an identity rather than a derivative: `dF/dN = F/n - v dF/dV`. With `n = 1` that is
/// `f + Z - 1`, because `dF/dV` is `-(Z-1)/v` - the pressure's own relation - so
///
/// ```text
/// ln phi = f + Z - 1 - ln Z
/// ```
///
/// **The mixture branch is not this.** NeqSim sums three analytic per-component derivatives
/// there (`dF_HC_SAFTdN + dF_DISP_SAFTdN + dFCPAdN`), and the dispersion's is a pair sum
/// whose cross parameters depend on both compositions. That is not ported, and a mixture is
/// refused rather than given the pure formula - which is right only at `m_bar = 1` and
/// silently wrong otherwise.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if more than one component is given.
/// * [`AzothError::OutOfRange`] if the state's compressibility is not positive, where the
///   logarithm is not defined.
pub fn ln_phi_pure(components: &[MieComponent], x: &[f64], t: f64, v: f64) -> Result<f64> {
    if components.len() != 1 {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "this is the pure component's fugacity coefficient and {} components were \
                 given. The mixture's is NeqSim's `dF_DISP_SAFTdN`, a pair sum whose cross \
                 parameters depend on both compositions, and it is not derived",
                components.len()
            ),
        ));
    }
    let s = saft_vr_mie::state(components, x, t, v)?;
    let z = saft_vr_mie::pressure_over_rt(components, x, t, v)? * v;
    if !z.is_finite() || z <= 0.0 {
        return Err(AzothError::OutOfRange {
            field: "Z".to_string(),
            value: z,
            detail: "the fugacity coefficient is `F/n - v dF/dV - ln Z`, and a \
                     compressibility at or below zero is not a state"
                .to_string(),
        });
    }
    Ok(s.f() + z - 1.0 - z.ln())
}

/// Every component's fugacity coefficient at a state, as `ln phi`.
///
/// `ln phi_i = d(nF)/dn_i - ln Z` with the derivative of the extensive `A^R/(RT)` at fixed
/// temperature and volume, and **NeqSim branches on the component count**:
///
/// - **One component** is the identity `dF/dN = F/n - v dF/dV`, which is `f + Z - 1`, so
///   the result is `f + Z - 1 - ln Z`.
/// - **A mixture** sums two analytic per-component derivatives. The hard-sphere and chain
///   part is
///   `m_i a_hs + m_bar a_hs' eta_i - (m_i - 1) ln g_hs - m_min1 (g_hs'/g_hs) eta_i`,
///   which is the *same* expression PC-SAFT's composition derivative reduces to - the only
///   difference between the two models here is that the `g_hs` is the Mie-weighted one. The
///   dispersion is `m_i (2 row_i - a_disp) + m_bar (da_disp/d eta) eta_i`, where `row_i` is
///   the sum of the pair terms against `i` and the first half is the quadratic form's
///   derivative in the segment fractions.
///
/// `eta_i = (pi/6) N_A m_i d_i^3/v` is the packing fraction's own derivative, and the `eta`
/// derivatives of `g_hs` and of the dispersion are the central differences NeqSim takes.
///
/// No association term: SAFT-VR-Mie's association is not in this port's scope, and a
/// mixture whose components have sites would need `dFCPAdN` beside these.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if the compressibility is not positive, or the state is not
///   one.
/// * [`AzothError::InvalidInput`] as [`saft_vr_mie::state`].
pub fn ln_phi(components: &[MieComponent], x: &[f64], t: f64, v: f64) -> Result<Vec<f64>> {
    let state = saft_vr_mie::state(components, x, t, v)?;
    let z = saft_vr_mie::pressure_over_rt(components, x, t, v)? * v;
    if !z.is_finite() || z <= 0.0 {
        return Err(AzothError::OutOfRange {
            field: "Z".to_string(),
            value: z,
            detail: "the fugacity coefficients are `d(nF)/dn_i - ln Z`, and a \
                     compressibility at or below zero is not a state"
                .to_string(),
        });
    }
    if components.len() == 1 {
        return Ok(vec![state.f() + z - 1.0 - z.ln()]);
    }

    let eta = state.eta;
    // NeqSim's step, in the packing fraction.
    let step = (eta.abs() * 1.0e-5).max(1.0e-12);
    let high = eta + step;
    let low = (eta - step).max(1.0e-15);
    let step = (high - low) / 2.0;

    let one = 1.0 - eta;
    let a_hs_eta = (4.0 - 2.0 * eta) / one.powi(3);
    let g_eta = (saft_vr_mie::chain_contact_value(components, x, t, high, &state.d)
        - saft_vr_mie::chain_contact_value(components, x, t, low, &state.d))
        / (2.0 * step);
    let dispersion = |e: f64| -> Result<f64> {
        let (a1, a2, a3) = saft_vr_mie::dispersion_at(components, x, t, e, &state.d)?;
        Ok(a1 + a2 + a3)
    };
    let dispersion_eta = (dispersion(high)? - dispersion(low)?) / (2.0 * step);

    let rows = saft_vr_mie::dispersion_row_sums(components, x, t, eta)?;
    let a_disp = state.a1 + state.a2 + state.a3;

    let mut out = Vec::with_capacity(components.len());
    for (i, component) in components.iter().enumerate() {
        let eta_i = std::f64::consts::PI / 6.0 * AVOGADRO * component.m * state.d[i].powi(3) / v;
        let hard_chain = component.m * state.a_hs + state.m_bar * a_hs_eta * eta_i
            - (component.m - 1.0) * state.g_hs.ln()
            - state.m_minus_1 * g_eta / state.g_hs * eta_i;
        let dispersion_i =
            component.m * (2.0 * rows[i] - a_disp) + state.m_bar * dispersion_eta * eta_i;
        out.push(hard_chain + dispersion_i - z.ln());
    }
    Ok(out)
}
