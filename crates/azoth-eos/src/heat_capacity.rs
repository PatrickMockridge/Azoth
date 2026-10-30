//! One phase's heat capacities at constant pressure and volume.
//!
//! `cp - cv = -T (dP/dT)^2 / (dP/dV)`, the cubic's own relation, with `cp` the molar
//! departure capacity the enthalpy/entropy model already carries. `eos.sound_speed` needs
//! `gamma = cp/cv` for its speed; `process.leak_blowdown` needs the same capacities,
//! weighted over the phases, for the class's `getGamma()`. This is that one arithmetic.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result};

use crate::mixture::{Mixture, ReducedParameters};
use crate::molar_enthalpy_entropy::{IdealGasModel, molar_enthalpy_entropy};
use crate::pr_molar_volume::MOLAR_GAS_CONSTANT;

/// One phase's capacities and the pressure derivative they rest on.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhaseCapacity {
    /// Molar isobaric heat capacity, J/(mol K).
    pub cp: f64,
    /// Molar isochoric heat capacity, J/(mol K).
    pub cv: f64,
    /// The cubic's `(dP/dV)_T`, Pa m3/mol; negative for a stable phase.
    pub d_p_d_v: f64,
}

/// One phase's `cp` and `cv` from the cubic's closed-form derivatives.
///
/// # Errors
/// `InvalidInput` when the state gives no admissible derivative or the capacities are not
/// finite and positive - a sound speed and a blowdown both need a stable phase.
pub fn phase_capacity(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    reduced: &ReducedParameters,
    t: ThermodynamicTemperature,
    p: Pressure,
    composition: &[f64],
    root: f64,
) -> Result<PhaseCapacity> {
    let state = mixture.phase_state_at(reduced, composition, root)?;
    let constant = MOLAR_GAS_CONSTANT;
    let v = root * constant * t.value / p.value;
    let (delta1, delta2) = (mixture.cubic().delta1(), mixture.cubic().delta2());

    // `A` is dimensionless (`aP/(RT)^2`), so the dimensional `a` and `b` are these.
    let a = state.a_mix * (constant * t.value) * (constant * t.value) / p.value;
    let b = state.b_mix * constant * t.value / p.value;
    let d_a_d_t = a * state.psi_bar / t.value;
    let f = (v + delta1 * b) * (v + delta2 * b);
    let d_f = 2.0 * v + (delta1 + delta2) * b;
    let d_p_d_v = -constant * t.value / ((v - b) * (v - b)) + a * d_f / (f * f);
    let d_p_d_t = constant / (v - b) - d_a_d_t / f;
    if d_p_d_v.is_nan() || d_p_d_v >= 0.0 || !d_p_d_v.is_finite() || !d_p_d_t.is_finite() {
        return Err(AzothError::invalid_input(
            "state",
            format!(
                "the cubic's pressure-volume derivative is {d_p_d_v} Pa m3/mol at {:.3} K and \
                 {:.3} Pa, and a sound speed needs a stable phase",
                t.value, p.value
            ),
        ));
    }

    let cp = molar_enthalpy_entropy(mixture, ideal_gas, t, p, composition, root)?
        .cp
        .value;
    let cv = cp + t.value * d_p_d_t * d_p_d_t / d_p_d_v;
    if cv.is_nan() || cv <= 0.0 || !cp.is_finite() {
        return Err(AzothError::invalid_input(
            "state",
            format!("the phase's heat capacities are cp={cp}, cv={cv}, and both must be positive"),
        ));
    }
    Ok(PhaseCapacity { cp, cv, d_p_d_v })
}
