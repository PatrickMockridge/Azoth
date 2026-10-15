//! `eos.iapws_saturation_pressure` - the IAPWS-IF97 Region 4 saturation pressure of water.
//!
//! ```text
//! p_sat = (2C / (-B + sqrt(B**2 - 4*A*C)))**4,  theta = T + n9/(T - n10)
//! ```
//!
//! Spec: `specs/calcs/eos/iapws_saturation_pressure.toml`. Oracle:
//! `validation/neqsim/captures/iapws_if97_probe.tsv`.
//!
//! **It is the boundary `eos.water_phase` selects its region by**, and the release states the
//! two directions as two fits - `p4_T` forward and `T4_p` backward are not one equation
//! inverted - so the inverse is a second id's question and not this one's.
//!
//! **A different answer from `eos.pure_saturation` to a question that sounds the same.** That
//! one is the pressure where a Peng-Robinson fluid's two roots cross, which is a property of
//! the cubic; this is the reference curve for water.

use azoth_core::units::{ThermodynamicTemperature, pascals};
use azoth_core::{Result, apply_checks};

use crate::iapws_if97;
use crate::results::IapwsSaturationPressureResult;
use crate::spec_gen;

/// The saturation pressure of water at a temperature, by the Region 4 equation.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `T` is outside `[273.15, 647.096)` K - the
///   triple point to the critical point, which is the curve's own domain. A value outside it
///   is a different question rather than a bad number.
pub fn iapws_saturation_pressure(
    t: ThermodynamicTemperature,
) -> Result<IapwsSaturationPressureResult> {
    let spec = &spec_gen::IAPWS_SATURATION_PRESSURE_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // The module returns the release's own unit, megapascals; the model is SI.
    let p_sat = iapws_if97::saturation_pressure(t.value) * 1.0e6;

    Ok(IapwsSaturationPressureResult {
        p_sat: pascals(p_sat),
        warnings,
    })
}
