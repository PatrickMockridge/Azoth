//! `eos.ideal_gas_cp` - ideal-gas heat capacity from a polynomial.
//!
//! ```text
//! cp = cp_a + cp_b*T + cp_c*T**2 + cp_d*T**3 + cp_e*T**4
//! ```
//!
//! A port of `neqsim.thermo.component.Component.getCp0(double)`, and dimensional
//! throughout: the coefficients carry the powers of temperature in their units, so
//! `cp_a` is a heat capacity and `cp_b` is one per kelvin. That is what lets the five
//! `CPA`-`CPE` columns of NeqSim's `COMP.csv` be read as they are stored.
//!
//! Spec: `specs/calcs/eos/ideal_gas_cp.yaml`.

use azoth_core::units::{ThermodynamicTemperature, joules_per_mole_kelvin};
use azoth_core::{Result, apply_checks};

use crate::results::IdealGasCpResult;
use crate::spec_gen;

/// The ideal-gas heat capacity at a temperature, from a caller-supplied polynomial.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `T` is not positive.
///
/// A non-positive `cp` is **returned** carrying `OutOfValidRange` rather than refused:
/// it means the polynomial has been evaluated outside the range it was fitted over, the
/// arithmetic is well defined, and inspecting the limit is a legitimate thing for a
/// caller to be doing. See the spec.
///
/// # Example
/// ```
/// use azoth_core::units::kelvins;
/// use azoth_eos::ideal_gas_cp;
///
/// // Methane's coefficients as NeqSim ships them, at 300 K.
/// let r = ideal_gas_cp(
///     37.978352,
///     -0.07461815,
///     0.000301881,
///     -2.83e-07,
///     9.070574e-11,
///     kelvins(300.0),
/// )?;
/// assert!((r.cp.value - 35.855913494).abs() < 1e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn ideal_gas_cp(
    cp_a: f64,
    cp_b: f64,
    cp_c: f64,
    cp_d: f64,
    cp_e: f64,
    t: ThermodynamicTemperature,
) -> Result<IdealGasCpResult> {
    let spec = &spec_gen::IDEAL_GAS_CP_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "cp_a" => Some(cp_a),
            "cp_b" => Some(cp_b),
            "cp_c" => Some(cp_c),
            "cp_d" => Some(cp_d),
            "cp_e" => Some(cp_e),
            "T" => Some(t.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // Abbreviated multiplication rather than `powi`, because these are the terms as the
    // source writes them and a reader checking one against the other should not have to
    // expand a call.
    let temperature = t.value;
    let cp = cp_a
        + cp_b * temperature
        + cp_c * temperature * temperature
        + cp_d * temperature * temperature * temperature
        + cp_e * temperature * temperature * temperature * temperature;

    apply_checks(
        spec.derived_checks(),
        |quantity| (quantity == "cp").then_some(cp),
        &mut warnings,
    )?;

    Ok(IdealGasCpResult {
        cp: joules_per_mole_kelvin(cp),
        warnings,
    })
}
