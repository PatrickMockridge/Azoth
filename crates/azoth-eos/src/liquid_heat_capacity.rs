//! `eos.liquid_heat_capacity` - NeqSim's pure-component liquid heat-capacity
//! polynomial.
//!
//! ```text
//! cp = 1e-3*(c0 + c1*T + c2*T**2 + c3*T**3 + c4*T**4)
//! ```
//!
//! Spec: `specs/calcs/eos/liquid_heat_capacity.toml`, which records the internal-unit
//! coefficients and the `1e-3` that recovers `J/(mol*K)`.

use azoth_core::units::{ThermodynamicTemperature, joules_per_mole_kelvin};
use azoth_core::{Result, apply_checks};

use crate::results::LiquidHeatCapacityResult;
use crate::spec_gen;

/// The pure-component liquid heat capacity at a temperature, from NeqSim's polynomial.
///
/// `c0..c4` are the raw `cpliquid1..5`, caller-supplied.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `T` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::kelvins;
/// use azoth_eos::liquid_heat_capacity;
///
/// let r = liquid_heat_capacity(276370.0, -2090.1, 8.125, -0.014116, 9.37e-6, kelvins(300.0))?;
/// assert!((r.cp.value - 75.355).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T` is the symbol in the published equation
pub fn liquid_heat_capacity(
    c0: f64,
    c1: f64,
    c2: f64,
    c3: f64,
    c4: f64,
    T: ThermodynamicTemperature,
) -> Result<LiquidHeatCapacityResult> {
    let spec = &spec_gen::LIQUID_HEAT_CAPACITY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "c0" => Some(c0),
            "c1" => Some(c1),
            "c2" => Some(c2),
            "c3" => Some(c3),
            "c4" => Some(c4),
            "T" => Some(T.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let t = T.value;
    let cp = 1e-3 * (c0 + c1 * t + c2 * t * t + c3 * t * t * t + c4 * t * t * t * t);

    apply_checks(
        spec.derived_checks(),
        |name| (name == "cp").then_some(cp),
        &mut warnings,
    )?;

    Ok(LiquidHeatCapacityResult {
        cp: joules_per_mole_kelvin(cp),
        warnings,
    })
}
