//! `eos.ideal_gas_cp` - ideal-gas heat capacity from a four-term polynomial.
//!
//! ```text
//! Cp/R = a + b*theta + c*theta**2 + d*theta**3,   theta = T / (1000 K)
//! ```
//!
//! Spec: `specs/calcs/eos/ideal_gas_cp.yaml`
//!
//! # The coefficients are the caller's, and that is the point
//!
//! This is the first calc in this namespace whose constants are *fitted data* rather
//! than published coefficients of an equation. A table of them is exactly the
//! component databank this library has decided not to have, so the four arrive as
//! arguments and none ship. The same decision keeps `Cd`, `Cv`, `f_T`, `k_ij` and
//! `kappa_1` out.
//!
//! # Why `T/(1000 K)` rather than `T`
//!
//! Because that is the substitution which makes the polynomial dimensionless. A
//! published table gives `Cp/R = A + B T + C T**2 + D T**3` with `T` in kelvin and
//! prints `B`, `C`, `D` already multiplied by `10**3`, `10**6` and `10**9` - exactly
//! so that all four numbers come out dimensionless - and against `theta = T/(1000 K)`
//! those printed numbers are the coefficients. A caller copies them across unchanged,
//! and this calc needs no unit for a polynomial coefficient at all.

use azoth_core::units::{ThermodynamicTemperature, joules_per_mole_kelvin};
use azoth_core::{Result, apply_checks};

use crate::pr_molar_volume::MOLAR_GAS_CONSTANT;
use crate::results::IdealGasCpResult;
use crate::spec_gen;

/// The reference temperature the polynomial is written against, in kelvin.
///
/// A stated constant of the correlation's form rather than a fitted quantity: it is
/// what makes a table's four printed numbers dimensionless, and a caller whose
/// coefficients are quoted against another reference rescales them once.
pub const REFERENCE_TEMPERATURE: f64 = 1000.0;

/// The ideal-gas heat capacity at a temperature, from a caller-supplied polynomial.
///
/// The coefficients are dimensionless and describe `Cp/R`; the result is dimensioned.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `T` is not positive.
///
/// A non-positive `cp` is **returned** carrying `OutOfValidRange` rather than
/// refused: it means the polynomial has been evaluated outside the range it was
/// fitted over, the arithmetic is well defined, and inspecting the limit is a
/// legitimate thing for a caller to be doing. See the spec.
///
/// # Example
/// ```
/// use azoth_core::units::kelvins;
/// use azoth_eos::ideal_gas_cp;
///
/// let r = ideal_gas_cp(4.0, 1.0, -0.5, 0.1, kelvins(500.0))?;
/// assert!((r.cp_over_r - 4.3875).abs() < 1e-15);
/// assert!((r.cp.value - 36.47970473714734).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn ideal_gas_cp(
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    t: ThermodynamicTemperature,
) -> Result<IdealGasCpResult> {
    let spec = &spec_gen::IDEAL_GAS_CP_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "a" => Some(a),
            "b" => Some(b),
            "c" => Some(c),
            "d" => Some(d),
            "T" => Some(t.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let theta = t.value / REFERENCE_TEMPERATURE;
    let cp_over_r = a + b * theta + c * theta * theta + d * theta * theta * theta;
    let cp = cp_over_r * MOLAR_GAS_CONSTANT;

    apply_checks(
        spec.derived_checks(),
        |quantity| (quantity == "cp").then_some(cp),
        &mut warnings,
    )?;

    Ok(IdealGasCpResult {
        cp_over_r,
        cp: joules_per_mole_kelvin(cp),
        warnings,
    })
}
