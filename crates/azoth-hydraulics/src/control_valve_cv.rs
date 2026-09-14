//! `hydraulics.control_valve_cv` - liquid flow through a control valve.
//!
//! ```text
//! q = CV_TO_SI * Cv * sqrt(dP / SG)
//! ```
//!
//! Spec: `specs/calcs/hydraulics/control_valve_cv.toml`
//!
//! # The coefficient is an input, and its units are the hard part
//!
//! `Cv` is supplied by the caller: IEC 60534's published coefficient values are a
//! table, and this library does not reproduce tables - the same position as `Cd` in
//! [`crate::orifice_flow`].
//!
//! What makes this calc different is that `Cv` is *not* dimensionless. It is
//! defined as the flow of water in US gallons per minute at a pressure drop of one
//! pound per square inch, so it carries `gpm / sqrt(psi)`. The spec has to declare
//! it `dimensionless` because the unit vocabulary cannot express a fractional power
//! of a non-SI unit, so undoing that declaration correctly is this calculation's
//! whole job - which is why [`CV_TO_SI`] is named, derived and tested rather than
//! folded into an input description.
//!
//! # Which convention
//!
//! This takes the US `Cv`. The metric `Kv` is defined against `m**3/h` and `bar`
//! and they differ by about 1.156, so a caller holding `Kv` should convert once and
//! deliberately.

use azoth_core::units::{Pressure, cubic_meters_per_second};
use azoth_core::{Result, apply_checks};

use crate::results::ControlValveCvResult;
use crate::spec_gen;

/// One US gallon in cubic metres. Exact, by definition.
pub const GALLON_M3: f64 = 3.785_411_784e-3;

/// One pound-force in newtons. Exact, by definition.
pub const POUND_FORCE_N: f64 = 4.448_221_615_260_5;

/// One inch in metres. Exact, by definition.
pub const INCH_M: f64 = 0.025_4;

/// The conversion from the US `Cv` convention to SI.
///
/// `(1 gpm in m**3/s) / sqrt(1 psi in Pa)`, so that
/// `q[m**3/s] = CV_TO_SI * Cv * sqrt(dP[Pa] / SG)` is the conventional relation
/// `Q[gpm] = Cv * sqrt(dP[psi] / SG)` written for SI inputs. Every factor in it is
/// a definition rather than a measurement, so the constant is arithmetic.
///
/// Written as a literal rather than computed, because `f64::sqrt` is not available
/// in a `const` context. `the_conversion_constant_is_derivable` below re-derives it
/// from [`GALLON_M3`], [`INCH_M`] and [`POUND_FORCE_N`] and fails if the two ever
/// disagree - so this is a stored value with a proof, not a number somebody typed.
/// The Python side computes it from the same definitions rather than storing it, and
/// a cross-language test checks the two agree.
pub const CV_TO_SI: f64 = 7.598_054_212_083_37e-7;

/// One pound per square inch in pascals, from the two definitions above.
///
/// Not a `const` because `f64` division is not `const`, so it is a function. It
/// exists so the derivation appears once, in code, rather than in prose that can
/// drift from the arithmetic.
#[must_use]
pub fn psi_in_pascals() -> f64 {
    POUND_FORCE_N / (INCH_M * INCH_M)
}

/// Volumetric flow through a control valve. `Cv` and `SG` are dimensionless.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Cv` or `SG` is not positive, or if
///   `dP` is negative.
///
/// # Example
/// ```
/// use azoth_core::units::pascals;
/// use azoth_hydraulics::control_valve_cv;
///
/// let r = control_valve_cv(10.0, pascals(6894.757293168361), 1.0)?;
/// assert!((r.q.value - 6.309019640000001e-4).abs() < 1e-18);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Cv` and `SG` are the symbols in the published relation
pub fn control_valve_cv(Cv: f64, dP: Pressure, SG: f64) -> Result<ControlValveCvResult> {
    let spec = &spec_gen::CONTROL_VALVE_CV_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            // Dimensionless by declaration, not by nature: `Cv` carries
            // gpm/sqrt(psi) and `SG` is a ratio, so neither has a unit to extract.
            "Cv" => Some(Cv),
            "dP" => Some(dP.value),
            "SG" => Some(SG),
            _ => None,
        },
        &mut warnings,
    )?;

    // Guarded by the checks above, so SG is positive and dP is not negative.
    let q = CV_TO_SI * Cv * (dP.value / SG).sqrt();

    apply_checks(
        spec.derived_checks(),
        |name| (name == "q").then_some(q),
        &mut warnings,
    )?;

    Ok(ControlValveCvResult {
        q: cubic_meters_per_second(q),
        warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_conversion_constant_is_derivable() {
        // The point of the literal above. `const` cannot call sqrt, so the constant
        // is stored - and this re-derives it from the definitions of the gallon,
        // the pound-force and the inch, the same way the spec's derivation does.
        // Without this, `CV_TO_SI` would be a typed number with a citation
        // attached, which is the thing this project refuses to ship.
        let derived = (GALLON_M3 / 60.0) / psi_in_pascals().sqrt();
        assert!(
            (CV_TO_SI - derived).abs() < 1e-21,
            "CV_TO_SI is {CV_TO_SI}, but the definitions give {derived}"
        );
    }

    #[test]
    fn the_definitions_are_the_conventional_ones() {
        // Pinned because everything above rests on them, and because a "correction"
        // to a definition is the kind of edit that looks harmless.
        assert_eq!(GALLON_M3, 3.785_411_784e-3);
        assert_eq!(POUND_FORCE_N, 4.448_221_615_260_5);
        assert_eq!(INCH_M, 0.025_4);
        assert!((psi_in_pascals() - 6_894.757_293_168_361).abs() < 1e-9);
    }
}
