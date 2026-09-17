//! `eos.nitric_sulfuric_acid_vapor_pressure` - the three pure-component vapour pressures
//! of the water-nitric-sulfuric acid system.
//!
//! Spec: `specs/calcs/eos/nitric_sulfuric_acid_vapor_pressure.toml`, which records the
//! three correlations, the nitric-acid adjustment NeqSim made to Pennington's Antoine
//! pair, and why the stated range is a warning rather than an error.
//!
//! Ported from `thermo/util/empiric/NitricSulfuricAcidVaporPressure`. The activity
//! coefficients of the same system are a separate model,
//! [`crate::van_laar_acid_activity_coefficients`]; a phase pairs the two.

use azoth_core::units::{ThermodynamicTemperature, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::NitricSulfuricAcidVaporPressureResult;
use crate::spec_gen;

/// Pascal per millibar, the unit the water correlation is written in.
const MBAR_TO_PA: f64 = 100.0;

/// Pascal per torr, the unit the nitric-acid Antoine is written in.
const TORR_TO_PA: f64 = 133.322368_421;

/// Pascal per atmosphere, the unit the sulfuric-acid correlation is written in.
const ATM_TO_PA: f64 = 101_325.0;

/// The nitric-acid Antoine coefficients, NeqSim's refit of Pennington's pair.
///
/// `A` and `B` are fitted to the 83 C normal boiling point and to a 6.9 % rise in `P0` at
/// 273.15 K against Vandoni (1944)'s salting-out data; `C` is Pennington's, unchanged.
/// See the spec's assumptions - the port carries NeqSim's values rather than the paper's.
const HNO3_ANTOINE_A: f64 = 7.57628;
const HNO3_ANTOINE_B: f64 = 1470.385;
const HNO3_ANTOINE_C: f64 = 43.0;

/// The pure-component vapour pressures of water, nitric acid and sulfuric acid.
///
/// The three come back together because a phase over this system needs all three at one
/// temperature. Each is its own correlation: water a `log10 P/mbar` polynomial in `1/T`,
/// nitric acid an Antoine in `log10 P/torr`, sulfuric acid a `ln P/atm` straight line.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `T` is not positive, or at or below the nitric-acid
///   form's 43 K pole.
///
/// An out-of-range `T` carries `OUT_OF_VALID_RANGE` rather than failing: the arithmetic
/// is well defined outside 190-298 K, and inspecting the extrapolation is legitimate.
///
/// # Example
/// ```
/// use azoth_core::units::kelvins;
/// use azoth_eos::nitric_sulfuric_acid_vapor_pressure;
///
/// let r = nitric_sulfuric_acid_vapor_pressure(kelvins(273.15))?;
/// assert!((r.p_water.value - 610.3592249571752).abs() < 1e-9);
/// assert!((r.p_nitric_acid.value - 2052.9169266045096).abs() < 1e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn nitric_sulfuric_acid_vapor_pressure(
    t: ThermodynamicTemperature,
) -> Result<NitricSulfuricAcidVaporPressureResult> {
    let spec = &spec_gen::NITRIC_SULFURIC_ACID_VAPOR_PRESSURE_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // The nitric-acid form is the one with a pole inside the domain: `B/(T - C)` is
    // infinite at `T = C`. The spec bounds `T` away from zero, which does not reach 43 K,
    // so it is refused here - the one arithmetic hazard the three correlations have.
    if t.value <= HNO3_ANTOINE_C {
        return Err(AzothError::out_of_range(
            "T",
            t.value,
            "the nitric-acid Antoine form is `B/(T - 43.0)`, so 43 K and below are a \
             pole rather than a pressure",
        ));
    }

    apply_checks(
        spec.derived_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // Water: `log10(P/mbar) = 8.42926609 - 1827.17843/T - 71208.271/T^2`.
    let log10_mbar = 8.42926609 - 1827.17843 / t.value - 71208.271 / (t.value * t.value);
    let p_water = 10.0_f64.powf(log10_mbar) * MBAR_TO_PA;

    // Nitric acid: `log10(P/torr) = A - B/(T - C)`.
    let log10_torr = HNO3_ANTOINE_A - HNO3_ANTOINE_B / (t.value - HNO3_ANTOINE_C);
    let p_nitric_acid = 10.0_f64.powf(log10_torr) * TORR_TO_PA;

    // Sulfuric acid: `ln(P/atm) = -10156.0/T + 16.259`.
    let ln_atm = -10156.0 / t.value + 16.259;
    let p_sulfuric_acid = ln_atm.exp() * ATM_TO_PA;

    Ok(NitricSulfuricAcidVaporPressureResult {
        p_water: pascals(p_water),
        p_nitric_acid: pascals(p_nitric_acid),
        p_sulfuric_acid: pascals(p_sulfuric_acid),
        warnings,
    })
}
