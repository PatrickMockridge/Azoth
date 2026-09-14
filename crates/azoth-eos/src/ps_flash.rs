//! `eos.ps_flash` - the temperature a mixture reaches at a given pressure and entropy.
//!
//! Spec: `specs/models/eos/ps_flash.yaml`
//!
//! An outer bisection on temperature, for an **isentropic** unit operation - a
//! compressor, an expander, a turbine or a nozzle, assumed ideal - whose outlet
//! temperature is not known because entropy is conserved and temperature is not. The
//! bracket, the loop, the phase branch and the warning handling live in
//! [`crate::flash_property`], shared with `eos.ph_flash`.
//!
//! Each phase's entropy is taken at *that phase's own composition*, so the weighted sum
//! carries the entropy of mixing; summing the phases at the feed composition instead
//! would conserve entropy across a phase change, which is wrong.

use azoth_core::units::{MolarHeatCapacity, Pressure, ThermodynamicTemperature, kelvins};
use azoth_core::{Result, apply_checks};

use crate::algorithm_of;
use crate::flash_property::{Property, property_at, solve_temperature};
use crate::mixture::Mixture;
use crate::model_gen;
use crate::molar_enthalpy_entropy::IdealGasModel;
use crate::results::{PsFlashResult, PtFlashResult};

/// The molar entropy of a mixture at a temperature and pressure, and its split.
///
/// # Errors
/// Propagates whatever [`crate::pt_flash`] and [`crate::molar_enthalpy_entropy`] raise.
pub fn entropy_at(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<(f64, PtFlashResult)> {
    property_at(mixture, ideal_gas, t, p, z, Property::Entropy)
}

/// The temperature at which a mixture has a given molar entropy at a pressure.
///
/// `s` is a *difference* from the datum `ideal_gas` carries, not an absolute quantity -
/// the same caveat `eos.ph_flash` and `eos.molar_enthalpy_entropy` carry.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `p` is not positive, or a range check on the answer
///   fails.
/// * [`azoth_core::AzothError::InvalidInput`] if the spec declares no bracket, which would be a
///   generator bug rather than a caller's.
/// * [`azoth_core::AzothError::SolverNotConverged`] if no temperature on the bracket covers the
///   requested entropy, or if the bisection reaches its cap.
///
/// # Example
/// ```
/// use azoth_core::units::{joules_per_mole_kelvin, pascals, kelvins};
/// use azoth_eos::mixture::{Component, Mixture};
/// use azoth_eos::molar_enthalpy_entropy::IdealGasModel;
/// use azoth_eos::ps_flash;
///
/// let mixture = Mixture::new(
///     vec![
///         Component::new(kelvins(190.56), pascals(4_599_000.0), 0.0115)?,
///         Component::new(kelvins(425.12), pascals(3_796_000.0), 0.2002)?,
///     ],
///     vec![0.0, 0.01289789, 0.01289789, 0.0],
/// )?;
/// let ideal_gas = IdealGasModel {
///     cp_a: vec![3.0, 5.0],
///     cp_b: vec![0.0, 0.0],
///     cp_c: vec![0.0, 0.0],
///     cp_d: vec![0.0, 0.0],
///     cp_e: vec![0.0, 0.0],
/// };
/// let r = ps_flash::ps_flash(
///     &mixture,
///     &ideal_gas,
///     pascals(2_000_000.0),
///     joules_per_mole_kelvin(-38.61276026788922),
///     &[0.6, 0.4],
/// )?;
/// assert!((r.temperature.value - 300.0).abs() < 1.0e-6);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn ps_flash(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    p: Pressure,
    s: MolarHeatCapacity,
    z: &[f64],
) -> Result<PsFlashResult> {
    let spec = &model_gen::PS_FLASH_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "P" => Some(p.value),
            "S" => Some(s.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let algorithm = algorithm_of(spec)?;
    let solved = solve_temperature(
        mixture,
        ideal_gas,
        p,
        s.value,
        z,
        Property::Entropy,
        algorithm,
    )?;

    apply_checks(
        spec.derived_checks(),
        |quantity| (quantity == "T").then_some(solved.temperature),
        &mut warnings,
    )?;

    Ok(PsFlashResult {
        temperature: kelvins(solved.temperature),
        beta: solved.flash.beta,
        x: solved.flash.x,
        y: solved.flash.y,
        k: solved.flash.k,
        phase: solved.flash.phase,
        z_liquid: solved.flash.z_liquid,
        z_vapour: solved.flash.z_vapour,
        iterations: solved.iterations,
        residual: solved.residual,
        warnings: solved.warnings,
    })
}
