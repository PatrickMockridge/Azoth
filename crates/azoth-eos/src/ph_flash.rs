//! `eos.ph_flash` - the temperature a mixture reaches at a given pressure and enthalpy.
//!
//! Spec: `specs/models/eos/ph_flash.toml`
//!
//! A damped quasi-Newton on temperature, in `1/T`, over the enthalpy assembled from two
//! things this crate already has: the phase split at a trial temperature, from
//! [`crate::pt_flash`], and each phase's enthalpy, from
//! [`crate::molar_enthalpy_entropy`]. Upstream's `PHflash.solveQ`, which is the scheme
//! NeqSim runs by default - see the spec's assumption on the second-order alternative it
//! carries.
//!
//! The damping, the step clamp, the trial-temperature recovery and the phase branch live
//! in [`crate::flash_property`], shared with `eos.ps_flash`.

use azoth_core::units::{MolarEnergy, Pressure, ThermodynamicTemperature, kelvins};
use azoth_core::{Result, apply_checks};

use crate::algorithm_of;
use crate::flash_property::{Property, property_at, solve_temperature};
use crate::mixture::Mixture;
use crate::model_gen;
use crate::molar_enthalpy_entropy::IdealGasModel;
use crate::results::{PhFlashResult, PtFlashResult};

/// The molar enthalpy of a mixture at a temperature and pressure, and its split.
///
/// The composition the flash settles on is the equilibrium one, so this is the enthalpy
/// of the *feed* at that state - which is what makes it comparable with a duty a caller
/// supplied.
///
/// # Errors
/// Propagates whatever [`crate::pt_flash`] and [`crate::molar_enthalpy_entropy`] raise.
pub fn enthalpy_at(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<(f64, PtFlashResult)> {
    property_at(mixture, ideal_gas, t, p, z, Property::Enthalpy)
}

/// The temperature at which a mixture has a given molar enthalpy at a pressure.
///
/// `h` is a *difference* from the datum `ideal_gas` carries, not an absolute quantity:
/// two calls with different reference values are not comparable, and their difference is
/// a plausible number rather than an error.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `p` is not positive, or a range check on the answer
///   fails.
/// * [`azoth_core::AzothError::InvalidInput`] if the spec declares no starting temperature,
///   which would be a generator bug rather than a caller's.
/// * [`azoth_core::AzothError::SolverNotConverged`] if the iteration reaches its cap.
///
/// # Example
/// ```
/// use azoth_core::units::{joules_per_mole, pascals};
/// use azoth_eos::{databank, ph_flash};
///
/// let (mixture, ideal_gas) = databank::mixture_of(&["methane", "n-butane"], None)
///     .expect("the pair resolves");
/// let r = ph_flash::ph_flash(
///     &mixture,
///     &ideal_gas,
///     pascals(2_000_000.0),
///     joules_per_mole(-5121.329517317879),
///     &[0.6, 0.4],
/// )?;
/// assert!((r.temperature.value - 300.0).abs() < 1.0e-6);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn ph_flash(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    p: Pressure,
    h: MolarEnergy,
    z: &[f64],
) -> Result<PhFlashResult> {
    let spec = &model_gen::PH_FLASH_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "P" => Some(p.value),
            "H" => Some(h.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let algorithm = algorithm_of(spec)?;
    let solved = solve_temperature(
        mixture,
        ideal_gas,
        p,
        h.value,
        z,
        Property::Enthalpy,
        algorithm,
    )?;

    apply_checks(
        spec.derived_checks(),
        |quantity| (quantity == "T").then_some(solved.temperature),
        &mut warnings,
    )?;

    Ok(PhFlashResult {
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
