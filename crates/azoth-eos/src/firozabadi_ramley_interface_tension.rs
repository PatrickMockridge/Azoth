//! `eos.firozabadi_ramley_interface_tension` - the interface tension from the
//! Firoozabadi-Ramey corresponding-states correlation.
//!
//! Spec: `specs/calcs/eos/firozabadi_ramley_interface_tension.toml`, which
//! records NeqSim's `FirozabadiRamleyInterfaceTension.calcSurfaceTension`.

use azoth_core::units::{MassDensity, ThermodynamicTemperature, newtons_per_meter};
use azoth_core::{Result, apply_checks};

use crate::results::FirozabadiRamleyInterfaceTensionResult;
use crate::spec_gen;

/// The interface tension between two phases from their mass densities and the
/// first phase's pseudo-critical temperature.
///
/// `pseudo_critical_temperature` is `sum(x_i Tc_i)` over the **first** phase's
/// composition; the class takes it from the first argument of
/// `calcSurfaceTension`, so the order of the two densities is part of the input.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `T`, the pseudo-critical
///   temperature, or either density is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_cubic_meter, kelvins};
/// use azoth_eos::firozabadi_ramley_interface_tension;
///
/// let r = firozabadi_ramley_interface_tension(kelvins(300.0), kelvins(540.2432308699056),
///     kilograms_per_cubic_meter(690.6031246582897),
///     kilograms_per_cubic_meter(1000.0024520207028))?;
/// assert!((r.sigma.value - 0.05083361859364053).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T` is the symbol in the published equation
pub fn firozabadi_ramley_interface_tension(
    T: ThermodynamicTemperature,
    pseudo_critical_temperature: ThermodynamicTemperature,
    rho_1: MassDensity,
    rho_2: MassDensity,
) -> Result<FirozabadiRamleyInterfaceTensionResult> {
    let spec = &spec_gen::FIROZABADI_RAMLEY_INTERFACE_TENSION_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T.value),
            "pseudo_critical_temperature" => Some(pseudo_critical_temperature.value),
            "rho_1" => Some(rho_1.value),
            "rho_2" => Some(rho_2.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let delta_over_1000 = (rho_2.value - rho_1.value).abs() / 1000.0;
    let (a1, b1) = if delta_over_1000 < 0.2 {
        (2.2062, -0.94716)
    } else if delta_over_1000 < 0.5 {
        (2.915, -0.76852)
    } else {
        (3.3858, -0.62590)
    };
    let reduced_temperature = T.value / pseudo_critical_temperature.value;
    let temp1 = a1 * delta_over_1000.powf(b1 + 1.0) / reduced_temperature.powf(0.3125);
    let sigma = temp1.powi(4) / 1000.0;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "sigma").then_some(sigma),
        &mut warnings,
    )?;

    Ok(FirozabadiRamleyInterfaceTensionResult {
        sigma: newtons_per_meter(sigma),
        warnings,
    })
}
