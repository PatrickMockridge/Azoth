//! `eos.parachor_mixture_surface_tension` - the interface surface tension from the parachor
//! (Weinaug-Katz) correlation.
//!
//! Spec: `specs/models/eos/parachor_mixture_surface_tension.toml`, which records the mixture
//! form, the sign convention and what the class answers where this refuses.

use azoth_core::units::{MassDensity, MolarMass, newtons_per_meter};
use azoth_core::{AzothError, Result, apply_checks};

use crate::model_gen;
use crate::results::ParachorMixtureSurfaceTensionResult;

/// The surface tension of the interface between a gas and a liquid, from the parachor.
///
/// The mixture form of [`crate::parachor_surface_tension`]: the pure-component relation's mole
/// fractions divide out, and this one sums a per-component **molar-density** difference,
/// `rho_liquid/M_liquid * x_liquid_i - rho_gas/M_gas * x_gas_i`, over the components, raising
/// the total to the fourth.
///
/// `parachors` are the components' `PARACHOR` values, in NeqSim's mixed unit
/// `(mN/m)**(1/4) * cm**3/mol`, and both compositions are in that same component order.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if the three vectors differ in length.
/// * [`azoth_core::AzothError::OutOfRange`] if a density or a molar mass is not positive -
///   where the class divides, catches the exception and answers `0.0`.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
/// use azoth_eos::parachor_mixture_surface_tension;
///
/// let r = parachor_mixture_surface_tension(
///     &[77.3, 191.7],
///     kilograms_per_cubic_meter(19.938328330315997),
///     kilograms_per_mole(0.022959902730870334),
///     &[0.8356249351028912, 0.1643750648971088],
///     kilograms_per_cubic_meter(549.3851668858096),
///     kilograms_per_mole(0.054107691623917306),
///     &[0.09542082642782022, 0.9045791735721797],
/// )?;
/// assert!((r.sigma.value - 0.009424889197268284).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `M_gas` and `M_liquid` are the symbols in the published equation
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn parachor_mixture_surface_tension(
    parachors: &[f64],
    rho_gas: MassDensity,
    M_gas: MolarMass,
    x_gas: &[f64],
    rho_liquid: MassDensity,
    M_liquid: MolarMass,
    x_liquid: &[f64],
) -> Result<ParachorMixtureSurfaceTensionResult> {
    let spec = &model_gen::PARACHOR_MIXTURE_SURFACE_TENSION_SPEC;
    let mut warnings = Vec::new();

    if parachors.len() != x_gas.len() || parachors.len() != x_liquid.len() {
        return Err(AzothError::invalid_input(
            "parachors",
            format!(
                "{} parachor(s), {} gas fraction(s) and {} liquid fraction(s): the sum is over \
                 the components, so the three vectors are one entry each",
                parachors.len(),
                x_gas.len(),
                x_liquid.len()
            ),
        ));
    }

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "rho_gas" => Some(rho_gas.value),
            "rho_liquid" => Some(rho_liquid.value),
            "M_gas" => Some(M_gas.value),
            "M_liquid" => Some(M_liquid.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // `rho/M` per phase, which is the molar density each component's fraction weights.
    let molar_density_gas = rho_gas.value / M_gas.value;
    let molar_density_liquid = rho_liquid.value / M_liquid.value;

    let mut total = 0.0;
    for ((parachor, gas_fraction), liquid_fraction) in
        parachors.iter().zip(x_gas.iter()).zip(x_liquid.iter())
    {
        total += parachor
            * 1.0e-6
            * (molar_density_liquid * liquid_fraction - molar_density_gas * gas_fraction);
    }
    let sigma = 1.0e-3 * total.powi(4);

    apply_checks(
        spec.derived_checks(),
        |name| (name == "sigma").then_some(sigma),
        &mut warnings,
    )?;

    Ok(ParachorMixtureSurfaceTensionResult {
        sigma: newtons_per_meter(sigma),
        warnings,
    })
}
