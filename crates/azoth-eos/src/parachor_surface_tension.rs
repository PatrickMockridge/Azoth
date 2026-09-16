//! `eos.parachor_surface_tension` - the surface tension from the parachor
//! (Macleod-Sugden) correlation.
//!
//! Spec: `specs/calcs/eos/parachor_surface_tension.toml`, which records the pure
//! component form of NeqSim's `ParachorSurfaceTension.calcPureComponentSurfaceTension`.

use azoth_core::units::{MassDensity, MolarMass, newtons_per_meter};
use azoth_core::{Result, apply_checks};

use crate::results::ParachorSurfaceTensionResult;
use crate::spec_gen;

/// The surface tension of a pure component, from the Macleod-Sugden parachor
/// correlation.
///
/// `parachor` is in the mixed unit `(mN/m)**(1/4) * cm**3/mol` NeqSim stores; the
/// `1e-6` in the formula converts its `cm**3/mol` to `m**3/mol`, and the `1e-3`
/// recovers `N/m` from `mN/m`. `rho_l` and `rho_v` are the liquid and vapour mass
/// densities at the interface state.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `M` or `rho_l` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_mole, kilograms_per_cubic_meter};
/// use azoth_eos::parachor_surface_tension;
///
/// let r = parachor_surface_tension(77.3, kilograms_per_cubic_meter(422.0),
///     kilograms_per_cubic_meter(1.82), kilograms_per_mole(0.016043))?;
/// assert!((r.sigma.value - 0.016800304320340388).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `M` is the symbol in the published equation
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn parachor_surface_tension(
    parachor: f64,
    rho_l: MassDensity,
    rho_v: MassDensity,
    M: MolarMass,
) -> Result<ParachorSurfaceTensionResult> {
    let spec = &spec_gen::PARACHOR_SURFACE_TENSION_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "parachor" => Some(parachor),
            "rho_l" => Some(rho_l.value),
            "rho_v" => Some(rho_v.value),
            "M" => Some(M.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let molar_density_l = rho_l.value / M.value;
    let molar_density_v = rho_v.value / M.value;
    let sigma = 1.0e-3 * (parachor * 1.0e-6 * (molar_density_l - molar_density_v)).powi(4);

    apply_checks(
        spec.derived_checks(),
        |name| (name == "sigma").then_some(sigma),
        &mut warnings,
    )?;

    Ok(ParachorSurfaceTensionResult {
        sigma: newtons_per_meter(sigma),
        warnings,
    })
}
