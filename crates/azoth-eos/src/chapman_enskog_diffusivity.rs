//! `eos.chapman_enskog_diffusivity` - the gas binary diffusivity from the Chapman-Enskog
//! theory.
//!
//! Spec: `specs/calcs/eos/chapman_enskog_diffusivity.toml`, which records the correlation, the
//! Neufeld collision integral and the two routes over the Lennard-Jones parameters.

use azoth_core::units::{MolarMass, Pressure, ThermodynamicTemperature, square_meters_per_second};
use azoth_core::{Result, apply_checks};

use crate::results::ChapmanEnskogDiffusivityResult;
use crate::spec_gen;

/// The Neufeld coefficients of the Lennard-Jones collision integral.
///
/// `Neufeld, Janzen & Aziz (1972)`, the same eight numbers the class writes inline.
const A: f64 = 1.06036;
const B: f64 = 0.15610;
const C: f64 = 0.19300;
const D: f64 = 0.47635;
const E: f64 = 1.03587;
const F: f64 = 1.52996;
const G: f64 = 1.76474;
const H: f64 = 3.89411;

/// The collision integral `Omega_D(T/eps)`, Neufeld's fit.
///
/// Exposed because it is a function of the pair alone, and the pair parameters are what a caller
/// combines - see [`crate::databank::lennard_jones_pair`].
#[must_use]
pub fn collision_integral(t: f64) -> f64 {
    A / t.powf(B) + C / (D * t).exp() + E / (F * t).exp() + G / (H * t).exp()
}

/// The binary diffusivity of a gas pair, from the Chapman-Enskog theory.
///
/// **The pair's parameters, not the components'**: `sigma` and `eps` are NeqSim's
/// `binaryMolecularDiameter` and `binaryEnergyParameter`, which
/// [`crate::databank::lennard_jones_pair`] combines from two components' values.
///
/// **This is a gas phase's default diffusivity.** `GasPhysicalProperties` assigns the base
/// `Diffusivity` class to a gas, and the Fuller correlation is reached only by
/// `setDiffusionCoefficientModel` - so a caller that selects nothing takes this, over the
/// database's Lennard-Jones parameters. The same arithmetic on Poling's textbook parameters is a
/// second route, and the capture records both: measured on methane/nitrogen at 298 K and one
/// atmosphere, the database route answers `3.555e-5` m**2/s where Marrero & Mason measure
/// `2.2e-5`, and the textbook route `2.185e-5`.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `T`, `sigma` or `eps` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, kilograms_per_mole, pascals};
/// use azoth_eos::chapman_enskog_diffusivity;
///
/// let r = chapman_enskog_diffusivity(
///     kilograms_per_mole(0.016043),
///     kilograms_per_mole(0.0280135),
///     2.826253374e-10,
///     kelvins(140.43522570075095),
///     kelvins(298.15),
///     pascals(101325.0),
/// )?;
/// assert!((r.d.value - 3.555070057344887e-5).abs() < 1e-18);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `MA`, `MB` and `T` are the symbols in the published equation
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn chapman_enskog_diffusivity(
    MA: MolarMass,
    MB: MolarMass,
    sigma: f64,
    eps: ThermodynamicTemperature,
    T: ThermodynamicTemperature,
    P: Pressure,
) -> Result<ChapmanEnskogDiffusivityResult> {
    let spec = &spec_gen::CHAPMAN_ENSKOG_DIFFUSIVITY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "MA" => Some(MA.value),
            "MB" => Some(MB.value),
            "sigma" => Some(sigma),
            "eps" => Some(eps.value),
            "T" => Some(T.value),
            "P" => Some(P.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // The published constants are tuned to g/mol, angstrom, bar and cm**2/s.
    let pair_mass_g = 2.0 / (1.0 / (MA.value * 1000.0) + 1.0 / (MB.value * 1000.0));
    let sigma_angstrom = sigma * 1.0e10;
    let p_bar = P.value * 1.0e-5;

    let omega = collision_integral(T.value / eps.value);
    let d_cm2s = 0.00266 * T.value.powf(1.5)
        / (p_bar * pair_mass_g.sqrt() * sigma_angstrom * sigma_angstrom * omega);
    let d = d_cm2s * 1.0e-4;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "d").then_some(d),
        &mut warnings,
    )?;

    Ok(ChapmanEnskogDiffusivityResult {
        d: square_meters_per_second(d),
        warnings,
    })
}
