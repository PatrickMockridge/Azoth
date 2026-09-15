//! `eos.wilke_viscosity` - the gas mixture viscosity from Wilke's rule over
//! Chung's pure-component viscosities.
//!
//! Spec: `specs/models/eos/wilke_viscosity.toml`, which records the mixing rule and
//! points the pure-component half at [`crate::chung_viscosity`]. It is a *direct*
//! model: no iteration, so no `algorithm` block.

use azoth_core::units::{cubic_meters_per_mole, kelvins, kilograms_per_mole, pascal_seconds};
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::WilkeViscosityResult;
use crate::{chung_viscosity, model_gen};

/// The gas dynamic viscosity of a mixture, by Wilke's mixing rule over the pure
/// component Chung viscosities.
///
/// Every per-component constant is an SI magnitude: `Tc` in K, `Vc` in m**3/mol,
/// `M` in kg/mol, and `omega`, `dipole` (debye) and `kappa` dimensionless. `z` is
/// checked (non-negative, sums to one) rather than renormalised. Each pure
/// viscosity is [`crate::chung_viscosity`]'s at the *mixture* molar volume `V`.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the vectors disagree in length, or `z` is not a
///   composition.
/// * [`AzothError::OutOfRange`] if `T` or `V` is not positive, or a component's
///   `Tc`/`Vc` is not positive.
///
/// # Example
/// ```
/// use azoth_eos::wilke_viscosity;
///
/// let r = wilke_viscosity(
///     &[190.56, 369.83], &[9.9e-5, 0.000203], &[0.016043, 0.044097],
///     &[0.0115, 0.1523], &[0.0, 0.0], &[0.0, 0.0],
///     300.0, 2.2987856717900145e-3, &[0.5, 0.5],
/// )?;
/// assert!((r.mu.value - 9.543074369236108e-6).abs() < 1e-18);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tc`, `Vc`, `M`, `T` and `V` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn wilke_viscosity(
    Tc: &[f64],
    Vc: &[f64],
    M: &[f64],
    omega: &[f64],
    dipole: &[f64],
    kappa: &[f64],
    T: f64,
    V: f64,
    z: &[f64],
) -> Result<WilkeViscosityResult> {
    let spec = &model_gen::WILKE_VISCOSITY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T),
            "V" => Some(V),
            _ => None,
        },
        &mut warnings,
    )?;

    let n = Tc.len();
    if [
        Vc.len(),
        M.len(),
        omega.len(),
        dipole.len(),
        kappa.len(),
        z.len(),
    ]
    .iter()
    .any(|&length| length != n)
    {
        return Err(AzothError::invalid_input(
            "Tc",
            format!(
                "the per-component vectors disagree in length: Tc has {n} entries, \
                 Vc {}, M {}, omega {}, dipole {}, kappa {}, z {}",
                Vc.len(),
                M.len(),
                omega.len(),
                dipole.len(),
                kappa.len(),
                z.len()
            ),
        ));
    }
    if n == 0 {
        return Err(AzothError::invalid_input(
            "Tc",
            "a mixture of zero components has no viscosity",
        ));
    }
    if let Some(bad) = z.iter().position(|&value| value < 0.0) {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "z[{bad}] is {} but a mole fraction cannot be negative",
                z[bad]
            ),
        ));
    }
    let sum: f64 = z.iter().sum();
    if (sum - 1.0).abs() > 1.0e-09 {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "the mole fractions sum to {sum}, not to one. Renormalising them here \
                 would make a composition error invisible in every number downstream, \
                 so it is refused instead"
            ),
        ));
    }

    let mut pure_mu = Vec::with_capacity(n);
    for i in 0..n {
        let r = chung_viscosity(
            omega[i],
            kelvins(Tc[i]),
            cubic_meters_per_mole(Vc[i]),
            kilograms_per_mole(M[i]),
            dipole[i],
            kappa[i],
            kelvins(T),
            cubic_meters_per_mole(V),
        )?;
        warnings.extend(r.warnings.iter().cloned());
        pure_mu.push(r.mu.value);
    }

    let mut mu = 0.0;
    for i in 0..n {
        let mut denominator = 0.0;
        for j in 0..n {
            let phi = (1.0 + (pure_mu[i] / pure_mu[j]).sqrt() * (M[j] / M[i]).powf(0.25)).powi(2)
                / (8.0 * (1.0 + M[i] / M[j])).sqrt();
            denominator += z[j] * phi;
        }
        mu += z[i] * pure_mu[i] / denominator;
    }

    Ok(WilkeViscosityResult {
        mu: pascal_seconds(mu),
        warnings,
    })
}
