//! `eos.mason_saxena_conductivity` - the gas mixture conductivity from
//! Mason-Saxena mixing over Chung's pure-component conductivities.
//!
//! Spec: `specs/models/eos/mason_saxena_conductivity.toml`, which records the mixing
//! rule and points the pure-component half at [`crate::chung_conductivity`]. A
//! *direct* model: no iteration, so no `algorithm` block.

use azoth_core::units::{
    cubic_meters_per_mole, joules_per_mole_kelvin, kelvins, kilograms_per_mole,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::MasonSaxenaConductivityResult;
use crate::{chung_conductivity, model_gen};

/// The gas thermal conductivity of a mixture, by Mason-Saxena mixing over the pure
/// component Chung conductivities.
///
/// Every per-component constant is an SI magnitude: `Cv0` in J/(mol*K), `M` in
/// kg/mol, `Tc` in K, `Vc` in m**3/mol, and `omega`, `dipole` (debye) and `kappa`
/// dimensionless. `z` is checked (non-negative, sums to one) rather than
/// renormalised. Each pure conductivity is [`crate::chung_conductivity`]'s at `T`.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the vectors disagree in length, or `z` is not a
///   composition.
/// * [`AzothError::OutOfRange`] if `T` is not positive, or a component's `Tc`/`Vc` is
///   not positive.
///
/// # Example
/// ```
/// use azoth_eos::mason_saxena_conductivity;
///
/// let r = mason_saxena_conductivity(
///     &[27.544151394, 65.809299386], &[0.016043, 0.044097],
///     &[0.0115, 0.1523], &[190.56, 369.83], &[9.9e-5, 0.000203],
///     &[0.0, 0.0], &[0.0, 0.0], 300.0, &[0.5, 0.5],
/// )?;
/// assert!((r.k.value - 0.025062059480046174).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Cv0`, `Tc`, `Vc`, `M` and `T` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn mason_saxena_conductivity(
    Cv0: &[f64],
    M: &[f64],
    omega: &[f64],
    Tc: &[f64],
    Vc: &[f64],
    dipole: &[f64],
    kappa: &[f64],
    T: f64,
    z: &[f64],
) -> Result<MasonSaxenaConductivityResult> {
    let spec = &model_gen::MASON_SAXENA_CONDUCTIVITY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| (quantity == "T").then_some(T),
        &mut warnings,
    )?;

    let n = Cv0.len();
    if [
        M.len(),
        omega.len(),
        Tc.len(),
        Vc.len(),
        dipole.len(),
        kappa.len(),
        z.len(),
    ]
    .iter()
    .any(|&length| length != n)
    {
        return Err(AzothError::invalid_input(
            "Cv0",
            format!(
                "the per-component vectors disagree in length: Cv0 has {n} entries, \
                 M {}, omega {}, Tc {}, Vc {}, dipole {}, kappa {}, z {}",
                M.len(),
                omega.len(),
                Tc.len(),
                Vc.len(),
                dipole.len(),
                kappa.len(),
                z.len()
            ),
        ));
    }
    if n == 0 {
        return Err(AzothError::invalid_input(
            "Cv0",
            "a mixture of zero components has no conductivity",
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

    let mut pure_k = Vec::with_capacity(n);
    for i in 0..n {
        let r = chung_conductivity(
            joules_per_mole_kelvin(Cv0[i]),
            kilograms_per_mole(M[i]),
            omega[i],
            kelvins(Tc[i]),
            cubic_meters_per_mole(Vc[i]),
            dipole[i],
            kappa[i],
            kelvins(T),
        )?;
        warnings.extend(r.warnings.iter().cloned());
        pure_k.push(r.k.value);
    }

    let mut k = 0.0;
    for i in 0..n {
        let mut denominator = 0.0;
        for j in 0..n {
            let aij = (1.0 + (pure_k[i] / pure_k[j]).sqrt() * (M[i] / M[j]).powf(0.25)).powi(2)
                / (8.0 * (1.0 + M[i] / M[j])).sqrt();
            denominator += z[j] * aij;
        }
        k += z[i] * pure_k[i] / denominator;
    }

    Ok(MasonSaxenaConductivityResult {
        k: azoth_core::units::watts_per_meter_kelvin(k),
        warnings,
    })
}
