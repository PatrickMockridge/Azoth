//! `eos.srk_departure` - the Soave-Redlich-Kwong fugacity coefficient and departures.
//!
//! ```text
//! psi      = -kappa*sqrt(Tr) / (1 + kappa*(1 - sqrt(Tr)))
//! I        = ln((z + B) / z)
//! ln_phi   = z - 1 - ln(z - B) - C*I
//! h_dep_rt = (z - 1) + C*(psi - 1)*I
//! s_dep_r  = ln(z - B) + C*psi*I          where C = A / B
//! ```
//!
//! Spec: `specs/calcs/eos/srk_departure.toml`, which carries the provenance and the
//! Gibbs identity `h_dep_rt - s_dep_r = ln_phi`.
//!
//! The geometry is SRK's `delta = (1, 0)`, so `I = ln((z + B)/z)` and `C = A/B`; the
//! alpha derivative `psi` is the same Soave term both cubics share.

use azoth_core::{Result, apply_checks};

use crate::alpha_term::{AlphaTerm, Soave};
use crate::cubic::Cubic;
use crate::results::SrkDepartureResult;
use crate::spec_gen;

/// The Soave-Redlich-Kwong fugacity coefficient and departure functions, for one state.
///
/// All five arguments and all four outputs are dimensionless. `kappa` comes from
/// [`crate::srk_kappa`].
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `b_reduced <= 0` or `z <= b_reduced`.
///
/// # Example
/// ```
/// use azoth_eos::{srk_departure, srk_kappa};
///
/// let kappa = srk_kappa(0.152)?.kappa;
/// let d = srk_departure(0.19315231739218255, 0.02707510936404929, 0.8019551972557889, kappa, 0.8)?;
/// assert!((d.ln_phi + 0.1798730824039273).abs() < 1e-12);
/// assert!((d.h_dep_rt - d.s_dep_r - d.ln_phi).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn srk_departure(
    a_reduced: f64,
    b_reduced: f64,
    z: f64,
    kappa: f64,
    Tr: f64,
) -> Result<SrkDepartureResult> {
    let spec = &spec_gen::SRK_DEPARTURE_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "a_reduced" => Some(a_reduced),
            "b_reduced" => Some(b_reduced),
            "z" => Some(z),
            "kappa" => Some(kappa),
            "Tr" => Some(Tr),
            _ => None,
        },
        &mut warnings,
    )?;

    // The logarithmic derivative of the alpha function, from the same Soave term the
    // mixture layer uses.
    let term = Soave { kappa };
    let psi = term.psi(Tr);

    let i_term = Cubic::Srk.i_term(z, b_reduced);
    let coefficient = Cubic::Srk.coefficient(a_reduced, b_reduced);
    let ln_z_minus_b = (z - b_reduced).ln();

    let ln_phi = z - 1.0 - ln_z_minus_b - coefficient * i_term;
    let h_dep_rt = (z - 1.0) + coefficient * (psi - 1.0) * i_term;
    let s_dep_r = ln_z_minus_b + coefficient * psi * i_term;

    // The heat-capacity departure, from `Cp^R/R = y + T*(dy/dT)_P` with `y = h_dep_rt`.
    let t_da = a_reduced * (psi - 2.0);
    let t_db = -b_reduced;
    let t_dpsi = term.psi_t(Tr);
    let t_dc = coefficient * (psi - 1.0);

    let d_f_dz = Cubic::Srk.df_dz(z, a_reduced, b_reduced);
    let t_dfdt = Cubic::Srk.t_dfdt(z, a_reduced, b_reduced, t_da, t_db);
    let t_dz = -t_dfdt / d_f_dz;

    let d1 = Cubic::Srk.delta1();
    let d2 = Cubic::Srk.delta2();
    let n_plus = z + d1 * b_reduced;
    let n_minus = z + d2 * b_reduced;
    let t_di = (t_dz + d1 * t_db) / n_plus - (t_dz + d2 * t_db) / n_minus;

    let cp_dep_r = h_dep_rt
        + t_dz
        + t_dc * (psi - 1.0) * i_term
        + coefficient * t_dpsi * i_term
        + coefficient * (psi - 1.0) * t_di;

    apply_checks(
        spec.derived_checks(),
        |quantity| match quantity {
            "z_minus_b_reduced" => Some(z - b_reduced),
            _ => None,
        },
        &mut warnings,
    )?;

    Ok(SrkDepartureResult {
        ln_phi,
        h_dep_rt,
        s_dep_r,
        cp_dep_r,
        warnings,
    })
}
