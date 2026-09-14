//! `eos.rk_departure` - the Redlich-Kwong fugacity coefficient and departures.
//!
//! ```text
//! psi      = -1/2
//! I        = ln((z + B) / z)
//! ln_phi   = z - 1 - ln(z - B) - C*I
//! h_dep_rt = (z - 1) + C*(psi - 1)*I
//! s_dep_r  = ln(z - B) + C*psi*I          where C = A / B
//! ```
//!
//! Spec: `specs/calcs/eos/rk_departure.toml`. The geometry is SRK's `delta = (1, 0)`,
//! and the alpha derivative is the constant `-1/2` - RK's `1/sqrt(Tr)` has no fitted
//! coefficient, so `psi` and `T*dpsi/dT` are constants rather than functions of `Tr`.

use azoth_core::{Result, apply_checks};

use crate::alpha_term::{AlphaTerm, RkAlpha};
use crate::cubic::Cubic;
use crate::results::RkDepartureResult;
use crate::spec_gen;

/// The Redlich-Kwong fugacity coefficient and departure functions, for one state.
///
/// Three dimensionless arguments and four dimensionless outputs. There is no `kappa`
/// and no `Tr`: RK's alpha derivative is the constant `-1/2`, so neither appears.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `b_reduced <= 0` or `z <= b_reduced`.
///
/// # Example
/// ```
/// use azoth_eos::rk_departure;
///
/// let d = rk_departure(0.1866943088347048, 0.02707510936404929, 0.8119920001727409)?;
/// assert!((d.ln_phi + 0.17200180559684855).abs() < 1e-12);
/// assert!((d.h_dep_rt - d.s_dep_r - d.ln_phi).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn rk_departure(a_reduced: f64, b_reduced: f64, z: f64) -> Result<RkDepartureResult> {
    let spec = &spec_gen::RK_DEPARTURE_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "a_reduced" => Some(a_reduced),
            "b_reduced" => Some(b_reduced),
            "z" => Some(z),
            _ => None,
        },
        &mut warnings,
    )?;

    // RK's logarithmic derivative is the constant -1/2 and its temperature
    // derivative is zero; the `RkAlpha` term states both, and the argument it
    // ignores is the reduced temperature no departure quantity here needs.
    let term = RkAlpha;
    let psi = term.psi(1.0);
    let t_dpsi = term.psi_t(1.0);

    let i_term = Cubic::Rk.i_term(z, b_reduced);
    let coefficient = Cubic::Rk.coefficient(a_reduced, b_reduced);
    let ln_z_minus_b = (z - b_reduced).ln();

    let ln_phi = z - 1.0 - ln_z_minus_b - coefficient * i_term;
    let h_dep_rt = (z - 1.0) + coefficient * (psi - 1.0) * i_term;
    let s_dep_r = ln_z_minus_b + coefficient * psi * i_term;

    // The heat-capacity departure, from `Cp^R/R = y + T*(dy/dT)_P` with `y = h_dep_rt`.
    let t_da = a_reduced * (psi - 2.0);
    let t_db = -b_reduced;
    let t_dc = coefficient * (psi - 1.0);

    let d_f_dz = Cubic::Rk.df_dz(z, a_reduced, b_reduced);
    let t_dfdt = Cubic::Rk.t_dfdt(z, a_reduced, b_reduced, t_da, t_db);
    let t_dz = -t_dfdt / d_f_dz;

    let d1 = Cubic::Rk.delta1();
    let d2 = Cubic::Rk.delta2();
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

    Ok(RkDepartureResult {
        ln_phi,
        h_dep_rt,
        s_dep_r,
        cp_dep_r,
        warnings,
    })
}
