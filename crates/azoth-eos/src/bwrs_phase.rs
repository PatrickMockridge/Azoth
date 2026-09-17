//! `eos.bwrs_phase` - the BWRS (MBWR-32) phase state at a temperature and pressure.
//!
//! ```text
//! solve P(rho) = P for the molar density, then
//! Z     = 1 + rho dF/drho
//! ln_phi_i = d(nF)/dn_i - ln Z
//! h_res, s_res, cp_res = the residual Helmholtz departures
//! ```
//!
//! Spec: `specs/models/eos/bwrs_phase.toml`, which carries the assembly and what the
//! caller is responsible for: the MBWR-32 coefficients resolve by name, and the
//! ideal-gas part - which this model does not add - is `eos.molar_enthalpy_entropy`'s.

use azoth_core::units::{
    Pressure, ThermodynamicTemperature, joules_per_mole, joules_per_mole_kelvin,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::bwrs::{self, BwrsCoefficients};
use crate::bwrs_mixture;
use crate::model_gen;
use crate::results::BwrsPhaseResult;

/// The BWRS phase state of a mixture of MBWR-32 substances.
///
/// `coeffs` is the per-component coefficient set, resolved by name against
/// [`crate::databank::bwrs_coefficients`]. The density solve is Newton from the
/// ideal-gas guess in the coefficients' native mol/L, MPa units; the returned
/// departures are the residual (real minus ideal gas) properties in SI.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
/// * [`AzothError::InvalidInput`] if `z` is not a composition of `coeffs`'s length.
pub fn bwrs_phase(
    coeffs: &[BwrsCoefficients],
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<BwrsPhaseResult> {
    let spec = &model_gen::BWRS_PHASE_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "P" => Some(p.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let n = coeffs.len();
    if z.len() != n {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "a mixture of {n} components needs {n} mole fractions, but z has {}",
                z.len()
            ),
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
    if (sum - 1.0).abs() > 1.0e-9 {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "the composition sums to {sum}, not to one; renormalising here would hide a caller's error"
            ),
        ));
    }

    let tk = t.value;
    let p_mpa = p.value / 1.0e6;

    let bps: Vec<[f64; 9]> = coeffs.iter().map(|c| bwrs::bp(tk, &c.a)).collect();
    let bes: Vec<[f64; 6]> = coeffs.iter().map(|c| bwrs::be(tk, &c.a)).collect();
    let bts: Vec<[f64; 9]> = coeffs.iter().map(|c| bwrs::bp_dt(tk, &c.a)).collect();
    let ets: Vec<[f64; 6]> = coeffs.iter().map(|c| bwrs::be_dt(tk, &c.a)).collect();
    let btts: Vec<[f64; 9]> = coeffs.iter().map(|c| bwrs::bp_dt_dt(tk, &c.a)).collect();
    let etts: Vec<[f64; 6]> = coeffs.iter().map(|c| bwrs::be_dt_dt(tk, &c.a)).collect();
    let rhocs: Vec<f64> = coeffs.iter().map(|c| c.rhoc).collect();

    let (mb, me, gamma) = bwrs_mixture::mix(z, &bps, &bes, &rhocs);
    let (mbt, met, _) = bwrs_mixture::mix(z, &bts, &ets, &rhocs);
    let (mbtt, mett, _) = bwrs_mixture::mix(z, &btts, &etts, &rhocs);

    let rho = bwrs::solve_density(tk, p_mpa, &mb, &me, gamma);
    let z_factor = 1.0 + rho * bwrs::d_helmholtz_drho(tk, rho, &mb, &me, gamma);
    let ln_phi = bwrs_mixture::ln_fugacity_coefficients(tk, rho, z, &bps, &bes, &rhocs);
    let departure = bwrs::departure(tk, rho, &mb, &mbt, &mbtt, &me, &met, &mett, gamma);

    Ok(BwrsPhaseResult {
        z_factor,
        ln_phi,
        h_res: joules_per_mole(departure.h_res),
        s_res: joules_per_mole_kelvin(departure.s_res),
        cp_res: joules_per_mole_kelvin(departure.cp_res),
        warnings,
    })
}
