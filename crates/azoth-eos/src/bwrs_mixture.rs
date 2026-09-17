//! Mixture-level BWRS: linear coefficient mixing and the fugacity coefficients.
//!
//! The MBWR-32 coefficients mix linearly in the mole fraction - `mixBP[k] = sum_j x_j
//! BP_j[k]`, `mixBE[k] = sum_j x_j BE_j[k]` - and the critical density mixes the same way,
//! `rhoc_mix = sum_j x_j rhoc_j`, before `gamma = 1/rhoc_mix^2`. That is NeqSim's
//! `PhaseBWRSEos.computeMixedParameters`, restated in the kernel's native mol/L, MPa units.
//!
//! The fugacity coefficient follows the Helmholtz definition
//!
//! ```text
//! ln phi_i = d(nF)/dn_i - ln Z,   Z = 1 + rho dF/drho
//! ```
//!
//! where `d(nF)/dn_i` is the mole-number derivative of the extensive Helmholtz at fixed
//! temperature and volume. NeqSim differentiates this numerically (`PhaseBWRSEos.getdFdN`
//! perturbs `n_i` and renormalises the composition and density through the fixed volume),
//! and this module does the same, so the derivative is oracle-checked against NeqSim rather
//! than carried as a hand-derived composition derivative. NeqSim's own `getLogFugacity`
//! is never computed for BWRS - it returns zero - so the only faithful oracle is
//! `getdFdN` itself, which the tests check against. Its `Z` is the physical one (only
//! `calcPressure2` carries the ten-fold error), and this module's `Z = 1 + rho dF/drho`
//! agrees with it.

use crate::bwrs::{BwrsCoefficients, be, bp, d_helmholtz_drho, helmholtz};

/// The linearly-mixed coefficients at one composition: `mixBP`, `mixBE`, and `gamma`.
#[must_use]
pub fn mix(
    x: &[f64],
    bps: &[[f64; 9]],
    bes: &[[f64; 6]],
    rhocs: &[f64],
) -> ([f64; 9], [f64; 6], f64) {
    let mut mb = [0.0; 9];
    let mut me = [0.0; 6];
    let mut rhoc_mix = 0.0;
    for (j, &xj) in x.iter().enumerate() {
        for (k, slot) in mb.iter_mut().enumerate() {
            *slot += xj * bps[j][k];
        }
        for (k, slot) in me.iter_mut().enumerate() {
            *slot += xj * bes[j][k];
        }
        rhoc_mix += xj * rhocs[j];
    }
    (mb, me, 1.0 / (rhoc_mix * rhoc_mix))
}

/// The extensive Helmholtz `n A^res/(R T)` of `n` moles at molar density `rho`.
#[must_use]
pub fn extensive_helmholtz(
    t: f64,
    rho: f64,
    n: f64,
    x: &[f64],
    bps: &[[f64; 9]],
    bes: &[[f64; 6]],
    rhocs: &[f64],
) -> f64 {
    let (mb, me, gamma) = mix(x, bps, bes, rhocs);
    n * helmholtz(t, rho, &mb, &me, gamma)
}

/// `d(nF)/dn_i` at constant temperature and volume, by NeqSim's central difference:
/// perturb `n_i` by `n/10000`, renormalising the composition and the density through the
/// fixed total volume.
#[must_use]
pub fn d_f_dn(
    t: f64,
    rho: f64,
    x: &[f64],
    bps: &[[f64; 9]],
    bes: &[[f64; 6]],
    rhocs: &[f64],
) -> Vec<f64> {
    let n = x.iter().sum::<f64>();
    let dn = n / 10_000.0;
    (0..x.len())
        .map(|i| {
            let eval = |sign: f64| {
                let n_new = n + sign * dn;
                let x_new: Vec<f64> = x
                    .iter()
                    .enumerate()
                    .map(|(j, &xj)| (xj * n + if j == i { sign * dn } else { 0.0 }) / n_new)
                    .collect();
                let rho_new = rho * n_new / n;
                extensive_helmholtz(t, rho_new, n_new, &x_new, bps, bes, rhocs)
            };
            (eval(1.0) - eval(-1.0)) / (2.0 * dn)
        })
        .collect()
}

/// `ln phi_i = d(nF)/dn_i - ln Z`, with `Z = 1 + rho dF/drho`.
#[must_use]
pub fn ln_fugacity_coefficients(
    t: f64,
    rho: f64,
    x: &[f64],
    bps: &[[f64; 9]],
    bes: &[[f64; 6]],
    rhocs: &[f64],
) -> Vec<f64> {
    let (mb, me, gamma) = mix(x, bps, bes, rhocs);
    let ln_z = (1.0 + rho * d_helmholtz_drho(t, rho, &mb, &me, gamma)).ln();
    d_f_dn(t, rho, x, bps, bes, rhocs)
        .into_iter()
        .map(|d| d - ln_z)
        .collect()
}

/// `ln phi_i` for a mixture of [`BwrsCoefficients`], computing the temperature-dependent
/// coefficients once.
#[must_use]
pub fn ln_fugacity(t: f64, rho: f64, x: &[f64], comps: &[BwrsCoefficients]) -> Vec<f64> {
    let bps: Vec<[f64; 9]> = comps.iter().map(|c| bp(t, &c.a)).collect();
    let bes: Vec<[f64; 6]> = comps.iter().map(|c| be(t, &c.a)).collect();
    let rhocs: Vec<f64> = comps.iter().map(|c| c.rhoc).collect();
    ln_fugacity_coefficients(t, rho, x, &bps, &bes, &rhocs)
}
