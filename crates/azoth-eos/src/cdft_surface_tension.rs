//! `eos.cdft_surface_tension` - NeqSim's `CDFTSurfaceTension`.
//!
//! Spec: `specs/models/eos/cdft_surface_tension.toml`. A tanh density profile on a
//! 2048-point grid, convolved with a step kernel whose width is one molecular diameter;
//! the interface width is chosen by a fixed scan and golden section, and the tension is
//! the integral of the excess grand potential.

use azoth_core::units::{ThermodynamicTemperature, newtons_per_meter, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::model_gen;
use crate::results::CdftSurfaceTensionResult;
use crate::{Cubic, Mixture};

/// NeqSim's gas constant.
const R: f64 = 8.314_462_1;
/// NeqSim's `avagadroNumber`, which is the older `6.023e23`.
const NA: f64 = 6.023e23;
/// The class's grid.
const NGRID: usize = 2048;
/// The domain half-width, in molecular diameters.
const DOMAIN_HALF_WIDTH_IN_D: f64 = 15.0;
/// The kernel half-width, in molecular diameters.
const ATTRACTIVE_RANGE: f64 = 0.5;
/// Scan points in the pure branch.
const SCAN_PURE: usize = 80;
/// Golden-section steps in the pure branch.
const GOLDEN_PURE: usize = 60;
/// Scan points in the mixture branch.
const SCAN_MIXTURE: usize = 60;
/// Golden-section steps in the mixture branch.
const GOLDEN_MIXTURE: usize = 50;

/// The classical density functional theory's interface tension between two states.
///
/// The class consumes a flashed system and reads only each phase's molar volume and mole
/// fractions, so this takes those four quantities: one component runs the pure branch,
/// whose potentials are analytic; more run the mixture branch, whose equilibrium
/// potentials are a central difference of the mixture's free-energy density and whose
/// attraction is `sqrt(a_i a_j)` with no interaction term.
///
/// # Errors
/// * [`AzothError::InvalidInput`] for a composition of the wrong length or one that does
///   not sum to one, a non-positive molar volume or temperature, or a variational minimum
///   that is not positive.
/// * [`AzothError::OutOfRange`] if an input is outside the spec's range.
///
/// # Example
/// ```
/// use azoth_core::units::kelvins;
/// use azoth_eos::Cubic;
/// use azoth_eos::cdft_surface_tension::cdft_surface_tension;
/// use azoth_eos::databank;
///
/// let (mixture, _) = databank::mixture_of(&["methane"], Cubic::Pr, None)?;
/// let r = cdft_surface_tension(
///     &mixture,
///     kelvins(150.0),
///     &[1.0],
///     &[1.0000000011004573],
///     4.128785631897625e-5,
///     9.709562218655069e-4,
/// )?;
/// assert!((r.sigma.value - 0.002875574857984931).abs() < 1e-10);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn cdft_surface_tension(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    x_liquid: &[f64],
    x_vapour: &[f64],
    molar_volume_liquid: f64,
    molar_volume_vapour: f64,
) -> Result<CdftSurfaceTensionResult> {
    let spec = &model_gen::CDFT_SURFACE_TENSION_SPEC;
    let mut warnings = Vec::new();

    let count = mixture.len();
    if x_liquid.len() != count || x_vapour.len() != count {
        return Err(AzothError::invalid_input(
            "x_liquid",
            format!(
                "{count} component(s) against {} liquid and {} vapour mole fraction(s): each \
                 is one entry per component",
                x_liquid.len(),
                x_vapour.len()
            ),
        ));
    }
    for (name, composition) in [("x_liquid", x_liquid), ("x_vapour", x_vapour)] {
        let total: f64 = composition.iter().sum();
        if (total - 1.0).abs() > 1.0e-6 {
            return Err(AzothError::invalid_input(
                name,
                format!("the mole fractions sum to {total}, not one"),
            ));
        }
    }
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "x_liquid" => x_liquid.iter().copied().reduce(f64::min),
            "x_vapour" => x_vapour.iter().copied().reduce(f64::min),
            "molar_volume_liquid" => Some(molar_volume_liquid),
            "molar_volume_vapour" => Some(molar_volume_vapour),
            _ => None,
        },
        &mut warnings,
    )?;

    // The cubic's a_i(T) and b_i are pressure-independent; the reduced parameters are read
    // at one bar and converted back, so the model needs no pressure.
    let reference = pascals(1.0e5);
    let reduced = mixture.reduced_parameters(t, reference)?;
    let rt = R * t.value;
    let a: Vec<f64> = reduced
        .a
        .iter()
        .map(|value| value * rt * rt / reference.value)
        .collect();
    let b: Vec<f64> = reduced
        .b
        .iter()
        .map(|value| value * rt / reference.value)
        .collect();
    // The class picks its constants by system class name: SRK takes (1, 0), everything
    // else the Peng-Robinson pair.
    let (delta1, delta2) = if mixture.cubic() == Cubic::Srk {
        (1.0, 0.0)
    } else {
        (1.0 + 2.0f64.sqrt(), 1.0 - 2.0f64.sqrt())
    };

    // The component densities the class derives: `x_i` over the phase's molar volume.
    let rho_liquid: Vec<f64> = x_liquid
        .iter()
        .map(|value| value / molar_volume_liquid)
        .collect();
    let rho_vapour: Vec<f64> = x_vapour
        .iter()
        .map(|value| value / molar_volume_vapour)
        .collect();

    let sigma = if count == 1 {
        let (rho_l, rho_v) = ordered_total(rho_liquid[0], rho_vapour[0]);
        pure_sigma(rho_l, rho_v, t.value, a[0], b[0], delta1, delta2)
    } else {
        let (liquid, vapour) = ordered_vectors(&rho_liquid, &rho_vapour);
        mixture_sigma(liquid, vapour, t.value, &a, &b, delta1, delta2)
    }?;

    if sigma <= 0.0 {
        return Err(AzothError::invalid_input(
            "molar_volume_liquid",
            format!(
                "the variational search's minimum is {sigma:.6e} N/m, so the excess grand \
                 potential never rises above the bulk; the class clamps the answer to 0.0, \
                 and this port refuses a state it cannot give a tension to"
            ),
        ));
    }
    Ok(CdftSurfaceTensionResult {
        sigma: newtons_per_meter(sigma),
        warnings,
    })
}

/// The two totals, liquid first: the class swaps them if the "liquid" is the less dense.
fn ordered_total(liquid: f64, vapour: f64) -> (f64, f64) {
    if liquid < vapour {
        (vapour, liquid)
    } else {
        (liquid, vapour)
    }
}

fn ordered_vectors<'a>(liquid: &'a [f64], vapour: &'a [f64]) -> (&'a [f64], &'a [f64]) {
    let total = |values: &[f64]| values.iter().sum::<f64>();
    if total(liquid) < total(vapour) {
        (vapour, liquid)
    } else {
        (liquid, vapour)
    }
}

/// The variational search of the pure branch.
#[allow(clippy::too_many_arguments)]
fn pure_sigma(
    rho_l: f64,
    rho_v: f64,
    t: f64,
    a: f64,
    b: f64,
    delta1: f64,
    delta2: f64,
) -> Result<f64> {
    let d_mol = (6.0 * b / (std::f64::consts::PI * NA)).cbrt();
    let mu_eq = R * t * rho_l.ln()
        + mu_repulsive_excess(rho_l, t, b)
        + mu_attractive(rho_l, a, b, delta1, delta2);
    let p_bulk = pressure(rho_l, t, a, b, delta1, delta2);
    let width = |delta: f64| {
        pure_sigma_for_width(
            rho_l, rho_v, t, a, b, delta1, delta2, d_mol, delta, mu_eq, p_bulk,
        )
    };
    Ok(variational_minimum(width, d_mol, SCAN_PURE, GOLDEN_PURE))
}

/// The variational search of the mixture branch.
fn mixture_sigma(
    rho_liq: &[f64],
    rho_vap: &[f64],
    t: f64,
    a: &[f64],
    b: &[f64],
    delta1: f64,
    delta2: f64,
) -> Result<f64> {
    let count = rho_liq.len();
    let aij: Vec<Vec<f64>> = (0..count)
        .map(|i| (0..count).map(|j| (a[i] * a[j]).sqrt()).collect())
        .collect();
    let total: f64 = rho_liq.iter().sum();
    let b_average: f64 = (0..count).map(|i| rho_liq[i] / total * b[i]).sum();
    let d_mol = (6.0 * b_average / (std::f64::consts::PI * NA)).cbrt();

    let eps = total * 1.0e-7;
    let mut mu_eq = Vec::with_capacity(count);
    for i in 0..count {
        let mut plus = rho_liq.to_vec();
        let mut minus = rho_liq.to_vec();
        plus[i] += eps;
        minus[i] -= eps;
        let f_plus = mixture_free_energy_density(&plus, &aij, b, t, delta1, delta2);
        let f_minus = mixture_free_energy_density(&minus, &aij, b, t, delta1, delta2);
        mu_eq.push((f_plus - f_minus) / (2.0 * eps));
    }
    let p_bulk = mixture_pressure(rho_liq, &aij, b, t, delta1, delta2);

    let width = |delta: f64| {
        mixture_sigma_for_width(
            rho_liq, rho_vap, t, &aij, b, delta1, delta2, d_mol, delta, &mu_eq, p_bulk,
        )
    };
    Ok(variational_minimum(
        width,
        d_mol,
        SCAN_MIXTURE,
        GOLDEN_MIXTURE,
    ))
}

/// The class's scan followed by a golden section, both over `[0.3, 30]` diameters.
fn variational_minimum(
    width: impl Fn(f64) -> f64,
    d_mol: f64,
    scan_points: usize,
    golden_steps: usize,
) -> f64 {
    let (lo, hi) = (0.3 * d_mol, 30.0 * d_mol);
    let mut best_delta = d_mol;
    let mut best = f64::INFINITY;
    for k in 0..scan_points {
        let delta = lo + (hi - lo) * k as f64 / (scan_points - 1) as f64;
        let sigma = width(delta);
        if sigma < best {
            best = sigma;
            best_delta = delta;
        }
    }
    let mut ref_lo = lo.max(best_delta - 3.0 * d_mol);
    let mut ref_hi = hi.min(best_delta + 3.0 * d_mol);
    let ratio = (5.0f64.sqrt() + 1.0) / 2.0;
    for _ in 0..golden_steps {
        let c = ref_hi - (ref_hi - ref_lo) / ratio;
        let d = ref_lo + (ref_hi - ref_lo) / ratio;
        if width(c) < width(d) {
            ref_hi = d;
        } else {
            ref_lo = c;
        }
    }
    width(0.5 * (ref_lo + ref_hi))
}

/// One width's tension in the pure branch.
#[allow(clippy::too_many_arguments)]
fn pure_sigma_for_width(
    rho_l: f64,
    rho_v: f64,
    t: f64,
    a: f64,
    b: f64,
    delta1: f64,
    delta2: f64,
    d_mol: f64,
    delta: f64,
    mu_eq: f64,
    p_bulk: f64,
) -> f64 {
    let half = DOMAIN_HALF_WIDTH_IN_D * d_mol;
    let dz = 2.0 * half / (NGRID - 1) as f64;
    let n_half = (ATTRACTIVE_RANGE * d_mol / dz).ceil() as isize;
    let n_points = (2 * n_half + 1) as f64;
    let w = -2.0 * a / (n_points * dz);

    let average = 0.5 * (rho_l + rho_v);
    let difference = rho_l - rho_v;
    let rho: Vec<f64> = (0..NGRID)
        .map(|i| {
            let z = -half + i as f64 * dz;
            average + 0.5 * difference * (-z / delta).tanh()
        })
        .collect();
    let convolution = step_convolve(&rho, n_half, w, dz);

    let mut sigma = 0.0;
    for i in 0..NGRID {
        let mut ri = rho[i].max(1.0e-10);
        if b * ri >= 0.999 {
            ri = 0.999 / b;
        }
        let f_repulsive = ri * R * t * (ri.ln() - 1.0 - (1.0 - b * ri).ln());
        let f_attractive_nonlocal = 0.5 * ri * convolution[i];
        let f_correction = f_attractive_local(ri, a, b, delta1, delta2) + a * ri * ri;
        sigma += (f_repulsive + f_attractive_nonlocal + f_correction - mu_eq * ri + p_bulk) * dz;
    }
    sigma
}

/// One width's tension in the mixture branch.
#[allow(clippy::too_many_arguments)]
fn mixture_sigma_for_width(
    rho_liq: &[f64],
    rho_vap: &[f64],
    t: f64,
    aij: &[Vec<f64>],
    b: &[f64],
    delta1: f64,
    delta2: f64,
    d_mol: f64,
    delta: f64,
    mu_eq: &[f64],
    p_bulk: f64,
) -> f64 {
    let count = rho_liq.len();
    let half = DOMAIN_HALF_WIDTH_IN_D * d_mol;
    let dz = 2.0 * half / (NGRID - 1) as f64;
    let n_half = (ATTRACTIVE_RANGE * d_mol / dz).ceil() as isize;
    let n_points = (2 * n_half + 1) as f64;

    let mut profile = vec![vec![0.0; NGRID]; count];
    let mut total = vec![0.0; NGRID];
    for ic in 0..count {
        let average = 0.5 * (rho_liq[ic] + rho_vap[ic]);
        let difference = rho_liq[ic] - rho_vap[ic];
        for i in 0..NGRID {
            let z = -half + i as f64 * dz;
            let mut value = average + 0.5 * difference * (-z / delta).tanh();
            if value < 1.0e-30 {
                value = 1.0e-30;
            }
            profile[ic][i] = value;
            total[i] += value;
        }
    }
    let mut averaged = vec![vec![0.0; NGRID]; count];
    for (jc, row) in profile.iter().enumerate() {
        for (i, value) in averaged[jc].iter_mut().enumerate() {
            let mut sum = 0.0;
            for k in -n_half..=n_half {
                sum += row[clamp_index(i as isize + k)];
            }
            *value = sum / n_points;
        }
    }

    let mut sigma = 0.0;
    for i in 0..NGRID {
        let mut rho_total = total[i].max(1.0e-10);
        let mut b_mix = 0.0;
        let mut a_mix = 0.0;
        for ic in 0..count {
            let xi = profile[ic][i] / rho_total;
            b_mix += xi * b[ic];
            for jc in 0..count {
                let xj = profile[jc][i] / rho_total;
                a_mix += xi * xj * aij[ic][jc];
            }
        }
        if b_mix * rho_total >= 0.999 {
            rho_total = 0.999 / b_mix;
        }
        let mut f_ideal = 0.0;
        for row in &profile {
            let ri = row[i];
            if ri > 1.0e-30 {
                f_ideal += ri * R * t * (ri.ln() - 1.0);
            }
        }
        let f_repulsive = -rho_total * R * t * (1.0 - b_mix * rho_total).ln();
        let mut f_attractive_nonlocal = 0.0;
        for ic in 0..count {
            for jc in 0..count {
                f_attractive_nonlocal -= aij[ic][jc] * profile[ic][i] * averaged[jc][i];
            }
        }
        let f_correction = f_attractive_local(rho_total, a_mix, b_mix, delta1, delta2)
            + a_mix * rho_total * rho_total;
        let mut mu_rho = 0.0;
        for ic in 0..count {
            mu_rho += mu_eq[ic] * profile[ic][i];
        }
        sigma +=
            (f_ideal + f_repulsive + f_attractive_nonlocal + f_correction - mu_rho + p_bulk) * dz;
    }
    sigma
}

/// A grid index clamped to the array, as the class clamps it at both ends.
fn clamp_index(index: isize) -> usize {
    index.clamp(0, NGRID as isize - 1) as usize
}

/// The class's step convolution, with its index clamping at both ends.
fn step_convolve(rho: &[f64], n_half: isize, w: f64, dz: f64) -> Vec<f64> {
    rho.iter()
        .enumerate()
        .map(|(i, _)| {
            let mut sum = 0.0;
            for k in -n_half..=n_half {
                sum += rho[clamp_index(i as isize + k)];
            }
            w * sum * dz
        })
        .collect()
}

/// `RT (-ln(1 - b rho) + b rho/(1 - b rho))`.
fn mu_repulsive_excess(rho: f64, t: f64, b: f64) -> f64 {
    let b_rho = (b * rho).min(0.999);
    R * t * (-(1.0 - b_rho).ln() + b_rho / (1.0 - b_rho))
}

/// The attractive chemical potential of a cubic EOS at a molar density.
fn mu_attractive(rho: f64, a: f64, b: f64, delta1: f64, delta2: f64) -> f64 {
    let u = b * rho;
    let denom1 = 1.0 + delta1 * u;
    let denom2 = 1.0 + delta2 * u;
    a / ((delta1 - delta2) * b) * (denom2 / denom1).ln() - a * rho / (denom1 * denom2)
}

/// The attractive Helmholtz free-energy density of a cubic EOS.
fn f_attractive_local(rho: f64, a: f64, b: f64, delta1: f64, delta2: f64) -> f64 {
    let u = b * rho;
    rho * a / ((delta1 - delta2) * b) * ((1.0 + delta2 * u) / (1.0 + delta1 * u)).ln()
}

/// The cubic EOS pressure at a molar density.
fn pressure(rho: f64, t: f64, a: f64, b: f64, delta1: f64, delta2: f64) -> f64 {
    let b_rho = b * rho;
    let p_repulsive = R * t * rho / (1.0 - b_rho);
    let denominator = (1.0 + delta1 * b_rho) * (1.0 + delta2 * b_rho);
    p_repulsive - a * rho * rho / denominator
}

/// The mixture's Helmholtz free-energy density at component densities.
fn mixture_free_energy_density(
    rho: &[f64],
    aij: &[Vec<f64>],
    b: &[f64],
    t: f64,
    delta1: f64,
    delta2: f64,
) -> f64 {
    let count = rho.len();
    let mut rho_total: f64 = rho.iter().sum();
    if rho_total < 1.0e-20 {
        return 0.0;
    }
    let mut b_mix = 0.0;
    let mut a_mix = 0.0;
    for i in 0..count {
        let xi = rho[i] / rho_total;
        b_mix += xi * b[i];
        for j in 0..count {
            let xj = rho[j] / rho_total;
            a_mix += xi * xj * aij[i][j];
        }
    }
    if b_mix * rho_total >= 0.999 {
        rho_total = 0.999 / b_mix;
    }
    let mut f_ideal = 0.0;
    for &value in rho {
        if value > 1.0e-30 {
            f_ideal += value * R * t * (value.ln() - 1.0);
        }
    }
    let f_repulsive = -rho_total * R * t * (1.0 - b_mix * rho_total).ln();
    let f_attractive = f_attractive_local(rho_total, a_mix, b_mix, delta1, delta2);
    f_ideal + f_repulsive + f_attractive
}

/// The mixture's EOS pressure at component densities.
fn mixture_pressure(
    rho: &[f64],
    aij: &[Vec<f64>],
    b: &[f64],
    t: f64,
    delta1: f64,
    delta2: f64,
) -> f64 {
    let count = rho.len();
    let rho_total: f64 = rho.iter().sum();
    let mut b_mix = 0.0;
    let mut a_mix = 0.0;
    for i in 0..count {
        let xi = rho[i] / rho_total;
        b_mix += xi * b[i];
        for j in 0..count {
            let xj = rho[j] / rho_total;
            a_mix += xi * xj * aij[i][j];
        }
    }
    pressure(rho_total, t, a_mix, b_mix, delta1, delta2)
}
