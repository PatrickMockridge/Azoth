//! `eos.gt_surface_tension` - NeqSim's `GTSurfaceTensionFullGT`.
//!
//! Spec: `specs/models/eos/gt_surface_tension.toml`. The class solves the full gradient
//! theory on a finite-difference grid: a damped Newton on the Euler-Lagrange equations,
//! a block-tridiagonal Jacobian, deferred correction, and an outer remeshing loop that
//! doubles the grid until the integral stops moving.

use azoth_core::units::{Pressure, ThermodynamicTemperature, newtons_per_meter};
use azoth_core::{AzothError, Result, apply_checks};

use crate::model_gen;
use crate::results::GtSurfaceTensionResult;
use crate::tp_multiflash::tp_multiflash;
use crate::{Mixture, ReducedParameters};

/// NeqSim's gas constant.
const R: f64 = 8.314_462_1;
/// The class's initial level: `Ngrid = 2^3 + 1`.
const NLEVEL: usize = 3;
/// The class's refinement cap.
const NREFINEMENTS: usize = 7;
/// Newton iterations per refinement.
const NEWTON_STEPS: usize = 20;
/// The damping bound, `maxRelChange`.
const MAX_REL_CHANGE: f64 = 0.5;
/// The Newton stopping test, `urel_Newton * maxrelchange < 1e-5`.
const NEWTON_TOLERANCE: f64 = 1.0e-5;
/// The outer stopping test, `|sigma - sigma_old| < 0.01 * sigma`.
const OUTER_TOLERANCE: f64 = 0.01;
/// The domain half-width, in nm.
const DOMAIN_HALF_WIDTH: f64 = 3.0;
/// The influence matrix's nm scaling, `1e18`.
const NM_SCALE: f64 = 1.0e18;

/// The state at one density vector: chemical potentials and their physical-density
/// derivative.
struct ProfileState {
    mu: Vec<f64>,
    dmu: Vec<Vec<f64>>,
}

/// The full gradient theory's interface tension of a flashed fluid.
///
/// The class approximates the density profile on an equi-spaced grid, solves the
/// Euler-Lagrange equations by a damped Newton whose Jacobian is block tridiagonal, and
/// remeshes until the integrated tension stops changing. The residual is
/// `Delta mu(rho) - C * (rho[i-1] - 2 rho[i] + rho[i+1]) / h^2`, and the tension is
/// `integral drho/dz^T C drho/dz dz` in the grid's nanometre scale.
///
/// # Errors
/// * [`AzothError::InvalidInput`] for a single-phase feed, an absent or negative
///   component, a bulk state out of equilibrium, a singular block solve, or a profile
///   that does not stay finite.
/// * [`AzothError::OutOfRange`] if an input is outside the spec's range.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, pascals};
/// use azoth_eos::Cubic;
/// use azoth_eos::databank;
/// use azoth_eos::gt_surface_tension::gt_surface_tension;
///
/// let (mixture, _) = databank::mixture_of(&["methane", "n-butane"], Cubic::Srk, None)?;
/// let r = gt_surface_tension(
///     &mixture,
///     &[9.634812734820117e-21, 2.7571945977326497e-19],
///     kelvins(300.0),
///     pascals(2.0e6),
///     &[0.5, 0.5],
/// )?;
/// assert!((r.sigma.value - 0.00950471504464287).abs() < 1e-6);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn gt_surface_tension(
    mixture: &Mixture,
    influence_parameters: &[f64],
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<GtSurfaceTensionResult> {
    let spec = &model_gen::GT_SURFACE_TENSION_SPEC;
    let mut warnings = Vec::new();

    let count = mixture.len();
    if influence_parameters.len() != count || z.len() != count {
        return Err(AzothError::invalid_input(
            "influence_parameters",
            format!(
                "{count} component(s) against {} influence parameter(s) and {} mole \
                 fraction(s): every one is one entry per component",
                influence_parameters.len(),
                z.len()
            ),
        ));
    }

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "P" => Some(p.value),
            "influence_parameters" => influence_parameters.iter().copied().reduce(f64::min),
            "z" => z.iter().copied().reduce(f64::min),
            _ => None,
        },
        &mut warnings,
    )?;

    let flash = tp_multiflash(mixture, t, p, z)?;
    if flash.phase_count != 2 {
        let message = if flash.phase_count == 1 {
            "the feed is single phase at this state, so there is no interface to have a \
             tension; the full gradient theory solves a profile between two phases"
                .to_string()
        } else {
            format!(
                "the feed splits into {} phases at this state, and the class profiles its \
                 first two; this port refuses a state that is not two phases rather than \
                 dropping one",
                flash.phase_count
            )
        };
        return Err(AzothError::invalid_input("z", message));
    }
    if flash.z_factor[0] <= 0.0 || flash.z_factor[1] <= 0.0 {
        return Err(AzothError::invalid_input(
            "z",
            "a phase's compressibility root is not positive, so it has no density",
        ));
    }
    for fraction in flash.x[0].iter().chain(&flash.x[1]) {
        if *fraction <= 0.0 {
            return Err(AzothError::invalid_input(
                "z",
                "a component is absent from one of the phases. The class adds a trace amount \
                 and reflashes; this port refuses the state rather than perturbing the feed",
            ));
        }
    }

    let rt = R * t.value;
    let rho_ph1: Vec<f64> = flash.x[0]
        .iter()
        .map(|fraction| fraction * p.value / (flash.z_factor[0] * rt))
        .collect();
    let rho_ph2: Vec<f64> = flash.x[1]
        .iter()
        .map(|fraction| fraction * p.value / (flash.z_factor[1] * rt))
        .collect();
    let reduced = mixture.reduced_parameters(t, p)?;
    let a_abs: Vec<f64> = reduced
        .a
        .iter()
        .map(|value| value * rt * rt / p.value)
        .collect();
    let b_abs: Vec<f64> = reduced.b.iter().map(|value| value * rt / p.value).collect();
    let cubic = mixture.cubic();
    let (delta1, delta2) = (cubic.delta1(), cubic.delta2());

    let equilibrium_1 = state_at(
        mixture, &reduced, &a_abs, &b_abs, delta1, delta2, &rho_ph1, rt, p.value,
    )?;
    let equilibrium_2 = state_at(
        mixture, &reduced, &a_abs, &b_abs, delta1, delta2, &rho_ph2, rt, p.value,
    )?;
    // The class's `initmu` refuses a flash whose two phases do not carry the same
    // chemical potentials; the test is on the molar potentials, whose common constant
    // cancels in the difference.
    for i in 0..count {
        let difference = (equilibrium_1.mu[i] - equilibrium_2.mu[i]).abs();
        let scale = equilibrium_1.mu[i].abs().max(equilibrium_2.mu[i].abs());
        if difference > 1.0e-6 * scale.max(rt) {
            return Err(AzothError::invalid_input(
                "z",
                format!(
                    "the flash's two phases do not have the same chemical potential for \
                     component {i}: {:.6e} J/mol apart, and the class refuses a bulk state \
                     that is not solved",
                    equilibrium_1.mu[i] - equilibrium_2.mu[i]
                ),
            ));
        }
    }

    let cij: Vec<Vec<f64>> = (0..count)
        .map(|i| {
            (0..count)
                .map(|j| NM_SCALE * (influence_parameters[i] * influence_parameters[j]).sqrt())
                .collect()
        })
        .collect();

    let mut ngrid = (1usize << NLEVEL) + 1;
    let mut half_width = DOMAIN_HALF_WIDTH;
    let half = ngrid / 2;
    let mut rhomat: Vec<Vec<f64>> = (0..ngrid)
        .map(|index| {
            if index < half {
                rho_ph1.clone()
            } else {
                rho_ph2.clone()
            }
        })
        .collect();

    let mut sigma = 0.0;
    let mut sigma_old = 9.9e99;
    for refinement in 0..NREFINEMENTS {
        newton(
            &cij,
            half_width,
            &mut rhomat,
            mixture,
            &reduced,
            &a_abs,
            &b_abs,
            delta1,
            delta2,
            rt,
            p.value,
            &equilibrium_1.mu,
            count,
        )?;
        let h = 2.0 * half_width / (ngrid - 1) as f64;
        let drhodz = density_gradient(&rhomat, h);
        sigma = sigma_calc(h, &drhodz, &cij);
        let std = std_integral(h, &drhodz, &cij);
        if refinement > 0 && (sigma - sigma_old).abs() < OUTER_TOLERANCE * sigma {
            break;
        }
        sigma_old = sigma;

        let h = 2.0 * half_width / (ngrid - 1) as f64;
        let n_half = (6.0 * std / h).ceil() as usize;
        let n_ideal = 2 * n_half + 1;
        let new_half_width = n_half as f64 * h;
        let new_ngrid = 2 * n_ideal - 1;
        let old_grid = linspace(-half_width, half_width, ngrid);
        let new_grid = linspace(-new_half_width, new_half_width, new_ngrid);
        let mut next: Vec<Vec<f64>> = Vec::with_capacity(new_ngrid);
        let mut kk = 0usize;
        for &x in &new_grid {
            if x < old_grid[0] {
                next.push(rho_ph1.clone());
            } else if x > old_grid[ngrid - 1] {
                next.push(rho_ph2.clone());
            } else {
                while kk + 1 < ngrid && x > old_grid[kk + 1] {
                    kk += 1;
                }
                let alpha = (old_grid[kk + 1] - x) / (old_grid[kk + 1] - old_grid[kk]);
                let mut row = vec![0.0; count];
                for j in 0..count {
                    row[j] = alpha * rhomat[kk][j] + (1.0 - alpha) * rhomat[kk + 1][j];
                }
                next.push(row);
            }
        }
        if let Some(first) = next.first_mut() {
            *first = rho_ph1.clone();
        }
        if let Some(last) = next.last_mut() {
            *last = rho_ph2.clone();
        }
        half_width = new_half_width;
        ngrid = new_ngrid;
        rhomat = next;
    }

    Ok(GtSurfaceTensionResult {
        sigma: newtons_per_meter(sigma.abs()),
        warnings,
    })
}

/// One Newton solve for a grid, modifying the density matrix in place.
#[allow(clippy::too_many_arguments)]
fn newton(
    cij: &[Vec<f64>],
    half_width: f64,
    rhomat: &mut [Vec<f64>],
    mixture: &Mixture,
    reduced: &ReducedParameters,
    a_abs: &[f64],
    b_abs: &[f64],
    delta1: f64,
    delta2: f64,
    rt: f64,
    reference_pressure: f64,
    mu_equilibrium: &[f64],
    count: usize,
) -> Result<()> {
    let ngrid = rhomat.len();
    let h = 2.0 * half_width / (ngrid - 1) as f64;
    let mut jacobian = vec![vec![vec![0.0; count]; count]; ngrid];
    let mut dmu = vec![vec![0.0; count]; ngrid];
    let mut residual = vec![vec![0.0; count]; ngrid];
    for i in 0..ngrid {
        let state = state_at(
            mixture,
            reduced,
            a_abs,
            b_abs,
            delta1,
            delta2,
            &rhomat[i],
            rt,
            reference_pressure,
        )?;
        for j in 0..count {
            dmu[i][j] = state.mu[j] - mu_equilibrium[j];
            jacobian[i][j].copy_from_slice(&state.dmu[j]);
        }
        for j in 0..count {
            residual[i][j] = dmu[i][j];
            for k in 0..count {
                residual[i][j] -= jacobian[i][j][k] * rhomat[i][k];
            }
        }
    }
    add_deferred_correction(&mut residual, &dmu, ngrid, count);

    let xgrid = linspace(-half_width, half_width, ngrid);
    // The class's centre-of-gravity shift assigns row *references*, so the row it shifts
    // from and the boundary row become one array that the Newton update moves together.
    let mut alias: Option<(usize, usize)> = None;
    for _ in 0..NEWTON_STEPS {
        let previous: Vec<Vec<f64>> = rhomat.to_vec();
        direct_solve(&residual, &jacobian, cij, h, rhomat, count)?;
        follow_alias(rhomat, alias);

        let mut max_change = 0.0f64;
        for i in 1..ngrid - 1 {
            for j in 0..count {
                let change = (rhomat[i][j] - previous[i][j]).abs();
                if change > max_change * previous[i][j].abs() {
                    max_change = change / previous[i][j].abs();
                }
            }
        }
        let damping = if max_change > 0.0 {
            (MAX_REL_CHANGE / max_change).min(1.0)
        } else {
            1.0
        };
        for i in 1..ngrid - 1 {
            for j in 0..count {
                rhomat[i][j] = damping * rhomat[i][j] + (1.0 - damping) * previous[i][j];
            }
        }
        follow_alias(rhomat, alias);
        if damping * max_change < NEWTON_TOLERANCE {
            break;
        }

        // Re-centre the profile on its centre of gravity, as the class does.
        let drhodz = density_gradient(rhomat, h);
        let mut sum_z = 0.0;
        let mut sum = 0.0;
        for i in 0..drhodz.len() {
            let mut tmp = 0.0;
            for j in 0..count {
                for k in 0..count {
                    tmp += drhodz[i][j] * cij[j][k] * drhodz[i][k];
                }
            }
            sum_z += 0.5 * (xgrid[i + 1] + xgrid[i]) * tmp;
            sum += tmp;
        }
        if sum != 0.0 {
            // `Math.round` is `floor(x + 0.5)`, which differs from Rust's round at -0.5.
            let correction = (sum_z / sum / h + 0.5).floor() as isize;
            if correction > 0 {
                let shift = correction as usize;
                for i in 1..ngrid - shift {
                    rhomat[i] = rhomat[i + shift].clone();
                }
                for i in ngrid - shift..ngrid - 1 {
                    rhomat[i] = rhomat[ngrid - 1].clone();
                }
                alias = Some((ngrid - 1 - shift, ngrid - 1));
            } else if correction < 0 {
                let shift = (-correction) as usize;
                for i in (shift..ngrid - 1).rev() {
                    rhomat[i] = rhomat[i - shift].clone();
                }
                for i in 1..shift {
                    rhomat[i] = rhomat[0].clone();
                }
                alias = Some((shift, 0));
            }
        }

        for i in 1..ngrid - 1 {
            let state = state_at(
                mixture,
                reduced,
                a_abs,
                b_abs,
                delta1,
                delta2,
                &rhomat[i],
                rt,
                reference_pressure,
            )?;
            for j in 0..count {
                dmu[i][j] = state.mu[j] - mu_equilibrium[j];
                jacobian[i][j].copy_from_slice(&state.dmu[j]);
            }
            for j in 0..count {
                residual[i][j] = dmu[i][j];
                for k in 0..count {
                    residual[i][j] -= jacobian[i][j][k] * rhomat[i][k];
                }
            }
        }
        add_deferred_correction(&mut residual, &dmu, ngrid, count);
    }
    Ok(())
}

/// Carries an aliased row's value into the row it shares its array with.
fn follow_alias(rhomat: &mut [Vec<f64>], alias: Option<(usize, usize)>) {
    if let Some((source, shared)) = alias {
        let value = rhomat[source].clone();
        rhomat[shared] = value;
    }
}

/// `rres[i] += (dmu[i-1] - 2 dmu[i] + dmu[i+1]) / 12` for the inner grid points.
fn add_deferred_correction(
    residual: &mut [Vec<f64>],
    dmu: &[Vec<f64>],
    ngrid: usize,
    count: usize,
) {
    for i in 1..ngrid - 1 {
        for j in 0..count {
            residual[i][j] += (dmu[i - 1][j] - 2.0 * dmu[i][j] + dmu[i + 1][j]) / 12.0;
        }
    }
}

/// The block-tridiagonal Newton step: `A x_{i-1} + B_i x_i + C x_{i+1} = rhs_i`, solved
/// by the block Thomas algorithm and added to `rhomat` in place.
fn direct_solve(
    residual: &[Vec<f64>],
    jacobian: &[Vec<Vec<f64>>],
    cij: &[Vec<f64>],
    h: f64,
    rhomat: &mut [Vec<f64>],
    count: usize,
) -> Result<()> {
    let ngrid = rhomat.len();
    let h2 = h * h;
    let blocks = ngrid - 2;
    let mut rhs: Vec<Vec<f64>> = vec![vec![0.0; count]; blocks];
    let mut diagonal: Vec<Vec<Vec<f64>>> = vec![vec![vec![0.0; count]; count]; blocks];
    for (block, i) in (1..ngrid - 1).enumerate() {
        for j in 0..count {
            let mut value = residual[i][j];
            for k in 0..count {
                value += jacobian[i][j][k] * rhomat[i][k]
                    - cij[j][k] * (rhomat[i - 1][k] - 2.0 * rhomat[i][k] + rhomat[i + 1][k]) / h2;
            }
            rhs[block][j] = -value;
            for k in 0..count {
                diagonal[block][j][k] = jacobian[i][j][k] + 2.0 / h2 * cij[j][k];
            }
        }
    }
    let off_diagonal: Vec<Vec<f64>> = cij
        .iter()
        .map(|row| row.iter().map(|v| -v / h2).collect())
        .collect();

    // Forward elimination on the block tridiagonal system.
    let mut modified_diagonal = diagonal.clone();
    let mut modified_rhs = rhs.clone();
    for block in 1..blocks {
        let inverse = invert(&modified_diagonal[block - 1])?;
        let product = mat_mul(&off_diagonal, &inverse);
        let correction = mat_mul(&product, &off_diagonal);
        for j in 0..count {
            for k in 0..count {
                modified_diagonal[block][j][k] -= correction[j][k];
            }
        }
        let rhs_correction = mat_vec(&product, &modified_rhs[block - 1]);
        for j in 0..count {
            modified_rhs[block][j] -= rhs_correction[j];
        }
    }
    // Back substitution.
    let mut solution = vec![vec![0.0; count]; blocks];
    let inverse = invert(&modified_diagonal[blocks - 1])?;
    solution[blocks - 1] = mat_vec(&inverse, &modified_rhs[blocks - 1]);
    for block in (0..blocks - 1).rev() {
        let mut rhs_block = modified_rhs[block].clone();
        let sub = mat_vec(&off_diagonal, &solution[block + 1]);
        for j in 0..count {
            rhs_block[j] -= sub[j];
        }
        let inverse = invert(&modified_diagonal[block])?;
        solution[block] = mat_vec(&inverse, &rhs_block);
    }

    for (block, i) in (1..ngrid - 1).enumerate() {
        for j in 0..count {
            rhomat[i][j] += solution[block][j];
            if !rhomat[i][j].is_finite() {
                return Err(AzothError::invalid_input(
                    "z",
                    "the full gradient theory's Newton step left a non-finite density",
                ));
            }
        }
    }
    Ok(())
}

/// `drho/dz` on each interval.
fn density_gradient(rhomat: &[Vec<f64>], h: f64) -> Vec<Vec<f64>> {
    let count = rhomat[0].len();
    (0..rhomat.len() - 1)
        .map(|i| {
            (0..count)
                .map(|j| (rhomat[i + 1][j] - rhomat[i][j]) / h)
                .collect()
        })
        .collect()
}

/// `sigma = h * 1e-9 * sum_jk C_jk sum_i drho_i[j] drho_i[k]`.
fn sigma_calc(h: f64, drhodz: &[Vec<f64>], cij: &[Vec<f64>]) -> f64 {
    let count = cij.len();
    let mut sigma = 0.0;
    for j in 0..count {
        for k in 0..count {
            let mut sum = 0.0;
            for row in drhodz {
                sum += row[j] * row[k];
            }
            sigma += cij[j][k] * sum;
        }
    }
    sigma * h * 1.0e-9
}

/// The interface width: the second moment of the tension integrand.
///
/// **The class's accumulator is a single assignment, not a sum** -
/// `dum = drhodz[j][m] * C[m][n] * drhodz[j][n]` inside the `m`/`n` loops - so its width
/// is the last column's contribution. The port reproduces the assignment: it sets the
/// grid the outer loop re-solves on, and the answer carries that grid.
fn std_integral(h: f64, drhodz: &[Vec<f64>], cij: &[Vec<f64>]) -> f64 {
    let count = cij.len();
    let mut sum = 0.0;
    let mut z_sum = 0.0;
    let mut z2_sum = 0.0;
    for (index, row) in drhodz.iter().enumerate() {
        let z = index as f64 * h;
        let mut value = 0.0;
        for m in 0..count {
            for n in 0..count {
                // The class's own accumulator: an assignment, not a sum. Reproduced here so
                // the remesh width is its width; the docstring above the function records it.
                value = row[m] * cij[m][n] * row[n];
            }
        }
        sum += value * h;
        z_sum += value * ((z + h) * (z + h) - z * z) / 2.0;
        z2_sum += value * ((z + h).powi(3) - z.powi(3)) / 3.0;
    }
    if sum == 0.0 {
        return 0.0;
    }
    let mean = z_sum / sum;
    (z2_sum / sum - mean * mean).max(0.0).sqrt()
}

fn linspace(a: f64, b: f64, count: usize) -> Vec<f64> {
    let dx = (b - a) / (count - 1) as f64;
    (0..count).map(|index| a + index as f64 * dx).collect()
}

fn mat_mul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let count = a.len();
    (0..count)
        .map(|i| {
            (0..count)
                .map(|j| (0..count).map(|k| a[i][k] * b[k][j]).sum())
                .collect()
        })
        .collect()
}

fn mat_vec(a: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    a.iter()
        .map(|row| row.iter().zip(v).map(|(a, v)| a * v).sum())
        .collect()
}

/// The inverse of a small dense matrix, by Gaussian elimination with partial pivoting.
fn invert(matrix: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
    let count = matrix.len();
    let mut augmented: Vec<Vec<f64>> = (0..count)
        .map(|i| {
            let mut row = matrix[i].clone();
            row.extend((0..count).map(|j| if i == j { 1.0 } else { 0.0 }));
            row
        })
        .collect();
    for column in 0..count {
        let mut pivot = column;
        for row in column + 1..count {
            if augmented[row][column].abs() > augmented[pivot][column].abs() {
                pivot = row;
            }
        }
        if augmented[pivot][column].abs() < 1.0e-300 {
            return Err(AzothError::invalid_input(
                "z",
                "the full gradient theory's block is singular; the class's band solver \
                 would fail on it",
            ));
        }
        augmented.swap(column, pivot);
        let divisor = augmented[column][column];
        for value in &mut augmented[column] {
            *value /= divisor;
        }
        for row in 0..count {
            if row == column {
                continue;
            }
            let factor = augmented[row][column];
            if factor == 0.0 {
                continue;
            }
            let pivot_row = augmented[column].clone();
            for (value, pivot) in augmented[row].iter_mut().zip(&pivot_row) {
                *value -= factor * pivot;
            }
        }
    }
    Ok((0..count).map(|i| augmented[i][count..].to_vec()).collect())
}

/// The T,V state at a density vector: `mu_i = RT (ln(rho_i RT) + dA_R/dn_i)` and
/// `dmu_i/drho_k = RT (delta_ik / rho_i + d2A_R/dn_i dn_k / N)`.
#[allow(clippy::too_many_arguments)]
fn state_at(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    a_abs: &[f64],
    b_abs: &[f64],
    delta1: f64,
    delta2: f64,
    density: &[f64],
    rt: f64,
    reference_pressure: f64,
) -> Result<ProfileState> {
    let count = density.len();
    if density
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(AzothError::invalid_input(
            "z",
            "the full gradient theory's density vector is not a state: a component density \
             is non-positive or not finite",
        ));
    }
    let total: f64 = density.iter().sum();
    let x: Vec<f64> = density.iter().map(|value| value / total).collect();
    let volume = 1.0 / total;
    let mut b_mix = 0.0;
    let mut a_mix = 0.0;
    for i in 0..count {
        b_mix += x[i] * b_abs[i];
        for j in 0..count {
            let kij = reduced.kij[i * count + j];
            a_mix += x[i] * x[j] * (1.0 - kij) * (a_abs[i] * a_abs[j]).sqrt();
        }
    }
    if !volume.is_finite() || volume <= b_mix {
        return Err(AzothError::invalid_input(
            "z",
            "the full gradient theory's density vector is not a state: its molar volume is \
             at or below the mixture covolume",
        ));
    }
    let pressure =
        rt / (volume - b_mix) - a_mix / ((volume + delta1 * b_mix) * (volume + delta2 * b_mix));
    let compressibility = pressure * volume / rt;
    let ratio = pressure / reference_pressure;
    let mut reduced_step = reduced.clone();
    for value in &mut reduced_step.a {
        *value *= ratio;
    }
    for value in &mut reduced_step.b {
        *value *= ratio;
    }
    let gradient = mixture.helmholtz_gradient(&reduced_step, &x, compressibility)?;
    let hessian = mixture.helmholtz_hessian(&reduced_step, &x, compressibility)?;
    let mu: Vec<f64> = (0..count)
        .map(|i| rt * ((density[i] * rt).ln() + gradient[i]))
        .collect();
    let mut dmu = hessian;
    for i in 0..count {
        for value in &mut dmu[i] {
            *value *= rt / total;
        }
        dmu[i][i] += rt / density[i];
    }
    Ok(ProfileState { mu, dmu })
}
