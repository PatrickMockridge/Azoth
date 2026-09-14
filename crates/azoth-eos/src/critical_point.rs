//! `eos.critical_point` - the mixture critical point.
//!
//! Spec: `specs/models/eos/critical_point.yaml`
//!
//! # What a critical point is, and why it is not where the spinodal closes
//!
//! For a pure component the critical point is where the cubic's three roots merge, and
//! `dP/dV = d2P/dV2 = 0` finds it exactly. That route does **not** generalise, and the
//! way it fails is worth stating because it is silent: in reduced variables those two
//! conditions depend only on `(a_tilde, v)` and have a *single universal root*, so they
//! predict the same `Z_c = 0.3074` for every mixture - exact for one component, wrong
//! for all the rest, and wrong in a way its own name conceals. This model therefore
//! implements the two conditions that do generalise, from Heidemann & Khalil (1980).
//!
//! # The two conditions
//!
//! At a fixed composition, with `Q` the scaled Helmholtz Hessian of
//! [`crate::mixture::Mixture::criticality_matrix`]:
//!
//! 1. the **algebraically smallest** eigenvalue of `Q` vanishes;
//! 2. the **cubic form** - the third derivative of the Helmholtz energy along the
//!    eigenvector of that eigenvalue - vanishes.
//!
//! Both are evaluated at constant temperature and volume, so the state is `(T, V)` and
//! this model works in `(T, V)` throughout rather than in `(T, P)`. That is not a
//! detail: at a critical point the cubic is degenerate and its root is ill-conditioned,
//! and solving for the pressure *from the volume* through the explicit equation of
//! state avoids ever asking the cubic a question it cannot answer there.
//!
//! # Why the smallest eigenvalue and not `det(Q)`
//!
//! The determinant is the product of every eigenvalue, so it vanishes when *any* of
//! them does - including ones whose vanishing is not criticality - and being a product
//! it is badly scaled for a Newton step. The eigenvalue is the quantity whose vanishing
//! is the condition. Note also "algebraically smallest" rather than "smallest in
//! magnitude": the eigenvalue that crosses zero at a critical point is the one that
//! goes from negative to positive, and below the critical point it is not the one
//! nearest zero.
//!
//! # The cubic form, and where the `-1` comes from
//!
//! The cubic form is a central difference of the Hessian along the eigenvector, which
//! needs no third-derivative tensor:
//!
//! ```text
//! C(u) = sum_ijk u_i u_j u_k d3(A^R/RT)/dn_i dn_j dn_k  -  sum_i u_i**3 / n_i**2
//! ```
//!
//! **The second sum is not optional.** `A^R` is the *residual* Helmholtz energy, and
//! the ideal part's third derivative at constant volume is `delta_ijk / n_i**2`.
//! Omitting it leaves a constant offset that is exactly `1` at unit composition - so a
//! pure component's cubic form reads `1.000000000` at its critical point instead of
//! zero, and the outer Newton converges on that offset rather than on the critical
//! point.
//!
//! # The iteration
//!
//! Nested, as the published method describes: an inner Newton on `T` drives the
//! eigenvalue to zero at fixed `V`, and an outer Newton on `V` drives the cubic form to
//! zero. The eigenvector is re-derived at each state rather than carried, because the
//! outer derivative is taken with the eigenvector that belongs to each perturbed state.

use azoth_core::units::{ThermodynamicTemperature, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::mixture::{Mixture, ReducedParameters};
use crate::model_gen;
use crate::pr_molar_volume::MOLAR_GAS_CONSTANT;
use crate::results::CriticalPointResult;

/// The reference pressure the reduced parameters are un-reduced at.
///
/// `eos.pr_alpha_ab` returns `A_i = a_i P/(R^2 T^2)` and `B_i = b_i P/(R T)`, so the
/// dimensioned pair is a division by any pressure. A reference rather than a
/// re-derivation of the alpha function, which would be a second expression of it.
const REFERENCE_PRESSURE: f64 = 1.0e5;

/// The composition step the cubic form is a central difference over.
///
/// Fixed rather than adaptive, which puts a floor of about `1e-10` on the outer
/// residual: the difference quotient's error falls as the step squared and its
/// round-off rises as the step squared inverted, and `1e-4` is near the sum's minimum
/// for a third derivative in `f64`. See the spec's notes.
const CUBIC_FORM_STEP: f64 = 1.0e-4;

/// Guard for the Jacobi sweeps: an off-diagonal below this is treated as zero.
const EIGEN_TINY: f64 = 1.0e-300;

/// A symmetric matrix's eigenvalues and their unit eigenvectors.
pub struct Eigen {
    /// Eigenvalues, ascending.
    pub values: Vec<f64>,
    /// `vectors[i]` is the unit eigenvector belonging to `values[i]`.
    pub vectors: Vec<Vec<f64>>,
}

/// Eigen-decomposition of a symmetric matrix, by cyclic Jacobi rotations.
///
/// Jacobi rather than a tridiagonal reduction because the matrices here are small -
/// one row per component - and because it is unconditionally stable and needs no
/// pivoting: it converges for any symmetric input, and the rotations are orthogonal by
/// construction.
///
/// The eigenvalues come back **ascending**, so the criticality condition can name
/// `values[0]` rather than searching for it.
///
/// # Why not the characteristic polynomial
///
/// The cubic solver in this namespace finds roots of a degree-three polynomial in
/// closed form. Doing the same here would mean a polynomial of degree `N`, whose
/// coefficients are sums of minors - a computation whose conditioning collapses long
/// before the matrix does. The rotations touch only the matrix itself.
#[must_use]
pub fn symmetric_eigen(matrix: &[Vec<f64>], tolerance: f64, max_iterations: u32) -> Eigen {
    let n = matrix.len();
    if n == 0 {
        return Eigen {
            values: Vec::new(),
            vectors: Vec::new(),
        };
    }
    let mut a: Vec<Vec<f64>> = matrix.to_vec();
    // `vectors[k][i]` is component `k` of eigenvector `i`; the identity to start.
    let mut vectors: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();

    for _ in 0..max_iterations {
        let off: f64 = a
            .iter()
            .enumerate()
            .flat_map(|(p, row)| {
                row.iter()
                    .enumerate()
                    .filter(move |(q, _)| *q != p)
                    .map(|(_, value)| value * value)
            })
            .sum();
        if off.sqrt() <= tolerance {
            break;
        }
        for p in 0..n - 1 {
            for q in p + 1..n {
                if a[p][q] == 0.0 || a[p][q].abs() < EIGEN_TINY {
                    continue;
                }
                // The rotation that annihilates `a[p][q]`. `t` is the smaller root of
                // `t^2 + 2 theta t - 1 = 0`, which is the choice that keeps the
                // rotation angle under 45 degrees and the iteration stable.
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let sign = if theta >= 0.0 { 1.0 } else { -1.0 };
                let t = sign / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                // Columns `p` and `q` of every row.
                for row in a.iter_mut() {
                    let (akp, akq) = (row[p], row[q]);
                    row[p] = c * akp - s * akq;
                    row[q] = s * akp + c * akq;
                }
                // Rows `p` and `q`, which are two rows of the same matrix and so need
                // a split rather than two mutable borrows of the whole thing.
                let (above, below) = a.split_at_mut(q);
                let row_p = &mut above[p];
                let row_q = &mut below[0];
                for (apk, aqk) in row_p.iter_mut().zip(row_q.iter_mut()) {
                    let (old_p, old_q) = (*apk, *aqk);
                    *apk = c * old_p - s * old_q;
                    *aqk = s * old_p + c * old_q;
                }
                for row in vectors.iter_mut() {
                    let (vkp, vkq) = (row[p], row[q]);
                    row[p] = c * vkp - s * vkq;
                    row[q] = s * vkp + c * vkq;
                }
            }
        }
    }

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| a[i][i].partial_cmp(&a[j][j]).expect("no NaN in a rotation"));
    Eigen {
        values: order.iter().map(|&i| a[i][i]).collect(),
        vectors: order
            .iter()
            .map(|&i| (0..n).map(|k| vectors[k][i]).collect())
            .collect(),
    }
}

/// `(a_i, b_i)` in SI, un-reduced from `eos.pr_alpha_ab`.
///
/// The reduced pair is `A_i = a_i P/(R^2 T^2)` and `B_i = b_i P/(R T)`, so any pressure
/// inverts the reduction. Going through the registered calc rather than writing the
/// alpha function a second time is deliberate: there is one expression of `alpha`, and
/// this is not it.
///
/// # Errors
/// Propagates [`Mixture::reduced_parameters`]'s errors.
pub fn dimensional_parameters(mixture: &Mixture, temperature: f64) -> Result<(Vec<f64>, Vec<f64>)> {
    let reduced = mixture.reduced_parameters(
        ThermodynamicTemperature::new::<azoth_core::units::kelvin>(temperature),
        pascals(REFERENCE_PRESSURE),
    )?;
    let scale_a =
        MOLAR_GAS_CONSTANT * MOLAR_GAS_CONSTANT * temperature * temperature / REFERENCE_PRESSURE;
    let scale_b = MOLAR_GAS_CONSTANT * temperature / REFERENCE_PRESSURE;
    Ok((
        reduced.a.iter().map(|v| v * scale_a).collect(),
        reduced.b.iter().map(|v| v * scale_b).collect(),
    ))
}

/// The pressure of the explicit equation of state at `(T, V, z)`.
///
/// `P = R T/(V - b) - a/(V^2 + 2 b V - b^2)`, evaluated directly. No cubic is solved,
/// which is the point: at a critical point the cubic is degenerate, and this is where
/// that degeneracy lives.
///
/// # Errors
/// Propagates [`dimensional_parameters`]'s errors.
pub fn pressure_at_volume(
    mixture: &Mixture,
    temperature: f64,
    volume: f64,
    composition: &[f64],
) -> Result<f64> {
    let (a, b) = dimensional_parameters(mixture, temperature)?;
    let count = composition.len();
    let b_mix: f64 = (0..count).map(|i| composition[i] * b[i]).sum();
    let mut a_mix = 0.0;
    for i in 0..count {
        for j in 0..count {
            a_mix +=
                composition[i] * composition[j] * (1.0 - mixture.kij(i, j)) * (a[i] * a[j]).sqrt();
        }
    }
    let repulsion = MOLAR_GAS_CONSTANT * temperature / (volume - b_mix);
    let attraction = a_mix / (volume * volume + 2.0 * b_mix * volume - b_mix * b_mix);
    Ok(repulsion - attraction)
}

/// The third derivative of the total `A/RT` along `vector`.
///
/// A central difference of the Hessian along the direction, which is the identity
///
/// ```text
/// u^T H(n + s u) u - u^T H(n - s u) u  =  2 s * (third derivative along u) + O(s^3)
/// ```
///
/// evaluated at unit total moles, minus the ideal part's third derivative
/// `sum_i u_i**3 / n_i**2`. See the module documentation for why that term is there.
fn cubic_form(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    compressibility: f64,
    composition: &[f64],
    vector: &[f64],
) -> Result<f64> {
    let count = composition.len();
    let mut values = Vec::with_capacity(2);
    for sign in [1.0, -1.0] {
        let shifted: Vec<f64> = (0..count)
            .map(|i| composition[i] + sign * CUBIC_FORM_STEP * vector[i])
            .collect();
        let hessian = mixture.helmholtz_hessian(reduced, &shifted, compressibility)?;
        let mut total = 0.0;
        for i in 0..count {
            for j in 0..count {
                total += vector[i] * hessian[i][j] * vector[j];
            }
        }
        values.push(total);
    }
    let residual = (values[0] - values[1]) / (2.0 * CUBIC_FORM_STEP);
    let ideal: f64 = (0..count)
        .map(|i| vector[i].powi(3) / composition[i].powi(2))
        .sum();
    Ok(residual - ideal)
}

/// The smallest eigenpair of `Q` at `(T, V)`, and the state it sits in.
fn state(
    mixture: &Mixture,
    temperature: f64,
    volume: f64,
    composition: &[f64],
) -> Result<(f64, Vec<f64>, f64, ReducedParameters)> {
    let pressure = pressure_at_volume(mixture, temperature, volume, composition)?;
    let compressibility = pressure * volume / (MOLAR_GAS_CONSTANT * temperature);
    let reduced = mixture.reduced_parameters(
        ThermodynamicTemperature::new::<azoth_core::units::kelvin>(temperature),
        pascals(pressure),
    )?;
    let matrix = mixture.criticality_matrix(&reduced, composition, compressibility)?;
    let eigen = symmetric_eigen(&matrix, 1.0e-15, 100);
    Ok((
        eigen.values[0],
        eigen.vectors[0].clone(),
        compressibility,
        reduced,
    ))
}

/// The critical point of a mixture of composition `z`.
///
/// `z` is the composition of the single phase whose critical point is wanted, not a
/// feed being split - so unlike the flash's, this model does not decide what happens to
/// it. Whether that composition could exist at the returned state is a different
/// question and is not asked.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `z` is the wrong length, has a negative entry, or
///   does not sum to one. Renormalising it here would make a caller's error invisible
///   in every number downstream, so it is refused instead.
/// * [`AzothError::OutOfRange`] if the state the iteration converged on is not a state.
///   A diverged iteration can converge on a negative temperature or a volume inside the
///   co-volume, and those are refused rather than reported.
/// * [`AzothError::SolverNotConverged`] if the two conditions do not reach the
///   tolerance in the spec's iteration cap. A critical point that is nearly one is not
///   one.
pub fn critical_point(mixture: &Mixture, z: &[f64]) -> Result<CriticalPointResult> {
    let spec = &model_gen::CRITICAL_POINT_SPEC;
    let mut warnings = Vec::new();

    let algorithm = spec.algorithm.expect("a procedure");
    let inner = algorithm.inner.expect("this scheme nests");
    let residual_tolerance = algorithm.tolerance;
    let max_iterations = algorithm.max_iterations;
    let eigenvalue_tolerance = inner.tolerance;
    let inner_iterations = inner.max_iterations;

    check_composition(z, mixture.len())?;

    // Kay's rule for the temperature and four times the co-volume, which are the
    // starting values the method is normally stated with. `b` is a per-component
    // constant - `omega_b R Tc/Pc`, with no temperature in it - so evaluating it at the
    // starting temperature is a convenience rather than a dependence.
    let mut temperature: f64 = (0..mixture.len())
        .map(|i| z[i] * mixture.components()[i].tc.value)
        .sum();
    let (_, b) = dimensional_parameters(mixture, temperature)?;
    let mut volume = 4.0 * (0..mixture.len()).map(|i| z[i] * b[i]).sum::<f64>();

    let mut iterations: u32 = 0;
    let mut residual = f64::INFINITY;
    while iterations < max_iterations {
        iterations += 1;

        // Inner: Newton on T drives the smallest eigenvalue to zero at this volume. It
        // stops on the eigenvalue rather than on the step, so `inner.tolerance` is a
        // bound on the quantity being solved for and not on how far the last step went.
        for _ in 0..inner_iterations {
            let (eigenvalue, _, _, _) = state(mixture, temperature, volume, z)?;
            if eigenvalue.abs() < eigenvalue_tolerance {
                break;
            }
            let probe = temperature * 1.0e-6;
            let (ahead, _, _, _) = state(mixture, temperature + probe, volume, z)?;
            let derivative = (ahead - eigenvalue) / probe;
            if derivative == 0.0 || !derivative.is_finite() {
                break;
            }
            // The step is capped at a tenth of the temperature, so a Newton step taken
            // from a poor starting point cannot leave the region where the state is
            // defined - where the equation of state has no critical point to find.
            let step = (-eigenvalue / derivative).clamp(-0.1 * temperature, 0.1 * temperature);
            temperature += step;
        }

        let (eigenvalue, eigenvector, compressibility, reduced) =
            state(mixture, temperature, volume, z)?;
        let form = cubic_form(mixture, &reduced, compressibility, z, &eigenvector)?;
        residual = eigenvalue.abs().max(form.abs());
        if residual < residual_tolerance {
            break;
        }

        let probe = volume * 1.0e-6;
        let (_, ahead_vector, ahead_compressibility, ahead_reduced) =
            state(mixture, temperature, volume + probe, z)?;

        let ahead_form = cubic_form(
            mixture,
            &ahead_reduced,
            ahead_compressibility,
            z,
            &ahead_vector,
        )?;
        let derivative = (ahead_form - form) / probe;
        if derivative == 0.0 || !derivative.is_finite() {
            break;
        }
        let volume_step = (-form / derivative).clamp(-0.1 * volume, 0.1 * volume);
        // Damped by a half: the cubic form is a finite difference of a Hessian, and a
        // full step overshoots on the first passes from the coarse starting volume.
        volume += 0.5 * volume_step;
    }

    if residual >= residual_tolerance {
        return Err(AzothError::SolverNotConverged {
            iterations,
            residual,
            tolerance: residual_tolerance,
        });
    }

    let pressure = pressure_at_volume(mixture, temperature, volume, z)?;
    let compressibility = pressure * volume / (MOLAR_GAS_CONSTANT * temperature);

    // The bounds are on the *result*, which is unusual for this tree: every input is a
    // vector or a matrix, so there is no scalar argument whose range a caller could
    // violate. What can go wrong is the iteration, and these are how a diverged one
    // announces itself.
    apply_checks(
        spec.derived_checks(),
        |quantity| match quantity {
            "tc" => Some(temperature),
            "pc" => Some(pressure),
            "z_c" => Some(compressibility),
            _ => None,
        },
        &mut warnings,
    )?;

    for (field, value) in [
        ("tc", temperature),
        ("pc", pressure),
        ("z_c", compressibility),
    ] {
        // Written as an explicit finiteness-and-positivity test rather than
        // `!(value > 0.0)`, which reads as a double negative and is the shape clippy
        // flags. Both reject NaN; only the explicit form also rejects an infinity.
        if !value.is_finite() || value <= 0.0 {
            return Err(AzothError::OutOfRange {
                field: field.to_string(),
                value,
                detail: "the critical point this composition was found to have is not a \
                         state; the iteration left the region where the equation of \
                         state is defined"
                    .to_string(),
            });
        }
    }

    Ok(CriticalPointResult {
        tc: ThermodynamicTemperature::new::<azoth_core::units::kelvin>(temperature),
        pc: pascals(pressure),
        vc: azoth_core::units::cubic_meters_per_mole(volume),
        z_c: compressibility,
        iterations,
        residual,
        warnings,
    })
}

/// A composition must be one entry per component, in `[0, 1]`, summing to one.
fn check_composition(values: &[f64], n: usize) -> Result<()> {
    if values.len() != n {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "a mixture of {n} components needs {n} mole fractions, but z has {}",
                values.len()
            ),
        ));
    }
    if let Some(bad) = values.iter().position(|&value| value < 0.0) {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "z[{bad}] is {} but a mole fraction cannot be negative",
                values[bad]
            ),
        ));
    }
    let sum: f64 = values.iter().sum();
    if (sum - 1.0).abs() > 1.0e-09 {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "the composition sums to {sum}, not to one. Renormalising it here would \
                 make a caller's error invisible in every number downstream, so it is \
                 refused instead"
            ),
        ));
    }
    Ok(())
}
