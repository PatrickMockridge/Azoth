//! `eos.gt_surface_tension_simple` - NeqSim's `GTSurfaceTensionSimple`.
//!
//! Spec: `specs/models/eos/gt_surface_tension_simple.toml`. The class walks a linear
//! density path like `LGTSurfaceTension`, but at each step it solves for the composition
//! gradients that equalise the chemical potentials - an inner Newton on the `dmu/dn`
//! matrix, which is where this solver is dearer than LGT.

use azoth_core::units::{Pressure, ThermodynamicTemperature, newtons_per_meter};
use azoth_core::{AzothError, Result, Warning, WarningCode, apply_checks};

use crate::model_gen;
use crate::pt_flash::pt_flash;
use crate::results::{GtSurfaceTensionSimpleResult, Phase};
use crate::{Mixture, ReducedParameters};

/// NeqSim's gas constant.
const R: f64 = 8.314_462_1;
/// The class's fixed number of steps along the reference component's density.
const STEPS: usize = 200;
/// The class's inner-loop cap, `iterations < 1200`.
const INNER_STEPS: usize = 1200;
/// The class's inner-loop tolerance, `err > 1e-15`.
const INNER_TOLERANCE: f64 = 1.0e-15;

/// The state at one density vector: pressure, chemical potentials and `dmu/dn` at fixed T,V.
struct ProfileState {
    pressure: f64,
    mu: Vec<f64>,
    dmu: Vec<Vec<f64>>,
}

/// The GT-simple interface tension of a flashed fluid.
///
/// At each of [`STEPS`] density steps the class builds the linear system
/// `sqrt(c_i) dmu_ref - sqrt(c_ref) dmu_i` over the non-reference components, solves it
/// for the density gradients, and iterates until the chemical-potential residuals are
/// below [`INNER_TOLERANCE`]. The accumulation is the same `sqrt(2 kappa mu_times_rho)`
/// integral LGT uses. **This class's interaction rule assigns zero in both arms of its
/// water condition**, so a water pair contributes no interaction term; the caller supplies
/// the matrix and the class's own answer for water is recorded rather than assumed.
///
/// # Errors
/// * [`AzothError::InvalidInput`] for a single-phase feed, an absent component, a
///   singular inner system, or a non-positive integrand.
/// * [`AzothError::OutOfRange`] if an input is outside the spec's range.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, pascals};
/// use azoth_eos::Cubic;
/// use azoth_eos::databank;
/// use azoth_eos::gt_surface_tension_simple::gt_surface_tension_simple;
///
/// let (mixture, _) = databank::mixture_of(&["methane", "n-butane"], Cubic::Srk, None)?;
/// let r = gt_surface_tension_simple(
///     &mixture,
///     &[9.634812734820117e-21, 2.7571945977326497e-19],
///     1.0,
///     &[vec![0.0, 0.0], vec![0.0, 0.0]],
///     kelvins(300.0),
///     pascals(2.0e6),
///     &[0.5, 0.5],
/// )?;
/// assert!((r.sigma.value - 0.00951973185903182).abs() < 1e-6);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn gt_surface_tension_simple(
    mixture: &Mixture,
    influence_parameters: &[f64],
    reference_component: f64,
    pair_interactions: &[Vec<f64>],
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<GtSurfaceTensionSimpleResult> {
    let spec = &model_gen::GT_SURFACE_TENSION_SIMPLE_SPEC;
    let mut warnings = Vec::new();

    let count = mixture.len();
    if influence_parameters.len() != count
        || pair_interactions.len() != count
        || pair_interactions.iter().any(|row| row.len() != count)
        || z.len() != count
    {
        return Err(AzothError::invalid_input(
            "influence_parameters",
            format!(
                "{count} component(s) against {} influence parameter(s), {} interaction \
                 row(s) and {} mole fraction(s): every one is one entry per component",
                influence_parameters.len(),
                pair_interactions.len(),
                z.len()
            ),
        ));
    }
    if !reference_component.is_finite()
        || reference_component < 0.0
        || reference_component.fract() != 0.0
        || reference_component as usize >= count
    {
        return Err(AzothError::invalid_input(
            "reference_component",
            format!(
                "`{reference_component}` is not a whole index in 0..{count}: it is the \
                 component the profile integrates along"
            ),
        ));
    }
    let reference = reference_component as usize;

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

    let flash = pt_flash(mixture, t, p, z)?;
    if flash.phase != Phase::TwoPhase {
        return Err(AzothError::invalid_input(
            "z",
            "the feed is single phase at this state, so there is no interface to have a \
             tension; a gradient-theory profile runs between two phases",
        ));
    }
    if flash.z_vapour <= 0.0 || flash.z_liquid <= 0.0 {
        return Err(AzothError::invalid_input(
            "z",
            "a phase's compressibility root is not positive, so it has no density",
        ));
    }
    for fraction in flash.y.iter().chain(&flash.x) {
        if *fraction <= 0.0 {
            return Err(AzothError::invalid_input(
                "z",
                "a component is absent from one of the phases, so its chemical potential \
                 would be `ln(0)`; this port refuses the partial integral the class would \
                 return at the resulting NaN",
            ));
        }
    }

    let rt = R * t.value;
    let rho_vapour: Vec<f64> = flash
        .y
        .iter()
        .map(|fraction| fraction * p.value / (flash.z_vapour * rt))
        .collect();
    let rho_liquid: Vec<f64> = flash
        .x
        .iter()
        .map(|fraction| fraction * p.value / (flash.z_liquid * rt))
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
    let sqrt_influence: Vec<f64> = influence_parameters
        .iter()
        .map(|value| value.sqrt())
        .collect();

    let equilibria = state_at(
        mixture,
        &reduced,
        &a_abs,
        &b_abs,
        delta1,
        delta2,
        &rho_vapour,
        rt,
        p.value,
    )?;
    let mu_equi = equilibria.mu;

    let delta_rho: Vec<f64> = (0..count)
        .map(|i| (rho_liquid[i] - rho_vapour[i]) / STEPS as f64)
        .collect();
    if delta_rho[reference] == 0.0 {
        return Err(AzothError::invalid_input(
            "reference_component",
            "the reference component's density is the same in both phases, so the profile \
             has no integration variable",
        ));
    }

    let mut density = rho_vapour.clone();
    let mut sigma = 0.0;
    for step in 1..STEPS {
        let state = state_at(
            mixture, &reduced, &a_abs, &b_abs, delta1, delta2, &density, rt, p.value,
        )?;
        let mut gradient =
            initial_gradient(&state.dmu, &sqrt_influence, reference, &delta_rho, &density)?;
        let mut converged = false;
        for _ in 0..INNER_STEPS {
            let candidate: Vec<f64> = (0..count).map(|i| density[i] + gradient[i]).collect();
            if step == 1 && candidate.iter().any(|v| !v.is_finite() || *v <= 0.0) {
                eprintln!(
                    "BAD step {step}: candidate {candidate:?} gradient {gradient:?} density {density:?}"
                );
            }
            let inner = state_at(
                mixture, &reduced, &a_abs, &b_abs, delta1, delta2, &candidate, rt, p.value,
            )?;
            let (b, f) =
                residual_system(&mu_equi, &inner.mu, &inner.dmu, &sqrt_influence, reference);
            let x = solve_linear(&f, &b).ok_or_else(|| {
                AzothError::invalid_input(
                    "influence_parameters",
                    format!(
                        "the GT-simple inner system is singular at step {step}; the class's \
                         LU solver throws there and the loop continues with the previous \
                         gradients"
                    ),
                )
            })?;
            let total: f64 = candidate.iter().sum();
            let error: f64 = x.iter().map(|value| value.abs()).sum::<f64>() / total;
            // `x` is over the non-reference components; the update is by component index.
            let others: Vec<usize> = (0..count).filter(|index| *index != reference).collect();
            for (position, &i) in others.iter().enumerate() {
                gradient[i] += x[position];
            }
            if error <= INNER_TOLERANCE {
                converged = true;
                break;
            }
        }
        if !converged && warnings.is_empty() {
            // **The class carries this rather than failing.** Its `err > 1e-15` test is
            // below the floor a double-precision Newton reaches on these states - the
            // loop stalls around `1e-11` - so it runs its 1200 steps and integrates the
            // last iterate. The port does the same and says so.
            warnings.push(Warning::new(
                WarningCode::SolverNotConverged,
                format!(
                    "the inner mu-equality loop reached {INNER_STEPS} steps at step {step}                      without meeting its {INNER_TOLERANCE} tolerance; the last iterate is                      integrated, as the class does"
                ),
            ));
        }

        let final_density: Vec<f64> = (0..count).map(|i| density[i] + gradient[i]).collect();
        let final_state = state_at(
            mixture,
            &reduced,
            &a_abs,
            &b_abs,
            delta1,
            delta2,
            &final_density,
            rt,
            p.value,
        )?;

        let mut kappa = 0.0;
        for i in 0..count {
            let ki = if i == reference {
                1.0
            } else {
                gradient[i] / delta_rho[reference]
            };
            for k in 0..count {
                let kk = if k == reference {
                    1.0
                } else {
                    gradient[k] / delta_rho[reference]
                };
                kappa += sqrt_influence[i]
                    * sqrt_influence[k]
                    * ki
                    * kk
                    * (1.0 - pair_interactions[i][k]);
            }
        }
        let mu_times_rho: f64 = (0..count)
            .map(|i| final_density[i] * (final_state.mu[i] - mu_equi[i]))
            .sum::<f64>()
            - (final_state.pressure - p.value);
        if kappa <= 0.0 || mu_times_rho <= 0.0 {
            return Err(AzothError::invalid_input(
                "z",
                format!(
                    "the gradient-theory integrand is not positive at step {step} \
                     (kappa = {kappa}, mu_times_rho = {mu_times_rho}); the class takes its \
                     square root and stops at the NaN, and this port refuses a partial \
                     integral"
                ),
            ));
        }
        sigma += (2.0 * kappa * mu_times_rho).sqrt() * delta_rho[reference];
        density = final_density;
    }

    Ok(GtSurfaceTensionSimpleResult {
        sigma: newtons_per_meter(sigma.abs()),
        warnings,
    })
}

/// The state at one density vector: pressure, chemical potentials and `dmu/dn` at fixed T,V.
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
            "the GT-simple profile's density vector is not a state: a component density is \
             non-positive or not finite. The class's local system returns a NaN potential \
             there, its loop breaks, and it returns the partial integral accumulated so far",
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
            format!(
                "the GT-simple profile's density vector is not a state: its molar volume \
                 {volume} is at or below the mixture covolume {b_mix}. The class evaluates \
                 its local system there anyway and the potentials become NaN"
            ),
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
    // The Helmholtz derivatives are taken at mole fractions (their dimensionless frame);
    // the chemical potential is the class's physical-density form, and the Jacobian is
    // `dmu/dn` at that density, which is the fraction-frame Hessian divided by the total.
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
    Ok(ProfileState { pressure, mu, dmu })
}

/// The class's initial linearisation: solve for the density gradients from the current
/// `dmu/dn`, then clamp a step that would move a component by more than ten per cent.
fn initial_gradient(
    dmu: &[Vec<f64>],
    sqrt_influence: &[f64],
    reference: usize,
    delta_rho: &[f64],
    density: &[f64],
) -> Result<Vec<f64>> {
    let count = dmu.len();
    // The reference component's step is the constant the class sets once, outside the
    // inner loop; only the others are solved for.
    let mut gradient = vec![0.0; count];
    gradient[reference] = delta_rho[reference];
    let others: Vec<usize> = (0..count).filter(|index| *index != reference).collect();
    if others.is_empty() {
        return Ok(gradient);
    }
    let mut f = vec![vec![0.0; others.len()]; others.len()];
    let mut b = vec![0.0; others.len()];
    for (row, &i) in others.iter().enumerate() {
        b[row] = sqrt_influence[reference] * dmu[i][reference]
            - sqrt_influence[i] * dmu[reference][reference];
        for (column, &k) in others.iter().enumerate() {
            f[row][column] =
                sqrt_influence[i] * dmu[reference][k] - sqrt_influence[reference] * dmu[i][k];
        }
    }
    // **The class's initial step solves the transposed system**: `solveTranspose(b^T)`
    // is `F^T x = b`, while the inner loop's `solve(b)` is `F x = b`.
    let transposed: Vec<Vec<f64>> = (0..others.len())
        .map(|i| (0..others.len()).map(|j| f[j][i]).collect())
        .collect();
    let solved = solve_linear(&transposed, &b).ok_or_else(|| {
        AzothError::invalid_input(
            "influence_parameters",
            "the GT-simple initial system is singular; the class's LU solver throws there",
        )
    })?;
    for (row, &i) in others.iter().enumerate() {
        let step = solved[row] * delta_rho[reference];
        gradient[i] = if (step / density[i]).abs() > 0.1 {
            solved[row].signum() * density[i]
        } else {
            step
        };
    }
    Ok(gradient)
}

/// The class's inner residual system: the fugacity-equality residual and its Jacobian.
fn residual_system(
    mu_equi: &[f64],
    mu: &[f64],
    dmu: &[Vec<f64>],
    sqrt_influence: &[f64],
    reference: usize,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let count = dmu.len();
    let others: Vec<usize> = (0..count).filter(|index| *index != reference).collect();
    let mut f = vec![vec![0.0; others.len()]; others.len()];
    let mut b = vec![0.0; others.len()];
    for (row, &i) in others.iter().enumerate() {
        b[row] = -sqrt_influence[i] * (mu_equi[reference] - mu[reference])
            + sqrt_influence[reference] * (mu_equi[i] - mu[i]);
        for (column, &k) in others.iter().enumerate() {
            f[row][column] =
                -sqrt_influence[i] * dmu[reference][k] + sqrt_influence[reference] * dmu[i][k];
        }
    }
    (b, f)
}

/// Gaussian elimination with partial pivoting, on the small systems the inner loop builds.
fn solve_linear(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let count = b.len();
    if count == 0 {
        return Some(Vec::new());
    }
    let mut matrix = a.to_vec();
    let mut rhs = b.to_vec();
    for column in 0..count {
        let mut pivot = column;
        for row in (column + 1)..count {
            if matrix[row][column].abs() > matrix[pivot][column].abs() {
                pivot = row;
            }
        }
        if matrix[pivot][column].abs() < 1.0e-300 {
            return None;
        }
        matrix.swap(column, pivot);
        rhs.swap(column, pivot);
        for row in (column + 1)..count {
            let factor = matrix[row][column] / matrix[column][column];
            let pivot_tail: Vec<f64> = matrix[column][column..].to_vec();
            for (offset, value) in matrix[row].iter_mut().skip(column).enumerate() {
                *value -= factor * pivot_tail[offset];
            }
            rhs[row] -= factor * rhs[column];
        }
    }
    let mut solution = vec![0.0; count];
    for row in (0..count).rev() {
        let mut sum = rhs[row];
        for k in (row + 1)..count {
            sum -= matrix[row][k] * solution[k];
        }
        solution[row] = sum / matrix[row][row];
    }
    Some(solution)
}
