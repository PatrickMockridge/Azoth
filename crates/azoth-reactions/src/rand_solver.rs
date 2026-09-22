//! The modified-RAND reactive equilibrium solve, from
//! `flashops/reactiveflash/ModifiedRANDSolver.java`.
//!
//! Simultaneous chemical and phase equilibrium by minimising the Gibbs energy subject to the
//! element balances: the unknowns are the moles `n[j][i]` and one Lagrange multiplier `λ_k`
//! per element row, and each pass corrects them together:
//!
//! ```text
//! e[j][i] = g0_i + ln x_j,i + ln phi_j,i - sum_k lambda_k A_ki     the potential error
//! C_kl    = sum_j sum_i A_ki A_li n_j,i                            the RAND matrix
//! rhs_k   = (b_k - sum A n) + sum A n e                            the inventory plus the error
//! n_j,i  *= exp(alpha (-e_j,i + sum_k A_ki dlambda_k))             the damped step
//! lambda += alpha dlambda
//! ```
//!
//! # What is reproduced and what is not, yet
//!
//! Reproduced: the potential error, the matrix and its right-hand side, the Tikhonov term,
//! the symmetric diagonal scaling, the damped step with its backtracking line search, the
//! damping recovery rule and the two convergence tests. **Not yet ported**: `DIISAccelerator`
//! (`diis/`, 225 lines) and the multiphase and ionic branches - a gas-phase ion is pinned to
//! `EPS` there, and the reference state of an electrolyte phase is corrected by
//! `getLogInfiniteDiluteFugacity`. Those need the phase's *type* and its ions, which is the
//! P8 seam, so this is the single-phase neutral path the captured benchmark states take.
//!
//! # `ln phi` is a closure, not a dependency
//!
//! The solver needs a fugacity coefficient at a *trial* composition, and NeqSim gets one by
//! writing the moles into its own phase and calling `init(1)`. Here the caller supplies it -
//! which keeps this crate free of an equation of state, and keeps the question of *which*
//! phase model answers a separate one.

use azoth_core::{AzothError, Result};

/// The residual both the potential error and the element balance must come in under, from
/// `ModifiedRANDSolver.TOL`.
pub const TOL: f64 = 1.0e-9;

/// The pass cap, from `MAX_ITER`.
pub const MAX_ITERATIONS: usize = 500;

/// The floor a mole number is kept above, from `EPS`.
pub const EPS: f64 = 1.0e-30;

/// The gas constant NeqSim's solver carries, from `R_GAS`. **Not `azoth-core`'s**: this is a
/// literal in the class, and the standard potentials are computed with it.
pub const R_GAS: f64 = 8.314462;

/// The reference pressure the `ln(P/P_ref)` term is taken against, from `P_REF`.
pub const P_REF: f64 = 1.0;

/// The temperature the Cp polynomials are integrated from, from `computeG0`'s `T0`.
pub const T0: f64 = 298.15;

/// The singular-pivot floor `solveLinear` refuses below, from `mx < 1e-30`.
pub const LINEAR_PIVOT_FLOOR: f64 = 1.0e-30;

/// A component's ideal-gas heat-capacity polynomial and its two formation properties.
///
/// The five coefficients are `Cp = A + B T + C T^2 + D T^3 + E T^4`, read from the component
/// databank by `eos.ideal_gas_cp`'s reader.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermoData {
    /// Ideal-gas enthalpy of formation at [`T0`], in J/mol.
    pub enthalpy_of_formation: f64,
    /// Ideal-gas absolute entropy at [`T0`], in J/(mol*K).
    pub absolute_entropy: f64,
    /// Gibbs energy of formation, in J/mol - the fallback when neither of the other two is
    /// stated.
    pub gibbs_energy_of_formation: f64,
    /// The heat-capacity polynomial's five coefficients, `A` through `E`.
    pub cp: [f64; 5],
}

/// The standard chemical potentials `g0`, from `computeG0`'s neutral branch.
///
/// `g0_i = (hT - T sT) / (R T) + ln(P/P_ref)`, where `hT` and `sT` are the formation
/// properties corrected from [`T0`] to `T` by integrating the Cp polynomial:
///
/// ```text
/// hT = dHf + integral(T0,T) Cp dT
/// sT = S0  + integral(T0,T) Cp/T dT
/// ```
///
/// The three fallbacks are the class's own, in its order: with `dHf` and `S0` but no Cp data
/// the polynomial is dropped, with only `dGf298` the Gibbs energy is used directly, and with
/// nothing at all the potential is `ln(P/P_ref)` alone - which is the ideal-gas potential of a
/// substance whose formation data the databank does not carry.
#[must_use]
pub fn standard_potentials(data: &[ThermoData], temperature: f64, pressure: f64) -> Vec<f64> {
    let rt = R_GAS * temperature;
    let ln_p = (pressure / P_REF).ln();
    data.iter()
        .map(|entry| {
            let [cp_a, cp_b, cp_c, cp_d, cp_e] = entry.cp;
            let has_cp_data = cp_a.abs() > 1.0e-10 || cp_b.abs() > 1.0e-10;
            let has_thermo = entry.enthalpy_of_formation.abs() > 1.0e-10
                || entry.absolute_entropy.abs() > 1.0e-10;

            if has_thermo && has_cp_data {
                let dt = temperature - T0;
                let dt2 = temperature.powi(2) - T0.powi(2);
                let dt3 = temperature.powi(3) - T0.powi(3);
                let dt4 = temperature.powi(4) - T0.powi(4);
                let dt5 = temperature.powi(5) - T0.powi(5);
                let ln_ratio = (temperature / T0).ln();
                let delta_h = cp_a * dt
                    + cp_b / 2.0 * dt2
                    + cp_c / 3.0 * dt3
                    + cp_d / 4.0 * dt4
                    + cp_e / 5.0 * dt5;
                let delta_s = cp_a * ln_ratio
                    + cp_b * dt
                    + cp_c / 2.0 * dt2
                    + cp_d / 3.0 * dt3
                    + cp_e / 4.0 * dt4;
                let h_t = entry.enthalpy_of_formation + delta_h;
                let s_t = entry.absolute_entropy + delta_s;
                (h_t - temperature * s_t) / rt + ln_p
            } else if has_thermo {
                (entry.enthalpy_of_formation - temperature * entry.absolute_entropy) / rt + ln_p
            } else if entry.gibbs_energy_of_formation.abs() > 1.0e-10 {
                entry.gibbs_energy_of_formation / rt + ln_p
            } else {
                ln_p
            }
        })
        .collect()
}

/// What the solve answers with.
#[derive(Debug, Clone, PartialEq)]
pub struct RandSolution {
    /// The equilibrium moles, one per component.
    pub moles: Vec<f64>,
    /// The element Lagrange multipliers at the answer, one per element row.
    pub lambda: Vec<f64>,
    /// Passes taken.
    pub iterations: u32,
    /// The largest absolute potential error over the species holding more than `1e-10` of the
    /// total.
    pub max_error: f64,
    /// `max |sum A n - b|`, the element balance's own residual.
    pub element_residual: f64,
    /// Whether both came in under [`TOL`].
    pub converged: bool,
}

/// The single-phase neutral RAND solve.
///
/// `ln_phi` is handed a *mole-fraction* composition and answers each component's
/// `ln(phi_i)` there - NeqSim writes the moles into its own phase and reads
/// `getFugacityCoefficient`, and this is the same question asked of whatever phase model the
/// caller has.
///
/// # Errors
/// [`AzothError::InvalidInput`] where the shapes disagree, and whatever `ln_phi` raises.
// Indexed rather than iterated, deliberately: every loop here is NeqSim's own loop, over a
// matrix and a phase's components, and the correspondence is what makes the port checkable
// against the class line by line.
#[allow(clippy::needless_range_loop)]
pub fn solve(
    a_matrix: &[Vec<f64>],
    g0: &[f64],
    b: &[f64],
    feed_moles: &[f64],
    ln_phi: &mut dyn FnMut(&[f64]) -> Result<Vec<f64>>,
) -> Result<RandSolution> {
    let ne = b.len();
    let nc = feed_moles.len();
    if g0.len() != nc || a_matrix.len() != ne {
        return Err(AzothError::InvalidInput {
            field: "a_matrix".to_string(),
            reason: format!(
                "{ne} element row(s), {} potential(s) and {} component(s)",
                g0.len(),
                nc
            ),
        });
    }
    for row in a_matrix {
        if row.len() != nc {
            return Err(AzothError::InvalidInput {
                field: "a_matrix".to_string(),
                reason: format!("a row has {} entries against {nc} component(s)", row.len()),
            });
        }
    }

    // `initialize`: one phase, so `beta` is 1 and the moles are the feed's, floored.
    let mut n: Vec<f64> = feed_moles.iter().map(|moles| moles.max(EPS)).collect();
    let mut total: f64 = n.iter().sum();
    let mut fractions: Vec<f64> = n.iter().map(|moles| moles / total).collect();

    let mut ln_phi_here = ln_phi(&fractions)?;
    let mut lambda = initial_lambda(a_matrix, g0, &fractions, &ln_phi_here);

    let mut damping = 1.0_f64;
    let mut previous_residual = f64::MAX;
    let mut stagnation = 0_u32;
    let mut final_error = f64::MAX;
    let mut final_element = f64::MAX;
    let mut converged = false;
    let mut iterations = 0_u32;

    for iteration in 0..MAX_ITERATIONS {
        iterations = iteration as u32 + 1;

        // The potential error, per component.
        let mut error = vec![0.0_f64; nc];
        for i in 0..nc {
            let x_i = fractions[i].max(EPS);
            let mut row_sum = 0.0;
            for k in 0..ne {
                row_sum += lambda[k] * a_matrix[k][i];
            }
            let value = g0[i] + x_i.ln() + ln_phi_here[i] - row_sum;
            error[i] = if value.is_finite() { value } else { 0.0 };
        }

        // The RAND matrix and its right-hand side.
        let mut c = vec![vec![0.0_f64; ne]; ne];
        let mut rhs = vec![0.0_f64; ne];
        for k in 0..ne {
            let mut element_sum = 0.0;
            let mut element_error = 0.0;
            for i in 0..nc {
                element_sum += a_matrix[k][i] * n[i];
                element_error += a_matrix[k][i] * n[i] * error[i];
            }
            rhs[k] = (b[k] - element_sum) + element_error;
            for l in 0..ne {
                let mut value = 0.0;
                for i in 0..nc {
                    value += a_matrix[k][i] * a_matrix[l][i] * n[i];
                }
                c[k][l] = value;
            }
            // Tikhonov regularization, as the class does it.
            c[k][k] += (1.0e-10 * c[k][k].abs()).max(1.0e-14);
        }

        // Symmetric diagonal scaling, for the conditioning the ions cause.
        let scale: Vec<f64> = (0..ne)
            .map(|k| {
                let diagonal = c[k][k].abs();
                if diagonal > 1.0e-30 {
                    1.0 / diagonal.sqrt()
                } else {
                    1.0
                }
            })
            .collect();
        for k in 0..ne {
            for l in 0..ne {
                c[k][l] *= scale[k] * scale[l];
            }
            rhs[k] *= scale[k];
        }

        let Some(mut delta) = solve_linear(&c, &rhs) else {
            break;
        };
        for k in 0..ne {
            delta[k] *= scale[k];
        }

        let n_old = n.clone();
        let lambda_old = lambda.clone();

        // The damped step, with up to five halvings while the element residual does not fall.
        let mut alpha = damping;
        let mut accepted = false;
        for _ in 0..5 {
            n.copy_from_slice(&n_old);
            lambda.copy_from_slice(&lambda_old);
            for i in 0..nc {
                let mut correction = alpha * -error[i];
                for k in 0..ne {
                    correction += alpha * a_matrix[k][i] * delta[k];
                }
                if !correction.is_finite() {
                    correction = 0.0;
                }
                let correction = correction.clamp(-3.0, 3.0);
                let stepped = n[i] * correction.exp();
                n[i] = if !stepped.is_finite() || stepped < EPS {
                    EPS
                } else {
                    stepped
                };
            }
            for k in 0..ne {
                lambda[k] += alpha * delta[k];
            }
            total = n.iter().sum::<f64>().max(EPS);
            let element_residual = element_residual(a_matrix, b, &n);
            if element_residual < previous_residual * 1.5 || alpha < 0.05 {
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }

        if !accepted {
            // The class's own fallback: restore and take a tenth of the step.
            n.copy_from_slice(&n_old);
            lambda.copy_from_slice(&lambda_old);
            let alpha = 0.1_f64;
            for i in 0..nc {
                let mut correction = alpha * -error[i];
                for k in 0..ne {
                    correction += alpha * a_matrix[k][i] * delta[k];
                }
                if !correction.is_finite() {
                    correction = 0.0;
                }
                let correction = correction.clamp(-3.0, 3.0);
                let stepped = n[i] * correction.exp();
                n[i] = if !stepped.is_finite() || stepped < EPS {
                    EPS
                } else {
                    stepped
                };
            }
            for k in 0..ne {
                lambda[k] += alpha * delta[k];
            }
            total = n.iter().sum::<f64>().max(EPS);
        }

        fractions = n.iter().map(|moles| moles / total).collect();
        ln_phi_here = ln_phi(&fractions)?;

        // The convergence test: the worst potential error among the species that matter, and
        // the element balance.
        let mut max_error = 0.0_f64;
        for i in 0..nc {
            if n[i] > 1.0e-10 * total {
                max_error = max_error.max(error[i].abs());
            }
        }
        let element = element_residual(a_matrix, b, &n);
        final_error = max_error;
        final_element = element;
        let residual = max_error.max(element);

        if max_error < TOL && element < TOL {
            converged = true;
            break;
        }

        // The damping recovery rule, single phase.
        if residual < previous_residual * 0.9 {
            damping = (damping * 1.5).min(1.0);
            stagnation = 0;
        } else if residual > previous_residual * 1.1 {
            damping = (damping * 0.5).max(0.01);
            stagnation += 1;
        } else {
            stagnation += 1;
        }
        previous_residual = residual;

        // A long stagnation at a small residual is accepted, as the class accepts it.
        if stagnation > 50 && max_error < 1.0e-3 && element < 1.0e-4 {
            converged = true;
            break;
        }
    }

    Ok(RandSolution {
        moles: n,
        lambda,
        iterations,
        max_error: final_error,
        element_residual: final_element,
        converged,
    })
}

/// `max |sum_i A_k,i n_i - b_k|` over the element rows, from `computeElementResidual`.
fn element_residual(a_matrix: &[Vec<f64>], b: &[f64], n: &[f64]) -> f64 {
    let mut worst = 0.0_f64;
    for (k, row) in a_matrix.iter().enumerate() {
        let sum: f64 = row.iter().zip(n).map(|(a, n)| a * n).sum();
        worst = worst.max((sum - b[k]).abs());
    }
    worst
}

/// `initializeLambda`: the least-squares multipliers for the initial potentials,
/// `lambda = (A A^T)^-1 A h` with `h_i = g0_i + ln x_i + ln phi_i`.
#[allow(clippy::needless_range_loop)] // the same index loops the class writes
fn initial_lambda(
    a_matrix: &[Vec<f64>],
    g0: &[f64],
    fractions: &[f64],
    ln_phi: &[f64],
) -> Vec<f64> {
    let ne = a_matrix.len();
    let nc = g0.len();
    let h: Vec<f64> = (0..nc)
        .map(|i| g0[i] + fractions[i].max(EPS).ln() + ln_phi[i])
        .collect();
    let mut ata = vec![vec![0.0_f64; ne]; ne];
    let mut ath = vec![0.0_f64; ne];
    for k in 0..ne {
        for l in 0..ne {
            for i in 0..nc {
                ata[k][l] += a_matrix[k][i] * a_matrix[l][i];
            }
        }
        for i in 0..nc {
            ath[k] += a_matrix[k][i] * h[i];
        }
    }
    solve_linear(&ata, &ath).unwrap_or_else(|| vec![0.0; ne])
}

/// `solveLinear`: Gaussian elimination with partial pivoting, refusing a pivot under
/// `1e-30`.
#[allow(clippy::needless_range_loop)] // the elimination's own loops, as the class writes them
fn solve_linear(matrix: &[Vec<f64>], rhs: &[f64]) -> Option<Vec<f64>> {
    let dim = rhs.len();
    let mut augment = vec![vec![0.0_f64; dim + 1]; dim];
    for i in 0..dim {
        augment[i][..dim].copy_from_slice(&matrix[i][..dim]);
        augment[i][dim] = rhs[i];
    }

    for column in 0..dim {
        let mut pivot = column;
        let mut largest = augment[column][column].abs();
        for row in (column + 1)..dim {
            if augment[row][column].abs() > largest {
                largest = augment[row][column].abs();
                pivot = row;
            }
        }
        if largest < LINEAR_PIVOT_FLOOR {
            return None;
        }
        augment.swap(column, pivot);
        for row in (column + 1)..dim {
            let factor = augment[row][column] / augment[column][column];
            for k in column..=dim {
                augment[row][k] -= factor * augment[column][k];
            }
        }
    }

    let mut solution = vec![0.0_f64; dim];
    for i in (0..dim).rev() {
        let mut sum = augment[i][dim];
        for k in (i + 1)..dim {
            sum -= augment[i][k] * solution[k];
        }
        solution[i] = sum / augment[i][i];
    }
    Some(solution)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_linear_solve_refuses_a_singular_pivot() {
        let matrix = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        assert_eq!(solve_linear(&matrix, &[1.0, 2.0]), None);
        let good = vec![vec![2.0, 0.0], vec![0.0, 4.0]];
        assert_eq!(solve_linear(&good, &[2.0, 8.0]), Some(vec![1.0, 2.0]));
    }

    /// The three fallbacks, in the class's own order.
    #[test]
    fn the_standard_potentials_fall_back_the_way_the_class_does() {
        let with_cp = ThermoData {
            enthalpy_of_formation: -110_525.0,
            absolute_entropy: 197.66,
            gibbs_energy_of_formation: -137_168.0,
            cp: [25.56759, 6.096130, 4.054656, -2.671301, 0.131021],
        };
        let without_cp = ThermoData {
            cp: [0.0; 5],
            ..with_cp
        };
        let only_gibbs = ThermoData {
            enthalpy_of_formation: 0.0,
            absolute_entropy: 0.0,
            ..with_cp
        };
        let nothing = ThermoData {
            enthalpy_of_formation: 0.0,
            absolute_entropy: 0.0,
            gibbs_energy_of_formation: 0.0,
            cp: [0.0; 5],
        };

        let potentials =
            standard_potentials(&[with_cp, without_cp, only_gibbs, nothing], 600.0, 1.0);
        // The Cp polynomial and the constant-Cp fallback are different numbers, which is the
        // point of integrating it; the last two are the `dGf` and `ln(P)` branches.
        assert!(
            (potentials[0] - potentials[1]).abs() > 1.0,
            "{potentials:?}"
        );
        assert!((potentials[2] - (-137_168.0 / (R_GAS * 600.0))).abs() < 1e-12);
        assert_eq!(potentials[3], 0.0, "ln(P/P_ref) at the reference pressure");
    }
}
