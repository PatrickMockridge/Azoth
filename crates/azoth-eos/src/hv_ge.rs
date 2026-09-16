//! The Huron-Vidal excess-Gibbs model: a co-volume-weighted NRTL.
//!
//! NeqSim's `PhaseGENRTLmodifiedHV` feeds the Huron-Vidal mixing rule. Its activity
//! coefficient is the NRTL formula with two twists. First, `G_ij = B_i exp(-alpha_ij
//! tau_ij)` - the weight is the row's co-volume, which does not cancel. Second, the
//! `tau` matrix is mixed: pairs the database marks `HV` carry the fitted NRTL parameters
//! (`tau = Dij / T + DijT`), and every other pair carries the cubic's own excess energy,
//! `tau = Lambda * (A_j/B_j - 2 sqrt(A_i A_j)/(B_i + B_j) (1 - kij))`, so the rule stays
//! consistent with the pure-component equation of state.

/// The `tau` and `G` matrices the NRTL sums read, flattened row-major.
#[allow(clippy::too_many_arguments)] // the signature is the GE model's inputs
fn matrices(
    x: &[f64],
    t_kelvin: f64,
    a: &[f64],
    b: &[f64],
    kij: &[f64],
    hv_gij: &[f64],
    hv_gij_t: &[f64],
    hv_alpha: &[f64],
    hv_pairs: &[bool],
    lambda: f64,
) -> (Vec<f64>, Vec<f64>) {
    let n = x.len();
    let mut tau = vec![0.0; n * n];
    let mut g = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let tau_ij = if hv_pairs[i * n + j] {
                hv_gij[i * n + j] / t_kelvin + hv_gij_t[i * n + j]
            } else {
                lambda
                    * (a[j] / b[j]
                        - 2.0 * (a[i] * a[j]).sqrt() / (b[i] + b[j]) * (1.0 - kij[i * n + j]))
            };
            let alpha = if hv_pairs[i * n + j] {
                hv_alpha[i * n + j]
            } else {
                0.0
            };
            tau[i * n + j] = tau_ij;
            g[i * n + j] = b[i] * (-alpha * tau_ij).exp();
        }
    }
    (tau, g)
}

/// The natural logarithms of the activity coefficients for the Huron-Vidal GE model.
///
/// `a` and `b` are the *reduced* attraction and repulsion (`ReducedParameters::a` and
/// `::b`), `kij` the interaction matrix, `hv_gij` the fitted `Dij` in Kelvin, `hv_gij_t`
/// its temperature coefficient `DijT`, `hv_alpha` the fitted non-randomness, and `lambda`
/// the cubic's Huron-Vidal constant `ln((1+delta1)/(1+delta2)) / (delta1 - delta2)`.
/// `hv_gij`, `hv_gij_t` and `hv_alpha` are zero for pairs not marked `hv_pairs`.
#[must_use]
#[allow(clippy::too_many_arguments)] // the signature is the GE model's inputs
pub fn hv_ln_gamma(
    x: &[f64],
    t_kelvin: f64,
    a: &[f64],
    b: &[f64],
    kij: &[f64],
    hv_gij: &[f64],
    hv_gij_t: &[f64],
    hv_alpha: &[f64],
    hv_pairs: &[bool],
    lambda: f64,
) -> Vec<f64> {
    let n = x.len();
    let (tau, g) = matrices(
        x, t_kelvin, a, b, kij, hv_gij, hv_gij_t, hv_alpha, hv_pairs, lambda,
    );
    (0..n)
        .map(|i| {
            let mut first_num = 0.0;
            let mut first_den = 0.0;
            for j in 0..n {
                first_num += tau[j * n + i] * g[j * n + i] * x[j];
                first_den += g[j * n + i] * x[j];
            }
            let mut second = 0.0;
            for j in 0..n {
                let mut den = 0.0;
                let mut num = 0.0;
                for l in 0..n {
                    den += g[l * n + j] * x[l];
                    num += g[l * n + j] * tau[l * n + j] * x[l];
                }
                second += g[i * n + j] * x[j] / den * (tau[i * n + j] - num / den);
            }
            first_num / first_den + second
        })
        .collect()
}

/// The composition derivative `d ln gamma_i / dn_p`, flattened row-major over `(i, p)`.
///
/// NeqSim's `getLnGammadn`, divided by the total moles exactly as NeqSim's
/// `dlngammadn` is. The mole-number derivative is the form the Huron-Vidal fugacity
/// needs, because the attraction parameter mixes `ln gamma` with the composition.
#[must_use]
#[allow(clippy::too_many_arguments)] // the signature is the GE model's inputs
pub fn hv_d_ln_gamma_dn(
    x: &[f64],
    t_kelvin: f64,
    a: &[f64],
    b: &[f64],
    kij: &[f64],
    hv_gij: &[f64],
    hv_gij_t: &[f64],
    hv_alpha: &[f64],
    hv_pairs: &[bool],
    lambda: f64,
) -> Vec<f64> {
    let n = x.len();
    let (tau, g) = matrices(
        x, t_kelvin, a, b, kij, hv_gij, hv_gij_t, hv_alpha, hv_pairs, lambda,
    );

    // The first NRTL term's numerator and denominator, per component.
    let mut a_term = vec![0.0; n];
    let mut b_term = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            a_term[i] += tau[j * n + i] * g[j * n + i] * x[j];
            b_term[i] += g[j * n + i] * x[j];
        }
    }
    // The second NRTL term's `C_j` and `D_j`, per component.
    let mut c_term = vec![0.0; n];
    let mut d_term = vec![0.0; n];
    for j in 0..n {
        for l in 0..n {
            c_term[j] += g[l * n + j] * x[l];
            d_term[j] += g[l * n + j] * tau[l * n + j] * x[l];
        }
    }

    let total: f64 = x.iter().sum();
    let mut result = vec![0.0; n * n];
    for i in 0..n {
        for p in 0..n {
            let d_a_dn = tau[p * n + i] * g[p * n + i];
            let d_b_dn = g[p * n + i];
            let d_e_dn = g[i * n + p] * tau[i * n + p];
            let mut d_sum = 0.0;
            let mut f_sum = 0.0;
            let mut g_sum = 0.0;
            for f in 0..n {
                d_sum +=
                    g[p * n + f] * g[i * n + f] * tau[i * n + f] * x[f] / (c_term[f] * c_term[f]);
                f_sum += x[f] * g[p * n + f] * d_term[f] * g[i * n + f]
                    / (c_term[f] * c_term[f] * c_term[f]);
                g_sum +=
                    x[f] * g[p * n + f] * tau[p * n + f] * g[i * n + f] / (c_term[f] * c_term[f]);
            }
            result[i * n + p] = (d_a_dn / b_term[i] - a_term[i] * d_b_dn / (b_term[i] * b_term[i])
                + d_e_dn / c_term[p]
                - d_sum
                - d_term[p] * g[i * n + p] / (c_term[p] * c_term[p])
                + 2.0 * f_sum
                - g_sum)
                / total;
        }
    }
    result
}
