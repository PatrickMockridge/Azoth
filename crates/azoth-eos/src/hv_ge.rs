//! The Huron-Vidal excess-Gibbs model: a co-volume-weighted NRTL.
//!
//! NeqSim's `PhaseGENRTLmodifiedHV` feeds the Huron-Vidal mixing rule. Its activity
//! coefficient is the NRTL formula with two twists. First, `G_ij = B_i exp(-alpha_ij
//! tau_ij)` - the weight is the row's co-volume, which does not cancel. Second, the
//! `tau` matrix is mixed: pairs the database marks `HV` carry the fitted NRTL parameters
//! (`tau = Dij / T`), and every other pair carries the cubic's own excess energy,
//! `tau = Lambda * (A_j/B_j - 2 sqrt(A_i A_j)/(B_i + B_j) (1 - kij))`, so the rule stays
//! consistent with the pure-component equation of state.

/// The natural logarithms of the activity coefficients for the Huron-Vidal GE model.
///
/// `a` and `b` are the *reduced* attraction and repulsion (`ReducedParameters::a` and
/// `::b`), `kij` the interaction matrix, `hv_gij` the fitted `Dij` in Kelvin (zero for
/// pairs not marked `hv_pairs`), `hv_alpha` the fitted non-randomness, and `lambda` the
/// cubic's Huron-Vidal constant `ln((1+delta1)/(1+delta2)) / (delta1 - delta2)`.
#[must_use]
pub fn hv_ln_gamma(
    x: &[f64],
    t_kelvin: f64,
    a: &[f64],
    b: &[f64],
    kij: &[f64],
    hv_gij: &[f64],
    hv_alpha: &[f64],
    hv_pairs: &[bool],
    lambda: f64,
) -> Vec<f64> {
    let n = x.len();
    let mut tau = vec![0.0; n * n];
    let mut g = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let tau_ij = if hv_pairs[i * n + j] {
                hv_gij[i * n + j] / t_kelvin
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
