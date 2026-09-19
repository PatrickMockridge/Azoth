//! `eos.unifac_activity_coefficients` - the activity coefficients from the UNIFAC
//! group-contribution model.
//!
//! Spec: `specs/models/eos/unifac_activity_coefficients.toml`. A *direct* model: no
//! iteration, so no `algorithm` block.

use azoth_core::{AzothError, Result, apply_checks};

use crate::databank::UnifacParameters;
use crate::model_gen;
use crate::results::UnifacActivityCoefficientsResult;

/// The residual `ln Gamma_k` for one group, at a surface-fraction distribution `theta`,
/// and `T d ln Gamma_k / dT`.
///
/// `aij` is `G x G` row-major, so `aij[m][n]` is `aij[m * g + n]`. `daij_dt` is the same
/// shape and is the interaction's own temperature derivative, which a model whose tables
/// fit `a + b (T - 298.15) + c (T - 298.15)^2` has and the plain tables do not - an empty
/// slice means the matrix does not move, and the derivative is then zero.
///
/// **Every `T` here enters through `exp(-a_mn/T)`.** With `E_mn = exp(-a_mn/T)` that makes
/// `T dE_mn/dT = E_mn (a_mn/T - da_mn/dT)`, and the two sums are differentiated by the
/// quotient and chain rules - `d/dT ln s1 = (T d s1/dT)/s1/T`, and the `s3` term is a
/// quotient with `s2` in the denominator.
fn ln_gamma_group(
    k: usize,
    theta: &[f64],
    group_q: &[f64],
    aij: &[f64],
    daij_dt: &[f64],
    t: f64,
    g: usize,
) -> (f64, f64) {
    let moves = daij_dt.len() == aij.len();
    let d = |m: usize, n: usize| -> f64 { if moves { daij_dt[m * g + n] } else { 0.0 } };
    // `T dE_mn/dT` for each pair, which is what every sum below carries.
    let t_de = |m: usize, n: usize| -> f64 {
        let e = (-aij[m * g + n] / t).exp();
        e * (aij[m * g + n] / t - d(m, n))
    };

    let mut s1 = 0.0;
    let mut t_ds1 = 0.0;
    for m in 0..g {
        s1 += theta[m] * (-aij[m * g + k] / t).exp();
        t_ds1 += theta[m] * t_de(m, k);
    }
    let mut s3 = 0.0;
    let mut t_ds3 = 0.0;
    for m in 0..g {
        let mut s2 = 0.0;
        let mut t_ds2 = 0.0;
        for n in 0..g {
            s2 += theta[n] * (-aij[n * g + m] / t).exp();
            t_ds2 += theta[n] * t_de(n, m);
        }
        let e_km = (-aij[k * g + m] / t).exp();
        s3 += theta[m] * e_km / s2;
        t_ds3 += theta[m] * (t_de(k, m) * s2 - e_km * t_ds2) / (s2 * s2);
    }

    let ln_gamma = group_q[k] * (1.0 - s1.ln() - s3);
    // `T d/dT` of the same expression: the constant drops, `ln s1` becomes `T ds1/s1` and
    // the quotient term is the one accumulated above.
    let t_d_ln_gamma = group_q[k] * (-t_ds1 / s1 - t_ds3);
    (ln_gamma, t_d_ln_gamma)
}

/// The activity coefficients of a mixture, from UNIFAC.
///
/// `params.groups[i][k]` is the count of group `k` in component `i`, over the union of
/// the named components' groups (`N x G`). `group_r`/`group_q` are the per-group volume
/// and surface area, and `aij[m][n] = a_{main(m), main(n)}` in Kelvin. The component
/// volume and area are `R_i = sum_k groups[i][k] group_r[k]` and the same for `Q_i`;
/// the combinatorial term uses `Z = 10`, and the residual is the standard
/// `ln gamma^R_i = sum_k nu_k^i (ln Gamma_k^mix - ln Gamma_k^pure)`. `params` is
/// resolved by name through [`crate::databank::unifac_parameters`]; `x` is checked
/// rather than renormalised.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the shapes disagree (`groups` is `N x G`, `aij` is
///   `G x G`), a group count is negative, or `x` is not a composition.
/// * [`AzothError::OutOfRange`] if `T` is not positive.
///
/// # Example
/// ```
/// use azoth_eos::databank::UnifacParameters;
/// use azoth_eos::unifac_activity_coefficients;
///
/// let params = UnifacParameters {
///     groups: vec![1.0, 0.0, 0.0, 1.0],
///     group_r: vec![1.4311, 0.92],
///     group_q: vec![1.432, 1.4],
///     aij: vec![0.0, -181.0, 289.6, 0.0],
/// };
/// let r = unifac_activity_coefficients(&params, 298.15, &[0.5, 0.5])?;
/// assert!((r.gamma[0] - 1.1156815062468024).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T` and `x` are the symbols in the chemistry
pub fn unifac_activity_coefficients(
    params: &UnifacParameters,
    T: f64,
    x: &[f64],
) -> Result<UnifacActivityCoefficientsResult> {
    let spec = &model_gen::UNIFAC_ACTIVITY_COEFFICIENTS_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T),
            _ => None,
        },
        &mut warnings,
    )?;

    let n = x.len();
    if n == 0 {
        return Err(AzothError::invalid_input(
            "x",
            "a mixture of zero components has no activity coefficient",
        ));
    }
    let g = params.group_r.len();
    if g == 0 || params.group_q.len() != g {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "`group_r` has {g} entries and `group_q` has {}; both must be the same \
                 non-empty length",
                params.group_q.len()
            ),
        ));
    }
    if params.groups.len() != n * g {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "`groups` has {} entries but `x` has {n} components and there are {g} \
                 groups, which needs {}",
                params.groups.len(),
                n * g
            ),
        ));
    }
    if let Some(bad) = params.groups.iter().position(|&count| count < 0.0) {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "group count {bad} is {} but a count cannot be negative",
                bad
            ),
        ));
    }
    if params.aij.len() != g * g {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "`aij` has {} entries but must be {g} x {g}",
                params.aij.len()
            ),
        ));
    }
    if let Some(bad) = x.iter().position(|&value| value < 0.0) {
        return Err(AzothError::invalid_input(
            "x",
            format!(
                "x[{bad}] is {} but a mole fraction cannot be negative",
                x[bad]
            ),
        ));
    }
    let sum: f64 = x.iter().sum();
    if (sum - 1.0).abs() > 1.0e-09 {
        return Err(AzothError::invalid_input(
            "x",
            format!(
                "the mole fractions sum to {sum}, not to one. Renormalising them here \
                 would make a composition error invisible in every number downstream, \
                 so it is refused instead"
            ),
        ));
    }

    // The plain tables fit `a + b T + c T^2` about zero in the PSRK variant and carry no
    // offset at all here, so the interaction matrix does not move with the temperature and
    // the derivative is empty.
    let (ln_gamma, gamma, _) = unifac_ln_gamma(
        &UnifacBasis {
            groups: &params.groups,
            group_r: &params.group_r,
            group_q: &params.group_q,
            aij: &params.aij,
            daij_dt: &[],
        },
        T,
        x,
        Combinatorial::StavermanGuggenheim,
    );

    Ok(UnifacActivityCoefficientsResult {
        ln_gamma,
        gamma,
        warnings,
    })
}

/// A resolved group basis and its interaction matrix, as the kernel reads them.
///
/// `daij_dt` is the interaction's own temperature derivative, the same shape as `aij`; an
/// empty slice means the matrix does not move with the temperature, which is every model
/// whose tables do not fit `a + b (T - 298.15) + c (T - 298.15)^2`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct UnifacBasis<'a> {
    /// Per-component group counts, `N x G` row-major.
    pub groups: &'a [f64],
    /// The volume `R` of each group, length `G`.
    pub group_r: &'a [f64],
    /// The surface area `Q` of each group, length `G`.
    pub group_q: &'a [f64],
    /// The interaction matrix `G x G` row-major, in Kelvin.
    pub aij: &'a [f64],
    /// `d aij / dT`, the same shape, or empty.
    pub daij_dt: &'a [f64],
}

/// Which combinatorial term a UNIFAC `ln gamma` is built from.
///
/// **NeqSim's two UNIFAC components do not share one.** `ComponentGEUnifac` uses the
/// Staverman-Guggenheim form and `ComponentGEUnifacUMRPRU` uses the Flory-Huggins term
/// alone - no `ln(V/x)` and no `l_i`. Measured at 298.15 K, methane/water 0.98/0.02:
/// the two differ by `6.887827e-06` and `0.018913275` in `ln gamma`, which is the whole
/// of the difference between this library and NeqSim's UMR-CPA at that state.
pub enum Combinatorial {
    /// `ComponentGEUnifac`'s, and PSRK's.
    StavermanGuggenheim,
    /// `ComponentGEUnifacUMRPRU`'s.
    FloryHuggins,
}

/// The UNIFAC activity coefficients from a resolved group basis and interaction matrix.
///
/// The whole of the physics, shared by `eos.unifac_activity_coefficients` and
/// `eos.unifac_psrk_activity_coefficients`, which differ only in where the interaction
/// matrix comes from: the PSRK one evaluates `a + b T + c T^2` at the state's
/// temperature before calling this. `eos.unifac_umrpru_activity_coefficients` calls it
/// too, with a different interaction matrix *and* a different combinatorial term - the
/// one thing about that model the shared kernel cannot assume. Everything here is
/// checked by the caller, so this assumes `groups` is `N x G`, `aij` is `G x G` and `x`
/// is a composition of `N`.
pub(crate) fn unifac_ln_gamma(
    basis: &UnifacBasis<'_>,
    t: f64,
    x: &[f64],
    combinatorial: Combinatorial,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let UnifacBasis {
        groups,
        group_r,
        group_q,
        aij,
        daij_dt,
    } = *basis;
    let n = x.len();
    let g = group_r.len();
    let mut ri = vec![0.0; n];
    let mut qi = vec![0.0; n];
    for i in 0..n {
        for k in 0..g {
            ri[i] += groups[i * g + k] * group_r[k];
            qi[i] += groups[i * g + k] * group_q[k];
        }
    }

    let mut ln_gamma = vec![0.0; n];
    let mut gamma = vec![0.0; n];
    // `T d ln gamma_i/dT`, which only the residual carries: the combinatorial term is
    // built from `r_i` and `q_i`, which do not move with the temperature.
    let mut t_d_ln_gamma = vec![0.0; n];
    for i in 0..n {
        let mut t1 = 0.0;
        let mut t2 = 0.0;
        let mut suml = 0.0;
        for j in 0..n {
            t1 += x[j] * ri[j];
            t2 += x[j] * qi[j];
            suml += x[j] * (5.0 * (ri[j] - qi[j]) - (ri[j] - 1.0));
        }
        let v = x[i] * ri[i] / t1;
        let f = x[i] * qi[i] / t2;
        let lng_c = match combinatorial {
            Combinatorial::StavermanGuggenheim => {
                let li = 5.0 * (ri[i] - qi[i]) - (ri[i] - 1.0);
                (v / x[i]).ln() + 5.0 * qi[i] * (f / v).ln() + li - (v / x[i]) * suml
            }
            Combinatorial::FloryHuggins => -5.0 * qi[i] * ((v / f).ln() + 1.0 - v / f),
        };

        let denom: f64 = (0..n).map(|j| x[j] * qi[j]).sum();
        let mut qmix = vec![0.0; g];
        for l in 0..g {
            let num: f64 = (0..n).map(|j| x[j] * groups[j * g + l]).sum();
            qmix[l] = group_q[l] * num / denom;
        }
        let mut qcomp = vec![0.0; g];
        for l in 0..g {
            qcomp[l] = group_q[l] * groups[i * g + l] / qi[i];
        }
        let mut lng_r = 0.0;
        let mut t_d_lng_r = 0.0;
        for k in 0..g {
            let (mix, t_d_mix) = ln_gamma_group(k, &qmix, group_q, aij, daij_dt, t, g);
            let (comp, t_d_comp) = ln_gamma_group(k, &qcomp, group_q, aij, daij_dt, t, g);
            lng_r += groups[i * g + k] * (mix - comp);
            t_d_lng_r += groups[i * g + k] * (t_d_mix - t_d_comp);
        }

        let lng = lng_c + lng_r;
        ln_gamma[i] = lng;
        gamma[i] = lng.exp();
        t_d_ln_gamma[i] = t_d_lng_r;
    }

    (ln_gamma, gamma, t_d_ln_gamma)
}
