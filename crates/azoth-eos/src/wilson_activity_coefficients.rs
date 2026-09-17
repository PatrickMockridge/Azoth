//! `eos.wilson_activity_coefficients` - the activity coefficients from the
//! paraffin-wax Wilson model.
//!
//! Spec: `specs/models/eos/wilson_activity_coefficients.toml`. A *direct* model: no
//! iteration, so no `algorithm` block.

use azoth_core::{AzothError, Result, apply_checks};

use crate::mixture::Mixture;
use crate::model_gen;
use crate::results::WilsonActivityCoefficientsResult;

/// The gas constant, in J/(mol*K), the value NeqSim's `ThermodynamicConstantsInterface`
/// declares rather than the exact SI definition.
const R: f64 = 8.3144621;

/// The per-component "interaction energy" `lambda_i = -(2/z)(DeltaH_sub - R T)`, built
/// from a vaporization-enthalpy correlation, a total-transition correlation and a
/// fusion-temperature correlation over carbon number.
fn interaction_energy(m: f64, tc: f64, t: f64) -> f64 {
    let coordination = 6.0;
    let carbon = m / 0.014;
    let x = 1.0 - t / tc;
    let d0 = 5.2804 * x.powf(0.3333) + 12.865 * x.powf(0.8333) + 1.171 * x.powf(1.2083)
        - 13.166 * x
        + 0.4858 * x.powi(2)
        - 1.088 * x.powi(3);
    let d1 = 0.80022 * x.powf(0.3333) + 273.23 * x.powf(0.8333) + 465.08 * x.powf(1.2083)
        - 638.51 * x
        - 145.12 * x.powi(2)
        - 74.049 * x.powi(3);
    let d2 = 7.2543 * x.powf(0.3333) - 346.45 * x.powf(0.8333) - 610.48 * x.powf(1.2083)
        + 839.89 * x
        + 160.05 * x.powi(2)
        - 50.711 * x.powi(3);
    let omega = 0.052_075 + 0.044_894_6 * carbon - 0.000_185_397 * carbon * carbon;
    let dh_vap = R * tc * (d0 + omega * d1 + omega * omega * d2) * 4.1868;
    let dh_tot = (3.7791 * carbon - 12.654) * 1000.0;
    let tf = 374.5 + 0.2617 * m - 20.172 / m;
    let dh_f = 0.1426 * m * tf * 4.1868;
    let dh_trans = dh_tot - dh_f;
    let dh_sub = dh_vap + dh_f + dh_trans;
    -2.0 / coordination * (dh_sub - R * t)
}

/// The Wilson energy parameter `Lambda_ij`. The heavier-first branch collapses the
/// exponent to zero - the quirk in NeqSim's `getCharEnergyParamter` - so it is `1.0`.
fn char_energy(m: &[f64], tc: &[f64], t: f64, i: usize, j: usize) -> f64 {
    if i == j || m[i] > m[j] {
        return 1.0;
    }
    let li = interaction_energy(m[i], tc[i], t);
    let lj = interaction_energy(m[j], tc[j], t);
    (-(lj - li) / (R * t)).exp()
}

/// The activity coefficients of a mixture, from the paraffin-wax Wilson model.
///
/// `mixture` carries each component's molar mass (kg/mol) and critical temperature (K),
/// both of which feed the `lambda_i` correlation; `Lambda_ij` is `1.0` for `i == j` or
/// the heavier first component, else `exp(-(lambda_j - lambda_i) / (R T))`. The activity
/// coefficient is `ln gamma_c = 1 - ln(sum_i x_i Lambda_c_i) - sum_i x_i Lambda_i_c /
/// sum_j x_j Lambda_i_j`. `x` is checked rather than renormalised; a supercritical
/// component yields NaN, as in NeqSim.
///
/// # Errors
/// * [`AzothError::PropertyUnavailable`] if a component carries no molar mass.
/// * [`AzothError::InvalidInput`] if `x` is not a composition of `mixture`.
/// * [`AzothError::OutOfRange`] if `T` is not positive.
///
/// # Example
/// ```
/// use azoth_eos::databank::mixture_of;
/// use azoth_eos::wilson_activity_coefficients;
///
/// let (mixture, _) = mixture_of(&["n-butane", "nc12"], None)?;
/// let r = wilson_activity_coefficients(&mixture, 298.15, &[0.5, 0.5])?;
/// assert!((r.gamma[0] - 1.213061319425015).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T` and `x` are the symbols in the chemistry
pub fn wilson_activity_coefficients(
    mixture: &Mixture,
    T: f64,
    x: &[f64],
) -> Result<WilsonActivityCoefficientsResult> {
    let spec = &model_gen::WILSON_ACTIVITY_COEFFICIENTS_SPEC;
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
    if mixture.len() != n {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "the mixture has {} components but `x` has {n} entries",
                mixture.len()
            ),
        ));
    }
    // The correlation is over carbon number and the fusion temperature, both of which
    // are masses; a component without one is refused rather than given a zero, which
    // would put a zero on the wrong side of the `M_i > M_j` branch as well as in the
    // energy.
    let components = mixture.components();
    let m: Vec<f64> = components
        .iter()
        .map(|c| {
            c.molar_mass.ok_or_else(|| {
                AzothError::property_unavailable(
                    "component",
                    "molar mass",
                    "the paraffin-wax Wilson correlation needs a molar mass for the \
                     carbon number, and a card-added component carries none",
                )
            })
        })
        .collect::<Result<_>>()?;
    let tc: Vec<f64> = components.iter().map(|c| c.tc.value).collect();
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

    let mut ln_gamma = vec![0.0; n];
    let mut gamma = vec![0.0; n];
    for c in 0..n {
        let mut s1 = 0.0;
        for (i, &xi) in x.iter().enumerate() {
            s1 += xi * char_energy(&m, &tc, T, c, i);
        }
        let mut s2 = 0.0;
        for (i, &xi) in x.iter().enumerate() {
            let mut temp = 0.0;
            for (j, &xj) in x.iter().enumerate() {
                temp += xj * char_energy(&m, &tc, T, i, j);
            }
            s2 += xi * char_energy(&m, &tc, T, i, c) / temp;
        }
        let lng = 1.0 - s1.ln() - s2;
        ln_gamma[c] = lng;
        gamma[c] = lng.exp();
    }

    Ok(WilsonActivityCoefficientsResult {
        ln_gamma,
        gamma,
        warnings,
    })
}
