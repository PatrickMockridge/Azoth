//! `eos.van_laar_acid_activity_coefficients` - the activity coefficients from the
//! Taleb-Ponche-Mirabel Van Laar model for the water/nitric/sulfuric acid system.
//!
//! Spec: `specs/models/eos/van_laar_acid_activity_coefficients.toml`. A *direct*
//! model: no iteration, so no `algorithm` block.

use azoth_core::{AzothError, Result, apply_checks};

use crate::databank::VanLaarAcidParameters;
use crate::model_gen;
use crate::results::VanLaarAcidActivityCoefficientsResult;

/// The activity a species the model does not cover is given, so a liquid phase rejects
/// it. NeqSim's `ComponentGEVanLaarAcid.NON_MODELED_COMPONENT_PENALTY`.
const NON_MODELED_COMPONENT_PENALTY: f64 = 1.0e12;

/// System I (H2O-HNO3) Van Laar exponent `B` for water. `B_I_2 = 1/B_I_1`.
const B_I_1: f64 = 0.5695;
/// System II (H2O-H2SO4) Van Laar exponent `B` for water. `B_II_3 = 1/B_II_1`.
const B_II_1: f64 = 0.527;
/// System III (HNO3-H2SO4) Van Laar exponent `B` for nitric acid. `B_III_3 = 1/B_III_2`.
const B_III_2: f64 = 0.4;
/// System III (HNO3-H2SO4) Van Laar parameter `A` for nitric acid, fitted at 273 K.
const A_III_2: f64 = -250.52;
/// System III (HNO3-H2SO4) Van Laar parameter `A` for sulfuric acid, fitted at 273 K.
const A_III_3: f64 = -100.21;

/// System I `A_I,1(T)`.
fn a_i1(t: f64) -> f64 {
    -391.43 - 7.44e4 / t
}

/// System I `A_I,2(T)`.
fn a_i2(t: f64) -> f64 {
    -627.739 - 1.406e5 / t
}

/// System II `A_II,1(T)`.
fn a_ii1(t: f64) -> f64 {
    2.989e3 - 2.147e6 / t + 2.33e8 / (t * t)
}

/// System II `A_II,3(T)`.
fn a_ii3(t: f64) -> f64 {
    5.672e3 - 4.074e6 / t + 4.421e8 / (t * t)
}

/// `T * log10(gamma_1)` for water: equation (10a) of Taleb et al. (1996).
fn t_log10_gamma_water(x1: f64, x2: f64, x3: f64, t: f64) -> f64 {
    let b_iii_2 = B_III_2;
    let b_i_1 = B_I_1;
    let b_ii_1 = B_II_1;
    let num1 = a_i1(t) * (x2 * x2 + b_iii_2 * x2 * x3) - 0.5 * A_III_2 * b_i_1 * x2 * x3;
    let den1 = b_i_1 * x1 + x2 + b_iii_2 * x3;
    let num2 = a_ii1(t) * (x3 * x3 + b_iii_2 * x2 * x3) - 0.5 * A_III_3 * b_ii_1 * x2 * x3;
    let den2 = b_ii_1 * x1 + b_iii_2 * x2 + x3;
    num1 / (den1 * den1) + num2 / (den2 * den2)
}

/// `T * log10(gamma_2)` for nitric acid: equation (10b).
fn t_log10_gamma_nitric(x1: f64, x2: f64, x3: f64, t: f64) -> f64 {
    let b_i_2 = 1.0 / B_I_1;
    let b_ii_3 = 1.0 / B_II_1;
    let b_iii_2 = B_III_2;
    let b_ii_1 = B_II_1;
    let num1 = a_i2(t) * (x1 * x1 + b_ii_3 * x1 * x3) - 0.5 * a_ii3(t) * b_i_2 * x1 * x3;
    let den1 = x1 + x2 * b_i_2 + b_ii_3 * x3;
    let num2 = A_III_3 * (x3 * x3 + b_ii_1 * x1 * x3) - 0.5 * a_ii1(t) * b_iii_2 * x1 * x3;
    let den2 = b_ii_1 * x1 + b_iii_2 * x2 + x3;
    num1 / (den1 * den1) + num2 / (den2 * den2)
}

/// `T * log10(gamma_3)` for sulfuric acid: equation (10c).
fn t_log10_gamma_sulfuric(x1: f64, x2: f64, x3: f64, t: f64) -> f64 {
    let b_i_2 = 1.0 / B_I_1;
    let b_ii_3 = 1.0 / B_II_1;
    let b_iii_3 = 1.0 / B_III_2;
    let b_i_1 = B_I_1;
    let num1 = a_ii3(t) * (x1 * x1 + b_i_2 * x1 * x2) - 0.5 * a_i2(t) * b_ii_3 * x1 * x2;
    let den1 = x1 + x2 * b_i_2 + b_ii_3 * x3;
    let num2 = A_III_2 * (x2 * x2 + b_i_1 * x1 * x2) - 0.5 * a_i1(t) * b_iii_3 * x1 * x2;
    let den2 = b_i_1 * x1 + x2 + b_iii_3 * x3;
    num1 / (den1 * den1) + num2 / (den2 * den2)
}

/// The activity coefficients of a mixture containing water, nitric acid and sulfuric
/// acid, from the Van Laar model of Taleb, Ponche and Mirabel (1996).
///
/// The three modelled species are evaluated on their own mole-fraction basis: the
/// fractions of the acids in `x` are renormalised to sum to one, so a dissolved carrier
/// gas does not enter the ternary expression. A component the model does not cover is
/// given `NON_MODELED_COMPONENT_PENALTY` (`1.0e12`) rather than being refused, which is NeqSim's
/// behaviour and pushes it out of the liquid phase.
///
/// All the model's logarithms are base-10, and `gamma = 10^(T log10 gamma / T)`.
/// `params` is resolved by name through [`crate::databank::van_laar_acid_parameters`];
/// `x` is checked rather than renormalised.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `params` and `x` disagree in length, or `x` is not
///   a composition.
/// * [`AzothError::OutOfRange`] if `T` is not positive.
///
/// # Example
/// ```
/// use azoth_eos::databank::van_laar_acid_parameters;
/// use azoth_eos::van_laar_acid_activity_coefficients;
///
/// let params = van_laar_acid_parameters(&["water", "hno3", "h2so4"], None)?;
/// let r = van_laar_acid_activity_coefficients(&params, 250.0, &[0.5, 0.3, 0.2])?;
/// assert!((r.gamma[0] - 0.008700227753081819).abs() < 1e-18);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T` and `x` are the symbols in the chemistry
pub fn van_laar_acid_activity_coefficients(
    params: &VanLaarAcidParameters,
    T: f64,
    x: &[f64],
) -> Result<VanLaarAcidActivityCoefficientsResult> {
    let spec = &model_gen::VAN_LAAR_ACID_ACTIVITY_COEFFICIENTS_SPEC;
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
    if params.acid_index.len() != n {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "the mixture has {} components but `x` has {n} entries",
                params.acid_index.len()
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

    // The ternary expression is defined on the three acids' own basis.
    //
    // Dividing through by the total is NeqSim's form rather than a step that moves the
    // answer: each of the three expressions below is homogeneous of degree zero - a
    // quadratic numerator over a squared linear denominator - so it takes the same value
    // on any scaling of the basis. Measured, a water/nitric mixture at `0.5/0.3` and one
    // at `0.625/0.375` give the same gammas, and `x` scaled by any factor does too. The
    // renormalisation is kept because it is what the port source does, and because it
    // makes the empty basis below a case rather than a division by zero.
    let mut acids = [0.0_f64; 3];
    for (&index, &fraction) in params.acid_index.iter().zip(x) {
        if (1..=3).contains(&index) {
            acids[usize::from(index) - 1] += fraction;
        }
    }
    let total: f64 = acids.iter().sum();
    // No acid at all leaves the basis undefined; NeqSim reads that as pure water, which
    // makes every acid's gamma its own infinite-dilution value in water.
    let [x1, x2, x3] = if total <= 0.0 {
        [1.0, 0.0, 0.0]
    } else {
        [acids[0] / total, acids[1] / total, acids[2] / total]
    };

    let ternary = [
        10.0_f64.powf(t_log10_gamma_water(x1, x2, x3, T) / T),
        10.0_f64.powf(t_log10_gamma_nitric(x1, x2, x3, T) / T),
        10.0_f64.powf(t_log10_gamma_sulfuric(x1, x2, x3, T) / T),
    ];

    let mut ln_gamma = vec![0.0; n];
    let mut gamma = vec![0.0; n];
    for (i, &index) in params.acid_index.iter().enumerate() {
        let value = match index {
            1 => ternary[0],
            2 => ternary[1],
            3 => ternary[2],
            _ => NON_MODELED_COMPONENT_PENALTY,
        };
        gamma[i] = value;
        ln_gamma[i] = value.ln();
    }

    Ok(VanLaarAcidActivityCoefficientsResult {
        ln_gamma,
        gamma,
        warnings,
    })
}
