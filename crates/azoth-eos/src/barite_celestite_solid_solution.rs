//! `eos.barite_celestite_solid_solution` - NeqSim's `BariteCelestiteSolidSolution`.
//!
//! A Margules one-parameter solid between BaSO4 and SrSO4, with the solid composition found
//! by the class's bisection and the reported index taken from the Ba arm alone.

use azoth_core::{Result, apply_checks};

use crate::model_gen;
use crate::results::BariteCelestiteSolidSolutionResult;

/// The solid's composition and the class's saturation index.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's range.
///
/// # Example
/// ```
/// use azoth_eos::barite_celestite_solid_solution::barite_celestite_solid_solution;
///
/// let r = barite_celestite_solid_solution(2.3, 0.01, 0.0, 0.05, 1.08e-10, 3.44e-7)?;
/// assert_eq!(r.x_baso4, 1.0);
/// assert!((r.total_saturation_index.expect("the pure end member has an index")
///     - 6.665546248849069).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn barite_celestite_solid_solution(
    margules: f64,
    a_ba: f64,
    a_sr: f64,
    a_so4: f64,
    ksp_baso4: f64,
    ksp_srso4: f64,
) -> Result<BariteCelestiteSolidSolutionResult> {
    let spec = &model_gen::BARITE_CELESTITE_SOLID_SOLUTION_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "margules" => Some(margules),
            "a_ba" => Some(a_ba),
            "a_sr" => Some(a_sr),
            "a_so4" => Some(a_so4),
            "ksp_baso4" => Some(ksp_baso4),
            "ksp_srso4" => Some(ksp_srso4),
            _ => None,
        },
        &mut warnings,
    )?;

    let pure_index = |iap: f64, ksp: f64| {
        if iap > 0.0 && ksp > 0.0 {
            Some((iap / ksp).log10())
        } else {
            None
        }
    };

    let (x_ba, total) = if a_ba <= 0.0 && a_sr <= 0.0 {
        (0.5, None)
    } else if a_sr <= 1.0e-30 {
        (1.0, pure_index(a_ba * a_so4, ksp_baso4))
    } else if a_ba <= 1.0e-30 {
        (0.0, pure_index(a_sr * a_so4, ksp_srso4))
    } else {
        let objective = |x: f64| {
            let x_sr = 1.0 - x;
            let gamma_ba = solid_activity_coefficient(margules, x);
            let gamma_sr = solid_activity_coefficient(margules, x_sr);
            let r_ba = a_ba * a_so4 / (ksp_baso4 * x * gamma_ba);
            let r_sr = a_sr * a_so4 / (ksp_srso4 * x_sr * gamma_sr);
            r_ba - r_sr
        };
        let mut x_lo = 0.001;
        let mut x_hi = 0.999;
        let mut x = 0.5;
        for _ in 0..100 {
            let x_mid = 0.5 * (x_lo + x_hi);
            let f_mid = objective(x_mid);
            if f_mid.abs() < 1.0e-12 || (x_hi - x_lo) < 1.0e-10 {
                x = x_mid;
                break;
            }
            // The class re-evaluates the lower end every iteration, as written.
            if f_mid * objective(x_lo) < 0.0 {
                x_hi = x_mid;
            } else {
                x_lo = x_mid;
            }
            x = x_mid;
        }
        let gamma_ba = solid_activity_coefficient(margules, x);
        let si_ba = (a_ba * a_so4) / (ksp_baso4 * x * gamma_ba);
        (x, Some(si_ba.max(1.0e-30).log10()))
    };

    let supersaturated = total.is_some_and(|value| value > 0.0);
    Ok(BariteCelestiteSolidSolutionResult {
        x_baso4: x_ba,
        x_srso4: 1.0 - x_ba,
        total_saturation_index: total,
        supersaturated,
        warnings,
    })
}

/// `exp(W (1 - x)^2)`, the Margules one-parameter coefficient.
fn solid_activity_coefficient(margules: f64, x: f64) -> f64 {
    let other = 1.0 - x;
    (margules * other * other).exp()
}
