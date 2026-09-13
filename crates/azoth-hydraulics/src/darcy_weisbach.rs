//! `hydraulics.darcy_weisbach` - pressure drop in a straight pipe.
//!
//! ```text
//! dP = f * (L / D) * (rho * v**2 / 2)
//! ```
//!
//! Spec: `specs/calcs/hydraulics/darcy_weisbach.yaml`
//!
//! # A note on the sign convention
//!
//! `dP` is a *drop* and is returned positive. The direction of flow is not
//! modelled: velocity is taken as a speed, and a negative input is rejected
//! rather than interpreted as reversed flow.
//!
//! # Fittings are not included
//!
//! This is straight pipe only. Fitting losses are a separate calc
//! (`crane_k_factors`) computed by a different method, and combining them is a
//! modelling decision the caller should make deliberately rather than have
//! applied invisibly. The `azoth pipe` CLI performs that composition and shows
//! both contributions separately.

use crate::results::DarcyWeisbachResult;
use crate::reynolds_number::regime_for;
use crate::spec_gen;
use azoth_core::units::{DynamicViscosity, Length, MassDensity, Velocity, pascals};
use azoth_core::{Result, Warning, WarningCode, apply_checks};

/// Pressure drop over a length of straight pipe.
///
/// `mu` is optional, and the reason is worth stating. It is not needed for the
/// pressure drop itself: `f` already carries the flow information. It is needed
/// to *check* the flow regime, because the Reynolds number requires viscosity.
///
/// When `mu` is supplied the result reports `re` and `regime`, and a
/// transitional flow produces a `TransitionalFlow` warning. When it is omitted
/// the result reports neither, and carries a
/// [`WarningCode::RangeCheckSkipped`] warning saying the regime went unchecked.
///
/// That last part is the point. "Checked and fine" and "never checked" must not
/// look the same to a caller, and an omitted optional input is exactly how those
/// two states would otherwise become indistinguishable.
///
/// # Errors
/// [`azoth_core::AzothError::OutOfRange`] if `f`, `L`, `D`, `rho` or `v`
/// violates a hard bound in the spec - a non-positive friction factor, diameter
/// or density, a negative length or velocity.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_cubic_meter, meters, meters_per_second};
/// use azoth_hydraulics::darcy_weisbach;
/// use azoth_core::CalcResult; // for has_warning
///
/// // Without viscosity: pressure drop only, regime unchecked.
/// let r = darcy_weisbach(
///     0.02,
///     meters(100.0),
///     meters(0.1),
///     kilograms_per_cubic_meter(998.0),
///     meters_per_second(1.5),
///     None,
/// )?;
/// assert!((r.dp.value - 22455.0).abs() < 1e-9);
/// assert!(r.re.is_none());
/// assert!(r.has_warning(azoth_core::WarningCode::RangeCheckSkipped));
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `L` and `D` are the symbols in the published equation
pub fn darcy_weisbach(
    f: f64,
    L: Length,
    D: Length,
    rho: MassDensity,
    v: Velocity,
    mu: Option<DynamicViscosity>,
) -> Result<DarcyWeisbachResult> {
    let spec = &spec_gen::DARCY_WEISBACH_SPEC;
    let mut warnings = Vec::new();

    let (l_v, d_v, rho_v, v_v) = (L.value, D.value, rho.value, v.value);

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "f" => Some(f),
            "L" => Some(l_v),
            "D" => Some(d_v),
            "rho" => Some(rho_v),
            "v" => Some(v_v),
            _ => None,
        },
        &mut warnings,
    )?;

    // Guarded by the checks above, so D is non-zero here.
    let dp = f * (l_v / d_v) * (rho_v * v_v * v_v / 2.0);

    // Reynolds number only when viscosity was supplied. Everything downstream
    // of this is `Option`, which is what makes "unchecked" visible in the type
    // rather than only in the warnings.
    let (re, regime) = match mu {
        Some(mu) => {
            let re = rho_v * v_v * d_v / mu.value;
            (Some(re), Some(regime_for(re)))
        }
        None => (None, None),
    };

    // `re` resolves to None when mu was omitted, so the regime check emits
    // RangeCheckSkipped rather than quietly passing.
    apply_checks(
        spec.derived_checks(),
        |quantity| match quantity {
            "re" => re,
            "dp" => Some(dp),
            _ => None,
        },
        &mut warnings,
    )?;

    // The spec's transitional band check carries TRANSITIONAL_FLOW, so it has
    // already warned if applicable. Nothing extra is needed here - adding a
    // second warning for the same condition would just train callers to ignore
    // them.

    Ok(DarcyWeisbachResult {
        dp: pascals(dp),
        f,
        re,
        regime,
        warnings,
    })
}

/// Add a fitting loss to a straight-pipe pressure drop.
///
/// `dP_fittings = K * (rho * v**2 / 2)`, the velocity-head form of the same
/// minor loss the equivalent-length method expresses as added pipe length.
///
/// Provided as a named function rather than left to callers because the two
/// terms are computed by different methods and combined in a specific way: the
/// straight-pipe term uses `f * L/D` and the fitting term uses `K`. Presenting
/// the composition explicitly is what lets the CLI show both contributions
/// instead of a single unexplained total.
#[must_use]
pub fn add_fitting_loss(dp_straight: f64, k_total: f64, rho: f64, v: f64) -> f64 {
    dp_straight + k_total * (rho * v * v / 2.0)
}

/// A warning to attach when a composed result includes fitting losses based on
/// estimated data.
///
/// Takes the warning the `crane_k_factors` result already carries, so the
/// provenance is propagated rather than re-derived. Returns `None` when there is
/// nothing to say.
#[must_use]
pub fn propagate_estimated_data(warnings: &[Warning]) -> Option<Warning> {
    warnings
        .iter()
        .find(|w| w.code == WarningCode::EstimatedData)
        .cloned()
}
