//! `characterization.tbp_cut_boundaries` - a cumulative TBP curve cut into an assay.
//!
//! Spec: `specs/models/characterization/tbp_cut_boundaries.toml`. Oracle:
//! `validation/neqsim/captures/tbp_cut_boundaries_probe.tsv`.
//!
//! # The curve is stated as boundaries
//!
//! `addTBPCutBoundariesKelvin` takes N cumulative yields, N boiling points and N-1 gravities, and
//! makes N-1 cuts of them. Each interval's share is `(cumulative[i + 1] - cumulative[i]) / 100`,
//! its gravity is the one stated beside it, and `withBoilingRangeKelvin(lower, upper)` keeps both
//! boundaries and the midpoint between them - so a cut carries a range and an average rather than
//! only an average.
//!
//! # The Celsius entry point is the same arithmetic
//!
//! `addTBPCutBoundariesCelsius` adds 273.15 to its temperatures and calls this one. That is the
//! crossing `to_si` makes at the boundary, so there is one id here and a caller states its
//! temperatures in whatever unit it holds them in.

use azoth_core::units::{ThermodynamicTemperature, kelvins};
use azoth_core::{AzothError, Result, apply_checks};

use crate::model_gen;
use crate::results::TbpCutBoundariesResult;

/// `OilAssayCharacterisation.PERCENT_TOLERANCE`, on the curve's two ends and on each yield.
const PERCENT_TOLERANCE: f64 = 1.0e-8;

/// A cumulative TBP curve cut into an assay, one row per interval.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the prefix is empty, if the three vectors are not N, N and
///   N-1 long, if the yields do not span zero to a hundred, if either the yields or the
///   temperatures fail to increase strictly, or if a gravity is not finite and positive.
///
/// # Example
/// ```
/// use azoth_characterization::tbp_cut_boundaries;
///
/// let cuts = tbp_cut_boundaries(
///     "TBP",
///     &[0.0, 50.0, 100.0],
///     &[320.0, 420.0, 520.0],
///     &[0.70, 0.80],
/// )?;
/// assert_eq!(cuts.cut_name, vec!["TBP1", "TBP2"]);
/// assert_eq!(cuts.volume_fraction, vec![0.5, 0.5]);
/// // The midpoint of the pair that bounds it, not a boundary.
/// assert!((cuts.average_boiling_point[1].value - 470.0).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn tbp_cut_boundaries(
    name_prefix: &str,
    cumulative_volume_percent: &[f64],
    boiling_point: &[f64],
    specific_gravity: &[f64],
) -> Result<TbpCutBoundariesResult> {
    let spec = &model_gen::TBP_CUT_BOUNDARIES_SPEC;
    let mut warnings = Vec::new();
    apply_checks(spec.input_checks(), |_quantity| None, &mut warnings)?;

    if name_prefix.trim().is_empty() {
        return Err(AzothError::invalid_input(
            "name_prefix",
            "a cut's name is its prefix and its position, so an empty prefix names nothing",
        ));
    }
    let count = cumulative_volume_percent.len();
    if count < 2
        || boiling_point.len() != count
        || specific_gravity.len() != count - 1
    {
        return Err(AzothError::invalid_input(
            "cumulative_volume_percent",
            format!(
                "a TBP curve is N cumulative yields, N temperatures and N-1 gravities: {} yields, \
                 {} temperatures and {} gravities",
                count,
                boiling_point.len(),
                specific_gravity.len()
            ),
        ));
    }
    if cumulative_volume_percent[0].abs() > PERCENT_TOLERANCE
        || (cumulative_volume_percent[count - 1] - 100.0).abs() > PERCENT_TOLERANCE
    {
        return Err(AzothError::invalid_input(
            "cumulative_volume_percent",
            format!(
                "a cumulative volume yield spans zero to a hundred percent, and this curve runs \
                 from {} to {}",
                cumulative_volume_percent[0],
                cumulative_volume_percent[count - 1]
            ),
        ));
    }
    for (index, yield_percent) in cumulative_volume_percent.iter().enumerate() {
        if !yield_percent.is_finite()
            || *yield_percent < -PERCENT_TOLERANCE
            || *yield_percent > 100.0 + PERCENT_TOLERANCE
        {
            return Err(AzothError::invalid_input(
                format!("cumulative_volume_percent[{index}]"),
                format!("a cumulative yield is finite and between 0 and 100, and this is {yield_percent}"),
            ));
        }
        // The class's own `!(yield > previous)`, written out: a NaN and a non-increasing pair are
        // both refused, which is what "must strictly increase" means.
        if index > 0
            && (yield_percent.is_nan() || *yield_percent <= cumulative_volume_percent[index - 1])
        {
            return Err(AzothError::invalid_input(
                format!("cumulative_volume_percent[{index}]"),
                format!("the cumulative yields must strictly increase, and {yield_percent} follows {}", cumulative_volume_percent[index - 1]),
            ));
        }
    }
    for (index, temperature) in boiling_point.iter().enumerate() {
        if !temperature.is_finite() || *temperature <= 0.0 {
            return Err(AzothError::invalid_input(
                format!("boiling_point[{index}]"),
                format!("a boiling-point boundary is finite and positive, and this is {temperature}"),
            ));
        }
        if index > 0 && *temperature <= boiling_point[index - 1] {
            return Err(AzothError::invalid_input(
                format!("boiling_point[{index}]"),
                format!("the boiling-point boundaries must strictly increase, and {temperature} follows {}", boiling_point[index - 1]),
            ));
        }
    }
    for (index, gravity) in specific_gravity.iter().enumerate() {
        if !gravity.is_finite() || *gravity <= 0.0 {
            return Err(AzothError::invalid_input(
                format!("specific_gravity[{index}]"),
                format!("a specific gravity is finite and positive, and this is {gravity}"),
            ));
        }
    }

    let cut_name = (0..specific_gravity.len())
        .map(|index| format!("{name_prefix}{}", index + 1))
        .collect();
    let volume_fraction = cumulative_volume_percent
        .windows(2)
        .map(|pair| (pair[1] - pair[0]) / 100.0)
        .collect();
    let lower_boiling_point: Vec<ThermodynamicTemperature> =
        boiling_point[..boiling_point.len() - 1]
        .iter()
        .map(|value| kelvins(*value))
        .collect();
    let upper_boiling_point: Vec<ThermodynamicTemperature> = boiling_point[1..]
        .iter()
        .map(|value| kelvins(*value))
        .collect();
    let average_boiling_point = lower_boiling_point
        .iter()
        .zip(&upper_boiling_point)
        .map(|(lower, upper)| kelvins(0.5 * (lower.value + upper.value)))
        .collect();

    Ok(TbpCutBoundariesResult {
        cut_name,
        volume_fraction,
        specific_gravity: specific_gravity.to_vec(),
        lower_boiling_point,
        upper_boiling_point,
        average_boiling_point,
        warnings,
    })
}
