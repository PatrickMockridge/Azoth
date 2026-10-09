//! `characterization.tbp_grouping` - a phase's components binned by normal boiling point.
//!
//! Spec: `specs/models/characterization/tbp_grouping.toml`. Oracle:
//! `validation/neqsim/captures/tbp_grouping_probe.tsv`.
//!
//! # One name, three implementations
//!
//! `groupTBPfractions` is declared once and implemented three times. `PlusCharacterize`'s
//! returns `true` and computes nothing. `TBPCharacterize`'s does the binning and is reachable
//! only from inside its own package, because `Characterise.TBPCharacterise` has no accessor.
//! `Phase`'s is the one ported: public on a public class, though **not** on the
//! `PhaseInterface` that `SystemInterface.getPhase` returns, so a caller needs a cast.
//!
//! Nothing in NeqSim calls any of the three, and nothing in its tests either.
//!
//! # The bins are boiling points, not carbon numbers
//!
//! Fourteen thresholds in degrees Celsius, the lowest at 69.2, each opening a bin indexed from
//! six. Bins 0 to 5 exist in the returned array and no component can reach them, a component
//! below 69.2 reaches no bin at all, and two components between one pair of thresholds sum into
//! one bin - so the twenty entries do not add up to one and are not meant to.
//!
//! The comparison is on `getNormalBoilingPoint("C")`, which is the stored value less 273.15. This
//! takes the stored kelvin and subtracts the same 273.15, so the comparison sees the same number
//! rather than one round-tripped through a display unit.

use azoth_core::units::ThermodynamicTemperature;
use azoth_core::{AzothError, Result, apply_checks};

use crate::model_gen;
use crate::results::TbpGroupingResult;

/// The bins the answer always carries, of which only these last fourteen can be filled.
const BIN_COUNT: usize = 20;

/// The fourteen thresholds, in degrees Celsius, paired with the bin each opens.
///
/// Descending, because that is the order the class tests them in and a component above every
/// threshold would otherwise fall through to the last branch.
const THRESHOLDS: [(f64, usize); 14] = [
    (331.0, 19),
    (317.0, 18),
    (303.0, 17),
    (287.0, 16),
    (271.1, 15),
    (253.9, 14),
    (235.9, 13),
    (216.8, 12),
    (196.4, 11),
    (174.6, 10),
    (151.3, 9),
    (126.1, 8),
    (98.9, 7),
    (69.2, 6),
];

/// The kelvin-to-Celsius step `getNormalBoilingPoint("C")` makes.
const ABSOLUTE_ZERO_CELSIUS: f64 = 273.15;

/// A phase's components grouped into boiling-point bins.
///
/// # Errors
/// [`AzothError::InvalidInput`] if the two vectors are not the same length.
///
/// # Example
/// ```
/// use azoth_core::units::kelvins;
/// use azoth_characterization::tbp_grouping;
///
/// // n-heptane and n-octane, whose boiling points fall either side of the 98.9 C threshold.
/// let grouped = tbp_grouping(&[kelvins(371.6), kelvins(398.8)], &[0.4, 0.6])?;
/// assert_eq!(grouped.group_fraction.len(), 20);
/// assert!((grouped.group_fraction[6] - 0.4).abs() < 1e-15);
/// assert!((grouped.group_fraction[7] - 0.6).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn tbp_grouping(
    boiling_point: &[ThermodynamicTemperature],
    mole_fraction: &[f64],
) -> Result<TbpGroupingResult> {
    let spec = &model_gen::TBP_GROUPING_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        // The spec declares no range on either input, which is the house rule for a vector: a
        // bound on a vector is a bound on each of its entries, and `apply_checks` takes one value
        // per quantity. The length agreement below is the check this pair actually needs.
        |_quantity| None,
        &mut warnings,
    )?;

    if boiling_point.len() != mole_fraction.len() {
        return Err(AzothError::invalid_input(
            "mole_fraction",
            format!(
                "the two vectors are one component list and must be the same length: {} boiling \
                 points against {} mole fractions",
                boiling_point.len(),
                mole_fraction.len()
            ),
        ));
    }

    let mut group_fraction = vec![0.0; BIN_COUNT];
    for (temperature, fraction) in boiling_point.iter().zip(mole_fraction) {
        let celsius = temperature.value - ABSOLUTE_ZERO_CELSIUS;
        for (threshold, bin) in THRESHOLDS {
            if celsius >= threshold {
                group_fraction[bin] += fraction;
                break;
            }
        }
    }

    Ok(TbpGroupingResult {
        group_fraction,
        warnings,
    })
}