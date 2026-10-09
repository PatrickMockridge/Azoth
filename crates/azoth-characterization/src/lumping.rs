//! `characterization.lumping` - a cut table grouped into equal-mass lumps.
//!
//! Spec: `specs/models/characterization/lumping.toml`. Oracle:
//! `validation/neqsim/captures/plus_fraction_probe.tsv`.
//!
//! # Equal mass, not equal cut counts
//!
//! `StandardLumpingModel.generateLumpedComposition` walks the cuts accumulating `z*M` and
//! closes a lump when the running sum reaches a target. The target is `W/N` to begin with and
//! `(W - accumulated)/(N - k - 1)` after each lump, so every lump but the last is aimed at the
//! same mass and the last takes whatever is left.
//!
//! # The total mass is the fluid's, not the table's
//!
//! `W` and the mole-fraction total are accumulated over the *system's* TBP and plus rows, which
//! for a fluid with one plus end is `z_plus` and `z_plus * M_plus`. The cut table's own
//! `sum(z*M)` carries the abundance solve's residual instead and differs by about `1e-11`, and on
//! the heavy rows that is enough to move a partition boundary - the last lump reads `0.041863`
//! from the fluid's numbers against `0.048901` from the table's. That is why `molar_mass` and
//! `mole_fraction` are inputs here and why they are the *plus fraction's* rather than the
//! table's.
//!
//! # What is not here
//!
//! The class ends by removing the plus component and adding one `addTBPfraction` per lump, which
//! is a mutation of a `SystemInterface`. This id stops at the table: a pseudo-component still
//! cannot enter a `Mixture`, so the caller is handed the vectors and nothing is added to a fluid.

use azoth_core::units::{MassDensity, MolarMass, kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_core::{AzothError, Result, apply_checks};

use crate::model_gen;
use crate::results::LumpingResult;

/// `LumpingModel.numberOfPseudocomponents`' own default.
const DEFAULT_LUMPS: usize = 7;

/// The denominator floor the class adds to the lump count, so a zero count would not divide.
const COUNT_FLOOR: f64 = 1.0e-10;

/// A cut table grouped into equal-mass lumps.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if an input is outside its declared range.
/// * [`AzothError::InvalidInput`] if the three vectors are not the same length, if that length is
///   zero, or if more lumps are asked for than there are cuts.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
/// use azoth_characterization::lumping;
///
/// let lumps = lumping(
///     kilograms_per_mole(0.4),
///     0.1,
///     &[0.04, 0.03, 0.02, 0.01],
///     &[
///         kilograms_per_mole(0.5),
///         kilograms_per_mole(0.3),
///         kilograms_per_mole(0.2),
///         kilograms_per_mole(0.1),
///     ],
///     &[
///         kilograms_per_cubic_meter(500.0),
///         kilograms_per_cubic_meter(300.0),
///         kilograms_per_cubic_meter(200.0),
///         kilograms_per_cubic_meter(100.0),
///     ],
///     Some(2),
/// )?;
/// assert_eq!(lumps.fraction_of_heavy_end.len(), 2);
/// assert!((lumps.fraction_of_heavy_end[0] - 0.4).abs() < 1e-15);
/// assert!((lumps.fraction_of_heavy_end[1] - 0.6).abs() < 1e-15);
/// assert!((lumps.lump_molar_mass[0].value - 0.5).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn lumping(
    molar_mass: MolarMass,
    mole_fraction: f64,
    cut_z: &[f64],
    cut_molar_mass: &[MolarMass],
    cut_density: &[MassDensity],
    number_of_lumps: Option<usize>,
) -> Result<LumpingResult> {
    let spec = &model_gen::LUMPING_SPEC;
    let mut warnings = Vec::new();

    let lumps = number_of_lumps.unwrap_or(DEFAULT_LUMPS);
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "molar_mass" => Some(molar_mass.value),
            "mole_fraction" => Some(mole_fraction),
            "number_of_lumps" => number_of_lumps.map(|n| n as f64),
            _ => None,
        },
        &mut warnings,
    )?;

    let count = cut_z.len();
    if count == 0 {
        return Err(AzothError::invalid_input(
            "cut_z",
            "a cut table to group is at least one cut; there is nothing to lump here",
        ));
    }
    if cut_molar_mass.len() != count || cut_density.len() != count {
        return Err(AzothError::invalid_input(
            "cut_molar_mass",
            format!(
                "the three vectors are one table and must be the same length: {} cuts against {} \
                 molar masses and {} densities",
                count,
                cut_molar_mass.len(),
                cut_density.len()
            ),
        ));
    }
    if lumps > count {
        return Err(AzothError::invalid_input(
            "number_of_lumps",
            format!(
                "{lumps} lumps from {count} cuts: the partition closes a lump only on a cut, so \
                 asking for more groups than there are cuts leaves a count the loop cannot reach"
            ),
        ));
    }

    // The two totals, which the class reads off the fluid rather than off the table.
    let weight_total = mole_fraction * molar_mass.value;
    let mole_fraction_total = mole_fraction;
    let mut target = weight_total / (lumps as f64 + COUNT_FLOOR);

    let weight = |index: usize| cut_z[index] * cut_molar_mass[index].value;
    let mut fraction_of_heavy_end = Vec::with_capacity(lumps);
    let mut lump_mole_fraction = Vec::with_capacity(lumps);
    let mut lump_molar_mass = Vec::with_capacity(lumps);
    let mut lump_density = Vec::with_capacity(lumps);

    let mut accumulated = 0.0;
    let mut running = 0.0;
    let mut lump_z = 0.0;
    let mut lump_weight = 0.0;
    let mut denominator = 0.0;
    let mut opened = 1_usize;
    for index in 0..count {
        let w = weight(index);
        running += w;
        accumulated += w;
        lump_z += cut_z[index];
        lump_weight += w;
        denominator += w / cut_density[index].value;

        // **The target moves before the counter does**, so the k-th lump divides by `lumps - k - 1`.
        // On the last one that denominator is zero and the division is an infinity the class never
        // reads, because the loop ends on the same cut.
        let closes = (running >= target && lumps != opened) || index == count - 1;
        if !closes {
            continue;
        }
        let remaining = lumps - opened;
        // `remaining` is zero only on the last of `lumps` lumps, which can only be closed by the
        // final cut - so the target set here is never compared against anything. The class's own
        // division by zero produces an infinity it never reads; this states one instead, and both
        // languages state the same one.
        target = if remaining == 0 {
            f64::INFINITY
        } else {
            (weight_total - accumulated) / remaining as f64
        };
        fraction_of_heavy_end.push(lump_z / mole_fraction_total);
        lump_mole_fraction.push(lump_z);
        lump_molar_mass.push(kilograms_per_mole(lump_weight / lump_z));
        // The gravity is the mass-weighted harmonic mean, in the same scale the table uses: a
        // pure rescaling of the densities leaves the ratio unchanged, so no unit boundary is
        // crossed here.
        lump_density.push(kilograms_per_cubic_meter(lump_weight / denominator));
        running = 0.0;
        lump_z = 0.0;
        lump_weight = 0.0;
        denominator = 0.0;
        opened += 1;
    }

    Ok(LumpingResult {
        fraction_of_heavy_end,
        lump_mole_fraction,
        lump_molar_mass,
        lump_density,
        warnings,
    })
}
