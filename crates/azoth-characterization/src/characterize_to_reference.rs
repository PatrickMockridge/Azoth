//! `characterization.characterize_to_reference` - one fluid re-cut onto another's slate.
//!
//! Spec: `specs/models/characterization/characterize_to_reference.toml`. Oracle:
//! `validation/neqsim/captures/pseudo_component_combiner_probe.tsv`.
//!
//! # Read out of the middle of `PseudoComponentCombiner`
//!
//! The class is 1841 lines and every public method takes or returns a `SystemInterface`. This is
//! the part of it that is a function of two tables: `extractComponents` lifts each fluid's
//! pseudo-components into rows keyed by a boiling point that falls back to the molar mass,
//! `determineReferenceBoundaries` puts a cut at the midpoint of each adjacent pair of the
//! *reference's* keys, and `distributeToProfiles` walks the source binning each row into the group
//! its key falls in.
//!
//! A group's molar mass is its accumulated mass over its accumulated moles; its density is that
//! mass over the volume the rows would occupy, `sum(m) / sum(m/rho)` - the mass-weighted harmonic
//! mean of the gravities, which is neither their mean nor a mole-weighted one.
//!
//! # The empty groups are dropped
//!
//! A group whose accumulated mass or moles is at or below `1e-12` produces no profile and the
//! caller skips it, so a source entirely below the reference's first cut resolves to one component
//! against a reference of many. [`CharacterizeToReferenceResult::reference_index`] is what says
//! which cut that one is.
//!
//! # What is not here
//!
//! `combineReservoirFluids`, `characterizeToCommonSlate`, the delumping, `inheritReferenceProperties`,
//! `transferBinaryInteractionParameters`, `normalizeComposition` and `generateValidationReport`.
//! They are blocked on the card width - a pseudo-component still cannot enter a `Mixture` - and on
//! the `record_list` input kind their twenty-field rows would need.

use azoth_core::units::{MassDensity, MolarMass, kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_core::{AzothError, Result, apply_checks};

use crate::model_gen;
use crate::results::CharacterizeToReferenceResult;

/// `PseudoComponentCombiner.MASS_TOLERANCE`, on the accumulated mass and on a group's validity.
const MASS_TOLERANCE: f64 = 1.0e-12;

/// One row of a fluid's pseudo-component table, as `extractComponents` lifts it.
struct Row {
    /// Where it sits on the grid.
    key: f64,
    /// Its mole amount, and its mass, which are the two the group accumulates.
    moles: f64,
    /// Its gravity, which the group's volume is taken through.
    density: f64,
    mass: f64,
}

/// `!(value > MASS_TOLERANCE)`, NaN included: the class's own comparison drops an accumulation it
/// cannot order rather than reporting it, where `value <= MASS_TOLERANCE` would keep the NaN.
fn below_tolerance(value: f64) -> bool {
    value.is_nan() || value <= MASS_TOLERANCE
}

/// A group's running totals, as `PseudoComponentGroupBuilder` keeps them.
#[derive(Default)]
struct Group {
    mass: f64,
    moles: f64,
    volume: f64,
    density_mass: f64,
}

/// `PseudoComponentContribution.sortingKey`: the boiling point, or the molar mass without one.
fn sorting_key(boiling_point: f64, molar_mass: f64) -> f64 {
    if boiling_point.is_finite() && boiling_point > 0.0 {
        boiling_point
    } else {
        molar_mass
    }
}

/// A fluid's rows, in ascending sorting key, as `extractComponents` returns them.
///
/// A row with a non-positive mole amount is skipped rather than carried: `extractComponents` drops
/// it, so a zero here is a row that is not on the grid at all.
fn rows(
    moles: &[f64],
    molar_mass: &[MolarMass],
    density: &[MassDensity],
    boiling_point: &[f64],
    field: &str,
) -> Result<Vec<Row>> {
    let count = moles.len();
    if molar_mass.len() != count || density.len() != count || boiling_point.len() != count {
        return Err(AzothError::invalid_input(
            field,
            format!(
                "a fluid's table is four aligned vectors: {count} mole amounts against {} molar \
                 masses, {} densities and {} boiling points",
                molar_mass.len(),
                density.len(),
                boiling_point.len()
            ),
        ));
    }
    let mut rows: Vec<Row> = moles
        .iter()
        .zip(molar_mass)
        .zip(density)
        .zip(boiling_point)
        .filter(|(((amount, _), _), _)| **amount > 0.0)
        .map(|(((amount, mass), rho), tb)| Row {
            key: sorting_key(*tb, mass.value),
            moles: *amount,
            density: rho.value,
            mass: amount * mass.value,
        })
        .collect();
    rows.sort_by(|left, right| left.key.total_cmp(&right.key));
    Ok(rows)
}

/// A fluid's pseudo-components re-cut onto another fluid's slate.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if any of the seven vectors is not aligned with its own table.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
/// use azoth_characterization::characterize_to_reference;
///
/// // Two reference cuts, so one boundary at 550 K, and one source cut either side of it.
/// let grouped = characterize_to_reference(
///     &[1.0, 1.0],
///     &[kilograms_per_mole(0.1), kilograms_per_mole(0.2)],
///     &[
///         kilograms_per_cubic_meter(700.0),
///         kilograms_per_cubic_meter(820.0),
///     ],
///     &[380.0, 628.0],
///     &[kilograms_per_mole(0.15), kilograms_per_mole(0.25)],
///     &[400.0, 700.0],
/// )?;
/// assert_eq!(grouped.reference_index, vec![0.0, 1.0]);
/// assert_eq!(grouped.group_moles, vec![1.0, 1.0]);
/// assert!((grouped.group_molar_mass[0].value - 0.1).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)]
pub fn characterize_to_reference(
    source_moles: &[f64],
    source_molar_mass: &[MolarMass],
    source_density: &[MassDensity],
    source_boiling_point: &[f64],
    reference_molar_mass: &[MolarMass],
    reference_boiling_point: &[f64],
) -> Result<CharacterizeToReferenceResult> {
    let spec = &model_gen::CHARACTERIZE_TO_REFERENCE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(spec.input_checks(), |_quantity| None, &mut warnings)?;

    let source = rows(
        source_moles,
        source_molar_mass,
        source_density,
        source_boiling_point,
        "source_moles",
    )?;
    if reference_molar_mass.len() != reference_boiling_point.len() {
        return Err(AzothError::invalid_input(
            "reference_molar_mass",
            format!(
                "the reference is two aligned vectors: {} molar masses against {} boiling points",
                reference_molar_mass.len(),
                reference_boiling_point.len()
            ),
        ));
    }
    let mut reference: Vec<f64> = reference_molar_mass
        .iter()
        .zip(reference_boiling_point)
        .map(|(mass, tb)| sorting_key(*tb, mass.value))
        .collect();
    reference.sort_by(f64::total_cmp);

    // `determineReferenceBoundaries`: the midpoint of each adjacent pair, or the larger of the two
    // when either is not finite.
    let boundaries: Vec<f64> = reference
        .windows(2)
        .map(|pair| {
            let (low, high) = (pair[0], pair[1]);
            if low.is_finite() && high.is_finite() {
                0.5 * (low + high)
            } else {
                low.max(high)
            }
        })
        .collect();

    // `distributeToProfiles`: a forward walk, advancing while the row's key is **greater than** the
    // boundary, so a key exactly on one stays in the lower group.
    let mut groups: Vec<Group> = (0..reference.len()).map(|_| Group::default()).collect();
    let mut index = 0_usize;
    let mut boundary = boundaries.first().copied().unwrap_or(f64::INFINITY);
    for row in &source {
        while index + 1 < reference.len() && row.key > boundary {
            index += 1;
            boundary = boundaries.get(index).copied().unwrap_or(f64::INFINITY);
        }
        let group = &mut groups[index];
        group.mass += row.mass;
        group.moles += row.moles;
        if row.density > 0.0 {
            group.volume += row.mass / row.density;
            group.density_mass += row.mass * row.density;
        }
    }

    // `buildWithOverrides`: a group below the tolerance is empty rather than reported empty.
    let mut reference_index = Vec::new();
    let mut group_moles = Vec::new();
    let mut group_molar_mass = Vec::new();
    let mut group_density = Vec::new();
    for (index, group) in groups.iter().enumerate() {
        if below_tolerance(group.mass) || below_tolerance(group.moles) {
            continue;
        }
        reference_index.push(index as f64);
        group_moles.push(group.moles);
        group_molar_mass.push(kilograms_per_mole(group.mass / group.moles));
        group_density.push(kilograms_per_cubic_meter(if group.volume > MASS_TOLERANCE {
            group.mass / group.volume
        } else {
            group.density_mass / group.mass
        }));
    }

    Ok(CharacterizeToReferenceResult {
        reference_index,
        group_moles,
        group_molar_mass,
        group_density,
        warnings,
    })
}
