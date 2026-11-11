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
//! # The rows are records
//!
//! `PseudoComponentContribution` carries twenty-one fields and this id reads four of them. Two are
//! optional *per row* - a boiling point, without which the row sorts by its molar mass, and a
//! density, without which it contributes to no volume - and a parallel vector can say "absent" per
//! vector and not per entry. The row this id takes is the fields it reads; the other seventeen are
//! owed with the methods below.
//!
//! # The class's own options, less the two that need a fluid
//!
//! [`characterize_to_reference`] carries four of `CharacterizationOptions`' switches:
//! `inheritReferenceProperties`, `delumpBeforeRecharacterization` with its `delumpResolution`, and
//! `sharedImaginaryBoundaries`. The first reads the reference's densities, the last its mole
//! amounts, and the middle is [`crate::delumping`]'s arithmetic run over the source first.
//!
//! # The class's other methods are their own ids
//!
//! `combineReservoirFluids` and `characterizeToCommonSlate` weight several fluids before cutting,
//! and `transferBinaryInteractionParameters` moves their interaction tables;
//! [`crate::combine_reservoir_fluids`], [`crate::characterize_to_common_slate`] and
//! [`crate::transfer_binary_interaction_parameters`] are those, over rows rather than a
//! `SystemInterface`. What is not here is the two report-like methods: `normalizeComposition`,
//! which computes a total and never uses it, and `generateValidationReport`, which reports on
//! data.

use azoth_core::units::{MassDensity, MolarMass, kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_core::{AzothError, Result, apply_checks};

use crate::delumping::{Lump, equal_mass_boundaries, sorting_key, split};
use crate::model_gen;
use crate::results::CharacterizeToReferenceResult;

/// `PseudoComponentCombiner.MASS_TOLERANCE`, on the accumulated mass and on a group's validity.
const MASS_TOLERANCE: f64 = 1.0e-12;

/// `CharacterizationOptions.Builder`'s own `delumpResolution`, which an absent input takes.
const DELUMP_RESOLUTION: i64 = 12;

/// One pseudo-component of the source fluid, as this id takes it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SourceRow {
    /// Its mole amount. Only a positive one is read: `extractComponents` skips the rest, so a zero
    /// is a row that is not on the grid at all.
    pub moles: f64,
    /// Its molar mass, which is its sorting key when it states no boiling point.
    pub molar_mass: MolarMass,
    /// Its normal liquid density, which the group's volume is the mass over. A row that states
    /// none contributes to neither the volume nor the density-mass sum.
    pub density: Option<MassDensity>,
    /// Its normal boiling point. A row that states none is sorted by its molar mass.
    pub boiling_point: Option<f64>,
}

/// One cut of the reference fluid's slate.
///
/// The reference decides where the cuts are and nothing else: its own properties are never a
/// group's. Two of its fields are read by one variant each - `density` where
/// `inherit_reference_properties` is set, `moles` where `shared_imaginary_boundaries` is - and by
/// neither otherwise.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReferenceRow {
    /// Its molar mass, the fallback sorting key - and the group's own molar mass where the
    /// inheriting variant is asked for.
    pub molar_mass: MolarMass,
    /// Its normal boiling point, and with it the boundaries the source is cut at.
    pub boiling_point: Option<f64>,
    /// Its normal liquid density, read only by the inheriting variant.
    pub density: Option<MassDensity>,
    /// Its mole amount, read only by the shared-boundaries variant - the equal-mass cuts are
    /// placed on the masses the slate carries, so a row that states none is refused there.
    pub moles: Option<f64>,
}

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

/// The source table split by `characterization.delumping`'s arithmetic, its names dropped.
///
/// **`delumpContributions` names each sub-fraction after its parent, and this walk reads no name.**
/// The row it takes has none, so the names are the one thing the split cannot carry across - which
/// the spec states, and which `characterization.delumping` answers for a caller that wants them.
fn expand(source: &[SourceRow], resolution: i64) -> Vec<SourceRow> {
    let lumps: Vec<Lump> = source
        .iter()
        .map(|row| Lump {
            name: String::new(),
            moles: row.moles,
            molar_mass: row.molar_mass,
            density: row.density,
            boiling_point: row.boiling_point,
        })
        .collect();
    split(&lumps, resolution)
        .into_iter()
        .map(|fraction| SourceRow {
            moles: fraction.moles,
            molar_mass: fraction.molar_mass,
            density: fraction.density,
            boiling_point: fraction.boiling_point,
        })
        .collect()
}

/// The source's rows, in ascending sorting key, as `extractComponents` returns them.
fn rows(source: &[SourceRow]) -> Vec<Row> {
    let mut rows: Vec<Row> = source
        .iter()
        .filter(|row| row.moles > 0.0)
        .map(|row| Row {
            key: sorting_key(row.boiling_point, row.molar_mass.value),
            moles: row.moles,
            density: row.density.map_or(0.0, |rho| rho.value),
            mass: row.moles * row.molar_mass.value,
        })
        .collect();
    rows.sort_by(|left, right| left.key.total_cmp(&right.key));
    rows
}

/// A fluid's pseudo-components re-cut onto another fluid's slate.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if a declared input falls outside its range. There is
///   no alignment to refuse: the rows are records, and a row is complete by being one.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
/// use azoth_characterization::{
///     ReferenceRow, SourceRow, characterize_to_reference,
/// };
///
/// // Two reference cuts, so one boundary at 550 K, and one source cut either side of it.
/// let grouped = characterize_to_reference(
///     &[
///         SourceRow {
///             moles: 1.0,
///             molar_mass: kilograms_per_mole(0.1),
///             density: Some(kilograms_per_cubic_meter(700.0)),
///             boiling_point: Some(380.0),
///         },
///         SourceRow {
///             moles: 1.0,
///             molar_mass: kilograms_per_mole(0.2),
///             density: Some(kilograms_per_cubic_meter(820.0)),
///             boiling_point: Some(628.0),
///         },
///     ],
///     &[
///         ReferenceRow {
///             molar_mass: kilograms_per_mole(0.15),
///             boiling_point: Some(400.0),
///             density: None,
///             moles: None,
///         },
///         ReferenceRow {
///             molar_mass: kilograms_per_mole(0.25),
///             boiling_point: Some(700.0),
///             density: None,
///             moles: None,
///         },
///     ],
///     None,
///     None,
///     None,
///     None,
/// )?;
/// assert_eq!(grouped.reference_index, vec![0.0, 1.0]);
/// assert_eq!(grouped.group_moles, vec![1.0, 1.0]);
/// assert!((grouped.group_molar_mass[0].value - 0.1).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn characterize_to_reference(
    source: &[SourceRow],
    reference: &[ReferenceRow],
    inherit_reference_properties: Option<bool>,
    delump_before_recharacterization: Option<bool>,
    delump_resolution: Option<i64>,
    shared_imaginary_boundaries: Option<bool>,
) -> Result<CharacterizeToReferenceResult> {
    let spec = &model_gen::CHARACTERIZE_TO_REFERENCE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(spec.input_checks(), |_quantity| None, &mut warnings)?;
    // The class's own default entry point passes the flag false, so an absent one is false.
    let inherit = inherit_reference_properties.unwrap_or(false);
    let delump = delump_before_recharacterization.unwrap_or(false);
    let shared = shared_imaginary_boundaries.unwrap_or(false);
    let resolution = delump_resolution.unwrap_or(DELUMP_RESOLUTION);

    // **The class's own delumping stage**: the source table expanded before the walk bins it.
    let expanded;
    let source: &[SourceRow] = if delump {
        expanded = expand(source, resolution);
        &expanded
    } else {
        source
    };
    let source = rows(source);
    // **The rows are sorted, and not only their keys.** The walk reads the keys; the inheriting
    // variant also reads the row that sits at a group's own index, so the order has to be the
    // rows' and not a parallel vector's - and the shared-boundaries variant reads the same rows
    // for their mole amounts.
    let mut ordered: Vec<&ReferenceRow> = reference.iter().collect();
    ordered.sort_by(|left, right| {
        sorting_key(left.boiling_point, left.molar_mass.value)
            .total_cmp(&sorting_key(right.boiling_point, right.molar_mass.value))
    });
    let reference: Vec<f64> = ordered
        .iter()
        .map(|row| sorting_key(row.boiling_point, row.molar_mass.value))
        .collect();

    let boundaries: Vec<f64> = if shared {
        // The equal-mass cuts are placed on the masses the reference carries, so this is the one
        // variant that reads its mole amounts - and it refuses a slate that states none rather
        // than placing the cuts on nothing.
        let mut slate = Vec::with_capacity(ordered.len());
        for (index, row) in ordered.iter().enumerate() {
            let Some(moles) = row.moles else {
                return Err(AzothError::invalid_input(
                    format!("reference[{index}].moles"),
                    "the shared-boundaries variant places the cuts on the reference's own masses, \
                     and this row states no mole amount",
                ));
            };
            slate.push(Lump {
                name: String::new(),
                moles,
                molar_mass: row.molar_mass,
                density: row.density,
                boiling_point: row.boiling_point,
            });
        }
        equal_mass_boundaries(&slate, resolution)
    } else {
        // `determineReferenceBoundaries`: the midpoint of each adjacent pair, or the larger of the
        // two when either is not finite.
        reference
            .windows(2)
            .map(|pair| {
                let (low, high) = (pair[0], pair[1]);
                if low.is_finite() && high.is_finite() {
                    0.5 * (low + high)
                } else {
                    low.max(high)
                }
            })
            .collect()
    };

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
        // **The check happens even where the flag is off**, because the class reads the reference's
        // density only on that path - a reference slate that states none is a slate the inheriting
        // variant cannot answer for, and the refusal is what says so.
        let (molar_mass, density) = if inherit {
            let row = ordered[index];
            let Some(density) = row.density else {
                return Err(AzothError::invalid_input(
                    format!("reference[{index}].density"),
                    "the inheriting variant takes each group's density from the reference row, and \
                     this one states none",
                ));
            };
            (row.molar_mass, density)
        } else {
            (
                kilograms_per_mole(group.mass / group.moles),
                kilograms_per_cubic_meter(if group.volume > MASS_TOLERANCE {
                    group.mass / group.volume
                } else {
                    group.density_mass / group.mass
                }),
            )
        };
        group_molar_mass.push(molar_mass);
        group_density.push(density);
    }

    Ok(CharacterizeToReferenceResult {
        reference_index,
        group_moles,
        group_molar_mass,
        group_density,
        warnings,
    })
}
