//! `characterization.delumping` - a lumped pseudo-component table split into a finer grid.
//!
//! Spec: `specs/models/characterization/delumping.toml`. Oracle:
//! `validation/neqsim/captures/pseudo_component_delumping_probe.tsv`.
//!
//! # Two conservations, and they are different ones
//!
//! `delumpContributions` splits each lump into `resolution` sub-fractions whose mole *weights* sum
//! to one, so the parent's moles are conserved by the weights themselves - and then rescales the
//! molar-mass grid by a single factor so that `sum(n_k * M_k)` is the parent's mass exactly. The
//! moles are what the light-end bias was stated in, so it is the molar masses that move.
//!
//! The sub-fractions' normal boiling points are spread with those molar masses,
//! `Tb_k = Tb * M_k / M`, which is what lets them cross a reference slate's cut boundaries - and
//! the expanded list is re-sorted by sorting key so it can be fed straight to a cut.
//!
//! # One read of two fields, and sixteen copies
//!
//! The class's row carries twenty-one numbers. This stage *reads* five of them and *copies* the
//! rest to every sub-fraction unchanged, which is why the row here is five fields wide: the
//! helper is shared with the walk that redistributes the result, and that walk is where the
//! copies would be read.

use azoth_core::units::{MassDensity, MolarMass, kilograms_per_mole};
use azoth_core::{Result, apply_checks};

use crate::model_gen;
use crate::results::{DelumpingResult, SubFraction};

/// `PseudoComponentCombiner.MASS_TOLERANCE`, on the rescale's denominator.
pub(crate) const MASS_TOLERANCE: f64 = 1.0e-12;

/// The class's own `window` and `decay`: the grid's half-width, and how fast the bias falls.
const WINDOW: f64 = 0.4;
const DECAY: f64 = 1.5;

/// One lump of a pseudo-component table, as this stage takes it.
///
/// The five fields `delumpContributions` reads. `name` is not one of them - it is the string each
/// sub-fraction's name is built from, which is the only use the method makes of it.
#[derive(Debug, Clone, PartialEq)]
pub struct Lump {
    /// The lump's name, which each sub-fraction's is built from.
    pub name: String,
    /// Its mole amount, conserved across its sub-fractions.
    pub moles: f64,
    /// Its molar mass, which the grid is a window around.
    pub molar_mass: MolarMass,
    /// Its normal liquid density, copied to every sub-fraction rather than read.
    pub density: Option<MassDensity>,
    /// Its normal boiling point, which its sub-fractions' are spread from.
    pub boiling_point: Option<f64>,
}

/// `PseudoComponentContribution.sortingKey`: the boiling point, or the molar mass without one.
///
/// **The absence has two spellings and they are the same thing.** A row that states no boiling
/// point, and a row that states one which is not finite and positive, both fall back - the first
/// because the record has no field for it and the second because the class's own test is on the
/// value. A sentinel a caller wrote into a vector is now a field the caller leaves out.
///
/// Shared with [`crate::characterize_to_reference`], which sorts the same table by the same key.
pub(crate) fn sorting_key(boiling_point: Option<f64>, molar_mass: f64) -> f64 {
    match boiling_point {
        Some(tb) if tb.is_finite() && tb > 0.0 => tb,
        _ => molar_mass,
    }
}

/// `!(value > 0.0)`, NaN included: the class's own comparison lets a value it cannot order
/// through rather than splitting on one.
fn not_positive(value: f64) -> bool {
    value.is_nan() || value <= 0.0
}

/// One lump, as the sub-fraction that is the lump itself - what an unsplit one becomes.
fn unsplit(lump: &Lump) -> SubFraction {
    SubFraction {
        name: lump.name.clone(),
        moles: lump.moles,
        molar_mass: lump.molar_mass,
        density: lump.density,
        boiling_point: lump.boiling_point,
    }
}

/// `delumpContributions`, without the sort: the expanded table in the input's own order.
///
/// **A `resolution` of one or less returns the input unchanged, order included** - which is the
/// class's early return, and the reason this is a separate function: the sort below is part of the
/// splitting and not of the identity.
pub(crate) fn split(lumps: &[Lump], resolution: i64) -> Vec<SubFraction> {
    if resolution <= 1 {
        return lumps.iter().map(unsplit).collect();
    }
    let resolution = resolution as usize;
    let mut expanded = Vec::with_capacity(lumps.len() * resolution);
    for lump in lumps {
        // The class's own guard: a lump it cannot split is passed through rather than refused,
        // because the method defines what it does with one. Its `!(x > 0.0)` written the way
        // these crates write it, so a NaN and a non-positive both pass through.
        if not_positive(lump.moles) || not_positive(lump.molar_mass.value) {
            expanded.push(unsplit(lump));
            continue;
        }
        let position = |k: usize| k as f64 / (resolution - 1) as f64;
        let weights: Vec<f64> = (0..resolution).map(|k| (-DECAY * position(k)).exp()).collect();
        let weight_sum: f64 = weights.iter().sum();
        let grid: Vec<f64> = (0..resolution)
            .map(|k| lump.molar_mass.value * (1.0 - WINDOW + 2.0 * WINDOW * position(k)))
            .collect();

        let raw_mass: f64 = (0..resolution)
            .map(|k| lump.moles * (weights[k] / weight_sum) * grid[k])
            .sum();
        let parent_mass = lump.moles * lump.molar_mass.value;
        let scale = if raw_mass > MASS_TOLERANCE {
            parent_mass / raw_mass
        } else {
            1.0
        };

        for k in 0..resolution {
            let moles = lump.moles * (weights[k] / weight_sum);
            let molar_mass = grid[k] * scale;
            let ratio = molar_mass / lump.molar_mass.value;
            expanded.push(SubFraction {
                name: format!("{}_d{k}", lump.name),
                moles,
                molar_mass: kilograms_per_mole(molar_mass),
                density: lump.density,
                // The class's own `parentHasTb ? Tb * ratio : Tb`, with the absent spelling the
                // record uses in place of the non-finite one it would carry.
                boiling_point: match lump.boiling_point {
                    Some(tb) if tb.is_finite() && tb > 0.0 => Some(tb * ratio),
                    other => other,
                },
            });
        }
    }
    expanded.sort_by(|left, right| {
        sorting_key(left.boiling_point, left.molar_mass.value)
            .total_cmp(&sorting_key(right.boiling_point, right.molar_mass.value))
    });
    expanded
}

/// A lumped pseudo-component table split into a finer grid.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if a declared input falls outside its range. The
///   method itself refuses nothing: every lump table has an answer, and one it cannot split is
///   passed through.
///
/// # Example
/// ```
/// use azoth_core::units::kilograms_per_mole;
/// use azoth_characterization::{Lump, delumping};
///
/// let split = delumping(
///     &[Lump {
///         name: "C10_PC".to_string(),
///         moles: 2.0,
///         molar_mass: kilograms_per_mole(0.2),
///         density: None,
///         boiling_point: Some(500.0),
///     }],
///     2,
/// )?;
/// assert_eq!(split.fractions.len(), 2);
/// // The weights conserve the moles, and the one rescale conserves the mass.
/// let moles: f64 = split.fractions.iter().map(|f| f.moles).sum();
/// let mass: f64 = split.fractions.iter().map(|f| f.moles * f.molar_mass.value).sum();
/// assert!((moles - 2.0).abs() < 1e-12);
/// assert!((mass - 0.4).abs() < 1e-12);
/// // Spread with the molar mass, not held at the parent's.
/// assert!(split.fractions[0].boiling_point.unwrap() < 500.0);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn delumping(lumps: &[Lump], resolution: i64) -> Result<DelumpingResult> {
    let spec = &model_gen::DELUMPING_SPEC;
    let mut warnings = Vec::new();
    apply_checks(spec.input_checks(), |_quantity| None, &mut warnings)?;

    let fractions = split(lumps, resolution);

    apply_checks(spec.derived_checks(), |_quantity| None, &mut warnings)?;
    Ok(DelumpingResult { fractions, warnings })
}
