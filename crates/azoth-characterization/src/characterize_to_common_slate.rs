//! `characterization.characterize_to_common_slate` - several fluids onto one shared slate.
//!
//! Spec: `specs/models/characterization/characterize_to_common_slate.toml`. Oracle:
//! `validation/neqsim/captures/characterization_fluid_rebuild_probe.tsv`.
//!
//! # What the slate is
//!
//! Unlike [`crate::combine_reservoir_fluids`], the fluids are not blended: each comes back
//! separate, with its **own** mole amounts but an **identical** set of pseudo-component
//! properties. That is Pedersen's "common EoS" slate - several fluids, one equation of state -
//! and it is what makes a field's fluids comparable when a report and a model disagree about
//! how they were cut.
//!
//! The shared properties are weighted by `Wgt(j) * z_ij`, the fluid's slate weight times the
//! cut's mole fraction of that fluid's total. The molar mass is that mean of the cuts' own; the
//! density is the mean molar mass over the mean molar volume, so a cut's mass and its volume
//! enter separately.
//!
//! # The answer's mass can move
//!
//! A fluid keeps its own moles and receives the common molar mass and density, so its physical
//! mass changes wherever its own cuts differed from the slate. That is the slate's definition
//! rather than an error, and the spec says so.
//!
//! # The cut names are positions in the slate
//!
//! A cut that is invalid for any fluid is dropped from **every** fluid, and the names keep
//! their slate positions - `PC1_PC`, `PC3_PC` - so a cut can be matched across the fluids by
//! name even where one fluid does not carry it.

use azoth_core::{AzothError, Result};

use crate::combiner::{
    MASS_TOLERANCE, Profile, base_row, below, compute_common_slate, distribute, extract,
    profile_row, quantile_boundaries,
};
use crate::component_row::{CharacterizationComponent, FluidComponent};
use crate::results::CharacterizeToCommonSlateResult;

/// Several reservoir fluids re-characterized onto one shared pseudo-component slate.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `target_pseudo_components` is given and is not
///   positive.
/// * [`azoth_core::AzothError::InvalidInput`] if `weights` is not one per fluid, is negative
///   anywhere, does not have one positive entry, if a pseudo-component row states no molar mass
///   or density, if a row's fluid index is not a whole number from zero, or if no fluid is given.
///
/// # Example
/// ```
/// use azoth_characterization::{FluidComponent, characterize_to_common_slate};
///
/// let row = |fluid: f64, name: &str, moles: f64, mass: f64, tb: f64| FluidComponent {
///     fluid,
///     name: name.to_owned(),
///     pseudo_component: true,
///     moles,
///     molar_mass: Some(mass),
///     density: Some(750.0),
///     normal_boiling_point: Some(tb),
///     critical_temperature: Some(500.0),
///     critical_pressure: Some(3.0e6),
///     acentric_factor: Some(0.4),
///     critical_volume: None,
///     racket_z: None,
///     racket_z_cpa: None,
///     parachor: None,
///     critical_viscosity: None,
///     triple_point_temperature: None,
///     heat_of_fusion: None,
///     ideal_gas_enthalpy_of_formation: None,
///     cp_a: None,
///     cp_b: None,
///     cp_c: None,
///     cp_d: None,
///     attractive_m: None,
/// };
/// let slates = characterize_to_common_slate(
///     &[row(0.0, "C7_PC", 1.0, 0.100, 370.0), row(1.0, "C10_PC", 1.0, 0.140, 440.0)],
///     None,
///     None,
/// )?;
/// assert_eq!(slates.fluid_index, vec![0.0, 1.0]);
/// assert_eq!(slates.components.len(), 2);
/// assert_eq!(slates.components[0].name, "PC1_PC");
/// // The two fluids share the slate's molar mass, and each keeps its own moles.
/// assert_eq!(
///     slates.components[0].molar_mass,
///     slates.components[1].molar_mass
/// );
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn characterize_to_common_slate(
    components: &[FluidComponent],
    weights: Option<&[f64]>,
    target_pseudo_components: Option<i64>,
) -> Result<CharacterizeToCommonSlateResult> {
    let extractions = extract(components)?;
    if extractions.is_empty() {
        return Err(AzothError::invalid_input(
            "components",
            "a common slate needs at least one fluid",
        ));
    }
    let warnings = Vec::new();

    let resolved_weights = resolve_weights(weights, extractions.len())?;
    let target = resolve_target(target_pseudo_components, &extractions)?;

    // The class's own totals: every base component's amount and every pseudo-component's.
    let fluid_mole_totals: Vec<f64> = extractions
        .iter()
        .map(|extraction| {
            extraction.base.iter().map(|(_, moles)| moles).sum::<f64>()
                + extraction.pseudo.iter().map(|row| row.moles).sum::<f64>()
        })
        .collect();

    let any_pseudo = extractions
        .iter()
        .any(|extraction| !extraction.pseudo.is_empty());
    let mut rows: Vec<CharacterizationComponent> = Vec::new();
    let mut fluid_index: Vec<f64> = Vec::new();
    let mut fluid_moles: Vec<f64> = Vec::new();

    if any_pseudo {
        let mut pool = Vec::new();
        for (fluid, extraction) in extractions.iter().enumerate() {
            if extraction.pseudo.is_empty() || below(resolved_weights[fluid], 0.0) {
                continue;
            }
            let pseudo_mass: f64 = extraction.pseudo.iter().map(|row| row.mass).sum();
            if below(pseudo_mass, MASS_TOLERANCE) {
                continue;
            }
            let scale = resolved_weights[fluid] / pseudo_mass;
            pool.extend(
                extraction
                    .pseudo
                    .iter()
                    .map(|row| row.scaled_by_mass(scale)),
            );
        }
        pool.sort_by(|left, right| left.sorting_key().total_cmp(&right.sorting_key()));

        let boundaries = quantile_boundaries(&pool, target);
        let per_fluid: Vec<Vec<Profile>> = extractions
            .iter()
            .map(|extraction| distribute(&extraction.pseudo, &boundaries, target))
            .collect();
        let common =
            compute_common_slate(&per_fluid, &fluid_mole_totals, &resolved_weights, target);

        for (fluid, extraction) in extractions.iter().enumerate() {
            let mut total = 0.0;
            for (name, moles) in &extraction.base {
                rows.push(base_row(name.clone(), *moles));
                fluid_index.push(fluid as f64);
                total += moles;
            }
            for (cut, (shared, own)) in common.iter().zip(&per_fluid[fluid]).enumerate() {
                if !shared.is_valid() || !own.is_valid() {
                    continue;
                }
                let row = Profile {
                    moles: own.moles,
                    ..*shared
                };
                rows.push(profile_row(format!("PC{}_PC", cut + 1), &row));
                fluid_index.push(fluid as f64);
                total += own.moles;
            }
            fluid_moles.push(total);
        }
    } else {
        // The class's own early return: the fluids come back as themselves.
        for (fluid, extraction) in extractions.iter().enumerate() {
            let mut total = 0.0;
            for (name, moles) in &extraction.base {
                rows.push(base_row(name.clone(), *moles));
                fluid_index.push(fluid as f64);
                total += moles;
            }
            fluid_moles.push(total);
        }
    }

    Ok(CharacterizeToCommonSlateResult {
        components: rows,
        fluid_index,
        fluid_moles,
        resolved_pseudo_components: target as u32,
        warnings,
    })
}

/// `resolveWeights`: one weight per fluid, all ones when none is given.
fn resolve_weights(weights: Option<&[f64]>, count: usize) -> Result<Vec<f64>> {
    let Some(weights) = weights else {
        return Ok(vec![1.0; count]);
    };
    if weights.len() != count {
        return Err(AzothError::invalid_input(
            "weights",
            format!(
                "weights length ({}) must match number of fluids ({count})",
                weights.len()
            ),
        ));
    }
    let mut sum = 0.0;
    for weight in weights {
        if *weight < 0.0 {
            return Err(AzothError::invalid_input(
                "weights",
                "weights must be non-negative",
            ));
        }
        sum += weight;
    }
    if below(sum, 0.0) {
        return Err(AzothError::invalid_input(
            "weights",
            "at least one weight must be positive",
        ));
    }
    Ok(weights.to_vec())
}

/// The slate size: supplied, or the largest pseudo-component count among the fluids, at least
/// one.
fn resolve_target(
    target: Option<i64>,
    extractions: &[crate::combiner::Extraction],
) -> Result<usize> {
    let Some(target) = target else {
        let inferred = extractions
            .iter()
            .map(|extraction| extraction.pseudo.len())
            .max()
            .unwrap_or(0);
        return Ok(inferred.max(1));
    };
    usize::try_from(target)
        .ok()
        .filter(|target| *target > 0)
        .ok_or_else(|| AzothError::OutOfRange {
            field: "target_pseudo_components".to_string(),
            value: target as f64,
            detail: "a slate needs at least one pseudo-component".to_string(),
        })
}
