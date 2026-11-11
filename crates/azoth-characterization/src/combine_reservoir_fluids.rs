//! `characterization.combine_reservoir_fluids` - several fluids blended onto one slate.
//!
//! Spec: `specs/models/characterization/combine_reservoir_fluids.toml`. Oracle:
//! `validation/neqsim/captures/characterization_fluid_rebuild_probe.tsv`.
//!
//! # What the blend is
//!
//! Each fluid's base components merge by name, its pseudo-components are pooled, and the pooled
//! mass is cut into `target_pseudo_components` equal-mass groups. Each fluid's own list is then
//! distributed onto **those** groups, and the per-fluid groups blended: so a cut's properties
//! are a mixture of the fluids' own cuts on either side of it, weighted by how much of the
//! blend each fluid is.
//!
//! **The cut points are keys, not interpolations.** A boundary is the sorting key of the row
//! whose running mass first reaches the target, so a group's edge is always a boiling point (or
//! molar mass) some row carries, and a target nothing reaches is an infinity.
//!
//! # A row is a record
//!
//! `PseudoComponentCombination` carries twenty-one numbers per pseudo-component; the row here
//! carries all of them, because a rebuilt component needs every one and a blend has no fluid to
//! read them back from. A base component carries only its name and amount - its constants are
//! the databank's, and [`crate::assembly`] is what resolves them.

use azoth_core::{AzothError, Result, apply_checks};

use crate::combiner::{
    base_row, combine_profiles, distribute, extract, merge_bases, profile_row, quantile_boundaries,
};
use crate::component_row::{CharacterizationComponent, FluidComponent};
use crate::model_gen;
use crate::results::CombineReservoirFluidsResult;

/// Several reservoir fluids, blended into one fluid with a common pseudo-component slate.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `target_pseudo_components` is not positive.
/// * [`azoth_core::AzothError::InvalidInput`] if a pseudo-component row states no molar mass or
///   density, or a row's fluid index is not a whole number from zero.
///
/// # Example
/// ```
/// use azoth_characterization::{FluidComponent, combine_reservoir_fluids};
///
/// let row = |fluid: f64, name: &str, moles: f64, mass: f64, tb: f64| FluidComponent {
///     fluid,
///     name: name.to_owned(),
///     pseudo_component: true,
///     moles,
///     molar_mass: Some(mass),
///     density: Some(750.0),
///     normal_boiling_point: Some(tb),
///     critical_temperature: None,
///     critical_pressure: None,
///     acentric_factor: None,
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
/// let blended = combine_reservoir_fluids(
///     2,
///     &[row(0.0, "C7_PC", 1.0, 0.100, 370.0), row(1.0, "C10_PC", 1.0, 0.100, 440.0)],
/// )?;
/// assert_eq!(blended.components.len(), 2);
/// assert_eq!(blended.components[0].name, "PC1_PC");
/// assert_eq!(blended.components[1].name, "PC2_PC");
/// assert!((blended.total_moles - 2.0).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn combine_reservoir_fluids(
    target_pseudo_components: i64,
    components: &[FluidComponent],
) -> Result<CombineReservoirFluidsResult> {
    let spec = &model_gen::COMBINE_RESERVOIR_FLUIDS_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "target_pseudo_components" => Some(target_pseudo_components as f64),
            _ => None,
        },
        &mut warnings,
    )?;
    let target: usize = target_pseudo_components.try_into().map_err(|_| {
        AzothError::invalid_input(
            "target_pseudo_components",
            "a combined fluid needs at least one pseudo-component",
        )
    })?;

    let extractions = extract(components)?;
    let base = merge_bases(&extractions);
    let mut total_moles: f64 = base.iter().map(|(_, moles)| moles).sum();
    let mut rows: Vec<CharacterizationComponent> = base
        .iter()
        .map(|(name, moles)| base_row(name.clone(), *moles))
        .collect();

    // **The class's own early return**: a pool with nothing in it is a fluid made of its base
    // components, and the requested slate never enters the answer.
    let pooled: Vec<_> = extractions
        .iter()
        .flat_map(|extraction| extraction.pseudo.iter().cloned())
        .collect();
    if pooled.is_empty() {
        return Ok(CombineReservoirFluidsResult {
            components: rows,
            total_moles,
            warnings,
        });
    }

    let boundaries = quantile_boundaries(&pooled, target);
    let per_fluid = extractions
        .iter()
        .map(|extraction| distribute(&extraction.pseudo, &boundaries, target))
        .collect::<Vec<_>>();
    let fluid_mass_totals: Vec<f64> = per_fluid
        .iter()
        .map(|profiles| profiles.iter().map(|profile| profile.mass()).sum())
        .collect();
    let fluid_mole_totals: Vec<f64> = per_fluid
        .iter()
        .map(|profiles| profiles.iter().map(|profile| profile.moles).sum())
        .collect();

    let mut counter = 1;
    for profile in combine_profiles(&per_fluid, &fluid_mass_totals, &fluid_mole_totals) {
        if !profile.is_valid() {
            continue;
        }
        total_moles += profile.moles;
        rows.push(profile_row(format!("PC{counter}_PC"), &profile));
        counter += 1;
    }

    Ok(CombineReservoirFluidsResult {
        components: rows,
        total_moles,
        warnings,
    })
}
