//! `oil_assay_to_mixture` - an assay turned into a fluid a flash can take.
//!
//! This is `OilAssayCharacterisation.apply`, less the class: the derivation
//! ([`crate::oil_assay_derivation`]) resolves each cut's molar mass and mole amount, the cut's card
//! comes from [`crate::tbp_fraction_component`], and the mixture is built by
//! [`crate::mixture_from_characterization`] - constructed rather than mutated, which is the one
//! place this differs from the class's `addTBPfraction`/`addComponent` sequence.
//!
//! **The density the card is built from is the one the derivation resolved**, stated or implied by
//! the cut's API gravity. A cut the assay gives as a standard substance is a base component: its
//! constants are the databank's, and no card is built for it.
//!
//! It is a helper rather than an id because its answer is a fluid, and a fluid is not a result.

use azoth_core::units::kilograms_per_cubic_meter;
use azoth_core::{AzothError, Result};
use azoth_eos::{Cubic, IdealGasModel, Mixture};

use crate::assay_mass_fractions::AssayBasis;
use crate::assembly::mixture_from_characterization;
use crate::component_row::FluidComponent;
use crate::oil_assay_derivation::{AssayCutDerivation, oil_assay_derivation, specific_gravity};
use crate::tbp_fraction_component::tbp_fraction_component;

/// An oil assay as the mixture its cuts describe.
///
/// # Errors
/// * Everything [`crate::oil_assay_derivation`] refuses: a duplicate name or substance, a
///   standard component on a volume basis, a cut with neither a gravity nor a boiling point and
///   Watson factor, or a substance the databank does not carry.
/// * [`AzothError::InvalidInput`] for a cut whose card has no liquid root, or a mixture
///   [`crate::mixture_from_characterization`] cannot build.
pub fn oil_assay_to_mixture(
    basis: AssayBasis,
    cuts: &[AssayCutDerivation],
    total_assay_mass: f64,
    cubic: Cubic,
) -> Result<(Mixture, IdealGasModel)> {
    let derived = oil_assay_derivation(basis, cuts, total_assay_mass)?;

    // A derived row is named after its cut - `<cut>_PC` - or after the substance it is. The names
    // are matched case-insensitively, because the derivation matches substances that way, and two
    // cuts resolving to one name are refused here rather than silently sharing a row.
    let mut keys: Vec<String> = Vec::with_capacity(cuts.len());
    for cut in cuts {
        let key = cut
            .standard_component_name
            .clone()
            .unwrap_or_else(|| format!("{}_PC", cut.name))
            .to_lowercase();
        if keys.contains(&key) {
            return Err(AzothError::invalid_input(
                "cuts",
                format!("two cuts name one component `{key}` once the suffix is added"),
            ));
        }
        keys.push(key);
    }

    let mut rows = Vec::with_capacity(derived.component_name.len());
    for (index, name) in derived.component_name.iter().enumerate() {
        let key = name.to_lowercase();
        let position = keys.iter().position(|known| *known == key).ok_or_else(|| {
            AzothError::invalid_input(
                "cuts",
                format!("the derivation answered `{name}`, which no cut describes"),
            )
        })?;
        let cut = &cuts[position];
        let moles = derived.moles[index];
        if cut.standard_component_name.is_some() {
            rows.push(base_row(name.clone(), moles));
            continue;
        }
        let gravity = specific_gravity(cut).ok_or_else(|| {
            AzothError::invalid_input(
                "cuts",
                format!("`{}` states neither a density nor an API gravity", cut.name),
            )
        })?;
        let card = tbp_fraction_component(
            derived.molar_mass[index],
            kilograms_per_cubic_meter(gravity * 1000.0),
        )?;
        rows.push(FluidComponent {
            fluid: 0.0,
            name: name.clone(),
            pseudo_component: true,
            moles,
            molar_mass: Some(derived.molar_mass[index].value),
            density: Some(gravity * 1000.0),
            normal_boiling_point: Some(card.boiling_temperature.value),
            critical_temperature: Some(card.critical_temperature.value),
            critical_pressure: Some(card.critical_pressure.value),
            acentric_factor: Some(card.acentric_factor),
            critical_volume: Some(card.critical_volume.value),
            racket_z: Some(card.racket_z),
            racket_z_cpa: Some(card.racket_z_cpa),
            parachor: Some(card.parachor),
            critical_viscosity: Some(card.critical_viscosity.value),
            triple_point_temperature: Some(card.triple_point_temperature.value),
            heat_of_fusion: Some(card.heat_of_fusion.value),
            ideal_gas_enthalpy_of_formation: Some(card.ideal_gas_enthalpy_of_formation.value),
            cp_a: Some(card.cp_a),
            cp_b: Some(card.cp_b),
            cp_c: Some(card.cp_c),
            cp_d: Some(card.cp_d),
            attractive_m: Some(card.attraction_exponent),
        });
    }

    mixture_from_characterization(&rows, cubic)
}

/// A standard substance's row: its name and amount, every property absent because the databank
/// has them.
fn base_row(name: String, moles: f64) -> FluidComponent {
    FluidComponent {
        fluid: 0.0,
        name,
        pseudo_component: false,
        moles,
        molar_mass: None,
        density: None,
        normal_boiling_point: None,
        critical_temperature: None,
        critical_pressure: None,
        acentric_factor: None,
        critical_volume: None,
        racket_z: None,
        racket_z_cpa: None,
        parachor: None,
        critical_viscosity: None,
        triple_point_temperature: None,
        heat_of_fusion: None,
        ideal_gas_enthalpy_of_formation: None,
        cp_a: None,
        cp_b: None,
        cp_c: None,
        cp_d: None,
        attractive_m: None,
    }
}
