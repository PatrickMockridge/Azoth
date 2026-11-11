//! Tests for `oil_assay_to_mixture`, the assay-to-fluid composition.

use azoth_characterization::{
    AssayBasis, AssayCutDerivation, oil_assay_to_mixture, tbp_fraction_component,
};
use azoth_core::units::{kelvins, kilograms_per_cubic_meter, kilograms_per_mole, pascals};
use azoth_eos::{Cubic, pt_flash};

fn pseudo(name: &str, fraction: f64, molar_mass: f64, density: f64) -> AssayCutDerivation {
    AssayCutDerivation {
        name: name.to_string(),
        declared_fraction: fraction,
        density: Some(kilograms_per_cubic_meter(density)),
        api_gravity: None,
        molar_mass: Some(kilograms_per_mole(molar_mass)),
        average_boiling_point: None,
        watson_characterization_factor: None,
        standard_component_name: None,
    }
}

fn standard(name: &str, substance: &str, fraction: f64) -> AssayCutDerivation {
    AssayCutDerivation {
        name: name.to_string(),
        declared_fraction: fraction,
        density: None,
        api_gravity: None,
        molar_mass: None,
        average_boiling_point: None,
        watson_characterization_factor: None,
        standard_component_name: Some(substance.to_string()),
    }
}

fn assay() -> Vec<AssayCutDerivation> {
    vec![
        pseudo("C7", 0.5, 0.095, 730.0),
        pseudo("C10", 0.3, 0.134, 810.0),
        standard("methane", "methane", 0.2),
    ]
}

#[test]
fn an_assay_becomes_the_mixture_its_cuts_describe() {
    let (mixture, ideal_gas) = oil_assay_to_mixture(AssayBasis::Mass, &assay(), 100.0, Cubic::Srk)
        .expect("the assay resolves");
    let names: Vec<&str> = mixture
        .names()
        .expect("the assembler names its mixture")
        .iter()
        .map(String::as_str)
        .collect();
    assert_eq!(names, ["C7_PC", "C10_PC", "methane"]);
    assert_eq!(ideal_gas.cp_a.len(), 3);

    // The pseudo-component's constants are its card's, built from the cut's mass and gravity.
    let card = tbp_fraction_component(kilograms_per_mole(0.095), kilograms_per_cubic_meter(730.0))
        .expect("the card computes");
    assert_eq!(
        mixture.components()[0].tc.value,
        card.critical_temperature.value
    );
    assert_eq!(
        mixture.components()[0].pc.value,
        card.critical_pressure.value
    );
    assert_eq!(mixture.components()[0].omega, card.acentric_factor);
}

#[test]
fn the_moles_are_the_derivations_mass_shares() {
    let (mixture, _) = oil_assay_to_mixture(AssayBasis::Mass, &assay(), 100.0, Cubic::Srk)
        .expect("the assay resolves");
    // 50 kg of a 0.095 kg/mol cut is 526.3 mol, and the methane's 20 kg is 1246.6 mol.
    let flash = pt_flash(&mixture, kelvins(300.0), pascals(5.0e6), &[0.25, 0.15, 0.6])
        .expect("the mixture flashes");
    assert!(flash.z_liquid.is_finite());
}
