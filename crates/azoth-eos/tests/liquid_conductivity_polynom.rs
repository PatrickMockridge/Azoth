//! Spec-driven tests for the `eos.liquid_conductivity_polynom` model.

use azoth_core::units::{kelvins, kilograms_per_mole};
use azoth_eos::liquid_conductivity_polynom::liquid_conductivity_polynom;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.liquid_conductivity_polynom";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::LiquidConductivityPolynomResult {
    let flat = case
        .matrix("liquid_conductivity")
        .expect("the case states the coefficient matrix");
    assert_eq!(flat.len() % 3, 0, "three coefficients per component");
    let rows: Vec<[f64; 3]> = flat.chunks(3).map(|row| [row[0], row[1], row[2]]).collect();
    let molar_mass: Vec<_> = case
        .vector("molar_mass")
        .expect("the case states the molar masses")
        .iter()
        .map(|value| kilograms_per_mole(*value))
        .collect();
    liquid_conductivity_polynom(
        &rows,
        &molar_mass,
        case.vector("z")
            .expect("the case states the mole fractions"),
        kelvins(case.input("T").expect("the case states T")),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = azoth_eos::model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        common::assert_close(
            result.k.value,
            case.expected_value("k").expect("the case states k"),
            case.tolerance,
            &format!("{context} (k)"),
        );
        common::assert_consistent(&result, context);
    }
}

/// **The mean is weighted by mass, not by moles** - which is the one thing a port can get wrong
/// here and still look plausible.
///
/// The two components' polynomials are made to differ by a factor of a hundred, so a mole-weighted
/// mean and a mass-weighted one land far apart: 0.5 kPa*m/K against 0.99 for equal mole fractions
/// whose molar masses differ by a factor of a hundred.
#[test]
fn the_mean_is_weighted_by_mass() {
    let coefficients = [[1.0, 0.0, 0.0], [100.0, 0.0, 0.0]];
    let k = liquid_conductivity_polynom(
        &coefficients,
        &[kilograms_per_mole(0.002), kilograms_per_mole(0.200)],
        &[0.5, 0.5],
        kelvins(300.0),
    )
    .expect("the mean runs")
    .k
    .value;
    // Equal moles at a hundredfold molar mass are 1:100 by mass, so the heavy component
    // carries 100/101 of the mean.
    let expected = 100.0 * 100.0 / 101.0 + 1.0 / 101.0;
    assert!(
        (k - expected).abs() < 1e-9,
        "k = {k} against the mass-weighted {expected}"
    );
    assert!(k > 50.0, "a mole-weighted mean would be 50.5, not {k}");
}

/// The floor: a polynomial that would go negative answers `1e-10` rather than a negative
/// conductivity - which NeqSim's own light-gas rows reach above a few hundred kelvin.
#[test]
fn a_negative_polynomial_is_floored() {
    let k = liquid_conductivity_polynom(
        &[[1.0, -1.0, 0.0]],
        &[kilograms_per_mole(0.016)],
        &[1.0],
        kelvins(400.0),
    )
    .expect("the mean runs")
    .k
    .value;
    assert!((k - 1.0e-10).abs() < 1e-18, "k = {k}");
}
