//! Spec-driven tests for the `characterization.set_heavy_tbp_to_plus` model.
//!
//! Every expected value is a row of `validation/neqsim/captures/set_heavy_tbp_probe.tsv`.

use azoth_characterization::{SetHeavyTbpToPlusResult, model_gen, set_heavy_tbp_to_plus};
use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_core::{AzothError, spec::TestCase};
use azoth_test_support as common;

const MODEL_ID: &str = "characterization.set_heavy_tbp_to_plus";

/// The rows, from the four columns the case states.
fn slate(name: &str, case: &TestCase) -> Vec<String> {
    common::list_input(case, name)
        .iter()
        .map(|value| (*value).to_string())
        .collect()
}

fn call(case: &TestCase) -> SetHeavyTbpToPlusResult {
    let names = slate("names", case);
    let names: Vec<&str> = names.iter().map(String::as_str).collect();
    let molar_mass: Vec<_> = case
        .vector("molar_mass")
        .expect("the case states a slate")
        .iter()
        .map(|value| kilograms_per_mole(*value))
        .collect();
    let density: Vec<_> = case
        .vector("normal_liquid_density")
        .expect("the case states a slate")
        .iter()
        .map(|value| kilograms_per_cubic_meter(*value))
        .collect();
    set_heavy_tbp_to_plus(
        &names,
        &molar_mass,
        case.vector("mole_fraction")
            .expect("the case states a slate"),
        &density,
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = format!("{}::{}", spec.id, case.id);
        for (name, got, want) in [
            (
                "plus_index",
                result.plus_index,
                case.expected_value("plus_index"),
            ),
            (
                "plus_molar_mass",
                result.plus_molar_mass.value,
                case.expected_value("plus_molar_mass"),
            ),
            (
                "plus_mole_fraction",
                result.plus_mole_fraction,
                case.expected_value("plus_mole_fraction"),
            ),
            (
                "plus_density",
                result.plus_density.value,
                case.expected_value("plus_density"),
            ),
            (
                "last_tbp_density",
                result.last_tbp_density.value,
                case.expected_value("last_tbp_density"),
            ),
            (
                "density_intercept",
                result.density_intercept.value,
                case.expected_value("density_intercept"),
            ),
        ] {
            let Some(want) = want else { continue };
            common::assert_close(got, want, case.tolerance, &format!("{context}::{name}"));
        }
        assert_eq!(
            result.plus_name,
            case.expected_string("plus_name")
                .expect("the case names the plus row"),
            "{context}: the plus row's name"
        );
        common::assert_consistent(&result, &context);
    }
}

/// One fluid's names and molar masses, and its densities in the unit the spec declares.
fn rows(names: &[&str], molar_mass: &[f64], density: &[f64]) -> (Vec<String>, Vec<f64>, Vec<f64>) {
    (
        names.iter().map(|name| (*name).to_string()).collect(),
        molar_mass.to_vec(),
        density.to_vec(),
    )
}

fn density_of(values: &[f64]) -> Vec<azoth_core::units::MassDensity> {
    values
        .iter()
        .map(|value| kilograms_per_cubic_meter(*value))
        .collect()
}

fn mass_of(values: &[f64]) -> Vec<azoth_core::units::MolarMass> {
    values
        .iter()
        .map(|value| kilograms_per_mole(*value))
        .collect()
}

#[test]
fn the_carbon_number_is_read_out_of_the_name_and_not_out_of_the_position() {
    // **The answer is a name and not a row.** The last row carries the *smallest* carbon number
    // here, so a walk that took the last row, or sorted by molar mass, answers differently - and
    // the capture's first fluid is this one.
    let (names, mass, density) = rows(
        &["C7_PC", "C20_PC", "C12_PC", "C9_PC"],
        &[0.1, 0.35, 0.2, 0.14],
        &[720.0, 900.0, 800.0, 750.0],
    );
    let names: Vec<&str> = names.iter().map(String::as_str).collect();
    let result = set_heavy_tbp_to_plus(&names, &mass_of(&mass), &[0.25; 4], &density_of(&density))
        .expect("computes");
    assert_eq!(result.plus_index, 1.0);
    assert_eq!(result.plus_name, "C20_PC");
    // `densLastTBP` is the row **below** the winner in index order, which is `C7_PC`'s 0.72 and
    // not `C9_PC`'s 0.75.
    assert!((result.last_tbp_density.value - 720.0).abs() < 1e-12);
}

#[test]
fn two_digit_names_are_read_as_two_digits() {
    // `C10_PC` reads 10 because its fourth character is the underscore; reading one character
    // would leave `C9_PC` the winner of this slate.
    let (names, mass, density) = rows(
        &["C6_PC", "C9_PC", "C10_PC", "C11_PC"],
        &[0.086, 0.13, 0.15, 0.17],
        &[680.0, 750.0, 780.0, 800.0],
    );
    let names: Vec<&str> = names.iter().map(String::as_str).collect();
    let result = set_heavy_tbp_to_plus(&names, &mass_of(&mass), &[0.25; 4], &density_of(&density))
        .expect("computes");
    assert_eq!(result.plus_index, 3.0);
    assert_eq!(result.plus_name, "C11_PC");
}

#[test]
fn a_name_that_carries_no_carbon_number_is_refused_by_name() {
    // **The class logs this and carries on**, choosing the plus end from whatever parsed. This
    // refuses instead, and names the row - the capture's third fluid is the class's answer to the
    // same slate.
    let (names, mass, density) = rows(
        &["C7_PC", "AB_PC", "C13_PC"],
        &[0.1, 0.2, 0.22],
        &[720.0, 800.0, 840.0],
    );
    let names: Vec<&str> = names.iter().map(String::as_str).collect();
    let refused = set_heavy_tbp_to_plus(
        &names,
        &mass_of(&mass),
        &[1.0 / 3.0; 3],
        &density_of(&density),
    );
    let error = refused.expect_err("`AB_PC` carries no carbon number");
    assert!(
        error.to_string().contains("AB_PC"),
        "the refusal names the row: {error}"
    );
    assert!(matches!(error, AzothError::InvalidInput { .. }));
}

#[test]
fn a_plus_end_at_the_first_row_is_refused_because_the_row_below_it_is_what_the_class_reads() {
    // The capture's fourth fluid: one row, so the winner is index zero and the class reads index
    // -1. A longer slate whose largest carbon number happens to be first is the same read.
    let one = set_heavy_tbp_to_plus(&["C9_PC"], &mass_of(&[0.15]), &[1.0], &density_of(&[780.0]));
    assert!(matches!(one, Err(AzothError::InvalidInput { .. })));

    let (names, mass, density) = rows(&["C20_PC", "C7_PC"], &[0.35, 0.1], &[900.0, 720.0]);
    let names: Vec<&str> = names.iter().map(String::as_str).collect();
    let first = set_heavy_tbp_to_plus(&names, &mass_of(&mass), &[0.5, 0.5], &density_of(&density));
    assert!(matches!(first, Err(AzothError::InvalidInput { .. })));
}

#[test]
fn components_that_are_not_the_slates_rows_are_never_handed_over() {
    // The class's first test is `isIsTBPfraction`, so the capture's sixth fluid - methane ahead
    // of three TBP rows - skips it. This id takes the slate, so the filter is the caller's; a
    // light end passed anyway carries no carbon number and is refused like any other name.
    let (names, mass, density) = rows(&["methane", "C7_PC"], &[0.016, 0.1], &[422.182, 720.0]);
    let names: Vec<&str> = names.iter().map(String::as_str).collect();
    let refused =
        set_heavy_tbp_to_plus(&names, &mass_of(&mass), &[0.5, 0.5], &density_of(&density));
    assert!(matches!(refused, Err(AzothError::InvalidInput { .. })));
}

#[test]
fn the_slate_is_four_aligned_columns() {
    let (names, mass, density) = rows(
        &["C7_PC", "C20_PC", "C9_PC"],
        &[0.1, 0.35, 0.14],
        &[720.0, 900.0, 750.0],
    );
    let names: Vec<&str> = names.iter().map(String::as_str).collect();
    assert!(matches!(
        set_heavy_tbp_to_plus(
            &names[..2],
            &mass_of(&mass),
            &[0.25; 3],
            &density_of(&density)
        ),
        Err(AzothError::InvalidInput { .. })
    ));
    assert!(matches!(
        set_heavy_tbp_to_plus(&[], &[], &[], &[]),
        Err(AzothError::InvalidInput { .. })
    ));
}
