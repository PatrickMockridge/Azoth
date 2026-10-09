//! Spec-driven tests for the `characterization.characterize_to_reference` model.
//!
//! Every expected value is a row of
//! `validation/neqsim/captures/pseudo_component_combiner_probe.tsv`.

use azoth_characterization::{
    CharacterizeToReferenceResult, characterize_to_reference, model_gen,
};
use azoth_core::units::{MassDensity, MolarMass, kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_core::{AzothError, spec::TestCase};
use azoth_test_support as common;

const MODEL_ID: &str = "characterization.characterize_to_reference";

fn call(case: &TestCase) -> CharacterizeToReferenceResult {
    let molar_masses = |name: &str| -> Vec<_> {
        case.vector(name)
            .expect("the case states a table")
            .iter()
            .map(|v| kilograms_per_mole(*v))
            .collect()
    };
    let densities: Vec<_> = case
        .vector("source_density")
        .expect("the case states a table")
        .iter()
        .map(|v| kilograms_per_cubic_meter(*v))
        .collect();
    characterize_to_reference(
        case.vector("source_moles").expect("a table"),
        &molar_masses("source_molar_mass"),
        &densities,
        case.vector("source_boiling_point").expect("a table"),
        &molar_masses("reference_molar_mass"),
        case.vector("reference_boiling_point")
            .expect("a table"),
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
        for (name, actual, expected) in [
            (
                "group_moles",
                result.group_moles.clone(),
                case.expected_vector("group_moles"),
            ),
            (
                "group_molar_mass",
                result.group_molar_mass.iter().map(|q| q.value).collect(),
                case.expected_vector("group_molar_mass"),
            ),
            (
                "group_density",
                result.group_density.iter().map(|q| q.value).collect(),
                case.expected_vector("group_density"),
            ),
        ] {
            let Some(expected) = expected else { continue };
            assert_eq!(actual.len(), expected.len(), "{context}::{name}: the group count");
            for (index, (got, want)) in actual.iter().zip(expected).enumerate() {
                common::assert_close(
                    *got,
                    *want,
                    case.tolerance,
                    &format!("{context}::{name}[{index}]"),
                );
            }
        }
        common::assert_consistent(&result, &context);
    }
}

/// One fluid's table, for the tests that build their own.
fn fluid(
    moles: &[f64],
    molar_mass: &[f64],
    density: &[f64],
    boiling_point: &[f64],
) -> (Vec<f64>, Vec<MolarMass>, Vec<MassDensity>, Vec<f64>) {
    (
        moles.to_vec(),
        molar_mass.iter().map(|v| kilograms_per_mole(*v)).collect(),
        density
            .iter()
            .map(|v| kilograms_per_cubic_meter(*v))
            .collect(),
        boiling_point.to_vec(),
    )
}

#[test]
fn a_key_exactly_on_a_boundary_stays_in_the_lower_group() {
    // The walk advances while the key is **greater than** the boundary, so a key of exactly 550.0
    // against a boundary of 550.0 stays in the lower group - and the upper group, having nothing
    // left, is then dropped rather than reported empty. A port that used `>=` would return two
    // groups here instead of one.
    let run = |key: f64| {
        let (moles, masses, densities, keys) = fluid(&[1.0, 1.0], &[0.1, 0.2], &[700.0, 820.0], &[380.0, key]);
        characterize_to_reference(
            &moles, &masses, &densities, &keys,
            &[kilograms_per_mole(0.15), kilograms_per_mole(0.25)],
            &[400.0, 700.0],
        )
        .expect("computes")
    };
    assert_eq!(run(550.0).reference_index, vec![0.0]);
    assert_eq!(run(550.001).reference_index, vec![0.0, 1.0]);
    assert_eq!(run(550.0).group_moles, vec![2.0]);
}

#[test]
fn the_rows_are_sorted_rather_than_taken_as_given() {
    // `extractComponents` sorts both tables by key, so a source stated in any order gives the same
    // answer - which is what makes the groups contiguous runs rather than whatever order the
    // caller happened to use.
    let sorted = |order: [usize; 3]| {
        let (moles, masses, densities, keys) = fluid(&[1.0, 1.0, 1.0], &[0.1, 0.2, 0.3], &[700.0, 820.0, 880.0], &[380.0, 500.0, 700.0]);
        let pick = |values: &[f64]| order.iter().map(|i| values[*i]).collect::<Vec<_>>();
        characterize_to_reference(
            &pick(&moles),
            &order
                .iter()
                .map(|i| kilograms_per_mole(masses[*i].value))
                .collect::<Vec<_>>(),
            &order
                .iter()
                .map(|i| kilograms_per_cubic_meter(densities[*i].value))
                .collect::<Vec<_>>(),
            &pick(&keys),
            &[kilograms_per_mole(0.15), kilograms_per_mole(0.28)],
            &[450.0, 620.0],
        )
        .expect("computes")
    };
    assert_eq!(sorted([0, 1, 2]), sorted([2, 1, 0]));
    assert_eq!(sorted([0, 1, 2]), sorted([1, 2, 0]));
}

#[test]
fn a_zero_amount_row_is_not_on_the_grid() {
    // `extractComponents` skips a component with no moles, so a zero here is not a row that lands
    // in a group and adds nothing - it is a row the walk never sees.
    let (moles, masses, densities, keys) = fluid(&[1.0, 0.0, 1.0], &[0.1, 0.2, 0.3], &[700.0, 820.0, 880.0], &[380.0, 500.0, 700.0]);
    let result = characterize_to_reference(
        &moles,
        &masses,
        &densities,
        &keys,
        &[kilograms_per_mole(0.15), kilograms_per_mole(0.28)],
        &[450.0, 620.0],
    )
    .expect("computes");
    assert_eq!(result.group_moles.iter().sum::<f64>(), 2.0);
}

#[test]
fn a_table_has_to_be_internally_aligned() {
    let (moles, masses, densities, keys) = fluid(&[1.0, 1.0], &[0.1], &[700.0, 820.0], &[380.0, 500.0]);
    assert!(matches!(
        characterize_to_reference(
            &moles, &masses, &densities, &keys,
            &[kilograms_per_mole(0.15)],
            &[450.0, 620.0],
        ),
        Err(AzothError::InvalidInput { .. })
    ));
}
