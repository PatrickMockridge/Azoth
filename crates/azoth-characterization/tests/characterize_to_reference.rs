//! Spec-driven tests for the `characterization.characterize_to_reference` model.
//!
//! Every expected value is a row of
//! `validation/neqsim/captures/pseudo_component_combiner_probe.tsv`.

use azoth_characterization::{
    CharacterizeToReferenceResult, ReferenceRow, SourceRow, characterize_to_reference, model_gen,
};
use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_core::{AzothError, spec::TestCase};
use azoth_test_support as common;

const MODEL_ID: &str = "characterization.characterize_to_reference";

/// The source table, as the kernel takes it.
///
/// A field the case omits is a row that states none, which is what the record exists for:
/// `record_number` returns `None` rather than failing and the `Option` carries it through.
fn source_of(case: &TestCase) -> Vec<SourceRow> {
    common::record_input(case, "source")
        .iter()
        .map(|row| SourceRow {
            moles: common::record_number(row, "moles").expect("every source row states its moles"),
            molar_mass: kilograms_per_mole(
                common::record_number(row, "molar_mass").expect("every source row states a mass"),
            ),
            density: common::record_number(row, "density").map(kilograms_per_cubic_meter),
            boiling_point: common::record_number(row, "boiling_point"),
        })
        .collect()
}

fn reference_of(case: &TestCase) -> Vec<ReferenceRow> {
    common::record_input(case, "reference")
        .iter()
        .map(|row| ReferenceRow {
            molar_mass: kilograms_per_mole(
                common::record_number(row, "molar_mass")
                    .expect("every reference row states a mass"),
            ),
            boiling_point: common::record_number(row, "boiling_point"),
            density: common::record_number(row, "density").map(kilograms_per_cubic_meter),
        })
        .collect()
}

fn call(case: &TestCase) -> CharacterizeToReferenceResult {
    characterize_to_reference(
        &source_of(case),
        &reference_of(case),
        case.flag("inherit_reference_properties"),
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

/// A source row from four numbers, the last two absent when `None`.
///
/// The kernel takes rows, so a test states one - and a test that wants a row without a density or
/// without a boiling point states `None`, which is the absence a parallel vector had no form for.
fn row(moles: f64, molar_mass: f64, density: Option<f64>, boiling_point: Option<f64>) -> SourceRow {
    SourceRow {
        moles,
        molar_mass: kilograms_per_mole(molar_mass),
        density: density.map(kilograms_per_cubic_meter),
        boiling_point,
    }
}

/// A reference row from two numbers, the boiling point absent when `None`.
fn ref_row(molar_mass: f64, boiling_point: Option<f64>) -> ReferenceRow {
    ReferenceRow {
        molar_mass: kilograms_per_mole(molar_mass),
        boiling_point,
        density: None,
    }
}

#[test]
fn a_key_exactly_on_a_boundary_stays_in_the_lower_group() {
    // The walk advances while the key is **greater than** the boundary, so a key of exactly 550.0
    // against a boundary of 550.0 stays in the lower group - and the upper group, having nothing
    // left, is then dropped rather than reported empty. A port that used `>=` would return two
    // groups here instead of one.
    let run = |key: f64| {
        characterize_to_reference(
            &[
                row(1.0, 0.1, Some(700.0), Some(380.0)),
                row(1.0, 0.2, Some(820.0), Some(key)),
            ],
            &[ref_row(0.15, Some(400.0)), ref_row(0.25, Some(700.0))],
            None,
        )
        .expect("computes")
    };
    assert_eq!(run(550.0).reference_index, vec![0.0]);
    assert_eq!(run(550.001).reference_index, vec![0.0, 1.0]);
    assert_eq!(run(550.0).group_moles, vec![2.0]);
}

#[test]
fn a_row_without_a_boiling_point_sorts_by_its_molar_mass() {
    // **The absence a vector could only spell as a sentinel.** `sortingKey` falls back to the molar
    // mass when the boiling point is not finite and positive, and a caller expressing that in a
    // vector had to write `0.0` or `NaN` into it. Here the field is simply left out, and the answer
    // has to be the one the sentinel gave: this row sorts between the two reference keys, 400 and
    // 700 K, so it lands in the lower group.
    let result = characterize_to_reference(
        &[
            row(1.0, 0.1, Some(700.0), Some(380.0)),
            row(1.0, 0.2, None, None),
        ],
        &[ref_row(0.15, Some(400.0)), ref_row(0.25, Some(700.0))],
        None,
    )
    .expect("computes");
    assert_eq!(result.reference_index, vec![0.0], "0.2 kg/mol sorts below 400 K's own key");
    assert_eq!(result.group_moles, vec![2.0]);

    // And the other half of the same field: a row that states no density adds its mass to the
    // group's numerator and no volume to its denominator, which is the class's own arithmetic and
    // not a correction of it - so the group's density is *higher* than the row that did state one.
    assert!((result.group_density[0].value - 0.3 / (0.1 / 700.0)).abs() < 1e-9);
}

#[test]
fn the_rows_are_sorted_rather_than_taken_as_given() {
    // `extractComponents` sorts both tables by key, so a source stated in any order gives the same
    // answer - which is what makes the groups contiguous runs rather than whatever order the
    // caller happened to use.
    let sorted = |order: [usize; 3]| {
        let table = [
            row(1.0, 0.1, Some(700.0), Some(380.0)),
            row(1.0, 0.2, Some(820.0), Some(500.0)),
            row(1.0, 0.3, Some(880.0), Some(700.0)),
        ];
        let source: Vec<SourceRow> = order.iter().map(|i| table[*i]).collect();
        characterize_to_reference(&source, &[ref_row(0.15, Some(450.0)), ref_row(0.28, Some(620.0))], None)
        .expect("computes")
    };
    assert_eq!(sorted([0, 1, 2]), sorted([2, 1, 0]));
    assert_eq!(sorted([0, 1, 2]), sorted([1, 2, 0]));
}

#[test]
fn a_zero_amount_row_is_not_on_the_grid() {
    // `extractComponents` skips a component with no moles, so a zero here is not a row that lands
    // in a group and adds nothing - it is a row the walk never sees.
    let result = characterize_to_reference(
        &[
            row(1.0, 0.1, Some(700.0), Some(380.0)),
            row(0.0, 0.2, Some(820.0), Some(500.0)),
            row(1.0, 0.3, Some(880.0), Some(700.0)),
        ],
        &[ref_row(0.15, Some(450.0)), ref_row(0.28, Some(620.0))],
        None,
    )
    .expect("computes");
    assert_eq!(result.group_moles.iter().sum::<f64>(), 2.0);
}

#[test]
fn the_inheriting_variant_changes_two_numbers_and_not_the_walk() {
    // The capture's first two blocks: the same fluid twice, once with the flag. The groups and
    // their moles are identical, and each molar mass and density is the reference row's own.
    let source = [
        row(1.0, 0.10, Some(700.0), Some(378.77)),
        row(1.5, 0.14, Some(760.0), Some(444.67)),
        row(2.0, 0.18, Some(800.0), Some(502.058)),
        row(1.0, 0.24, Some(840.0), Some(573.86)),
        row(0.5, 0.32, Some(880.0), Some(646.972)),
    ];
    let reference = vec![
        ReferenceRow {
            molar_mass: kilograms_per_mole(0.12),
            boiling_point: Some(412.832),
            density: Some(kilograms_per_cubic_meter(740.0)),
        },
        ReferenceRow {
            molar_mass: kilograms_per_mole(0.20),
            boiling_point: Some(527.8),
            density: Some(kilograms_per_cubic_meter(820.0)),
        },
        ReferenceRow {
            molar_mass: kilograms_per_mole(0.30),
            boiling_point: Some(630.83),
            density: Some(kilograms_per_cubic_meter(870.0)),
        },
    ];

    let own = characterize_to_reference(&source, &reference, None).expect("computes");
    let inherited = characterize_to_reference(&source, &reference, Some(true)).expect("computes");

    // The walk is the same walk: the groups and their moles are untouched.
    assert_eq!(own.reference_index, inherited.reference_index);
    assert_eq!(own.group_moles, inherited.group_moles);
    // And the two numbers the flag governs are the reference rows' own.
    let masses: Vec<f64> = inherited.group_molar_mass.iter().map(|q| q.value).collect();
    assert_eq!(masses, vec![0.12, 0.20, 0.30]);
    let got: Vec<f64> = inherited.group_density.iter().map(|q| q.value).collect();
    assert_eq!(got, vec![740.0, 820.0, 870.0]);
    // The un-flagged answer is the first capture block, so the two are not the same answer.
    let own_masses: Vec<f64> = own.group_molar_mass.iter().map(|q| q.value).collect();
    assert_ne!(own_masses, masses);
}

#[test]
fn the_inheriting_variant_refuses_a_reference_row_with_no_density() {
    // The reference's densities are read by that path and not by the other, so a slate that states
    // none is answerable one way and not the other - which is a refusal rather than a silent zero.
    let source = [row(1.0, 0.10, Some(700.0), Some(378.77))];
    let bare = vec![
        ReferenceRow {
            molar_mass: kilograms_per_mole(0.12),
            boiling_point: Some(412.832),
            density: None,
        },
        ReferenceRow {
            molar_mass: kilograms_per_mole(0.20),
            boiling_point: Some(527.8),
            density: None,
        },
    ];
    assert!(characterize_to_reference(&source, &bare, None).is_ok());
    assert!(matches!(
        characterize_to_reference(&source, &bare, Some(true)),
        Err(AzothError::InvalidInput { .. })
    ));
}
