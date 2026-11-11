//! Spec-driven tests for the `characterization.transfer_binary_interaction_parameters` model.
//!
//! The expected matrices are asserted by `python/tests/test_transfer_binary_interaction_parameters.py`,
//! which reads them from this model's raw spec: a case's `expected` block has no field for a
//! matrix, so the generator drops it and `test_cross_impl` carries the values to the other
//! kernel. What is here is the shape the answer must have and the two refusals - a missing
//! second base name, and a matrix that is not `N x N`.

use azoth_characterization::{
    InteractionComponent, model_gen, transfer_binary_interaction_parameters,
};
use azoth_core::spec::TestCase;
use azoth_test_support as common;

const MODEL_ID: &str = "characterization.transfer_binary_interaction_parameters";

/// One record as the kernel takes it.
fn component(row: &[azoth_core::spec::RecordField]) -> InteractionComponent {
    InteractionComponent {
        name: common::record_name(row, "name")
            .expect("every row is named")
            .to_string(),
        pseudo_component: row
            .iter()
            .find(|field| field.name == "pseudo_component")
            .and_then(|field| match field.value {
                azoth_core::spec::FieldValue::Flag(value) => Some(value),
                _ => None,
            })
            .expect("every row states its kind"),
    }
}

fn components(case: &TestCase, name: &str) -> Vec<InteractionComponent> {
    common::record_input(case, name)
        .iter()
        .map(|row| component(row))
        .collect()
}

fn reshape(flat: &[f64], n: usize) -> Vec<Vec<f64>> {
    flat.chunks(n).map(<[f64]>::to_vec).collect()
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let reference = components(case, "reference_components");
        let target = components(case, "target_components");
        let reference_kij = reshape(
            case.matrix("reference_kij")
                .expect("every case states the reference matrix"),
            reference.len(),
        );
        let target_kij = reshape(
            case.matrix("target_kij")
                .expect("every case states the target matrix"),
            target.len(),
        );
        let result = transfer_binary_interaction_parameters(
            &reference,
            &target,
            &reference_kij,
            &target_kij,
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));
        let context = format!("{}::{}", spec.id, case.id);

        assert_eq!(
            result.target_kij.len(),
            target.len(),
            "{context}: one row per target component"
        );
        for (row, values) in result.target_kij.iter().enumerate() {
            assert_eq!(
                values.len(),
                target.len(),
                "{context}: row {row} has one column per target component"
            );
            assert_eq!(values[row], 0.0, "{context}: the diagonal is never written");
        }
        for first in 0..target.len() {
            for second in (first + 1)..target.len() {
                assert_eq!(
                    result.target_kij[first][second], result.target_kij[second][first],
                    "{context}: the pair ({first}, {second}) is symmetric"
                );
            }
        }
        common::assert_consistent(&result, &context);
    }
}

#[test]
fn a_pair_whose_second_name_is_missing_is_refused() {
    // The capture's third transfer: the reference states methane-CO2 and the target has no
    // CO2. The class skips a missing *first* name and throws on the second, and the refusal
    // says which name it could not find.
    let reference = [
        InteractionComponent {
            name: "methane".to_string(),
            pseudo_component: false,
        },
        InteractionComponent {
            name: "CO2".to_string(),
            pseudo_component: false,
        },
    ];
    let target = [
        InteractionComponent {
            name: "methane".to_string(),
            pseudo_component: false,
        },
        InteractionComponent {
            name: "C7_PC".to_string(),
            pseudo_component: true,
        },
    ];
    let error = transfer_binary_interaction_parameters(
        &reference,
        &target,
        &[vec![0.0, 0.09], vec![0.09, 0.0]],
        &[vec![0.0, 0.0], vec![0.0, 0.0]],
    )
    .expect_err("the target has no CO2");
    assert!(
        error.to_string().contains("CO2") && error.to_string().contains("target_components"),
        "the refusal should name the component and the input: {error}"
    );
}

#[test]
fn a_reference_name_the_target_lacks_as_the_first_name_is_skipped() {
    // The class's own `continue`: a base component the target does not carry contributes
    // nothing, and that is not an error - provided nothing else states a pair with it, since
    // the reverse direction would look the same missing name up as the second one and refuse.
    let reference = [
        InteractionComponent {
            name: "propane".to_string(),
            pseudo_component: false,
        },
        InteractionComponent {
            name: "methane".to_string(),
            pseudo_component: false,
        },
    ];
    let target = [InteractionComponent {
        name: "methane".to_string(),
        pseudo_component: false,
    }];
    let result = transfer_binary_interaction_parameters(
        &reference,
        &target,
        &[vec![0.0, 0.0], vec![0.0, 0.0]],
        &[vec![0.0]],
    )
    .expect("the skipped name is not an error");
    assert_eq!(result.target_kij, vec![vec![0.0]]);
}

#[test]
fn a_matrix_that_is_not_square_is_refused() {
    let components = [InteractionComponent {
        name: "methane".to_string(),
        pseudo_component: false,
    }];
    let error = transfer_binary_interaction_parameters(
        &components,
        &components,
        &[vec![0.0, 0.0]],
        &[vec![0.0]],
    )
    .expect_err("a matrix is N x N");
    assert!(
        error.to_string().contains("reference_kij"),
        "the refusal should name the matrix: {error}"
    );
}
