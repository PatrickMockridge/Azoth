//! Spec-driven tests for the `eos.effective_diffusion` model.
//!
//! The oracle is `validation/neqsim/EffectiveDiffusionProbe.java` and its capture
//! `captures/effective_diffusion_probe.tsv`, which holds the matrix, the vector the class
//! computes from it and the vector recomputed from it in one run.

use azoth_core::AzothError;
use azoth_eos::effective_diffusion::effective_diffusion;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.effective_diffusion";

fn rows(case: &azoth_core::spec::TestCase, name: &str) -> Vec<Vec<f64>> {
    let flat = case.matrix(name).expect("the case states the matrix");
    let side = (flat.len() as f64).sqrt() as usize;
    assert_eq!(side * side, flat.len(), "the matrix is square");
    flat.chunks(side).map(<[f64]>::to_vec).collect()
}

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::EffectiveDiffusionResult {
    effective_diffusion(
        &rows(case, "binary_diffusion"),
        case.vector("x")
            .expect("the case states the mole fractions"),
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
        let expected = case
            .expected_vector("effective_diffusion")
            .unwrap_or_else(|| panic!("the case `{}` declares effective_diffusion", case.id));
        for (i, (&got, &want)) in result.effective_diffusion.iter().zip(expected).enumerate() {
            common::assert_close(
                got,
                want,
                case.tolerance,
                &format!("{context} (effective_diffusion[{i}])"),
            );
        }
        common::assert_consistent(&result, context);
    }
}

/// **Only the row of the component being asked about is read**, which is what makes an
/// asymmetric matrix the normal case: swapping `D_01` alone moves `D_eff_0` and leaves
/// `D_eff_1` where it was.
#[test]
fn an_asymmetric_matrix_is_read_one_row_at_a_time() {
    let matrix = vec![vec![1.0e-9, 2.0e-9], vec![3.0e-9, 4.0e-9]];
    let x = [0.4, 0.6];
    let base = effective_diffusion(&matrix, &x).expect("the assembly runs");
    let swapped = effective_diffusion(&[vec![1.0e-9, 9.0e-9], vec![3.0e-9, 4.0e-9]], &x)
        .expect("the assembly runs");
    assert!(
        base.effective_diffusion[0] != swapped.effective_diffusion[0],
        "row 0 decides D_eff_0"
    );
    assert_eq!(
        base.effective_diffusion[1], swapped.effective_diffusion[1],
        "and row 1 is untouched by row 0's change"
    );
}

/// The three shape refusals and the two the arithmetic forces.
#[test]
fn the_refusals() {
    let square = vec![vec![1.0e-9, 2.0e-9], vec![3.0e-9, 4.0e-9]];
    assert!(matches!(
        effective_diffusion(&square, &[1.0]),
        Err(AzothError::InvalidInput { .. })
    ));
    assert!(matches!(
        effective_diffusion(&[vec![1.0e-9, 2.0e-9]], &[0.5, 0.5]),
        Err(AzothError::InvalidInput { .. })
    ));
    assert!(matches!(
        effective_diffusion(&[vec![1.0e-9], vec![3.0e-9, 4.0e-9]], &[0.5, 0.5]),
        Err(AzothError::InvalidInput { .. })
    ));
    assert!(
        effective_diffusion(&square, &[-0.1, 1.1]).is_err(),
        "a negative mole fraction is out of range"
    );
    assert!(
        matches!(
            effective_diffusion(&[vec![1.0e-9, 0.0], vec![3.0e-9, 4.0e-9]], &[0.5, 0.5]),
            Err(AzothError::OutOfRange { .. })
        ),
        "a zero pair coefficient is divided by"
    );
}
