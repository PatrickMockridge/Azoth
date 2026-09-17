//! Spec-driven tests for the `eos.unifac_activity_coefficients` model.

use azoth_core::{AzothError, CalcResult};
use azoth_eos::databank::{self, UnifacParameters};
use azoth_eos::{model_gen, unifac_activity_coefficients};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.unifac_activity_coefficients";

/// The methanol/water group tables the databank carries, which are also the doc-test's
/// and the refusal tests' stand-ins.
fn methanol_water() -> UnifacParameters {
    UnifacParameters {
        groups: vec![1.0, 0.0, 0.0, 1.0],
        group_r: vec![1.4311, 0.92],
        group_q: vec![1.432, 1.4],
        aij: vec![0.0, -181.0, 289.6, 0.0],
    }
}

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::UnifacActivityCoefficientsResult {
    let params = databank::unifac_parameters(case.list("components").expect("components"))
        .unwrap_or_else(|e| panic!("case `{}` should resolve but failed: {e}", case.id));
    unifac_activity_coefficients(
        &params,
        common::input(case, "T"),
        case.vector("x").expect("x"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        for (name, actual) in [("ln_gamma", &result.ln_gamma), ("gamma", &result.gamma)] {
            let expected = case
                .expected_vector(name)
                .unwrap_or_else(|| panic!("the case declares {name}"));
            for (i, (&got, &want)) in actual.iter().zip(expected).enumerate() {
                common::assert_close(
                    got,
                    want,
                    case.tolerance,
                    &format!("{context} ({name}[{i}])"),
                );
            }
        }
        common::assert_consistent(&result, context);
    }
}

#[test]
fn the_result_is_clean_at_an_ordinary_state() {
    let r = unifac_activity_coefficients(&methanol_water(), 298.15, &[0.5, 0.5]).unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

/// The resolution itself, which is what the migration moved into the databank. Every
/// case above would still pass against a case file regenerated from the same mistake,
/// so this is the test that pins the group decomposition and the union's ordering: the
/// columns are the named components' subgroups sorted by subgroup number, and a
/// resolver that sorted them any other way would permute the matrices.
#[test]
fn the_databank_resolves_the_methanol_water_group_tables() {
    assert_eq!(
        databank::unifac_parameters(&["methanol", "water"]).unwrap(),
        methanol_water()
    );
}

/// Methanol is subgroup 15 (CH3OH) and water subgroup 16 (H2O), so the union is
/// `[15, 16]` in that order: each component's row is a one in its own column, and
/// methanol - listed first - takes the first column.
///
/// The subgroup numbers are worth pinning because they are easy to mistake for the
/// *main* groups the interaction matrix is keyed by: methanol's main group is 6 and
/// water's is 7, and subgroups 6 and 7 are a different pair of groups entirely -
/// `R = 1.1167` and `1.1173`, both in main group 2. The resolver reads the subgroup from
/// `UNIFACcomp.csv` and the main group from `UNIFACGroupParam.csv`, which is why the two
/// numberings appear together in one test.
#[test]
fn the_group_union_is_sorted_by_subgroup_and_rows_follow_it() {
    let named = databank::unifac_parameters(&["methanol", "water"]).unwrap();
    assert_eq!(named.groups, vec![1.0, 0.0, 0.0, 1.0]);
    assert_eq!(named.group_r, vec![1.4311, 0.92]);

    let reversed = databank::unifac_parameters(&["water", "methanol"]).unwrap();
    assert_eq!(reversed.groups, vec![0.0, 1.0, 1.0, 0.0]);
    assert_eq!(reversed.group_r, named.group_r);
    assert_eq!(reversed.group_q, named.group_q);
    assert_eq!(reversed.aij, named.aij);
}

/// The union is a property of the *set* of components, not of the order they were
/// listed in, so the two orderings must give byte-identical tables.
///
/// This is the invariant the sort exists for. The columns of `groups` and the entries
/// of `group_r`, `group_q` and `aij` permute together, so an unsorted union still
/// computes the *same activity coefficients* - a sabotage that reverses the insertion
/// order leaves every case in the spec passing. What it breaks is determinism: the same
/// mixture named two ways would describe itself with two different matrices, which is
/// the kind of difference that only shows up when someone compares two runs.
#[test]
fn the_union_does_not_depend_on_the_order_the_components_are_named() {
    let named = databank::unifac_parameters(&["methanol", "water", "acetone"]).unwrap();
    let reversed = databank::unifac_parameters(&["acetone", "water", "methanol"]).unwrap();
    assert_eq!(named.group_r, reversed.group_r);
    assert_eq!(named.group_q, reversed.group_q);
    assert_eq!(named.aij, reversed.aij);
    // Each component's row moves with the names, so the two `groups` differ; what must
    // not is the ordering of the columns, which is what the two calls above compare.
    assert_ne!(named.groups, reversed.groups);
}

/// A name with no group decomposition is refused rather than given zero groups, which
/// would make its `R` and `Q` zero and its activity coefficient NaN.
#[test]
fn a_name_without_a_group_decomposition_is_refused() {
    let err = databank::unifac_parameters(&["methanol", "unobtainium"]).unwrap_err();
    assert!(
        matches!(err, AzothError::PropertyUnavailable { .. }),
        "{err:?}"
    );
}

#[test]
fn a_shape_mismatch_is_refused() {
    let mut params = methanol_water();
    params.groups = vec![1.0, 0.0, 0.0];
    let err = unifac_activity_coefficients(&params, 298.15, &[0.5, 0.5]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("components"));
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let err = unifac_activity_coefficients(&methanol_water(), 298.15, &[0.6, 0.6]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}
