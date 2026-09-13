//! Spec-driven tests for `hydraulics.crane_k_factors`.
//!
//! # What these tests can and cannot check
//!
//! The coefficients in `data/fittings/crane_k_factors.csv` are estimated dummy
//! values, not engineering data. There is therefore no correct answer for this
//! calc to be checked against, and nothing here can detect a wrong coefficient.
//!
//! These tests validate the arithmetic, the registry lookup and the error
//! handling. They are honest about the gap: the numerical assertions compare
//! against values derived from the same dummy coefficients, so they would all
//! still pass if the coefficients were nonsense.

use azoth_test_support as common;

use azoth_core::AzothError;
use azoth_core::spec::TestCase;
use azoth_hydraulics::{crane_k_factors, fittings, known_fittings, spec_gen};

const CALC_ID: &str = "hydraulics.crane_k_factors";

fn call(case: &TestCase) -> azoth_hydraulics::KFactorsResult {
    crane_k_factors(
        common::list_input(case, "fittings"),
        common::input(case, "f_t"),
    )
    .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    common::assert_skips_are_explained(spec);

    let mut executed = 0;
    for case in spec.all_tests() {
        if !case.is_active() {
            continue;
        }
        match case.kind {
            "worked_example" | "reference" => {
                let result = call(case);
                common::assert_close(
                    result.k_total,
                    common::expected(case, "k_total"),
                    case.tolerance,
                    &format!("{}::{} (k_total)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));

                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |q| match q {
                        "f_t" => Some(common::input(case, "f_t")),
                        "n_fittings" => Some(result.components.len() as f64),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
            }
            "property" => match case.property {
                Some("symmetry") => symmetry(),
                other => panic!("{}::{}: unhandled property {other:?}", spec.id, case.id),
            },
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 3,
        "expected several active cases, ran {executed}"
    );
}

/// K is linear in `f_t` and additive over fittings.
#[test]
fn symmetry() {
    let fittings = ["90_elbow", "gate_valve_open"];

    // Linear in f_t: doubling it must exactly double K.
    let single = crane_k_factors(&fittings, 0.018).unwrap().k_total;
    let double = crane_k_factors(&fittings, 0.036).unwrap().k_total;
    common::assert_close(double, 2.0 * single, 1e-15, "K is linear in f_t");

    // Order-independent: reversing the list must not change the total.
    let reversed: Vec<&str> = fittings.iter().rev().copied().collect();
    let reversed_k = crane_k_factors(&reversed, 0.018).unwrap().k_total;
    common::assert_close(reversed_k, single, 1e-15, "K is order-independent");

    // Additive: a fitting's contribution is independent of its neighbours.
    let a = crane_k_factors(&["90_elbow"], 0.018).unwrap().k_total;
    let b = crane_k_factors(&["gate_valve_open"], 0.018)
        .unwrap()
        .k_total;
    common::assert_close(single, a + b, 1e-15, "K is additive over fittings");
}

#[test]
fn unknown_fitting_is_an_error_not_a_zero() {
    // Silently treating an unknown fitting as zero loss would under-report
    // pressure drop, which is the dangerous direction to be wrong in.
    let err = crane_k_factors(&["90_elbow", "no_such_fitting"], 0.018).unwrap_err();
    assert!(matches!(err, AzothError::UnknownFitting { .. }));
    assert!(err.to_string().contains("no_such_fitting"));
}

#[test]
fn empty_fitting_list_is_an_error() {
    // An empty list would return K = 0, which reads as "no fittings" when the
    // caller may have meant to supply some.
    let err = crane_k_factors(&[], 0.018).unwrap_err();
    assert_eq!(err.field(), Some("n_fittings"));
}

#[test]
fn non_positive_friction_factor_is_an_error() {
    let err = crane_k_factors(&["90_elbow"], 0.0).unwrap_err();
    assert_eq!(err.field(), Some("f_t"));
}

#[test]
fn components_echo_the_registry_coefficients() {
    // The per-fitting breakdown must come from the registry rather than being
    // recomputed, so a caller can see exactly which coefficients were used.
    let result = crane_k_factors(&["90_elbow", "gate_valve_open"], 0.018).unwrap();
    for component in &result.components {
        let row = fittings::find(&component.fitting_id).unwrap();
        assert_eq!(component.n_ld, row.n_ld, "n_ld must come from the registry");
        common::assert_close(
            component.k,
            0.018 * row.n_ld,
            1e-15,
            &format!("k for {}", component.fitting_id),
        );
    }
    let expected_total: f64 = result.components.iter().map(|c| c.k).sum();
    common::assert_close(result.k_total, expected_total, 1e-15, "k_total is the sum");
}

#[test]
fn the_registry_is_entirely_estimated_dummy_data_right_now() {
    // This test exists to fail loudly the day someone populates the registry
    // from a real source. At that point this assertion breaks, and whoever did
    // the work is forced to update the spec's verification notes and the docs
    // rather than leaving them claiming the data is placeholder. Deleting this
    // test is the correct response to that failure, not updating the expected
    // number.
    let rows = fittings::registry().unwrap();
    let estimated: Vec<&str> = rows
        .iter()
        .filter(|r| r.is_estimated())
        .map(|r| r.id.as_str())
        .collect();
    assert_eq!(
        estimated.len(),
        rows.len(),
        "some rows are no longer estimated_dummy ({:?} of {}). If the coefficients have \
         now been verified against the primary standard, delete this test and update the \
         spec's verification block.",
        estimated,
        rows.len()
    );
}

#[test]
fn known_fittings_lists_the_registry() {
    let ids = known_fittings().unwrap();
    assert!(ids.contains(&"90_elbow"));
    assert!(ids.contains(&"gate_valve_open"));
    assert_eq!(ids.len(), fittings::registry().unwrap().len());
}
