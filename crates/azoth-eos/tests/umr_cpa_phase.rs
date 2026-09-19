//! Spec-driven tests for the `eos.umr_cpa_phase` model.

use azoth_eos::{model_gen, umr_cpa_phase};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.umr_cpa_phase";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::UmrCpaPhaseResult {
    let names: Vec<String> = azoth_test_support::list_input(case, "components")
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    umr_cpa_phase(
        &names,
        azoth_core::units::kelvins(common::input(case, "T")),
        azoth_core::units::pascals(common::input(case, "P")),
        case.vector("z").expect("z"),
        case.string("compressed_phase").expect("compressed_phase"),
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
        let tolerance = case.tolerance;
        common::assert_close(
            result.z_factor,
            case.expected_value("z_factor").expect("z_factor"),
            tolerance,
            &format!("{context} (z_factor)"),
        );
        let ln_phi = case.expected_vector("ln_phi").expect("ln_phi");
        for (i, (&got, &want)) in result.ln_phi.iter().zip(ln_phi).enumerate() {
            common::assert_close(got, want, tolerance, &format!("{context} (ln_phi[{i}])"));
        }
        if let Some(want) = case.expected_value("h_res") {
            common::assert_close(
                result.h_res.value,
                want,
                tolerance,
                &format!("{context} (h_res)"),
            );
        }
        if let Some(want) = case.expected_value("s_res") {
            common::assert_close(
                result.s_res.value,
                want,
                tolerance,
                &format!("{context} (s_res)"),
            );
        }
        common::assert_consistent(&result, context);
    }
}

/// The model refuses a component that carries no `UMRCPA_MC1..5` set.
///
/// **A refusal rather than a fallback, and the reason is structural.** NeqSim's
/// `ComponentUMRCPA.setAttractiveTerm` chooses per component - term 22 where the set is,
/// term 19 seeded with `MCPR1..3` for a non-associating component without one and term 1
/// with `mCPA` for an associating one - and this library states the alpha once for the
/// mixture. A mixture of a component that carries the set and one that does not would need
/// two attractive terms at once, so the model reads the set it can express.
#[test]
fn a_component_without_the_umr_cpa_alpha_is_refused() {
    // Methane carries `UMRCPA_MC1..5`; n-pentane carries none, and neither associates.
    let err = umr_cpa_phase(
        &["methane".to_string(), "n-pentane".to_string()],
        azoth_core::units::kelvins(298.15),
        azoth_core::units::pascals(70.0e5),
        &[0.5, 0.5],
        "vapour",
    )
    .expect_err("a component without the set should be refused");
    assert!(
        matches!(err, azoth_core::AzothError::PropertyUnavailable { .. }),
        "{err:?}"
    );
}
