//! Spec-driven tests for the `eos.gt_surface_tension` model.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::gt_surface_tension::gt_surface_tension;
use azoth_eos::{Cubic, databank, model_gen};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.gt_surface_tension";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::GtSurfaceTensionResult {
    let names = case
        .list("components")
        .expect("the case declares components");
    let cubic = match case.string("eos").unwrap_or("srk") {
        "pr" => Cubic::Pr,
        "srk" => Cubic::Srk,
        "rk" => Cubic::Rk,
        other => panic!("case `{}` names an unhandled cubic `{other}`", case.id),
    };
    let (mixture, _) =
        databank::mixture_of(names, cubic, None).expect("the case's components resolve");
    gt_surface_tension(
        &mixture,
        case.vector("influence_parameters")
            .expect("the case states the influence parameters"),
        kelvins(case.input("T").expect("the case states T")),
        pascals(case.input("P").expect("the case states P")),
        case.vector("z").expect("the case states z"),
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
        common::assert_close(
            result.sigma.value,
            case.expected_value("sigma").expect("the case states sigma"),
            case.tolerance,
            &format!("{context} (sigma)"),
        );
        common::assert_consistent(&result, context);
    }
}

#[test]
fn a_single_phase_feed_is_refused() {
    let (mixture, _) =
        databank::mixture_of(&["methane"], Cubic::Srk, None).expect("methane resolves");
    let error = gt_surface_tension(
        &mixture,
        &[9.634812734820117e-21],
        kelvins(400.0),
        pascals(1.0e5),
        &[1.0],
    )
    .expect_err("a supercritical single phase has no interface");
    assert!(error.to_string().contains("single phase"), "{error}");
}
