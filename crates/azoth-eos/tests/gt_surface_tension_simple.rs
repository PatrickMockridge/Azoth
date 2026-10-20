//! Spec-driven tests for the `eos.gt_surface_tension_simple` model.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::gt_surface_tension_simple::gt_surface_tension_simple;
use azoth_eos::{Cubic, databank, model_gen};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.gt_surface_tension_simple";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::GtSurfaceTensionSimpleResult {
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
    let count = names.len();
    let flat = case
        .matrix("pair_interactions")
        .expect("the case states the interaction matrix");
    let interactions: Vec<Vec<f64>> = flat.chunks(count).map(<[f64]>::to_vec).collect();
    gt_surface_tension_simple(
        &mixture,
        case.vector("influence_parameters")
            .expect("the case states the influence parameters"),
        case.input("reference_component")
            .expect("the case states the reference component"),
        &interactions,
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

/// **The class's own CO2/water state is refused.** Its profile collapses at step 3 -
/// `den_interface` is zero from there - and it returns the two-step partial integral
/// `0.0145737553` without a word; a partial integral is not a surface tension.
#[test]
fn a_collapsing_profile_is_refused() {
    let (mixture, _) =
        databank::mixture_of(&["CO2", "water"], Cubic::Srk, None).expect("the pair resolves");
    let error = gt_surface_tension_simple(
        &mixture,
        &[3.453431599008472e-20, 1.9939776242117276e-20],
        1.0,
        &[vec![0.0, 0.0], vec![0.0, 0.0]],
        kelvins(313.15),
        pascals(5.0e6),
        &[0.02, 0.98],
    )
    .expect_err("a collapsed profile is not a state");
    assert!(
        error.to_string().contains("not a state")
            || error.to_string().contains("integrand")
            || error.to_string().contains("ln(0)"),
        "{error}"
    );
}

/// A single-phase feed has no interface; the class would walk a profile between a phase
/// and itself.
#[test]
fn a_single_phase_feed_is_refused() {
    let (mixture, _) =
        databank::mixture_of(&["methane"], Cubic::Srk, None).expect("methane resolves");
    let error = gt_surface_tension_simple(
        &mixture,
        &[9.634812734820117e-21],
        0.0,
        &[vec![0.0]],
        kelvins(400.0),
        pascals(1.0e5),
        &[1.0],
    )
    .expect_err("a supercritical single phase has no interface");
    assert!(error.to_string().contains("single phase"), "{error}");
}
