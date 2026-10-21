//! Spec-driven tests for the `eos.cdft_surface_tension` model.

use azoth_core::units::kelvins;
use azoth_eos::cdft_surface_tension::cdft_surface_tension;
use azoth_eos::{Cubic, databank, model_gen};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.cdft_surface_tension";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::CdftSurfaceTensionResult {
    let names = case
        .list("components")
        .expect("the case declares components");
    let cubic = match case.string("eos").unwrap_or("srk") {
        "pr" => Cubic::Pr,
        "srk" => Cubic::Srk,
        other => panic!("case `{}` names an unhandled cubic `{other}`", case.id),
    };
    let (mixture, _) =
        databank::mixture_of(names, cubic, None).expect("the case's components resolve");
    cdft_surface_tension(
        &mixture,
        kelvins(case.input("T").expect("the case states T")),
        case.vector("x_liquid").expect("the case states x_liquid"),
        case.vector("x_vapour").expect("the case states x_vapour"),
        case.input("molar_volume_liquid")
            .expect("the case states molar_volume_liquid"),
        case.input("molar_volume_vapour")
            .expect("the case states molar_volume_vapour"),
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
fn a_clamped_minimum_is_refused() {
    let (mixture, _) =
        databank::mixture_of(&["n-heptane", "water"], Cubic::Srk, None).expect("both resolve");
    let error = cdft_surface_tension(
        &mixture,
        kelvins(300.0),
        &[9.208893557273499e-26, 1.0],
        &[0.9995963504210494, 4.0364957895044424e-4],
        2.3850393536458033e-5,
        1.681108733850542e-4,
    )
    .expect_err("the class clamps this state's answer to zero, and the port refuses it");
    assert!(error.to_string().contains("clamps"), "{error}");
}

#[test]
fn a_density_vector_of_the_wrong_length_is_refused() {
    let (mixture, _) =
        databank::mixture_of(&["methane", "n-butane"], Cubic::Srk, None).expect("both resolve");
    let error = cdft_surface_tension(
        &mixture,
        kelvins(300.0),
        &[1.0],
        &[0.5, 0.5],
        1.0e-4,
        1.0e-2,
    )
    .expect_err("two components cannot take a one-entry liquid composition");
    assert!(
        error.to_string().contains("one entry per component"),
        "{error}"
    );
}
