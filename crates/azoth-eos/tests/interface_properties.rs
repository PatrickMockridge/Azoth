//! Spec-driven tests for the `eos.interface_properties` model.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::interface_properties::interface_properties;
use azoth_eos::{Cubic, databank, model_gen};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.interface_properties";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::InterfacePropertiesResult {
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
    interface_properties(
        &mixture,
        case.string("selection")
            .expect("the case states the selection"),
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
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
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
        assert_eq!(
            result.selected_model,
            case.expected_string("selected_model")
                .expect("the case states selected_model"),
            "{context}: the slot's class"
        );
        common::assert_consistent(&result, context);
    }
}

/// **The pair decides the slot.** Model 0 on a gas-oil feed takes Parachor; on an
/// oil-water feed it takes Firozabadi-Ramley, which is the liquid-liquid slot.
#[test]
fn the_pair_decides_the_slot() {
    let spec = model_gen::model(MODEL_ID).expect("the model");
    let case = |id: &str| {
        spec.cases
            .iter()
            .find(|case| case.id == id)
            .expect("the case is in the spec")
    };
    let gas_oil = call(case("model_0_on_a_gas_oil_pair"));
    assert_eq!(gas_oil.selected_model, "ParachorSurfaceTension");
    let liquid_liquid = call(case("model_0_on_a_liquid_liquid_pair"));
    assert_eq!(
        liquid_liquid.selected_model,
        "FirozabadiRamleyInterfaceTension"
    );
    assert_ne!(gas_oil.sigma.value, liquid_liquid.sigma.value);
}

/// The profile arms propagate their calculator's refusal rather than answering a
/// partial integral: model 1's gas-aqueous slot is GT-simple, whose profile collapses on
/// CO2/water at this state.
#[test]
fn a_profile_arms_refusal_is_the_calculators() {
    let (mixture, _) =
        databank::mixture_of(&["CO2", "water"], Cubic::Srk, None).expect("both resolve");
    let error = interface_properties(
        &mixture,
        "model_1",
        &[3.453431599008472e-20, 1.9939776242117276e-20],
        1.0,
        &[vec![0.0, 0.35], vec![0.35, 0.0]],
        kelvins(313.15),
        pascals(5.0e6),
        &[0.02, 0.98],
    )
    .expect_err("GT-simple's profile collapses on this state, and the arm carries that");
    assert!(
        error
            .to_string()
            .contains("GT-simple profile's density vector is not a state"),
        "{error}"
    );
}

#[test]
fn a_single_phase_feed_is_refused() {
    let (mixture, _) =
        databank::mixture_of(&["methane"], Cubic::Srk, None).expect("methane resolves");
    let error = interface_properties(
        &mixture,
        "model_0",
        &[9.634812734820117e-21],
        0.0,
        &[vec![0.0]],
        kelvins(400.0),
        pascals(1.0e5),
        &[1.0],
    )
    .expect_err("one supercritical phase has no pair to dispatch");
    assert!(error.to_string().contains("no pair to dispatch"), "{error}");
}
