//! Spec-driven tests for the `eos.interface_properties` model.

use azoth_core::units::{kelvins, kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_eos::interface_properties::interface_properties;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.interface_properties";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::InterfacePropertiesResult {
    interface_properties(
        case.string("interface")
            .expect("the case states the interface"),
        case.string("selection")
            .expect("the case states the selection"),
        case.vector("parachors")
            .expect("the case states the parachors"),
        &case
            .vector("tc")
            .expect("the case states the critical temperatures")
            .iter()
            .map(|value| kelvins(*value))
            .collect::<Vec<_>>(),
        kilograms_per_cubic_meter(case.input("rho_1").expect("the case states rho_1")),
        kilograms_per_cubic_meter(case.input("rho_2").expect("the case states rho_2")),
        kilograms_per_mole(
            case.input("molar_mass_1")
                .expect("the case states molar_mass_1"),
        ),
        kilograms_per_mole(
            case.input("molar_mass_2")
                .expect("the case states molar_mass_2"),
        ),
        case.vector("x_1").expect("the case states x_1"),
        case.vector("x_2").expect("the case states x_2"),
        kelvins(case.input("T").expect("the case states T")),
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

/// **The pair decides the slot.** The same selection on a liquid-liquid pair takes
/// Firozabadi-Ramley, and on a gas-oil pair Parachor; an oil/gas pair named as
/// `liquid_liquid` is the third slot, which is what NeqSim does with a reversed pair.
#[test]
fn the_pair_decides_the_slot() {
    let methane_butane = |interface: &str| {
        interface_properties(
            interface,
            "model_0",
            &[77.3, 191.7],
            &[kelvins(190.56), kelvins(425.12)],
            kilograms_per_cubic_meter(19.938328330315997),
            kilograms_per_cubic_meter(549.3851668858096),
            kilograms_per_mole(0.022959902730870334),
            kilograms_per_mole(0.054107691623917306),
            &[0.8356249351028912, 0.1643750648971088],
            &[0.09542082642782022, 0.9045791735721797],
            kelvins(300.0),
        )
        .expect("the state runs")
    };
    assert_eq!(
        methane_butane("gas_oil").selected_model,
        "ParachorSurfaceTension"
    );
    assert_eq!(
        methane_butane("liquid_liquid").selected_model,
        "FirozabadiRamleyInterfaceTension"
    );
    assert_ne!(
        methane_butane("gas_oil").sigma.value,
        methane_butane("liquid_liquid").sigma.value
    );
}

/// The three unported arms name their classes rather than approximating.
#[test]
fn the_solver_arms_are_refused_by_class() {
    for (selection, class) in [
        ("model_1", "GTSurfaceTension"),
        ("model_2", "LGTSurfaceTension"),
        ("model_4", "GTSurfaceTensionSimple"),
    ] {
        let error = interface_properties(
            "gas_oil",
            selection,
            &[77.3, 191.7],
            &[kelvins(190.56), kelvins(425.12)],
            kilograms_per_cubic_meter(19.938328330315997),
            kilograms_per_cubic_meter(549.3851668858096),
            kilograms_per_mole(0.022959902730870334),
            kilograms_per_mole(0.054107691623917306),
            &[0.8356249351028912, 0.1643750648971088],
            &[0.09542082642782022, 0.9045791735721797],
            kelvins(300.0),
        )
        .expect_err("the solver arms are refused");
        let message = error.to_string();
        assert!(
            message.contains(class),
            "the refusal names {class}, got {message}"
        );
    }
}
