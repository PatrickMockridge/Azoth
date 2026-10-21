//! Spec-driven tests for the `eos.friction_theory_conductivity` model.

use azoth_core::units::{
    cubic_meters_per_mole, joules_per_mole_kelvin, kelvins, kilograms_per_mole, pascal_seconds,
};
use azoth_eos::friction_theory_conductivity::friction_theory_conductivity;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.friction_theory_conductivity";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::FrictionTheoryConductivityResult {
    let vector = |name: &str| -> &'static [f64] {
        case.vector(name)
            .unwrap_or_else(|| panic!("case `{}` states {name}", case.id))
    };
    friction_theory_conductivity(
        &vector("cv0")
            .iter()
            .map(|v| joules_per_mole_kelvin(*v))
            .collect::<Vec<_>>(),
        &vector("tc").iter().map(|v| kelvins(*v)).collect::<Vec<_>>(),
        &vector("critical_volume")
            .iter()
            .map(|v| cubic_meters_per_mole(*v))
            .collect::<Vec<_>>(),
        vector("omega"),
        vector("dipole"),
        &vector("molar_mass")
            .iter()
            .map(|v| kilograms_per_mole(*v))
            .collect::<Vec<_>>(),
        vector("x"),
        kelvins(case.input("T").expect("the case states T")),
        pascal_seconds(case.input("viscosity").expect("the case states viscosity")),
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
        for (field, value) in [
            ("conductivity", result.conductivity.value),
            ("dilute_conductivity", result.dilute_conductivity.value),
            ("dense_conductivity", result.dense_conductivity.value),
        ] {
            common::assert_close(
                value,
                case.expected_value(field)
                    .unwrap_or_else(|| panic!("{context}: the case states {field}")),
                case.tolerance,
                &format!("{context} ({field})"),
            );
        }
        let expected = case
            .expected_vector("pure_component_conductivity")
            .expect("the case states the pure conductivities");
        assert_eq!(result.pure_component_conductivity.len(), expected.len());
        for (index, (value, want)) in result
            .pure_component_conductivity
            .iter()
            .zip(expected)
            .enumerate()
        {
            common::assert_close(
                *value,
                *want,
                case.tolerance,
                &format!("{context} (pure_component_conductivity[{index}])"),
            );
        }
        common::assert_consistent(&result, context);
    }
}

/// **The phase's viscosity is the caller's.** Halving it moves only the dense term, and
/// the dilute term - which Chung computes from the state, not the phase - does not move.
#[test]
fn the_viscosity_moves_only_the_dense_term() {
    let run = |viscosity: f64| {
        friction_theory_conductivity(
            &[joules_per_mole_kelvin(27.544151394)],
            &[kelvins(190.56)],
            &[cubic_meters_per_mole(9.9e-5)],
            &[0.0115],
            &[0.0],
            &[kilograms_per_mole(0.016043)],
            &[1.0],
            kelvins(300.0),
            pascal_seconds(viscosity),
        )
        .expect("the state runs")
    };
    let full = run(1.1403910448332615e-5);
    let half = run(0.5 * 1.1403910448332615e-5);
    assert_eq!(
        full.dilute_conductivity.value,
        half.dilute_conductivity.value
    );
    assert!(full.dense_conductivity.value > half.dense_conductivity.value);
}

#[test]
fn a_vector_of_the_wrong_length_is_refused() {
    let error = friction_theory_conductivity(
        &[joules_per_mole_kelvin(27.544151394)],
        &[kelvins(190.56)],
        &[cubic_meters_per_mole(9.9e-5)],
        &[0.0115],
        &[0.0],
        &[kilograms_per_mole(0.016043), kilograms_per_mole(0.058123)],
        &[1.0],
        kelvins(300.0),
        pascal_seconds(1.1403910448332615e-5),
    )
    .expect_err("one mole fraction cannot take two molar masses");
    assert!(
        error.to_string().contains("one entry per component"),
        "{error}"
    );
}
