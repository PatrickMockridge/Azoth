//! Spec-driven tests for the `eos.friction_theory_viscosity` model.

use azoth_core::units::{
    cubic_meters_per_mole, kelvins, kilograms_per_mole, pascal_seconds, pascals,
};
use azoth_eos::friction_theory_viscosity::friction_theory_viscosity;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.friction_theory_viscosity";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::FrictionTheoryViscosityResult {
    let vector = |name: &str| -> &'static [f64] {
        case.vector(name)
            .unwrap_or_else(|| panic!("case `{}` states {name}", case.id))
    };
    friction_theory_viscosity(
        &vector("tc").iter().map(|v| kelvins(*v)).collect::<Vec<_>>(),
        &vector("pc").iter().map(|v| pascals(*v)).collect::<Vec<_>>(),
        &vector("critical_volume")
            .iter()
            .map(|v| cubic_meters_per_mole(*v))
            .collect::<Vec<_>>(),
        vector("omega"),
        &vector("molar_mass")
            .iter()
            .map(|v| kilograms_per_mole(*v))
            .collect::<Vec<_>>(),
        &vector("critical_viscosity")
            .iter()
            .map(|v| pascal_seconds(*v))
            .collect::<Vec<_>>(),
        vector("viscosity_correction_factor"),
        vector("dipole"),
        vector("x"),
        kelvins(case.input("T").expect("the case states T")),
        pascals(
            case.input("p_repulsive")
                .expect("the case states p_repulsive"),
        ),
        pascals(
            case.input("p_attractive")
                .expect("the case states p_attractive"),
        ),
        case.flag("peng_robinson")
            .expect("the case states peng_robinson"),
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
            ("mu", result.mu.value),
            ("dilute_viscosity", result.dilute_viscosity.value),
            ("friction_viscosity", result.friction_viscosity.value),
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
            .expected_vector("pure_component_viscosity")
            .expect("the case states the pure viscosities");
        assert_eq!(result.pure_component_viscosity.len(), expected.len());
        for (index, (value, want)) in result
            .pure_component_viscosity
            .iter()
            .zip(expected)
            .enumerate()
        {
            common::assert_close(
                *value,
                *want,
                case.tolerance,
                &format!("{context} (pure_component_viscosity[{index}])"),
            );
        }
        common::assert_consistent(&result, context);
    }
}

/// **The two constant sets are not interchangeable**: the same fluid on PR and SRK
/// answers `1.3592e-4` and `1.3835e-4` Pa*s, and the PR composition differs too. A port
/// that ignored `peng_robinson` would still pass one case and fail the other.
#[test]
fn the_constant_set_moves_the_answer() {
    let run = |peng_robinson: bool, x_methane: f64, x_butane: f64, p_rep: f64, p_att: f64| {
        friction_theory_viscosity(
            &[kelvins(190.56), kelvins(425.12)],
            &[pascals(4.599e6), pascals(3.796e6)],
            &[
                cubic_meters_per_mole(9.9e-5),
                cubic_meters_per_mole(2.55e-4),
            ],
            &[0.0115, 0.2002],
            &[kilograms_per_mole(0.016043), kilograms_per_mole(0.058123)],
            &[pascal_seconds(1.5293e-5), pascal_seconds(2.5768e-5)],
            &[0.0, 0.0],
            &[0.0, 0.0],
            &[x_methane, x_butane],
            kelvins(300.0),
            pascals(p_rep),
            pascals(p_att),
            peng_robinson,
        )
        .expect("the state runs")
    };
    let srk = run(
        false,
        0.09542082642782022,
        0.9045791735721797,
        829.1753495002592,
        -809.1753495002592,
    );
    let pr = run(
        true,
        0.09718095876744977,
        0.9028190412325503,
        977.1223664955916,
        -957.1223664955916,
    );
    assert!((srk.mu.value - 1.383523806147524e-4).abs() < 1e-16);
    assert!((pr.mu.value - 1.3592326566884292e-4).abs() < 1e-16);
    assert!((srk.mu.value - pr.mu.value).abs() > 1e-6);
}

/// The floor: a state whose two terms would cancel below `1e-6` Pa*s answers the floor.
#[test]
fn the_answer_is_floored() {
    let result = friction_theory_viscosity(
        &[kelvins(190.56)],
        &[pascals(4.599e6)],
        &[cubic_meters_per_mole(9.9e-5)],
        &[0.0115],
        &[kilograms_per_mole(0.016043)],
        &[pascal_seconds(1.5293e-5)],
        &[0.0],
        &[0.0],
        &[1.0],
        kelvins(300.0),
        pascals(10.332089673913584),
        pascals(-0.33208967391358435),
        true,
    )
    .expect("the pure state runs");
    assert!(result.mu.value >= 1.0e-6);
}
