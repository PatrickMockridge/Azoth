//! Spec-driven tests for the `characterization.combine_reservoir_fluids` model.
//!
//! Every expected value is a row of
//! `validation/neqsim/captures/characterization_fluid_rebuild_probe.tsv`.

use azoth_characterization::{
    CharacterizationComponent, CombineReservoirFluidsResult, FluidComponent,
    combine_reservoir_fluids, model_gen,
};
use azoth_core::spec::{FieldValue, RecordField, TestCase};
use azoth_test_support as common;

const MODEL_ID: &str = "characterization.combine_reservoir_fluids";

/// One declared flag of a record row.
fn flag(row: &[RecordField], name: &str) -> Option<bool> {
    row.iter()
        .find(|field| field.name == name)
        .and_then(|field| match field.value {
            FieldValue::Flag(value) => Some(value),
            FieldValue::Name(_) | FieldValue::Number(_) => None,
        })
}

/// One record as the kernel takes it.
fn component(row: &[RecordField]) -> FluidComponent {
    let number = |name: &str| common::record_number(row, name);
    FluidComponent {
        fluid: number("fluid").expect("every row states its fluid"),
        name: common::record_name(row, "name")
            .expect("every row is named")
            .to_string(),
        pseudo_component: flag(row, "pseudo_component").expect("every row states its kind"),
        moles: number("moles").expect("every row states its moles"),
        molar_mass: number("molar_mass"),
        density: number("density"),
        normal_boiling_point: number("normal_boiling_point"),
        critical_temperature: number("critical_temperature"),
        critical_pressure: number("critical_pressure"),
        acentric_factor: number("acentric_factor"),
        critical_volume: number("critical_volume"),
        racket_z: number("racket_z"),
        racket_z_cpa: number("racket_z_cpa"),
        parachor: number("parachor"),
        critical_viscosity: number("critical_viscosity"),
        triple_point_temperature: number("triple_point_temperature"),
        heat_of_fusion: number("heat_of_fusion"),
        ideal_gas_enthalpy_of_formation: number("ideal_gas_enthalpy_of_formation"),
        cp_a: number("cp_a"),
        cp_b: number("cp_b"),
        cp_c: number("cp_c"),
        cp_d: number("cp_d"),
        attractive_m: number("attractive_m"),
    }
}

fn call(case: &TestCase) -> CombineReservoirFluidsResult {
    let target = case
        .input("target_pseudo_components")
        .expect("every case states a target");
    let components: Vec<FluidComponent> = common::record_input(case, "components")
        .iter()
        .map(|row| component(row))
        .collect();
    combine_reservoir_fluids(target as i64, &components)
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

/// One property, compared against the case's own, absence included.
fn assert_optional(actual: Option<f64>, expected: Option<f64>, tolerance: f64, at: &str) {
    match (actual, expected) {
        (None, None) => {}
        (Some(actual), Some(expected)) => common::assert_close(actual, expected, tolerance, at),
        (actual, expected) => panic!("{at}: {actual:?} against {expected:?}"),
    }
}

/// Every numeric field of one row, in the order the spec declares them.
fn assert_row(expected: &[RecordField], row: &CharacterizationComponent, at: &str, tolerance: f64) {
    assert_eq!(
        row.name,
        common::record_name(expected, "name").expect("every row is named"),
        "{at}: the name"
    );
    assert_eq!(
        flag(expected, "pseudo_component"),
        Some(row.pseudo_component),
        "{at}: the kind"
    );
    common::assert_close(
        row.moles,
        common::record_number(expected, "moles").expect("every row states moles"),
        tolerance,
        &format!("{at}::moles"),
    );
    let fields: [(&str, Option<f64>); 19] = [
        ("molar_mass", row.molar_mass),
        ("density", row.density),
        ("normal_boiling_point", row.normal_boiling_point),
        ("critical_temperature", row.critical_temperature),
        ("critical_pressure", row.critical_pressure),
        ("acentric_factor", row.acentric_factor),
        ("critical_volume", row.critical_volume),
        ("racket_z", row.racket_z),
        ("racket_z_cpa", row.racket_z_cpa),
        ("parachor", row.parachor),
        ("critical_viscosity", row.critical_viscosity),
        ("triple_point_temperature", row.triple_point_temperature),
        ("heat_of_fusion", row.heat_of_fusion),
        (
            "ideal_gas_enthalpy_of_formation",
            row.ideal_gas_enthalpy_of_formation,
        ),
        ("cp_a", row.cp_a),
        ("cp_b", row.cp_b),
        ("cp_c", row.cp_c),
        ("cp_d", row.cp_d),
        ("attractive_m", row.attractive_m),
    ];
    for (name, actual) in fields {
        assert_optional(
            actual,
            common::record_number(expected, name),
            tolerance,
            &format!("{at}::{name}"),
        );
    }
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = format!("{}::{}", spec.id, case.id);
        let expected = common::expected_record(case, "components");
        assert_eq!(
            result.components.len(),
            expected.len(),
            "{context}: the component count"
        );
        for (index, (row, claims)) in result.components.iter().zip(expected).enumerate() {
            assert_row(
                claims,
                row,
                &format!("{context}::components[{index}]"),
                case.tolerance,
            );
        }
        common::assert_close(
            result.total_moles,
            case.expected_value("total_moles")
                .expect("every case states its total"),
            case.tolerance,
            &format!("{context}::total_moles"),
        );
        common::assert_consistent(&result, &context);
    }
}

#[test]
fn a_target_of_zero_is_refused() {
    let error = combine_reservoir_fluids(0, &[]).expect_err("a slate needs a component");
    assert!(
        error.to_string().contains("target_pseudo_components"),
        "the refusal should name the input: {error}"
    );
}

#[test]
fn a_pseudo_row_without_a_molar_mass_is_refused() {
    let row = FluidComponent {
        fluid: 0.0,
        name: "C7_PC".to_string(),
        pseudo_component: true,
        moles: 1.0,
        molar_mass: None,
        density: Some(730.0),
        normal_boiling_point: Some(370.0),
        critical_temperature: None,
        critical_pressure: None,
        acentric_factor: None,
        critical_volume: None,
        racket_z: None,
        racket_z_cpa: None,
        parachor: None,
        critical_viscosity: None,
        triple_point_temperature: None,
        heat_of_fusion: None,
        ideal_gas_enthalpy_of_formation: None,
        cp_a: None,
        cp_b: None,
        cp_c: None,
        cp_d: None,
        attractive_m: None,
    };
    let error = combine_reservoir_fluids(2, &[row]).expect_err("the mass is the row's");
    assert!(
        error.to_string().contains("components[0].molar_mass"),
        "the refusal should name the row's field: {error}"
    );
}

#[test]
fn a_fluid_index_that_is_not_a_whole_number_is_refused() {
    let row = FluidComponent {
        fluid: 0.5,
        name: "methane".to_string(),
        pseudo_component: false,
        moles: 1.0,
        molar_mass: None,
        density: None,
        normal_boiling_point: None,
        critical_temperature: None,
        critical_pressure: None,
        acentric_factor: None,
        critical_volume: None,
        racket_z: None,
        racket_z_cpa: None,
        parachor: None,
        critical_viscosity: None,
        triple_point_temperature: None,
        heat_of_fusion: None,
        ideal_gas_enthalpy_of_formation: None,
        cp_a: None,
        cp_b: None,
        cp_c: None,
        cp_d: None,
        attractive_m: None,
    };
    let error = combine_reservoir_fluids(2, &[row]).expect_err("half a fluid is no fluid");
    assert!(
        error.to_string().contains("components[0].fluid"),
        "the refusal should name the row's field: {error}"
    );
}
