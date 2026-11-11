use azoth_core::provenance_gen::ALL_PROVENANCE;
use azoth_sheets::catalogue::{forms, to_json};
use azoth_sheets::form::CalcForm;

fn available() -> Vec<&'static CalcForm> {
    forms().iter().filter(|form| form.available).collect()
}

#[test]
fn every_calculation_is_listed_once_and_in_order() {
    let ids: Vec<&str> = forms().iter().map(|form| form.id).collect();
    assert_eq!(ids.len(), 98, "the registry is 98 calculations");
    let mut sorted = ids.clone();
    sorted.sort_unstable();
    assert_eq!(ids, sorted, "ids are in id order");
    let unique: std::collections::BTreeSet<&str> = ids.iter().copied().collect();
    assert_eq!(unique.len(), ids.len(), "no id appears twice");
}

#[test]
fn the_table_shaped_calculations_are_the_four_that_are_refused() {
    let refused: Vec<&str> = forms()
        .iter()
        .filter(|form| !form.available)
        .map(|form| form.id)
        .collect();
    assert_eq!(
        refused,
        [
            "hydraulics.crane_k_factors",
            "thermal.incremental_zone_analysis",
            "thermal.shell_and_tube_design",
            "thermal.zone_rating",
        ]
    );
    for form in forms().iter().filter(|form| !form.available) {
        assert!(
            form.refusal.is_some_and(|reason| !reason.is_empty()),
            "{} is refused without saying why",
            form.id
        );
    }
    assert_eq!(available().len(), 94);
}

#[test]
fn every_available_calculation_has_scalar_inputs_and_outputs() {
    for form in available() {
        for input in form.inputs {
            assert!(
                matches!(input.kind, "quantity" | "enum" | "boolean" | "string"),
                "{} has a {} input the field cannot draw",
                form.id,
                input.kind
            );
            if input.kind == "quantity" {
                assert!(
                    input.unit.is_some(),
                    "{}'s {} has no unit",
                    form.id,
                    input.name
                );
            }
        }
        for output in form.outputs {
            assert!(
                matches!(output.kind, "quantity" | "enum" | "boolean" | "string"),
                "{} has a {} output the result table cannot draw",
                form.id,
                output.kind
            );
        }
        assert!(!form.inputs.is_empty(), "{} has no inputs", form.id);
        assert!(!form.outputs.is_empty(), "{} has no outputs", form.id);
        assert!(
            !form.worked_example.is_empty(),
            "{} has no worked example for the sheet to load",
            form.id
        );
    }
}

#[test]
fn every_calculation_says_what_it_is_for_and_every_field_says_what_it_means() {
    for form in forms() {
        assert!(!form.summary.is_empty(), "{} has no summary", form.id);
        assert!(
            !form.summary.to_lowercase().contains("neqsim"),
            "{}'s summary names the upstream implementation",
            form.id
        );
        for input in form.inputs {
            assert!(
                !input.description.trim().is_empty(),
                "{}'s input {} has no description",
                form.id,
                input.name
            );
        }
        for output in form.outputs {
            assert!(
                !output.description.trim().is_empty(),
                "{}'s output {} has no description",
                form.id,
                output.name
            );
        }
        for assumption in form.assumptions {
            assert!(
                !assumption.to_lowercase().contains("neqsim"),
                "{}'s assumption names the upstream implementation",
                form.id
            );
        }
    }
}

#[test]
fn every_calculation_has_its_provenance() {
    for form in forms() {
        assert!(
            ALL_PROVENANCE
                .iter()
                .any(|provenance| provenance.calc_id == form.id),
            "{} has no provenance entry",
            form.id
        );
    }
}

#[test]
fn the_ui_projection_never_names_the_upstream_implementation() {
    let text = to_json().expect("the catalogue writes");
    assert!(
        !text.to_lowercase().contains("neqsim"),
        "the catalogue carries the upstream implementation's name"
    );
}

#[test]
fn provenance_and_the_worked_example_are_joined_into_the_wire_form() {
    let value: serde_json::Value = serde_json::from_str(&to_json().expect("writes")).expect("json");
    let pump = value
        .as_array()
        .expect("an array")
        .iter()
        .find(|entry| entry["id"] == "hydraulics.pump_power")
        .expect("pump power is registered");
    assert_eq!(pump["worked_example"]["rho"], 998.0);
    assert!(pump["provenance"]["verification"].is_string());
    assert!(pump["provenance"]["spec_sha256"].is_string());
    assert!(
        pump.get("source").is_none(),
        "no upstream source is on the wire"
    );
}
