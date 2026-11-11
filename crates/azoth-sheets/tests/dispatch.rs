use azoth_core::unit_vocab_gen::si_factor;
use azoth_sheets::catalogue::forms;
use azoth_sheets::dispatch::{self, Outcome};
use azoth_sheets::form::{CalcForm, ExampleValue};
use serde_json::{Map, Value, json};

/// The worked example, in SI, as a sheet would hold it.
fn inputs(form: &CalcForm) -> Value {
    let declarations = form
        .inputs
        .iter()
        .map(|input| (input.name, input.unit))
        .collect::<std::collections::BTreeMap<_, _>>();
    let mut object = Map::new();
    for entry in form.worked_example {
        let value = match entry.value {
            ExampleValue::Number(number) => {
                let unit = declarations.get(entry.name).copied().flatten();
                let factor = unit.and_then(si_factor).unwrap_or(1.0);
                json!(number * factor)
            }
            ExampleValue::Flag(flag) => json!(flag),
            ExampleValue::Name(name) => json!(name),
        };
        object.insert(entry.name.to_string(), value);
    }
    Value::Object(object)
}

fn close(actual: f64, expected: f64, tolerance: f64) -> bool {
    let scale = expected.abs().max(1.0);
    (actual - expected).abs() <= tolerance * scale
}

#[test]
fn every_available_calculation_runs_its_worked_example() {
    let mut failures = Vec::new();
    for form in forms().iter().filter(|form| form.available) {
        let inputs = inputs(form);
        let outcome = dispatch::run(form.id, &inputs);
        let Outcome::Ok(answer) = outcome else {
            failures.push(format!("{}: {:?}", form.id, outcome));
            continue;
        };
        for expected in form.expected {
            let Some(actual) = answer.outputs.get(expected.name) else {
                failures.push(format!("{}: no output `{}`", form.id, expected.name));
                continue;
            };
            let declared = form
                .outputs
                .iter()
                .find(|output| output.name == expected.name)
                .and_then(|output| output.unit);
            let factor = declared.and_then(si_factor).unwrap_or(1.0);
            match expected.value {
                ExampleValue::Number(number) => {
                    let parsed = actual
                        .as_f64()
                        .or_else(|| actual.get("magnitude_si").and_then(Value::as_f64));
                    let Some(parsed) = parsed else {
                        failures.push(format!(
                            "{}: output `{}` is not a number: {actual}",
                            form.id, expected.name
                        ));
                        continue;
                    };
                    if !close(parsed, number * factor, form.tolerance) {
                        failures.push(format!(
                            "{}: `{}` = {parsed}, expected {} (tolerance {})",
                            form.id,
                            expected.name,
                            number * factor,
                            form.tolerance
                        ));
                    }
                }
                ExampleValue::Flag(flag) => {
                    if actual.as_bool() != Some(flag) {
                        failures.push(format!(
                            "{}: `{}` = {actual}, expected {flag}",
                            form.id, expected.name
                        ));
                    }
                }
                ExampleValue::Name(name) => {
                    if actual.as_str() != Some(name) {
                        failures.push(format!(
                            "{}: `{}` = {actual}, expected {name}",
                            form.id, expected.name
                        ));
                    }
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} worked examples failed:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn a_hard_limit_is_an_outcome_with_the_field_that_violated_it() {
    let answer = dispatch::run(
        "hydraulics.pump_power",
        &json!({"rho": 998.0, "q": 0.01, "H": 30.0, "eta": 0.0}),
    );
    let Outcome::Error { message, field } = answer else {
        panic!("eta = 0 is a hard limit");
    };
    assert_eq!(field.as_deref(), Some("eta"));
    assert!(message.contains("eta"), "{message}");
}

#[test]
fn a_missing_input_is_an_outcome_and_not_a_panic() {
    let answer = dispatch::run("hydraulics.pump_power", &json!({"rho": 998.0}));
    let Outcome::Error { message, field } = answer else {
        panic!("a missing input must be refused");
    };
    assert_eq!(field.as_deref(), Some("q"));
    assert!(message.contains("required"), "{message}");
}

#[test]
fn the_four_table_shaped_calculations_are_refused_by_name() {
    for id in [
        "hydraulics.crane_k_factors",
        "thermal.incremental_zone_analysis",
        "thermal.shell_and_tube_design",
        "thermal.zone_rating",
    ] {
        let answer = dispatch::run(id, &json!({}));
        assert!(
            matches!(answer, Outcome::Error { .. }),
            "{id} must not run through a sheet"
        );
    }
}
