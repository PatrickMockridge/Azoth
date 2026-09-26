//! `azoth forms` and `azoth edit`, end to end.
//!
//! The two commands that are the middleware's doors on a command line: the palette as the forms a
//! front-end renders, and one command against a document. `tests/run.rs` covers the third, which
//! answers what a document computes; these cover what it *says* and what an edit makes of it.

use std::path::{Path, PathBuf};
use std::process::Command;

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn demo() -> String {
    root()
        .join("specs/flowsheets/demo.toml")
        .display()
        .to_string()
}

fn palette() -> String {
    root().join("specs/unit_ops").display().to_string()
}

/// The binary, run the way a user runs it.
fn azoth(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_azoth"))
        .args(args)
        .output()
        .expect("the built binary should run")
}

fn json(output: &std::process::Output) -> serde_json::Value {
    let text = String::from_utf8(output.stdout.clone()).expect("utf-8");
    serde_json::from_str(&text).unwrap_or_else(|error| panic!("{error}\n{text}"))
}

#[test]
fn the_palette_comes_back_as_a_form_per_entry() {
    let output = azoth(&["forms", "--palette", &palette()]);
    assert!(output.status.success(), "it should print");
    let document = json(&output);
    let entries = document["unit_ops"].as_array().expect("an array");
    assert_eq!(entries.len(), 29);
    assert_eq!(
        entries
            .iter()
            .filter(|entry| !entry["model"].is_null())
            .count(),
        28
    );
    assert_eq!(
        entries
            .iter()
            .filter(|entry| entry["runnable"] == true)
            .count(),
        28
    );

    // The pump, field by field: the unit comes from the vocabulary and the kind from the model's
    // own declaration, which is the pair a form needs to choose a widget.
    let pump = entries
        .iter()
        .find(|entry| entry["id"] == "unit_ops.pump")
        .expect("the pump is in the palette");
    assert_eq!(pump["name"], "Pump");
    assert_eq!(pump["model"], "process.pump");
    assert_eq!(pump["runnable"], true);
    assert_eq!(pump["ports"][0]["direction"], "in");
    assert_eq!(pump["ports"][0]["multiplicity"], "one");
    let pressure = pump["parameters"]
        .as_array()
        .expect("an array")
        .iter()
        .find(|parameter| parameter["name"] == "outlet_pressure")
        .expect("declared");
    assert_eq!(pressure["kind"], "quantity");
    assert_eq!(pressure["required"], true);
    assert_eq!(pressure["unit"], "Pa");
    assert_eq!(pressure["dimension"], "pressure");
    assert!(
        !pressure["description"]
            .as_str()
            .expect("a sentence")
            .is_empty()
    );

    // An enum carries its values, and a range its bounds - the model's own, with the rationale a
    // field shows on a violation.
    let column = entries
        .iter()
        .find(|entry| entry["id"] == "unit_ops.distillation_column")
        .expect("declared");
    let solver = column["parameters"]
        .as_array()
        .expect("an array")
        .iter()
        .find(|parameter| parameter["name"] == "solver_type")
        .expect("declared");
    assert_eq!(solver["kind"], "enum");
    assert!(
        solver["values"]
            .as_array()
            .expect("an array")
            .iter()
            .any(|value| value == "naphtali_sandholm"),
        "{solver}"
    );
    let efficiency = pump["parameters"]
        .as_array()
        .expect("an array")
        .iter()
        .find(|parameter| parameter["name"] == "isentropic_efficiency")
        .expect("declared");
    assert_eq!(efficiency["range"][0]["min"], 0.0);
    assert_eq!(efficiency["range"][0]["min_inclusive"], false);
    assert_eq!(efficiency["range"][0]["max"], 1.0);
    assert_eq!(efficiency["range"][0]["severity"], "error");
    assert!(
        !efficiency["range"][0]["rationale"]
            .as_str()
            .expect("a sentence")
            .is_empty()
    );
}

#[test]
fn the_tool_schema_is_the_command_model() {
    let output = azoth(&["forms", "--palette", &palette(), "--tools"]);
    assert!(output.status.success(), "it should print");
    let document = json(&output);
    let tools = document["tools"].as_array().expect("an array");
    assert_eq!(tools.len(), 15);

    let names: Vec<&str> = tools
        .iter()
        .map(|tool| tool["name"].as_str().expect("a name"))
        .collect();
    assert_eq!(
        names,
        [
            "add_instance",
            "remove_instance",
            "connect",
            "disconnect",
            "set_parameter",
            "unset_parameter",
            "add_input",
            "remove_input",
            "add_product",
            "remove_product",
            "add_recycle",
            "remove_recycle",
            "set_recycle",
            "unset_recycle_field",
            "set_position",
        ]
    );

    // The `unit` an `add_instance` takes is the palette's own list, which is the one place a
    // schema reaches past the command to the declaration.
    let add = &tools[0]["input_schema"];
    assert_eq!(add["additionalProperties"], false);
    assert_eq!(add["properties"]["command"]["const"], "add_instance");
    assert_eq!(
        add["properties"]["unit"]["enum"]
            .as_array()
            .expect("an array")
            .len(),
        29
    );
    assert_eq!(
        add["required"],
        serde_json::json!(["command", "id", "unit"])
    );
}

#[test]
fn an_edit_leaves_the_document_it_made_and_says_what_is_wrong_with_it() {
    // An instance the rest of the flowsheet was wired through: removing it leaves two ports unfed,
    // and the check that follows the edit is what reports them.
    let output = azoth(&[
        "edit",
        "--flowsheet",
        &demo(),
        "--palette",
        &palette(),
        "--command",
        r#"{"command":"remove_instance","id":"hx1"}"#,
    ]);
    assert_eq!(output.status.code(), Some(2), "it is refused");
    let text = String::from_utf8(output.stdout).expect("utf-8");
    assert!(text.contains("flowsheets.demo: refused"), "{text}");
    assert!(text.contains("p1.ports.outlet"), "{text}");
    assert!(text.contains("sep1.ports.feed"), "{text}");

    // And the same flowsheet with an edit that leaves it sound.
    let output = azoth(&[
        "edit",
        "--flowsheet",
        &demo(),
        "--palette",
        &palette(),
        "--command",
        r#"{"command":"set_parameter","instance":"p1","name":"outlet_pressure","value":4000000.0}"#,
    ]);
    assert_eq!(output.status.code(), Some(0), "it is sound");
    let text = String::from_utf8(output.stdout).expect("utf-8");
    assert!(text.contains("flowsheets.demo: OK"), "{text}");
}

/// **The envelope is one document**, and the CLI is a binding of it like any other: `--json`
/// prints exactly what a browser gets, with the edit applied and the run's values in it.
#[test]
fn the_envelope_carries_the_edited_document_and_the_run() {
    let output = azoth(&[
        "edit",
        "--flowsheet",
        &demo(),
        "--palette",
        &palette(),
        "--command",
        r#"{"command":"set_position","node":"instance:sep1","x":1234.0,"y":56.0}"#,
        "--run",
        "--json",
    ]);
    assert_eq!(output.status.code(), Some(0), "it runs");
    let document = json(&output);
    assert_eq!(document["ok"], true);
    assert_eq!(document["dirty"], false);
    assert_eq!(document["diagnostics"], serde_json::json!([]));

    // The edit is in the document *and* in the graph, because both are projected from one value.
    let text = document["flowsheet"]["document"]
        .as_str()
        .expect("the text");
    assert!(text.contains("[layout.instances]"), "{text}");
    let sep1 = document["flowsheet"]["graph"]["nodes"]
        .as_array()
        .expect("nodes")
        .iter()
        .find(|node| node["id"] == "instance:sep1")
        .expect("drawn");
    assert_eq!(
        sep1["position"],
        serde_json::json!({"x": 1234.0, "y": 56.0})
    );

    // And the values, through the codec the whole library shares.
    let outlet = &document["session"]["streams"]["p1.outlet"];
    assert_eq!(outlet["P"]["unit"], "Pa");
    assert_eq!(
        outlet["P"]["magnitude_si"], 2000000.0,
        "the document's own value"
    );
    assert_eq!(document["session"]["converged"], true);
    assert_eq!(document["session"]["tears"][0]["stream"], "recycle_1");

    // Without `--run` there is a document and no values, which is what `dirty` says.
    let output = azoth(&[
        "edit",
        "--flowsheet",
        &demo(),
        "--palette",
        &palette(),
        "--command",
        r#"{"command":"add_product","name":"purge"}"#,
        "--json",
    ]);
    let document = json(&output);
    assert_eq!(document["session"], serde_json::Value::Null);
    assert_eq!(document["dirty"], true);
    assert_eq!(document["run_error"], serde_json::Value::Null);
}

#[test]
fn a_command_that_is_not_one_is_refused() {
    let output = azoth(&[
        "edit",
        "--flowsheet",
        &demo(),
        "--palette",
        &palette(),
        "--command",
        r#"{"command":"delete_everything"}"#,
    ]);
    assert_eq!(output.status.code(), Some(2));
    let text = String::from_utf8(output.stderr).expect("utf-8");
    assert!(text.contains("the command is not one"), "{text}");
}
