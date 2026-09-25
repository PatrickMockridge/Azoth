//! The live document, and the one object a front-end reads it through.
//!
//! `docs/src/architecture/middleware.md`'s session/document layer, and the envelope that carries
//! it. What is measured here is the division the layer rests on: a check is cheap and happens on
//! every edit, a run is the physics and happens when asked, and between them the values are
//! marked stale rather than shown beside a document they no longer describe.

use std::path::{Path, PathBuf};

use azoth_process::executor::{Session, json};
use azoth_process::middleware::command::Command;
use azoth_process::middleware::envelope;
use azoth_process::middleware::session::Workspace;
use azoth_process::{ExecutionOrder, UnitOpSpec, load_palette};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn palette() -> Vec<UnitOpSpec> {
    load_palette(&root().join("specs/unit_ops")).expect("the shipped palette loads")
}

fn demo() -> String {
    std::fs::read_to_string(root().join("specs/flowsheets/demo.toml"))
        .expect("the shipped flowsheet is there")
}

fn opened() -> Workspace {
    Workspace::open(&demo(), palette()).expect("it reads")
}

#[test]
fn an_edit_rechecks_and_the_values_go_stale() {
    let mut workspace = opened();
    assert!(workspace.ok(), "the shipped document is sound");
    assert!(
        workspace.dirty(),
        "a document that has never run has no current values"
    );
    assert!(workspace.paths().is_empty(), "and no paths to show");
    assert!(workspace.value("p1.outlet.P").is_err());

    workspace.run().expect("the shipped document runs");
    assert!(!workspace.dirty(), "the values are the document's now");
    assert!(
        workspace.paths().contains(&"p1.outlet.P".to_string()),
        "{:?}",
        workspace.paths()
    );
    assert!(workspace.value("p1.outlet.P").expect("a value") > 0.0);

    // An edit re-checks, and marks the values stale rather than leaving them to be misread as
    // the document's.
    workspace
        .apply(&Command::SetParameter {
            instance: "p1".into(),
            name: "outlet_pressure".into(),
            value: serde_json::json!(4.0e6),
        })
        .expect("the edit lands");
    assert!(workspace.dirty(), "the values are older than the document");
    assert!(workspace.ok(), "and the document still runs");

    // A broken edit is reported by the check that follows it, in the same call.
    workspace
        .apply(&Command::RemoveInstance { id: "p1".into() })
        .expect("the edit lands");
    assert!(!workspace.ok());
    assert!(
        workspace
            .diagnostics()
            .iter()
            .any(|diagnostic| diagnostic.code() == "underfed_port"),
        "{:?}",
        workspace.diagnostics()
    );
}

/// **A failure the checker cannot see stays a failure.** Whether a named substance exists is the
/// databank's answer and not a document fact, so it arrives as a run error rather than as a
/// diagnostic - and the run that failed takes the previous values with it, because they were the
/// document as it was before the edit.
#[test]
fn a_run_that_fails_takes_the_values_with_it() {
    let mut workspace = opened();
    workspace.run().expect("it runs");
    assert!(workspace.value("p1.outlet.T").is_ok());

    // A *well-shaped* record naming a substance no databank carries: one component with one mole
    // fraction, so no shape rule fires and only the fluid's resolution can refuse it.
    let broken = demo()
        .replace(
            "components = [\"methane\", \"n-butane\"]",
            "components = [\"unobtainium\"]",
        )
        .replace("z = [0.9, 0.1]", "z = [1.0]");
    assert_ne!(broken, demo(), "the fixture changed something");
    let mut workspace = Workspace::open(&broken, palette()).expect("it reads");
    assert!(
        workspace.ok(),
        "the checker has no objection to a substance it cannot see: {:?}",
        workspace.diagnostics()
    );

    let error = workspace.run().expect_err("the run refuses it");
    assert!(error.to_string().contains("unobtainium"), "{error}");
    assert!(workspace.run_error().is_some());
    assert!(workspace.report().is_none(), "no report from a failed run");
    assert!(workspace.dirty());
    assert!(workspace.paths().is_empty());

    // And a successful run clears it.
    workspace = opened();
    workspace.run().expect("it runs");
    assert!(workspace.run_error().is_none());
}

/// **The envelope carries the codec's own report**, not a second rendering of it: the value
/// embedded in the envelope and the value `to_json` writes are the same JSON.
#[test]
fn the_envelope_carries_the_session_report_verbatim() {
    let mut workspace = opened();
    workspace.run().expect("it runs");
    let document = envelope::to_json(&workspace).expect("the envelope writes");
    let read: serde_json::Value = serde_json::from_str(&document).expect("it parses");

    // The same document, run directly, through the published codec.
    let flowsheet = workspace.flowsheet().clone();
    let session = Session::run(
        flowsheet,
        &palette(),
        &std::collections::BTreeMap::new(),
        ExecutionOrder::Insertion,
    )
    .expect("it runs");
    let published: serde_json::Value =
        serde_json::from_str(&json::to_json(&session).expect("it writes")).expect("it parses");

    assert_eq!(read["session"], published);
    assert_eq!(read["flowsheet"]["id"], "flowsheets.demo");
    assert_eq!(read["ok"], true);
    assert_eq!(read["dirty"], false);
    // Every path the session offers is in the envelope, in the session's own spelling.
    let paths: Vec<&str> = read["paths"]
        .as_array()
        .expect("an array")
        .iter()
        .filter_map(serde_json::Value::as_str)
        .collect();
    assert_eq!(paths, session.paths());
    assert!(paths.contains(&"recycle_1"), "{paths:?}");
}

/// The projection and the pass-throughs are the ones below, not a second copy.
#[test]
fn an_envelope_is_one_document_of_the_layers_below_it() {
    let mut workspace = opened();
    workspace.run().expect("it runs");
    let read: serde_json::Value =
        serde_json::from_str(&envelope::to_json(&workspace).expect("it writes")).expect("parses");

    assert_eq!(
        read["flowsheet"]["graph"]["nodes"].as_array().map(Vec::len),
        Some(6)
    );
    assert_eq!(
        read["flowsheet"]["graph"]["edges"].as_array().map(Vec::len),
        Some(6)
    );
    assert_eq!(read["diagnostics"], serde_json::json!([]));
    assert_eq!(read["run_error"], serde_json::Value::Null);
    // The order the next run takes is the class's default, and it is in the envelope so a control
    // reads it from the session rather than holding a copy of it.
    assert_eq!(read["execution_order"], "insertion");
    workspace.set_order(ExecutionOrder::Topological);
    let reordered: serde_json::Value =
        serde_json::from_str(&envelope::to_json(&workspace).expect("it writes")).expect("parses");
    assert_eq!(reordered["execution_order"], "topological");
    assert_eq!(
        reordered["dirty"], false,
        "an order is not a change to the document"
    );
    // The document is the text, which is what a save writes.
    let document = read["flowsheet"]["document"].as_str().expect("a string");
    assert!(
        document.starts_with("id = \"flowsheets.demo\""),
        "{document}"
    );
    assert_eq!(
        azoth_process::Flowsheet::from_toml(document).expect("it reads back"),
        *workspace.flowsheet()
    );
}
