//! `azoth run`, end to end.
//!
//! `tests/report.rs` covers the commands around it: the arguments a user types and the report
//! they read. This file covers the one that *computes* - the shipped flowsheet through the built
//! binary, in both renderings, and a document the checker refuses.
//!
//! **The two renderings are checked against each other rather than against a constant.** The
//! report and `--json` are two readings of one session, and a report that read a different field
//! than the codec would be the failure a constant would not catch. What is *also* pinned is the
//! NeqSim capture's `p1_outlet_T` band, which is the only external measurement of this run.

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

#[test]
fn the_shipped_flowsheet_runs_from_the_command_line() {
    let output = azoth(&["run", "--flowsheet", &demo(), "--palette", &palette()]);
    assert!(output.status.success(), "it should run clean");
    let text = String::from_utf8(output.stdout).expect("utf-8");

    // Every stream the run produced, by the endpoint that produced it.
    for endpoint in [
        "feed_1",
        "mix1.product",
        "p1.outlet",
        "hx1.outlet",
        "sep1.vapour",
        "sep1.liquid",
        "vapour_product",
    ] {
        assert!(text.contains(endpoint), "`{endpoint}` is missing:\n{text}");
    }
    // And the tear, with its convergence rather than only its destination.
    assert!(
        text.contains("recycle_1: iterations=2 solved=true"),
        "{text}"
    );
    assert!(text.contains("converged=true iterations=2"), "{text}");

    // The pump's outlet against `captures/process_flowsheet.tsv`, whose `p1_outlet_T` is
    // `416.011567832174`. The measured difference is `0.012 K` and is upstream of the execution
    // layer - `executor_calls_kernels.rs` holds the executor to `process.pump` at `5.7e-12`
    // relative - so the band is the capture's own and not the command's.
    let line = text
        .lines()
        .find(|line| line.trim_start().starts_with("p1.outlet:"))
        .expect("the pump's outlet is reported");
    let temperature: f64 = line
        .split("T=")
        .nth(1)
        .and_then(|rest| rest.split_whitespace().next())
        .expect("the line carries a temperature")
        .parse()
        .expect("it is a number");
    assert!(
        (temperature - 416.011567832174).abs() < 0.05,
        "the pump's outlet is at {temperature} K"
    );
}

/// **The codec's document, which is what a script reads.** `--json` is not a second rendering:
/// it is `executor::json`, the same document `azoth.process.run_flowsheet(...).document` hands a
/// Python caller, so a script that reads one reads the other.
#[test]
fn the_json_rendering_is_the_codec_s_document() {
    let output = azoth(&[
        "run",
        "--flowsheet",
        &demo(),
        "--palette",
        &palette(),
        "--json",
    ]);
    assert!(output.status.success(), "it should run clean");
    let text = String::from_utf8(output.stdout).expect("utf-8");
    let document: serde_json::Value =
        serde_json::from_str(&text).expect("`--json` prints valid JSON");

    assert_eq!(document["flowsheet"], "flowsheets.demo");
    assert_eq!(document["converged"], true);
    assert_eq!(document["iterations"], 2);
    assert_eq!(document["tears"][0]["stream"], "recycle_1");
    assert_eq!(document["tears"][0]["solved"], true);

    // The record's own field names, and a magnitude carrying the unit it is in - the shape the
    // Python bridge transports, which is why the document is the rendering a front-end wants.
    let outlet = &document["streams"]["p1.outlet"];
    assert_eq!(outlet["n"]["unit"], "mol/s");
    assert_eq!(outlet["P"]["unit"], "Pa");
    assert_eq!(outlet["T"]["unit"], "K");
    assert_eq!(outlet["h"]["unit"], "J/mol");
    assert_eq!(outlet["z"].as_array().map(Vec::len), Some(2));

    // **The two renderings agree about the same run.** This is the assertion a constant could not
    // make: the report reads the session's streams and the document reads the codec's, and a
    // report that had drifted onto another field would show up here and nowhere else.
    let report = azoth(&["run", "--flowsheet", &demo(), "--palette", &palette()]);
    let rendered = String::from_utf8(report.stdout).expect("utf-8");
    let line = rendered
        .lines()
        .find(|line| line.trim_start().starts_with("p1.outlet:"))
        .expect("the pump's outlet is reported");
    let reported_p: f64 = line
        .split("P=")
        .nth(1)
        .and_then(|rest| rest.split_whitespace().next())
        .expect("the line carries a pressure")
        .parse()
        .expect("it is a number");
    assert_eq!(
        Some(reported_p),
        outlet["P"]["magnitude_si"].as_f64(),
        "the report and the document disagree about the pump's outlet pressure"
    );
}

/// **A flowsheet that does not validate is not one to run**, and the checker's verdict is the
/// report rather than a stack trace or a wrong answer. The document here is the shipped one with
/// a connection reading a `hx1.vent` the entry does not declare, which is `UnknownPort` and an
/// `Error` — a `Warning` would not stop the run, and that is the line the two severities draw.
#[test]
fn a_document_the_checker_refuses_is_not_run() {
    let text = std::fs::read_to_string(root().join("specs/flowsheets/demo.toml"))
        .expect("the shipped flowsheet is there")
        .replace("from = \"hx1.outlet\"", "from = \"hx1.vent\"");
    let path = std::env::temp_dir().join("azoth-cli-run-refused.toml");
    std::fs::write(&path, text).expect("the fixture is writable");

    let output = azoth(&[
        "run",
        "--flowsheet",
        &path.display().to_string(),
        "--palette",
        &palette(),
    ]);
    let _ = std::fs::remove_file(&path);

    assert_eq!(output.status.code(), Some(2), "a refused document exits 2");
    let stderr = String::from_utf8(output.stderr).expect("utf-8");
    let stdout = String::from_utf8(output.stdout).expect("utf-8");
    assert!(
        stderr.contains("UnknownPort") || stdout.contains("UnknownPort"),
        "the checker's own diagnostic is reported:\n{stderr}{stdout}"
    );
    assert!(
        !stdout.contains("converged="),
        "a refused document was run anyway:\n{stdout}"
    );
}
