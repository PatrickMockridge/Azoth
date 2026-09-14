//! What `azoth` prints, and what it accepts.
//!
//! `tests/pipe.rs` covers the composition - velocity, friction factor, fitting loss.
//! This file covers the two things around it that a user meets first: the arguments
//! they type, and the report they read.
//!
//! The report has two rules of its own, both stated in `report.rs`, and both are
//! assertions here rather than intentions: assumptions the CLI supplied are printed
//! rather than implied, and warnings are printed in full rather than summarised. A
//! count is the failure mode the second rule exists against - it lets a reader take
//! the pressure drop and skip the caveats.

use std::process::Command;

use azoth_cli::cli::{Cli, Command as Subcommand};
use azoth_cli::pipe::{self, FlowUnit, FrictionMethod};
use azoth_cli::report;
use azoth_core::{Warning, WarningCode};
use clap::Parser;

fn headline_result() -> pipe::PipeResult {
    pipe::compute(
        "water",
        20.0,
        10.0,
        FlowUnit::CubicMetresPerHour,
        0.1,
        100.0,
        pipe::DEFAULT_ROUGHNESS_M,
        &["90_elbow".to_string(), "gate_valve_open".to_string()],
        FrictionMethod::Colebrook,
    )
    .expect("the headline example should compute")
}

fn headline_report(roughness_supplied: bool) -> String {
    report::render_pipe(
        &headline_result(),
        FlowUnit::CubicMetresPerHour,
        10.0,
        roughness_supplied,
    )
}

#[test]
fn the_report_says_which_values_the_cli_chose() {
    let out = headline_report(false);

    assert!(
        out.contains("default - override with --roughness"),
        "a roughness the user did not supply is printed as a default they can \
         override, not left to be inferred:\n{out}"
    );
    assert!(
        out.contains("derived from flow and bore"),
        "the velocity is derived from two inputs the user did supply, and says so:\n{out}"
    );
}

#[test]
fn a_supplied_roughness_is_not_called_a_default() {
    let out = headline_report(true);

    assert!(
        !out.contains("default - override with --roughness"),
        "the roughness was supplied, so nothing should call it a default:\n{out}"
    );
}

#[test]
fn the_report_carries_the_three_sections_and_the_total() {
    let result = headline_result();
    let out = headline_report(true);

    for section in ["inputs", "flow", "pressure drop"] {
        assert!(out.contains(section), "{section} is missing from:\n{out}");
    }
    assert!(
        out.contains(&format!("{:.4}", result.dp_total)),
        "the total pressure drop is not in the report:\n{out}"
    );
    assert!(out.contains("reynolds number"), "{out}");
    assert!(out.contains("90_elbow, gate_valve_open"), "{out}");
}

#[test]
fn a_clean_result_says_so() {
    let out = report::render_warnings(&[]);

    assert!(out.contains("no warnings"), "{out}");
    assert!(
        !out.contains("warning(s)"),
        "a clean result must not print a warning count:\n{out}"
    );
}

#[test]
fn warnings_are_printed_in_full_rather_than_counted() {
    let warnings = [
        Warning::for_field(
            WarningCode::OutOfValidRange,
            "velocity",
            "2.5 m/s is above the range this correlation was fitted over",
        ),
        Warning::new(
            WarningCode::TransitionalFlow,
            "the Reynolds number is between 2300 and 4000, where neither the laminar \
             nor the turbulent correlation applies",
        ),
    ];

    let out = report::render_warnings(&warnings);

    assert!(out.contains("2 warning(s)"), "{out}");
    assert!(
        out.contains("2.5 m/s is above the range this correlation was fitted over"),
        "the message itself must be printed; a count would let a reader take the \
         number and skip this:\n{out}"
    );
    assert!(out.contains("about: velocity"), "{out}");
    assert!(
        out.contains("the Reynolds number is between 2300 and 4000"),
        "a warning not attached to a field is printed too:\n{out}"
    );
    assert!(
        out.contains("it has not been\n  checked the way an in-range one has"),
        "the consequence of a range warning is stated, not left to the code:\n{out}"
    );
}

#[test]
fn the_fitting_list_marks_the_placeholders() {
    let out = report::render_fittings().expect("the embedded registry should parse");

    assert!(out.contains("90_elbow"), "{out}");
    assert!(
        out.contains("ESTIMATED DUMMY") || out.contains("placeholder"),
        "the fitting table is placeholder data and the listing must say so:\n{out}"
    );
}

#[test]
fn the_arguments_parse_and_the_defaults_are_the_documented_ones() {
    let cli = Cli::try_parse_from([
        "azoth",
        "pipe",
        "--flow",
        "10",
        "--diameter",
        "0.1",
        "--length",
        "100",
    ])
    .expect("the headline command should parse");

    let Subcommand::Pipe(args) = cli.command else {
        panic!("`pipe` parsed as another subcommand");
    };
    assert_eq!(args.fluid, "water");
    assert_eq!(args.temperature, 20.0);
    assert_eq!(args.flow_unit, "m3/h");
    assert_eq!(args.friction_method, "colebrook");
    assert!(args.roughness.is_none());
}

#[test]
fn a_missing_required_argument_is_refused() {
    let parsed = Cli::try_parse_from(["azoth", "pipe", "--flow", "10", "--length", "100"]);

    assert!(
        parsed.is_err(),
        "`--diameter` is required and a command without it must not parse"
    );
}

#[test]
fn fitting_ids_tolerate_the_spacing_a_caller_writes() {
    let cli = Cli::try_parse_from([
        "azoth",
        "pipe",
        "--flow",
        "10",
        "--diameter",
        "0.1",
        "--length",
        "100",
        "--fittings",
        " 90_elbow , gate_valve_open ,",
    ])
    .expect("the command should parse");
    let Subcommand::Pipe(args) = cli.command else {
        panic!("`pipe` parsed as another subcommand");
    };

    assert_eq!(
        args.fitting_ids(),
        vec!["90_elbow".to_string(), "gate_valve_open".to_string()],
        "a space around an id is the caller's, and a trailing comma names no fitting"
    );
}

/// The binary, run the way a user runs it.
fn azoth(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_azoth"))
        .args(args)
        .output()
        .expect("the built binary should run")
}

#[test]
fn the_headline_command_runs_end_to_end() {
    let output = azoth(&[
        "pipe",
        "--fluid",
        "water",
        "--flow",
        "10",
        "--diameter",
        "0.1",
        "--length",
        "100",
        "--fittings",
        "90_elbow,gate_valve_open",
    ]);

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        output.status.success(),
        "exit {:?}: {stdout}",
        output.status.code()
    );
    assert!(stdout.contains("azoth pipe"), "{stdout}");
    assert!(stdout.contains("no warnings"), "{stdout}");
}

#[test]
fn an_unknown_fitting_exits_two_and_names_the_input() {
    let output = azoth(&[
        "pipe",
        "--flow",
        "10",
        "--diameter",
        "0.1",
        "--length",
        "100",
        "--fittings",
        "nonsense",
    ]);

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(2), "stderr: {stderr}");
    assert!(
        stderr.contains("offending input"),
        "the error must name the input that was wrong, not just fail:\n{stderr}"
    );
}

#[test]
fn a_bad_argument_is_a_usage_error_with_a_readable_message() {
    let output = azoth(&[
        "pipe",
        "--flow",
        "not-a-number",
        "--diameter",
        "0.1",
        "--length",
        "100",
    ]);

    assert!(!output.status.success());
    assert_eq!(
        output.status.code(),
        Some(2),
        "clap exits 2 on a usage error, the same code as a refused input"
    );
}

#[test]
fn the_fittings_subcommand_runs() {
    let output = azoth(&["fittings"]);

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(output.status.success(), "exit {:?}", output.status.code());
    assert!(stdout.contains("known fitting ids"), "{stdout}");
    assert!(stdout.contains("90_elbow"), "{stdout}");
}
