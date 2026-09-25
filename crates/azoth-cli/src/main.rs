//! `azoth` - command-line interface for the azoth calculation library.
//!
//! ```text
//! azoth pipe --fluid water --flow 10 --diameter 0.1 --length 100 \
//!     --fittings "90_elbow,gate_valve_open"
//! ```
//!
//! The CLI exists to make the library runnable without writing code, and to
//! demonstrate the composition of several calcs in one boundary. It reports the
//! same warnings the library does, printed rather than buried, because a
//! command-line user is exactly the person most likely to take a number at face
//! value.
//!
//! The arguments are declared in [`azoth_cli::cli`] and the output is rendered in
//! [`azoth_cli::report`], both of which are testable. What is left here is the
//! wiring, which is the part an integration test against the built binary covers.

use azoth_cli::cli::{CheckArgs, Cli, Command, PipeArgs, RunArgs};
use azoth_cli::{check, pipe, report, run};
use clap::Parser;

fn main() -> std::process::ExitCode {
    let cli = Cli::parse();
    match cli.command {
        Command::Pipe(args) => run_pipe(args),
        Command::Fittings => match report::list_fittings() {
            Ok(()) => std::process::ExitCode::SUCCESS,
            Err(e) => fail(&e),
        },
        Command::Check(args) => run_check(args),
        Command::Run(args) => run_flowsheet(args),
    }
}

fn run_check(args: CheckArgs) -> std::process::ExitCode {
    match check::check(&args.flowsheet, &args.palette) {
        Ok(report) => {
            println!("{}", report.text);
            if report.defects == 0 {
                std::process::ExitCode::SUCCESS
            } else {
                std::process::ExitCode::from(2)
            }
        }
        Err(message) => {
            eprintln!("azoth: {message}");
            std::process::ExitCode::from(2)
        }
    }
}

fn run_flowsheet(args: RunArgs) -> std::process::ExitCode {
    match run::run(&args.flowsheet, &args.palette, args.json) {
        Ok(report) => {
            println!("{}", report.text);
            if report.ran {
                std::process::ExitCode::SUCCESS
            } else {
                std::process::ExitCode::from(2)
            }
        }
        Err(message) => {
            eprintln!("azoth: {message}");
            std::process::ExitCode::from(2)
        }
    }
}

fn run_pipe(args: PipeArgs) -> std::process::ExitCode {
    let flow_unit = match pipe::FlowUnit::parse(&args.flow_unit) {
        Ok(u) => u,
        Err(e) => return fail(&e),
    };
    let method = match pipe::FrictionMethod::parse(&args.friction_method) {
        Ok(m) => m,
        Err(e) => return fail(&e),
    };

    let roughness_supplied = args.roughness.is_some();
    let roughness = args.roughness.unwrap_or(pipe::DEFAULT_ROUGHNESS_M);

    let result = match pipe::compute(
        &args.fluid,
        args.temperature,
        args.flow,
        flow_unit,
        args.diameter,
        args.length,
        roughness,
        &args.fitting_ids(),
        method,
    ) {
        Ok(r) => r,
        Err(e) => return fail(&e),
    };

    report::print_pipe(&result, flow_unit, args.flow, roughness_supplied);
    std::process::ExitCode::SUCCESS
}

/// Print an error the way the library means it: what was wrong, and which input
/// was wrong. A bare "failed" would leave the user guessing at all of it.
fn fail(error: &azoth_core::AzothError) -> std::process::ExitCode {
    eprintln!("azoth: {error}");
    if let Some(field) = error.field() {
        eprintln!("  offending input: {field}");
    }
    std::process::ExitCode::from(2)
}
