//! `chemeng` - command-line interface for the chemeng calculation library.
//!
//! ```text
//! chemeng pipe --fluid water --flow 10 --diameter 0.1 --length 100 \
//!     --fittings "90_elbow,gate_valve_open"
//! ```
//!
//! The CLI exists to make the library runnable without writing code, and to
//! demonstrate the composition of several calcs in one boundary. It reports the
//! same warnings the library does, printed rather than buried, because a
//! command-line user is exactly the person most likely to take a number at face
//! value.

use chemeng_cli::{pipe, report};
use clap::{Parser, Subcommand};

/// Open, validated chemical engineering calculations.
#[derive(Debug, Parser)]
#[command(name = "chemeng", version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Pressure drop through a pipe with fittings.
    Pipe(PipeArgs),

    /// List the fitting ids the registry knows about.
    Fittings,
}

#[derive(Debug, clap::Args)]
struct PipeArgs {
    /// Fluid to use. Built in: water, air.
    #[arg(long, default_value = "water")]
    fluid: String,

    /// Fluid temperature, in degrees Celsius.
    #[arg(long, default_value_t = 20.0)]
    temperature: f64,

    /// Volumetric or mass flow rate, in the unit given by --flow-unit.
    #[arg(long)]
    flow: f64,

    /// Unit of --flow.
    #[arg(long, default_value = "m3/h")]
    flow_unit: String,

    /// Pipe internal diameter, in metres.
    #[arg(long)]
    diameter: f64,

    /// Pipe length, in metres.
    #[arg(long)]
    length: f64,

    /// Absolute pipe roughness, in metres. Defaults to 0.046 mm (commercial
    /// steel); the report says so when this default is used.
    #[arg(long)]
    roughness: Option<f64>,

    /// Comma-separated fitting ids, e.g. "90_elbow,gate_valve_open". Run
    /// `chemeng fittings` to see the registry.
    #[arg(long, default_value = "")]
    fittings: String,

    /// Which friction factor correlation to use.
    #[arg(long, default_value = "colebrook")]
    friction_method: String,
}

fn main() -> std::process::ExitCode {
    let cli = Cli::parse();
    match cli.command {
        Command::Pipe(args) => run_pipe(args),
        Command::Fittings => match report::list_fittings() {
            Ok(()) => std::process::ExitCode::SUCCESS,
            Err(e) => fail(&e),
        },
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
    let fittings: Vec<String> = args
        .fittings
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_string)
        .collect();

    let result = match pipe::compute(
        &args.fluid,
        args.temperature,
        args.flow,
        flow_unit,
        args.diameter,
        args.length,
        roughness,
        &fittings,
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
fn fail(error: &chemeng_core::ChemEngError) -> std::process::ExitCode {
    eprintln!("chemeng: {error}");
    if let Some(field) = error.field() {
        eprintln!("  offending input: {field}");
    }
    std::process::ExitCode::from(2)
}
