//! The command line itself: what `azoth` accepts.
//!
//! Here rather than in `main.rs` because a binary's items are not reachable from an
//! integration test, and what a user types is the part of a CLI that breaks first.
//!
//! There is no `--keycard`: neither subcommand resolves a substance, so a card passed to
//! one would change no answer. The sections that could feed them, `fluids` and
//! `fittings`, are `compiled`-stage and nothing reads them at run time.

use clap::{Parser, Subcommand};

/// Open, validated chemical engineering calculations.
#[derive(Debug, Parser)]
#[command(name = "azoth", version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    /// Pressure drop through a pipe with fittings.
    Pipe(PipeArgs),

    /// List the fitting ids the registry knows about.
    Fittings,

    /// Validate a flowsheet against the unit-operation palette.
    Check(CheckArgs),

    /// Run a flowsheet to its steady state.
    Run(RunArgs),

    /// Print the unit-operation palette as a form per entry.
    Forms(FormsArgs),

    /// Apply one command to a flowsheet and print what it made.
    Edit(EditArgs),

    /// Serve the flowsheet as MCP tools over stdio, for an agent.
    Mcp(McpArgs),
}

#[derive(Debug, clap::Args)]
pub struct McpArgs {
    /// The flowsheet to serve. It is read once and edited in place for the life of the process;
    /// the file itself is never written.
    #[arg(long)]
    pub flowsheet: std::path::PathBuf,

    /// The palette directory. Defaults to the shipped `specs/unit_ops`.
    #[arg(long, default_value = "specs/unit_ops")]
    pub palette: std::path::PathBuf,

    /// Do not run the document after each call.
    ///
    /// The default runs it, because the envelope is what every call answers with and the values
    /// are the half of it a caller cannot otherwise ask for — there is no `run` tool, since a tool
    /// that ran the session would be a command the command model does not have. The cost is that
    /// every call pays the physics; this is the flag for a document where that is too much.
    #[arg(long)]
    pub no_run: bool,
}

#[derive(Debug, clap::Args)]
pub struct FormsArgs {
    /// The palette directory. Defaults to the shipped `specs/unit_ops`.
    #[arg(long, default_value = "specs/unit_ops")]
    pub palette: std::path::PathBuf,

    /// Include the agent's tool schema, projected from the command model.
    #[arg(long)]
    pub tools: bool,
}

#[derive(Debug, clap::Args)]
pub struct EditArgs {
    /// The flowsheet to edit.
    #[arg(long)]
    pub flowsheet: std::path::PathBuf,

    /// The command, as the JSON object the middleware's command model reads.
    #[arg(long)]
    pub command: String,

    /// The palette directory. Defaults to the shipped `specs/unit_ops`.
    #[arg(long, default_value = "specs/unit_ops")]
    pub palette: std::path::PathBuf,

    /// Run the edited document as well, and report what it reached.
    #[arg(long)]
    pub run: bool,

    /// Print the whole envelope as JSON rather than as a report.
    #[arg(long)]
    pub json: bool,
}

#[derive(Debug, clap::Args)]
pub struct RunArgs {
    /// The flowsheet to run. It declares its own inputs.
    #[arg(long)]
    pub flowsheet: std::path::PathBuf,

    /// The palette directory. Defaults to the shipped `specs/unit_ops`.
    #[arg(long, default_value = "specs/unit_ops")]
    pub palette: std::path::PathBuf,

    /// Print the result as JSON rather than as a report.
    #[arg(long)]
    pub json: bool,
}

#[derive(Debug, clap::Args)]
pub struct CheckArgs {
    /// The flowsheet to validate.
    #[arg(long)]
    pub flowsheet: std::path::PathBuf,

    /// The palette directory. Defaults to the shipped `specs/unit_ops`.
    #[arg(long, default_value = "specs/unit_ops")]
    pub palette: std::path::PathBuf,
}

#[derive(Debug, clap::Args)]
pub struct PipeArgs {
    /// Fluid to use. Built in: water, air.
    #[arg(long, default_value = "water")]
    pub fluid: String,

    /// Fluid temperature, in degrees Celsius.
    #[arg(long, default_value_t = 20.0)]
    pub temperature: f64,

    /// Volumetric or mass flow rate, in the unit given by --flow-unit.
    #[arg(long)]
    pub flow: f64,

    /// Unit of --flow.
    #[arg(long, default_value = "m3/h")]
    pub flow_unit: String,

    /// Pipe internal diameter, in metres.
    #[arg(long)]
    pub diameter: f64,

    /// Pipe length, in metres.
    #[arg(long)]
    pub length: f64,

    /// Absolute pipe roughness, in metres. Defaults to 0.046 mm (commercial
    /// steel); the report says so when this default is used.
    #[arg(long)]
    pub roughness: Option<f64>,

    /// Comma-separated fitting ids, e.g. "90_elbow,gate_valve_open". Run
    /// `azoth fittings` to see the registry.
    #[arg(long, default_value = "")]
    pub fittings: String,

    /// Which friction factor correlation to use.
    #[arg(long, default_value = "colebrook")]
    pub friction_method: String,
}

impl PipeArgs {
    /// The fitting ids as a list, from the comma-separated argument.
    ///
    /// Whitespace around an id is the caller's, not part of the id: `"a, b"` and
    /// `"a,b"` name the same two fittings. Empty entries are dropped so a trailing
    /// comma is not an unknown fitting.
    #[must_use]
    pub fn fitting_ids(&self) -> Vec<String> {
        self.fittings
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(str::to_string)
            .collect()
    }
}
