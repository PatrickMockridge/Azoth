//! The command line itself: what `azoth` accepts.
//!
//! Here rather than in `main.rs` because a binary's items are not reachable from an
//! integration test, and what a user types is the part of a CLI that breaks first. A
//! mis-declared flag is not thermodynamics and no calculation test can see it.

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
