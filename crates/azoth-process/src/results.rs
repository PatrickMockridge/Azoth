//! Result types for the unit operations.
//!
//! The same contract as every other namespace's results: one struct per unit
//! operation, field names identical to the Python result dataclass and listed in
//! [`CalcResult::FIELDS`], with a test asserting the three agree. See
//! `crates/azoth-core/src/result.rs` for why the duplication is deliberate.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{CalcResult, Warning};
use azoth_eos::Phase;

/// Result of `process.separator`.
///
/// A feed split into a gas and a liquid. Both outlets are at the **same temperature and
/// pressure** - the flash's answer - because that is what a separator is: a vessel at
/// one state, with two streams leaving it. Reporting one `T` and one `P` rather than a
/// pair each is that statement rather than a saving.
#[derive(Debug, Clone)]
pub struct SeparatorResult {
    /// The temperature both outlets leave at. This is the flash's answer.
    pub temperature: ThermodynamicTemperature,
    /// The pressure both outlets leave at: the inlet pressure less `pressure_drop`.
    pub pressure: Pressure,
    /// The vapour fraction at that state, or `None` for a single-phase feed.
    ///
    /// `None` rather than a number outside `[0, 1]`, and the split is decided on
    /// [`Self::phase`] rather than on this - see the module note in `separator.rs`.
    pub beta: Option<f64>,
    /// Molar flow to the gas outlet, in mol/s. Zero when the feed is all liquid.
    pub gas_flow: f64,
    /// The gas outlet's mole fractions.
    ///
    /// The feed's composition when `gas_flow` is zero: a stream with no flow has no
    /// phase composition, and carrying the feed's through keeps the vector a
    /// composition rather than a row of zeros that no downstream unit could consume.
    pub gas_z: Vec<f64>,
    /// Molar flow to the liquid outlet, in mol/s. Zero when the feed is all vapour.
    pub liquid_flow: f64,
    /// The liquid outlet's mole fractions, by the same convention as [`Self::gas_z`].
    pub liquid_z: Vec<f64>,
    /// Which phase the feed was in.
    pub phase: Phase,
    /// Iterations the flash took. Zero is not reported: a single-phase feed still runs
    /// the flash, and its count is the evidence the answer is a converged one.
    pub iterations: u32,
    /// Caveats, deduplicated.
    pub warnings: Vec<Warning>,
}

impl CalcResult for SeparatorResult {
    const CALC_ID: &'static str = "process.separator";
    const FIELDS: &'static [&'static str] = &[
        "T",
        "P",
        "beta",
        "gas_flow",
        "gas_z",
        "liquid_flow",
        "liquid_z",
        "phase",
        "iterations",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}
