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

/// Result of `process.mixer`.
///
/// Several feeds blended into one. The outlet pressure is the lowest inlet pressure -
/// a mixer cannot deliver a stream at a pressure nobody supplied - and the temperature
/// is solved for rather than mixed, because the energy balance is at constant enthalpy
/// and the enthalpy of a mixture is not the average of its parts' temperatures.
#[derive(Debug, Clone)]
pub struct MixerResult {
    /// The outlet temperature, from the isenthalpic flash.
    pub temperature: ThermodynamicTemperature,
    /// The outlet pressure: the lowest of the inlet pressures.
    pub pressure: Pressure,
    /// The outlet molar flow, in mol/s. The sum of the inlets'.
    pub flow: f64,
    /// The outlet mole fractions - the flow-weighted blend of the inlets'.
    pub z_out: Vec<f64>,
    /// The vapour fraction at the outlet, or `None` for a single-phase feed.
    pub beta: Option<f64>,
    /// Which phase the blend is in.
    pub phase: Phase,
    /// Flash iterations taken.
    pub iterations: u32,
    /// Caveats, deduplicated.
    pub warnings: Vec<Warning>,
}

impl CalcResult for MixerResult {
    const CALC_ID: &'static str = "process.mixer";
    const FIELDS: &'static [&'static str] = &[
        "T",
        "P",
        "flow",
        "z_out",
        "beta",
        "phase",
        "iterations",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `process.throttling_valve`.
///
/// A stream's pressure dropped at constant enthalpy, so the temperature is an answer
/// rather than an input - and for a real gas it is not the inlet temperature, which is
/// the Joule-Thomson effect and the whole reason a valve is worth modelling.
///
/// **No flow and no composition are reported.** Neither is computed: a valve changes
/// nothing about how much is flowing or what it is, and reporting an input back as a
/// result invites a caller to treat it as one.
#[derive(Debug, Clone)]
pub struct ThrottlingValveResult {
    /// The outlet temperature, from the isenthalpic flash.
    pub temperature: ThermodynamicTemperature,
    /// The outlet pressure: the inlet pressure less the drop.
    pub pressure: Pressure,
    /// Which phase the stream is in at the outlet.
    pub phase: Phase,
    /// The vapour fraction at the outlet, or `None` for a single-phase stream.
    pub beta: Option<f64>,
    /// Flash iterations taken.
    pub iterations: u32,
    /// Caveats, deduplicated.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ThrottlingValveResult {
    const CALC_ID: &'static str = "process.throttling_valve";
    const FIELDS: &'static [&'static str] = &["T", "P", "phase", "beta", "iterations", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `process.heater`.
///
/// A duty applied at a fixed pressure, so the temperature is an answer. A cooler is the
/// same unit with a negative duty.
///
/// The shape is deliberately identical to [`ThrottlingValveResult`]'s, and the struct is
/// deliberately **not** shared with it: one result type per calculation is the rule
/// everywhere in this workspace (`crates/azoth-core/src/result.rs` records why), and a
/// shared struct would leave `result_fields` with no id to answer for one of the two.
#[derive(Debug, Clone)]
pub struct HeaterResult {
    /// The outlet temperature, from the isenthalpic flash.
    pub temperature: ThermodynamicTemperature,
    /// The outlet pressure: the inlet pressure less the drop.
    pub pressure: Pressure,
    /// Which phase the stream is in at the outlet.
    pub phase: Phase,
    /// The vapour fraction at the outlet, or `None` for a single-phase stream.
    pub beta: Option<f64>,
    /// Flash iterations taken.
    pub iterations: u32,
    /// Caveats, deduplicated.
    pub warnings: Vec<Warning>,
}

impl CalcResult for HeaterResult {
    const CALC_ID: &'static str = "process.heater";
    const FIELDS: &'static [&'static str] = &["T", "P", "phase", "beta", "iterations", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `process.splitter`.
///
/// One feed divided into branches at the same temperature, pressure and composition.
/// The only thing that differs between the branches is how much of the feed each
/// carries, which is why `flows` is a vector and `beta` is a single number.
#[derive(Debug, Clone)]
pub struct SplitterResult {
    /// The temperature every branch leaves at - the feed's.
    pub temperature: ThermodynamicTemperature,
    /// The pressure every branch leaves at - the feed's.
    pub pressure: Pressure,
    /// The phase the feed is in, and so the phase of every branch.
    pub phase: Phase,
    /// The vapour fraction, or `None` for a single-phase feed.
    pub beta: Option<f64>,
    /// Each branch's molar flow, in mol/s, in the order the fractions were given.
    pub flows: Vec<f64>,
    /// Flash iterations taken.
    pub iterations: u32,
    /// Caveats, deduplicated.
    pub warnings: Vec<Warning>,
}

impl CalcResult for SplitterResult {
    const CALC_ID: &'static str = "process.splitter";
    const FIELDS: &'static [&'static str] =
        &["T", "P", "phase", "beta", "flows", "iterations", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `process.compressor`.
///
/// A pressure rise at a stated isentropic efficiency. See `isentropic.rs` for the
/// procedure - a compressor, a pump and an expander run the same one and differ only in
/// which way the efficiency scales the ideal enthalpy change.
#[derive(Debug, Clone)]
pub struct CompressorResult {
    /// The outlet temperature, from the flash at the actual outlet enthalpy.
    pub temperature: ThermodynamicTemperature,
    /// The outlet pressure - the specification, not an answer.
    pub pressure: Pressure,
    /// The shaft power in watts, **signed**: positive when the machine puts energy into
    /// the fluid, negative when it takes it out. A compressor's is positive.
    pub power: f64,
    /// The vapour fraction at the outlet, or `None` for a single-phase stream.
    pub beta: Option<f64>,
    /// Which phase the stream is in at the outlet.
    pub phase: Phase,
    /// The temperature the fluid would have reached at **unit** efficiency.
    ///
    /// Reported because it is what an efficiency is measured against: it is the
    /// isentropic answer, it is where the loss is counted from, and a caller sizing an
    /// intercooler needs it and cannot recover it from the rest of the result.
    pub isentropic_temperature: ThermodynamicTemperature,
    /// Flash iterations taken, summed over this model's flashes.
    pub iterations: u32,
    /// Caveats, deduplicated.
    pub warnings: Vec<Warning>,
}

impl CalcResult for CompressorResult {
    const CALC_ID: &'static str = "process.compressor";
    const FIELDS: &'static [&'static str] = &[
        "T",
        "P",
        "power",
        "beta",
        "phase",
        "isentropic_temperature",
        "iterations",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `process.pump`.
///
/// The same procedure as [`CompressorResult`]'s and a separate struct for the same
/// reason [`HeaterResult`] gives. A pump is a compressor for a liquid, and NeqSim's
/// default path for one is literally the isentropic flash pair; the difference that
/// matters is that a liquid's temperature rise is small and its volume is nearly
/// constant, which is a property of the state and not of the model.
#[derive(Debug, Clone)]
pub struct PumpResult {
    /// The outlet temperature, from the flash at the actual outlet enthalpy.
    pub temperature: ThermodynamicTemperature,
    /// The outlet pressure - the specification, not an answer.
    pub pressure: Pressure,
    /// The shaft power in watts, signed: positive for work into the fluid.
    pub power: f64,
    /// The vapour fraction at the outlet, or `None` for a single-phase stream.
    pub beta: Option<f64>,
    /// Which phase the stream is in at the outlet.
    pub phase: Phase,
    /// The temperature the fluid would have reached at unit efficiency.
    pub isentropic_temperature: ThermodynamicTemperature,
    /// Flash iterations taken, summed over this model's flashes.
    pub iterations: u32,
    /// Caveats, deduplicated.
    pub warnings: Vec<Warning>,
}

impl CalcResult for PumpResult {
    const CALC_ID: &'static str = "process.pump";
    const FIELDS: &'static [&'static str] = &[
        "T",
        "P",
        "power",
        "beta",
        "phase",
        "isentropic_temperature",
        "iterations",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `process.expander`.
///
/// The same procedure again, with the efficiency multiplying rather than dividing -
/// which is what makes `power` negative. NeqSim's `Expander` inherits from its
/// `Compressor` and makes exactly that change, at `Expander.java:653`.
#[derive(Debug, Clone)]
pub struct ExpanderResult {
    /// The outlet temperature, from the flash at the actual outlet enthalpy.
    pub temperature: ThermodynamicTemperature,
    /// The outlet pressure - the specification, not an answer.
    pub pressure: Pressure,
    /// The shaft power in watts, signed: **negative**, because the fluid is doing the
    /// work. The sign is the only thing distinguishing this from a compressor's.
    pub power: f64,
    /// The vapour fraction at the outlet, or `None` for a single-phase stream.
    pub beta: Option<f64>,
    /// Which phase the stream is in at the outlet.
    pub phase: Phase,
    /// The temperature the fluid would have reached at unit efficiency - the *lowest*
    /// it could reach for this pressure drop.
    pub isentropic_temperature: ThermodynamicTemperature,
    /// Flash iterations taken, summed over this model's flashes.
    pub iterations: u32,
    /// Caveats, deduplicated.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ExpanderResult {
    const CALC_ID: &'static str = "process.expander";
    const FIELDS: &'static [&'static str] = &[
        "T",
        "P",
        "power",
        "beta",
        "phase",
        "isentropic_temperature",
        "iterations",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}
