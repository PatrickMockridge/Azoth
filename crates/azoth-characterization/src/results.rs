//! Result types for the characterisation calculations.
//!
//! The same contract as the other namespaces' results: one struct per calculation, field names
//! identical to the Python result dataclass and listed in [`CalcResult::FIELDS`], with a test
//! asserting the three agree. See `crates/azoth-core/src/result.rs` for why the duplication is
//! deliberate.

use azoth_core::{CalcResult, Warning};
use uom::si::f64::{MassDensity, MolarMass, Pressure, ThermodynamicTemperature};

/// Result of `characterization.tbp_cut_properties`.
#[derive(Debug, Clone, PartialEq)]
pub struct TbpCutPropertiesResult {
    /// Critical temperature.
    pub tc: ThermodynamicTemperature,
    /// Critical pressure. Every model's correlation is written in bar; this is pascals.
    pub pc: Pressure,
    /// Normal boiling point.
    pub boiling_temperature: ThermodynamicTemperature,
    /// Pitzer's acentric factor, by whichever correlation the model uses.
    pub acentric_factor: f64,
    /// The cubic alpha function's `m`, **absent** for the five models that have no exponent to
    /// offer: `RiaziDaubert` sets `calcm = false` and the other four inherit a method that throws.
    pub attraction_exponent: Option<f64>,
    /// Watson's characterization factor, `(1.8*Tb)^(1/3)/d`.
    pub watson_k: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for TbpCutPropertiesResult {
    const CALC_ID: &'static str = "characterization.tbp_cut_properties";
    const FIELDS: &'static [&'static str] = &[
        "tc",
        "pc",
        "boiling_temperature",
        "acentric_factor",
        "attraction_exponent",
        "watson_k",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `characterization.tbp_closure`.
#[derive(Debug, Clone, PartialEq)]
pub struct TbpClosureResult {
    /// The cut's molar mass, the number a cubic needs.
    pub molar_mass: MolarMass,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for TbpClosureResult {
    const CALC_ID: &'static str = "characterization.tbp_closure";
    const FIELDS: &'static [&'static str] = &["molar_mass", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `characterization.tbp_density`.
#[derive(Debug, Clone, PartialEq)]
pub struct TbpDensityResult {
    /// The cut's normal liquid density at 15 C.
    pub density: MassDensity,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for TbpDensityResult {
    const CALC_ID: &'static str = "characterization.tbp_density";
    const FIELDS: &'static [&'static str] = &["density", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}
