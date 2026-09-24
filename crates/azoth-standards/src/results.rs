//! Result types for the standard calculations.
//!
//! The same contract as the other crates' results: one struct per calculation, field names
//! identical to the Python result dataclass and listed in [`CalcResult::FIELDS`], with a
//! test asserting the three agree. See `crates/azoth-core/src/result.rs` for why the
//! duplication is deliberate.

use azoth_core::{CalcResult, Warning};
use uom::si::f64::{MassDensity, MolarEnergy, MolarMass, Ratio};

/// Result of `standards.iso6976`.
#[derive(Debug, Clone, PartialEq)]
pub struct Iso6976Result {
    /// The mixture's molar mass.
    pub molar_mass: MolarMass,
    /// The mixture's compression factor at the volumetric reference temperature.
    pub compression_factor: Ratio,
    /// Density relative to dry air at the volumetric reference temperature.
    pub relative_density: Ratio,
    /// The ideal-gas density at the reference pressure and volumetric reference temperature.
    pub density_ideal: MassDensity,
    /// The real-gas density there, which is the ideal one over the compression factor.
    pub density_real: MassDensity,
    /// Superior (gross) molar calorific value at the energy reference temperature.
    pub superior_calorific_value: MolarEnergy,
    /// Inferior (net) molar calorific value at the energy reference temperature.
    pub inferior_calorific_value: MolarEnergy,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for Iso6976Result {
    const CALC_ID: &'static str = "standards.iso6976";
    const FIELDS: &'static [&'static str] = &[
        "molar_mass",
        "compression_factor",
        "relative_density",
        "density_ideal",
        "density_real",
        "superior_calorific_value",
        "inferior_calorific_value",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}
