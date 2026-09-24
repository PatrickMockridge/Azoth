//! `process.three_phase_separator` - the separator's kernel as a registered id.
//!
//! Spec: `specs/models/process/three_phase_separator.toml`. The arithmetic is
//! [`crate::kernels::three_phase_separator`]; what is here is the boundary a case and a
//! cross-impl test address.
//!
//! **Three outlets, so fifteen fields.** A single-multiplicity port crosses as the record's
//! five fields per port, and this one has three of them, so a result is the triple written
//! out under the ports' own names rather than a list.

use azoth_core::units::{MolarEnergy, Power, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{CalcResult, Result, Warning, apply_checks};

use crate::kernels::three_phase_separator::{Entrainment, three_phase_separator as kernel};
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.three_phase_separator`.
#[derive(Debug, Clone, PartialEq)]
pub struct ThreePhaseSeparatorResult {
    /// Vapour outlet molar flow, mol/s.
    pub vapour_n: f64,
    /// Vapour outlet composition.
    pub vapour_z: Vec<f64>,
    /// Vapour outlet pressure.
    pub vapour_p: Pressure,
    /// Vapour outlet temperature.
    pub vapour_t: ThermodynamicTemperature,
    /// Vapour outlet molar enthalpy.
    pub vapour_h: MolarEnergy,
    /// Oil outlet molar flow, mol/s.
    pub light_liquid_n: f64,
    /// Oil outlet composition.
    pub light_liquid_z: Vec<f64>,
    /// Oil outlet pressure.
    pub light_liquid_p: Pressure,
    /// Oil outlet temperature.
    pub light_liquid_t: ThermodynamicTemperature,
    /// Oil outlet molar enthalpy.
    pub light_liquid_h: MolarEnergy,
    /// Aqueous outlet molar flow, mol/s.
    pub heavy_liquid_n: f64,
    /// Aqueous outlet composition.
    pub heavy_liquid_z: Vec<f64>,
    /// Aqueous outlet pressure.
    pub heavy_liquid_p: Pressure,
    /// Aqueous outlet temperature.
    pub heavy_liquid_t: ThermodynamicTemperature,
    /// Aqueous outlet molar enthalpy.
    pub heavy_liquid_h: MolarEnergy,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ThreePhaseSeparatorResult {
    const CALC_ID: &'static str = "process.three_phase_separator";
    const FIELDS: &'static [&'static str] = &[
        "vapour_n",
        "vapour_z",
        "vapour_p",
        "vapour_t",
        "vapour_h",
        "light_liquid_n",
        "light_liquid_z",
        "light_liquid_p",
        "light_liquid_t",
        "light_liquid_h",
        "heavy_liquid_n",
        "heavy_liquid_z",
        "heavy_liquid_p",
        "heavy_liquid_t",
        "heavy_liquid_h",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Flash a stream into vapour, oil and aqueous outlets.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the pressure drop takes the outlet below
/// zero, an entrainment fraction is outside `[0, 1]`, or the shapes disagree, and whatever
/// the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are twelve
pub fn three_phase_separator(
    components: &[String],
    feed_n: f64,
    feed_z: &[f64],
    feed_p: Pressure,
    feed_t: ThermodynamicTemperature,
    pressure_drop: Pressure,
    gas_in_aqueous: f64,
    gas_in_oil: f64,
    oil_in_aqueous: f64,
    oil_in_gas: f64,
    aqueous_in_gas: f64,
    aqueous_in_oil: f64,
    heat_input: Option<Power>,
) -> Result<ThreePhaseSeparatorResult> {
    let spec = &model_gen::THREE_PHASE_SEPARATOR_SPEC;
    let mut warnings = Vec::new();
    // One call over every declared input, for the reason `reactions.chemical_equilibrium`
    // gives. `heat_input` has no bound - a duty may heat or cool - so the `_` arm is the
    // whole of it, and an absent optional input cannot leave a check unrun because none
    // depends on it.
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "feed_t" => Some(feed_t.value),
            "pressure_drop" => Some(pressure_drop.value),
            "gas_in_aqueous" => Some(gas_in_aqueous),
            "gas_in_oil" => Some(gas_in_oil),
            "oil_in_aqueous" => Some(oil_in_aqueous),
            "oil_in_gas" => Some(oil_in_gas),
            "aqueous_in_gas" => Some(aqueous_in_gas),
            "aqueous_in_oil" => Some(aqueous_in_oil),
            _ => None,
        },
        &mut warnings,
    )?;

    let feed = Stream::from_pt(components.to_vec(), feed_z.to_vec(), feed_n, feed_p, feed_t)?;
    let entrainment = Entrainment {
        gas_in_aqueous,
        gas_in_oil,
        oil_in_aqueous,
        oil_in_gas,
        aqueous_in_gas,
        aqueous_in_oil,
    };
    let (vapour, oil, aqueous) = kernel(&feed, pressure_drop, heat_input, entrainment)?;

    Ok(ThreePhaseSeparatorResult {
        vapour_n: vapour.n,
        vapour_z: vapour.z,
        vapour_p: vapour.p,
        vapour_t: vapour.t,
        vapour_h: joules_per_mole(vapour.h.value),
        light_liquid_n: oil.n,
        light_liquid_z: oil.z,
        light_liquid_p: oil.p,
        light_liquid_t: oil.t,
        light_liquid_h: joules_per_mole(oil.h.value),
        heavy_liquid_n: aqueous.n,
        heavy_liquid_z: aqueous.z,
        heavy_liquid_p: aqueous.p,
        heavy_liquid_t: aqueous.t,
        heavy_liquid_h: joules_per_mole(aqueous.h.value),
        warnings,
    })
}
