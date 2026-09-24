//! `process.gas_scrubber` - `unit_ops.separator`'s kernel under the other entry.
//!
//! Spec: `specs/models/process/separator.toml`. The arithmetic is
//! [`crate::kernels::separator`]; what is here is the boundary a case and a cross-impl
//! test address.
//!
//! **Two outlets, so ten fields.** A single-multiplicity port crosses as the record's five
//! fields per port, and a separator has two of them, so a result is the pair written out
//! under the ports' own names rather than a list.

use azoth_core::units::{MolarEnergy, Power, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{CalcResult, Result, Warning, apply_checks};

use crate::kernels::gas_scrubber as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.separator`.
#[derive(Debug, Clone, PartialEq)]
pub struct GasScrubberResult {
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
    /// Liquid outlet molar flow, mol/s.
    pub liquid_n: f64,
    /// Liquid outlet composition.
    pub liquid_z: Vec<f64>,
    /// Liquid outlet pressure.
    pub liquid_p: Pressure,
    /// Liquid outlet temperature.
    pub liquid_t: ThermodynamicTemperature,
    /// Liquid outlet molar enthalpy.
    pub liquid_h: MolarEnergy,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for GasScrubberResult {
    const CALC_ID: &'static str = "process.gas_scrubber";
    const FIELDS: &'static [&'static str] = &[
        "vapour_n", "vapour_z", "vapour_p", "vapour_t", "vapour_h", "liquid_n", "liquid_z",
        "liquid_p", "liquid_t", "liquid_h", "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Flash a stream into vapour and liquid outlets.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the pressure drop takes the outlet below
/// zero, the entrainment fraction is outside `[0, 1]`, or the shapes disagree, and
/// whatever the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are eight
pub fn gas_scrubber(
    components: &[String],
    feed_n: f64,
    feed_z: &[f64],
    feed_p: Pressure,
    feed_t: ThermodynamicTemperature,
    pressure_drop: Pressure,
    gas_in_liquid: f64,
    heat_input: Option<Power>,
) -> Result<GasScrubberResult> {
    let spec = &model_gen::GAS_SCRUBBER_SPEC;
    let mut warnings = Vec::new();
    // One call over every declared input, for the reason `reactions.chemical_equilibrium`
    // gives. `heat_input` has no bound - a duty may heat or cool - so the `_` arm is the
    // whole of it, and an absent optional input cannot leave a check unrun because none
    // depends on it.
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "pressure_drop" => Some(pressure_drop.value),
            "gas_in_liquid" => Some(gas_in_liquid),
            _ => None,
        },
        &mut warnings,
    )?;

    let feed = Stream::from_pt(components.to_vec(), feed_z.to_vec(), feed_n, feed_p, feed_t)?;
    let (vapour, liquid) = kernel(&feed, pressure_drop, gas_in_liquid, heat_input)?;

    Ok(GasScrubberResult {
        vapour_n: vapour.n,
        vapour_z: vapour.z,
        vapour_p: vapour.p,
        vapour_t: vapour.t,
        vapour_h: joules_per_mole(vapour.h.value),
        liquid_n: liquid.n,
        liquid_z: liquid.z,
        liquid_p: liquid.p,
        liquid_t: liquid.t,
        liquid_h: joules_per_mole(liquid.h.value),
        warnings,
    })
}
