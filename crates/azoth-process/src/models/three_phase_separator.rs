//! `process.three_phase_separator` - the separator's kernel as a registered id.
//!
//! Spec: `specs/models/process/three_phase_separator.toml`. The arithmetic is
//! [`crate::kernels::three_phase_separator`]; what is here is the boundary a case and a
//! cross-impl test address.
//!
//! **Three outlets, so fifteen fields.** A single-multiplicity port crosses as the record's
//! five fields per port, and this one has three of them, so a result is the triple written
//! out under the ports' own names rather than a list.

use azoth_core::units::{
    Length, MolarEnergy, Power, Pressure, ThermodynamicTemperature, Time, joules_per_mole, meters,
    seconds,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use azoth_eos::phase_transport::{PhaseKind, phase_transport};
use azoth_eos::pt_flash;

use crate::kernels::three_phase_separator::{
    Entrainment, retention_times, three_phase_separator as kernel,
};
use crate::model_gen;
use crate::segment::phase::{mass_density, molar_mass_of};
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
    /// The stored water level, clamped at zero; present when a level was stated.
    pub water_level: Option<Length>,
    /// The stored oil level, clamped at zero; present as above.
    pub oil_level: Option<Length>,
    /// `getOilThickness`, present as above.
    pub oil_thickness: Option<Length>,
    /// `calcOilRetentionTime`, in seconds; absent with no vessel or no oil.
    pub oil_retention_time: Option<Time>,
    /// `calcWaterRetentionTime`, absent with no vessel or no aqueous phase.
    pub water_retention_time: Option<Time>,
    /// `calcInterfaceSettlingTime`: a 150 µm water droplet through the oil layer, in seconds.
    pub interface_settling_time: Option<Time>,
    /// Whether that time reaches the class's one-minute limit.
    pub settling_time_sufficient: Option<bool>,
    /// Caveats.
    pub warnings: Vec<Warning>,
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
    internal_diameter: Option<Length>,
    separator_length: Option<Length>,
    nll_fraction: Option<f64>,
    nil_fraction: Option<f64>,
    water_level: Option<Length>,
    oil_level: Option<Length>,
    gas_outlet_flow_fraction: Option<f64>,
    oil_outlet_flow_fraction: Option<f64>,
    water_outlet_flow_fraction: Option<f64>,
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
            "internal_diameter" => internal_diameter.map(|value| value.value),
            "separator_length" => separator_length.map(|value| value.value),
            "nll_fraction" => nll_fraction,
            "nil_fraction" => nil_fraction,
            _ => None,
        },
        &mut warnings,
    )?;

    // **The vessel block is both or neither.** A retention time needs a diameter and a
    // length, and the two level fractions are fractions of a vessel nobody stated.
    if internal_diameter.is_some() != separator_length.is_some() {
        let (named, missing) = if internal_diameter.is_some() {
            ("internal_diameter", "separator_length")
        } else {
            ("separator_length", "internal_diameter")
        };
        return Err(AzothError::invalid_input(
            missing,
            format!(
                "`{named}` was stated and `{missing}` was not, and a retention time needs both"
            ),
        ));
    }
    if internal_diameter.is_none() {
        if nll_fraction.is_some() {
            return Err(AzothError::invalid_input(
                "nll_fraction",
                "`nll_fraction` was stated without the vessel it is a fraction of",
            ));
        }
        if nil_fraction.is_some() {
            return Err(AzothError::invalid_input(
                "nil_fraction",
                "`nil_fraction` was stated without the vessel it is a fraction of",
            ));
        }
    }
    // **A level is stored state, and the setter clamps at zero.** A non-finite level is a
    // value no state carries, and `Math.max(0, NaN)` would store it.
    let water = match water_level {
        Some(value) => Some(clamp_level("water_level", value.value)?),
        None => None,
    };
    let oil = match oil_level {
        Some(value) => Some(clamp_level("oil_level", value.value)?),
        None => None,
    };

    // A zero flow is not a state, and the heat-per-mole divides by it.
    if feed_n <= 0.0 {
        return Err(azoth_core::AzothError::out_of_range(
            "feed_n",
            feed_n,
            "a unit operation needs a positive molar flow",
        ));
    }

    let feed = Stream::from_pt(components.to_vec(), feed_z.to_vec(), feed_n, feed_p, feed_t)?;
    let entrainment = Entrainment {
        gas_in_aqueous,
        gas_in_oil,
        oil_in_aqueous,
        oil_in_gas,
        aqueous_in_gas,
        aqueous_in_oil,
    };
    let (vapour, oil_out, aqueous_out) = kernel(&feed, pressure_drop, heat_input, entrainment)?;

    // The three valve fractions are stored state: the class's steady split never reads them.
    // They are declared because the class exposes and clamps them; `getMassBalance` is a
    // zero-by-construction residual the capture records, not a field.
    let _ = (
        gas_outlet_flow_fraction.map(|value| value.clamp(0.0, 1.0)),
        oil_outlet_flow_fraction.map(|value| value.clamp(0.0, 1.0)),
        water_outlet_flow_fraction.map(|value| value.clamp(0.0, 1.0)),
    );

    let retention = match (internal_diameter, separator_length) {
        (Some(diameter), Some(length)) => Some(retention_times(
            &oil_out,
            &aqueous_out,
            diameter,
            length,
            nll_fraction.unwrap_or(0.5),
            nil_fraction.unwrap_or(0.2),
        )?),
        _ => None,
    };
    let settling = match (internal_diameter, separator_length) {
        (Some(diameter), Some(length)) => interface_settling_minutes(
            &oil_out,
            &aqueous_out,
            diameter,
            length,
            nll_fraction.unwrap_or(0.5),
            nil_fraction.unwrap_or(0.2),
        )?,
        _ => None,
    };
    let has_levels = water.is_some() || oil.is_some();
    let water_value = water.unwrap_or(0.0);
    let oil_value = oil.unwrap_or(0.0);

    Ok(ThreePhaseSeparatorResult {
        vapour_n: vapour.n,
        vapour_z: vapour.z,
        vapour_p: vapour.p,
        vapour_t: vapour.t,
        vapour_h: joules_per_mole(vapour.h.value),
        light_liquid_n: oil_out.n,
        light_liquid_z: oil_out.z,
        light_liquid_p: oil_out.p,
        light_liquid_t: oil_out.t,
        light_liquid_h: joules_per_mole(oil_out.h.value),
        heavy_liquid_n: aqueous_out.n,
        heavy_liquid_z: aqueous_out.z,
        heavy_liquid_p: aqueous_out.p,
        heavy_liquid_t: aqueous_out.t,
        heavy_liquid_h: joules_per_mole(aqueous_out.h.value),
        water_level: has_levels.then(|| meters(water_value)),
        oil_level: has_levels.then(|| meters(oil_value)),
        oil_thickness: has_levels.then(|| meters((oil_value - water_value).max(0.0))),
        oil_retention_time: retention.and_then(|times| times.oil).map(seconds),
        water_retention_time: retention.and_then(|times| times.water).map(seconds),
        interface_settling_time: settling.map(|minutes| seconds(minutes * 60.0)),
        settling_time_sufficient: settling.map(|minutes| minutes >= 1.0),
        warnings,
    })
}

/// `ThreePhaseSeparator.calcInterfaceSettlingTime`, in **minutes** as the class answers it.
///
/// A **150 µm water droplet** settles through the oil layer by Stokes, on the oil phase's own
/// viscosity: `9.81 d^2 (rho_w - rho_o) / (18 mu)`, the layer over that velocity over 60.
/// The levels are the mechanical design's, with the class's `0.20`/`0.50` of the diameter
/// fallbacks. A non-positive velocity is the class's `Infinity`; `None` here, because a
/// non-finite time is not a value the record carries.
fn interface_settling_minutes(
    oil: &Stream,
    aqueous: &Stream,
    diameter: Length,
    _length: Length,
    nll_fraction: f64,
    nil_fraction: f64,
) -> Result<Option<f64>> {
    let (mixture, ideal_gas) = oil.mixture()?;
    let oil_flash = pt_flash(&mixture, oil.t, oil.p, &oil.z)?;
    let oil_density = mass_density(
        &mixture,
        oil.t,
        oil.p,
        &oil.z,
        oil_flash.z_liquid,
        molar_mass_of(&mixture, &oil.z)?,
    )?;
    let aqueous_flash = pt_flash(&mixture, aqueous.t, aqueous.p, &aqueous.z)?;
    let aqueous_density = mass_density(
        &mixture,
        aqueous.t,
        aqueous.p,
        &aqueous.z,
        aqueous_flash.z_liquid,
        molar_mass_of(&mixture, &aqueous.z)?,
    )?;
    let oil_viscosity =
        phase_transport(&mixture, &ideal_gas, PhaseKind::Oil, oil.t, oil.p, &oil.z)?
            .mu
            .value;
    if !(oil_density > 0.0 && aqueous_density > oil_density && oil_viscosity > 0.0) {
        return Ok(None);
    }
    let mut nil = nil_fraction * diameter.value;
    if nil <= 0.0 {
        nil = diameter.value * 0.20;
    }
    let mut nll = nll_fraction * diameter.value;
    if nll <= 0.0 {
        nll = diameter.value * 0.50;
    }
    let droplet = 150.0e-6_f64;
    let velocity =
        9.81 * droplet * droplet * (aqueous_density - oil_density) / (18.0 * oil_viscosity);
    if !velocity.is_finite() || velocity <= 0.0 {
        return Ok(None);
    }
    Ok(Some((nll - nil) / velocity.abs() / 60.0))
}

/// The class's setter clamp: a negative level is zero, and a non-finite one is refused.
fn clamp_level(field: &str, value: f64) -> Result<f64> {
    if !value.is_finite() {
        return Err(AzothError::invalid_input(
            field,
            format!("a level is a finite length, and {value} is not"),
        ));
    }
    Ok(value.max(0.0))
}

#[cfg(test)]
mod tests {
    use super::three_phase_separator;
    use azoth_core::units::{kelvins, meters, pascals};

    /// The settling capture's first vessel: the time is below the class's one-minute limit.
    /// **A boolean is not a case expectation** (`gen_registry` excludes flags), so this is
    /// the test that reads the flag back.
    #[test]
    fn the_settling_time_is_below_the_one_minute_limit() {
        let result = three_phase_separator(
            &[
                "methane".to_string(),
                "n-butane".to_string(),
                "water".to_string(),
            ],
            1.0,
            &[0.5, 0.3, 0.2],
            pascals(2_000_000.0),
            kelvins(300.0),
            pascals(0.0),
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            None,
            Some(meters(1.0)),
            Some(meters(3.0)),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .expect("the captured state settles");
        assert_eq!(result.settling_time_sufficient, Some(false));
        assert!(result.interface_settling_time.is_some());
    }
}
