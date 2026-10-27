//! `process.pipe` - the line's kernel as a registered id.
//!
//! Spec: `specs/models/process/pipe.toml`. The arithmetic is [`crate::kernels::pipe`]; what
//! is here is the boundary a case and a cross-impl test address.

use azoth_core::units::{
    Length, MolarEnergy, Pressure, ThermodynamicTemperature, joules_per_mole, meters, pascals,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{scalar, warnings as wire_warnings};
use crate::kernels::pipe as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// One row of `fittings`: an L/D ratio and how many of it.
#[derive(Debug, Clone, PartialEq)]
pub struct PipeFitting {
    /// The fitting's equivalent length over the pipe diameter.
    pub ld_ratio: f64,
    /// How many identical fittings.
    pub count: f64,
}

/// Result of `process.pipe`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PipeResult {
    /// Molar flow out, mol/s: the inlet's.
    pub outlet_n: f64,
    /// Outlet composition, the inlet's.
    pub outlet_z: Vec<f64>,
    /// Outlet pressure, which the line *solves* rather than takes.
    #[serde(serialize_with = "scalar")]
    pub outlet_p: Pressure,
    /// Outlet temperature, which is the inlet's.
    #[serde(serialize_with = "scalar")]
    pub outlet_t: ThermodynamicTemperature,
    /// Outlet molar enthalpy, which is not the inlet's - the drop is isothermal.
    #[serde(serialize_with = "scalar")]
    pub outlet_h: MolarEnergy,
    /// The length the pressure drop was computed over: the physical length plus the
    /// fittings' equivalent length.
    #[serde(serialize_with = "scalar")]
    pub effective_length: Length,
    /// `inlet_p - outlet_p`, which is `AdiabaticPipe.getPressureDrop()`.
    #[serde(serialize_with = "scalar")]
    pub pressure_drop: Pressure,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

impl PipeResult {
    /// The result of one kernel call.
    ///
    /// **The warnings are the caller's**: a flowsheet's are the *checker's* rather than
    /// [`apply_checks`]'.
    #[must_use]
    pub fn of(
        outlet: &Stream,
        pressure_drop: Pressure,
        effective_length: Length,
        warnings: Vec<Warning>,
    ) -> Self {
        Self {
            outlet_n: outlet.n,
            outlet_z: outlet.z.clone(),
            outlet_p: outlet.p,
            outlet_t: outlet.t,
            outlet_h: joules_per_mole(outlet.h.value),
            effective_length,
            pressure_drop,
            warnings,
        }
    }
}

/// Drop a stream's pressure along a line.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the geometry is not positive, and whatever the
/// flash, the density, the viscosity or the phase label refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are eight
pub fn pipe(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
    length: Length,
    diameter: Length,
    roughness: Length,
    fittings: Option<&[PipeFitting]>,
    outlet_pressure: Option<Pressure>,
    correlation: Option<&str>,
    number_of_increments: Option<usize>,
    elevation: Option<Length>,
) -> Result<PipeResult> {
    let spec = &model_gen::PIPE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "inlet_t" => Some(inlet_t.value),
            "length" => Some(length.value),
            "diameter" => Some(diameter.value),
            "roughness" => Some(roughness.value),
            "outlet_pressure" => outlet_pressure.map(|value| value.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let total_ld = fittings
        .unwrap_or(&[])
        .iter()
        .map(|fitting| {
            if !fitting.ld_ratio.is_finite() || fitting.ld_ratio < 0.0 {
                return Err(AzothError::invalid_input(
                    "fittings",
                    format!(
                        "a fitting's L/D is {}, and an equivalent length is non-negative",
                        fitting.ld_ratio
                    ),
                ));
            }
            if !fitting.count.is_finite() || fitting.count < 0.0 {
                return Err(AzothError::invalid_input(
                    "fittings",
                    format!(
                        "a fitting's count is {}, and a count is non-negative",
                        fitting.count
                    ),
                ));
            }
            Ok(fitting.ld_ratio * fitting.count)
        })
        .sum::<Result<f64>>()?;
    let effective_length = meters(length.value + total_ld * diameter.value);

    let feed = Stream::from_pt(
        components.to_vec(),
        inlet_z.to_vec(),
        inlet_n,
        inlet_p,
        inlet_t,
    )?;

    let correlation = correlation.unwrap_or("adiabatic");
    let result = match correlation {
        "adiabatic" => {
            if let Some(target) = outlet_pressure {
                let solved =
                    kernel::calc_flow(&feed, target, effective_length, diameter, roughness)?;
                let outlet = Stream::from_pt(
                    components.to_vec(),
                    inlet_z.to_vec(),
                    solved,
                    target,
                    inlet_t,
                )?;
                PipeResult::of(
                    &outlet,
                    pascals(inlet_p.value - target.value),
                    effective_length,
                    warnings,
                )
            } else {
                let outcome = kernel(&feed, effective_length, diameter, roughness)?;
                let pressure_drop = pascals(inlet_p.value - outcome.outlet.p.value);
                PipeResult::of(&outcome.outlet, pressure_drop, effective_length, warnings)
            }
        }
        "beggs_and_brills" => {
            if outlet_pressure.is_some() {
                return Err(AzothError::invalid_input(
                    "outlet_pressure",
                    "the flow solve is `AdiabaticPipe.calcFlow`'s, and the Beggs-and-Brills \
                     correlation has no such inverse here",
                ));
            }
            let march = azoth_eos::pipe_beggs_and_brills::pipe_beggs_and_brills(
                components,
                inlet_n,
                inlet_z,
                inlet_p,
                inlet_t,
                effective_length,
                diameter,
                elevation.unwrap_or_else(|| meters(0.0)),
                number_of_increments.unwrap_or(5),
                roughness,
                Some("pr"),
            )?;
            let outlet = Stream::from_pt(
                components.to_vec(),
                inlet_z.to_vec(),
                inlet_n,
                march.outlet_pressure,
                inlet_t,
            )?;
            PipeResult::of(
                &outlet,
                pascals(inlet_p.value - march.outlet_pressure.value),
                effective_length,
                warnings,
            )
        }
        other => {
            return Err(AzothError::invalid_input(
                "correlation",
                format!("{other} is not one of adiabatic or beggs_and_brills"),
            ));
        }
    };

    Ok(result)
}
