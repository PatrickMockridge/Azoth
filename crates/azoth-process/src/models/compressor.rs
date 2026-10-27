//! `process.compressor` - the Compressor's kernel as a registered id.
//!
//! Spec: `specs/models/process/compressor.toml`. The arithmetic is
//! [`crate::kernels::compressor`]; what is here is the boundary a case and a cross-impl
//! test address.
//!
//! **A delegation, because the route is `Pump.run`'s.** `Pump` takes it on `calculateAsCompressor` defaulting to true, and `Compressor.run`'s no-chart branch is the same three lines - so the kernel is `kernels::pump`'s and this adds nothing to it. **A stated curve selects the chart branch**, which is the compressor's own.

use azoth_core::units::{
    Length, MolarEnergy, Power, Pressure, ThermodynamicTemperature, VolumeRate,
    cubic_meter_per_hour, joules_per_mole, meter, meters,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};

use crate::compressor_chart::{CompressorChart, CurvePoint};
use crate::kernels::compressor as kernel;
use crate::kernels::compressor::compressor_with_chart;
use crate::model_gen;
use crate::stream::Stream;
use crate::unported;

/// One row of `curve`: a measured point, in the units the spec declares.
#[derive(Debug, Clone, PartialEq)]
pub struct CompressorCurvePoint {
    /// The speed the point was measured at, rpm.
    pub speed: f64,
    /// The measured actual flow.
    pub flow: VolumeRate,
    /// The measured polytropic head, in metres.
    pub head: Length,
    /// The measured polytropic efficiency, as a fraction.
    pub efficiency: f64,
}

/// Result of `process.compressor`.
#[derive(Debug, Clone, PartialEq)]
pub struct CompressorResult {
    /// Molar flow out, mol/s: the inlet's.
    pub outlet_n: f64,
    /// Outlet composition, the inlet's.
    pub outlet_z: Vec<f64>,
    /// Outlet pressure: the stated one, or the chart's solved ratio.
    pub outlet_p: Pressure,
    /// Outlet temperature, solved from the shifted enthalpy.
    pub outlet_t: ThermodynamicTemperature,
    /// Outlet molar enthalpy: the inlet's plus the step.
    pub outlet_h: MolarEnergy,
    /// The chart head the branch used, m; zero in the isentropic branch.
    pub polytropic_head: Length,
    /// The chart polytropic efficiency the branch used; zero in the isentropic branch.
    pub polytropic_efficiency: f64,
    /// The shaft power, W.
    pub power: Power,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

/// Raise a stream's pressure along an isentrope or the machine's chart.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the shapes disagree, the efficiency is outside
/// `(0, 1]`, a pump without a curve has no pressure, or the chart type is not the default;
/// and whatever the databank or the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input
pub fn compressor(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
    outlet_pressure: Option<Pressure>,
    isentropic_efficiency: f64,
    curve: Option<&[CompressorCurvePoint]>,
    speed: Option<f64>,
    chart_type: Option<&str>,
) -> Result<CompressorResult> {
    let spec = &model_gen::COMPRESSOR_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "inlet_t" => Some(inlet_t.value),
            "isentropic_efficiency" => Some(isentropic_efficiency),
            "speed" => speed,
            _ => None,
        },
        &mut warnings,
    )?;

    // The chart type's three non-default values are the declaration's `[[unported]]` rows.
    match chart_type.unwrap_or("fan_law") {
        "fan_law" => {}
        "interpolate" => return Err(unported::refuse("chart_type=interpolate")),
        "interpolate_and_extrapolate" => {
            return Err(unported::refuse("chart_type=interpolate_and_extrapolate"));
        }
        "khader_2015" => return Err(unported::refuse("chart_type=khader_2015")),
        other => {
            return Err(AzothError::invalid_input(
                "chart_type",
                format!("'{other}' is not one of the four declared chart types"),
            ));
        }
    }

    let feed = Stream::from_pt(
        components.to_vec(),
        inlet_z.to_vec(),
        inlet_n,
        inlet_p,
        inlet_t,
    )?;

    let points: Vec<CurvePoint> = curve
        .unwrap_or(&[])
        .iter()
        .map(|point| CurvePoint {
            speed: point.speed,
            flow: point.flow.get::<cubic_meter_per_hour>(),
            head: point.head.get::<meter>(),
            efficiency: point.efficiency,
        })
        .collect();

    let (out, head, efficiency, power) = if points.is_empty() {
        let pressure = outlet_pressure.ok_or_else(|| {
            AzothError::invalid_input(
                "outlet_pressure",
                "a compressor without a curve reaches only the pressure it is given, and this \
                 one has neither",
            )
        })?;
        let out = kernel(&feed, pressure, isentropic_efficiency)?;
        let power = azoth_core::units::watts(feed.n * (out.h.value - feed.h.value));
        (out, 0.0, 0.0, power)
    } else {
        let chart = CompressorChart::new(&points)?;
        let outcome = compressor_with_chart(&feed, &chart, speed.unwrap_or(3000.0))?;
        (
            outcome.outlet,
            outcome.polytropic_head,
            outcome.polytropic_efficiency,
            outcome.power,
        )
    };

    Ok(CompressorResult {
        outlet_n: out.n,
        outlet_z: out.z,
        outlet_p: out.p,
        outlet_t: out.t,
        outlet_h: joules_per_mole(out.h.value),
        polytropic_head: meters(head),
        polytropic_efficiency: efficiency,
        power,
        warnings,
    })
}
