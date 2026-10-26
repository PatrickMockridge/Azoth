//! `process.pump` - the pump's kernel as a registered id.
//!
//! Spec: `specs/models/process/pump.toml`. The arithmetic is
//! [`crate::kernels::pump`]; what is here is the boundary a case and a cross-impl test
//! address, which is the same split every other namespace makes.

use azoth_core::units::{
    DiffusionCoefficient, Length, MassDensity, MolarEnergy, Power, Pressure,
    ThermodynamicTemperature, VolumeRate, cubic_meter_per_hour, joules_per_mole, meter, meters,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};

use crate::kernels::pump as kernel;
use crate::kernels::pump::{
    cavitating, npsh_available, npsh_required, pump_with_chart, shaft_power,
};
use crate::model_gen;
use crate::pump_chart::{CurvePoint, PumpChart};
use crate::segment::phase::system_volumetric_flow;
use crate::stream::Stream;

/// One row of `curve`: a measured point, in the units the spec declares.
#[derive(Debug, Clone, PartialEq)]
pub struct PumpCurvePoint {
    /// The pump speed, rpm.
    pub speed: f64,
    /// The measured volumetric flow.
    pub flow: VolumeRate,
    /// The measured head.
    pub head: Length,
    /// The measured efficiency, as a fraction.
    pub efficiency: f64,
    /// The required NPSH at this point, where the chart states one.
    pub npsh_required: Option<Length>,
}

/// Result of `process.pump`.
#[derive(Debug, Clone, PartialEq)]
pub struct PumpResult {
    /// Molar flow out, mol/s. A pump moves the stream it is given, so this is the inlet's.
    pub outlet_n: f64,
    /// Outlet composition, one entry per component.
    pub outlet_z: Vec<f64>,
    /// Outlet pressure: the stated one, or the chart's head.
    pub outlet_p: Pressure,
    /// Outlet temperature, solved from the shifted enthalpy.
    pub outlet_t: ThermodynamicTemperature,
    /// Outlet molar enthalpy: the inlet's plus the shaft work per mole.
    pub outlet_h: MolarEnergy,
    /// The head the pump used, m: the chart's corrected head or the isentropic work over `M g`.
    pub head: Length,
    /// The efficiency the pump used: the chart's corrected value or the stated isentropic efficiency.
    pub efficiency: f64,
    /// The required NPSH at the operating flow, m.
    pub npsh_required: Length,
    /// The available NPSH at the inlet, m.
    pub npsh_available: Length,
    /// The cavitation assessment at the stated margin.
    pub cavitating: bool,
    /// The shaft power, W: the enthalpy rise times the molar flow.
    pub power: Power,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

/// Raise a stream's pressure, adding the pump's work as enthalpy.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the shapes disagree, the efficiency is outside
/// `(0, 1]`, the curve is not a curve, a pump with no chart has no outlet pressure, or the
/// flash refuses the state.
#[allow(clippy::too_many_arguments)] // one parameter per declared input
pub fn pump(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
    outlet_pressure: Option<Pressure>,
    isentropic_efficiency: f64,
    curve: Option<&[PumpCurvePoint]>,
    speed: Option<f64>,
    reference_density: Option<MassDensity>,
    viscosity: Option<DiffusionCoefficient>,
    use_viscosity_correction: Option<bool>,
    check_npsh: Option<bool>,
    npsh_margin: Option<f64>,
) -> Result<PumpResult> {
    let spec = &model_gen::PUMP_SPEC;
    let mut warnings = Vec::new();
    // One call over every declared input, for the reason `reactions.chemical_equilibrium`
    // gives: a closure that filters by name emits a warning for every check it does not
    // feed, which reads as a defect in the caller's state.
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "inlet_t" => Some(inlet_t.value),
            "isentropic_efficiency" => Some(isentropic_efficiency),
            "inlet_n" => Some(inlet_n),
            "speed" => speed,
            _ => None,
        },
        &mut warnings,
    )?;

    let feed = Stream::from_pt(
        components.to_vec(),
        inlet_z.to_vec(),
        inlet_n,
        inlet_p,
        inlet_t,
    )?;

    let speed = speed.unwrap_or(1000.0);
    let viscosity_cst = viscosity.map(|value| value.value * 1.0e6).unwrap_or(1.0);
    let check_npsh = check_npsh.unwrap_or(false);
    let npsh_margin = npsh_margin.unwrap_or(1.3);

    let points: Vec<CurvePoint> = curve
        .unwrap_or(&[])
        .iter()
        .map(|point| CurvePoint {
            speed: point.speed,
            flow: point.flow.get::<cubic_meter_per_hour>(),
            head: point.head.get::<meter>(),
            efficiency: point.efficiency,
            npsh: point.npsh_required.map(|npsh| npsh.get::<meter>()),
        })
        .collect();
    let chart = if points.is_empty() {
        None
    } else {
        Some(PumpChart::new(
            &points,
            reference_density.map(|density| density.value),
            use_viscosity_correction.unwrap_or(false),
        )?)
    };

    let (out, head, efficiency, required, available, cavitation) = match &chart {
        None => {
            let pressure = outlet_pressure.ok_or_else(|| {
                AzothError::invalid_input(
                    "outlet_pressure",
                    "a pump without a curve reaches only the pressure it is given, and this one \
                     has neither",
                )
            })?;
            let out = kernel::pump(&feed, pressure, isentropic_efficiency)?;
            // **The isentropic branch reports the head its work corresponds to.** The class's
            // chart getters have nothing to read without a curve, but the work per mole over
            // `M g` is the head the machine actually developed, and the efficiency is the one
            // the caller stated - so both outputs are a statement rather than a gap.
            let head = (out.h.value - feed.h.value)
                / (feed.molar_mass()?.value * crate::pump_chart::GRAVITY);
            let flow_m3hr = system_volumetric_flow(&feed)? * 3600.0;
            let required = npsh_required(None, flow_m3hr, speed);
            let available = npsh_available(&feed);
            let cavitation = cavitating(check_npsh, available, required, npsh_margin);
            (
                out,
                head,
                isentropic_efficiency,
                required,
                available,
                cavitation,
            )
        }
        Some(chart) => {
            let outcome =
                pump_with_chart(&feed, chart, speed, viscosity_cst, check_npsh, npsh_margin)?;
            (
                outcome.stream,
                outcome.head,
                outcome.efficiency,
                outcome.npsh_required,
                outcome.npsh_available,
                outcome.cavitating,
            )
        }
    };
    let power = shaft_power(&feed, &out);

    Ok(PumpResult {
        outlet_n: out.n,
        outlet_z: out.z,
        outlet_p: out.p,
        outlet_t: out.t,
        outlet_h: joules_per_mole(out.h.value),
        head: meters(head),
        efficiency,
        npsh_required: meters(required),
        npsh_available: meters(available),
        cavitating: cavitation,
        power,
        warnings,
    })
}
