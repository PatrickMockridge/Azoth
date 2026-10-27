//! `unit_ops.compressor` - the isentropic route, through the kernel the pump landed.

use azoth_core::units::{Power, Pressure, joules_per_mole, pascals, watts};
use azoth_core::{AzothError, Result};

use crate::compressor_chart::CompressorChart;
use crate::kernels::pump::pump;
use crate::segment::phase::system_volumetric_flow;
use crate::stream::Stream;

/// `ThermodynamicConstantsInterface.R`, J/(mol*K), as the class's own constant carries it.
const R_J_PER_MOL_K: f64 = 8.314;
/// `Compressor.MIN_EFFICIENCY`: the clamp `setPolytropicEfficiency` applies.
const MIN_EFFICIENCY: f64 = 1.0e-6;
/// `Compressor`'s own metre-to-specific-energy factor (`9.81`, not standard gravity).
const METRE_TO_JOULE_PER_KG: f64 = 9.81;

/// What the compressor's chart branch reports beside the outlet.
#[derive(Debug, Clone)]
pub struct ChartOutcome {
    /// The outlet stream.
    pub outlet: Stream,
    /// The chart head the branch used, m (the class's `polytropicHeadMeter`).
    pub polytropic_head: f64,
    /// The polytropic efficiency it used, as a fraction.
    pub polytropic_efficiency: f64,
    /// The shaft power, W.
    pub power: Power,
}

/// Raise a stream's pressure along an isentrope, dividing the step by the efficiency.
///
/// **This is `unit_ops.pump`'s arithmetic, and NeqSim's own code says so.** `Pump.run` takes
/// the compressor route because `calculateAsCompressor` defaults to `true`, and
/// `Compressor.run`'s no-chart branch is that route verbatim: set the outlet pressure,
/// `PSflash` at the *inlet's own entropy*, and move the enthalpy by
/// `(h_isentropic - h_in) / eta`. Writing the arithmetic out again here would be a second
/// implementation of one step, so this calls it.
///
/// **The entropy is derived, not carried.** The record has five fields and `s` is not one of
/// them: it is a function of `(T, P, z)`, and the flash the kernel already runs at the inlet
/// decides it. The captured isentropic row is the demonstration - at an efficiency of one the
/// outlet's entropy equals the inlet's to `1.3e-14`, and the step's entropy production is
/// `-1.4e-17` kJ/(mol K), which is zero.
///
/// **The shaft power is not an output**, for the reason `s` is not a field: `Compressor.getPower`
/// is `dH = h_out - h_in` over the machine, so the record this returns already carries it -
/// `n * (outlet_h - inlet_h)`, which the captured rows satisfy to their last digit.
///
/// Named out: the compressor chart and catalogue, the polytropic branch (the same rule with
/// the polytropic efficiency on the other side of the division), the gamma-based
/// pressure-ratio mode, the anti-surge machinery, the deposit model, and the
/// `useGERG2008`/`useLeachman`/`useVega` property overrides.
pub fn compressor(feed: &Stream, outlet_pressure: Pressure, efficiency: f64) -> Result<Stream> {
    pump(feed, outlet_pressure, efficiency)
}

/// Raise a stream along the compressor's chart, solving the pressure from the head.
///
/// **This is `Compressor.run`'s fixed-speed chart branch**, `Compressor.java:1466-1540`. The
/// inlet's actual volumetric flow reads the chart; the chart's polytropic efficiency is
/// clamped to `[1e-6, 1]` (the class's `normalizeEfficiency`); the polytropic exponent `n` is
/// the ideal-gas kappa divided by the efficiency; the pressure ratio is the head over the
/// polytropic flow work; and the enthalpy rise is `H M / eta` per mole, flashed to the outlet
/// pressure. **The class's `dH` is the total-joule form** `H M/eta * moles`, which is this
/// per-mole step times the flow.
///
/// **The pressure is an answer, not an input.** The chart branch sets
/// `thermoSystem.getPressure() * pressureRatio`, so a stated outlet pressure is ignored while
/// a curve is present - the same rule as the pump's chart.
///
/// # Errors
/// `InvalidInput` for a non-positive speed, and whatever the flash, the density or the heat
/// capacity refuses.
pub fn compressor_with_chart(
    feed: &Stream,
    chart: &CompressorChart,
    speed: f64,
) -> Result<ChartOutcome> {
    if speed.is_nan() || speed <= 0.0 {
        return Err(AzothError::invalid_input(
            "speed",
            format!("the chart's reduced flow divides by the speed, and {speed} rpm does not"),
        ));
    }

    let flow_m3hr = system_volumetric_flow(feed)? * 3600.0;
    let molar_mass = feed.molar_mass()?.value;
    // The cubic's own `Z`, which is `getZ()`: `p M/(rho_cubic R T)`.
    let z = feed.p.value * molar_mass / (feed.density()? * R_J_PER_MOL_K * feed.t.value);
    let cp = feed.molar_heat_capacity()?;
    let kappa = cp / (cp - R_J_PER_MOL_K);

    let efficiency = chart
        .polytropic_efficiency(flow_m3hr, speed)
        .clamp(MIN_EFFICIENCY, 1.0);
    // **The chart's head crosses in metres, and the class's metre branch converts it with
    // `9.81`, not standard gravity**: `polytropicFluidHead = head/1000*9.81` kJ/kg, so the
    // specific work is `head * 9.81` J/kg.
    let head_m = chart.polytropic_head(flow_m3hr, speed);
    let head_j_per_kg = head_m * METRE_TO_JOULE_PER_KG;
    let n = 1.0 / (1.0 - (kappa - 1.0) / kappa / efficiency);
    let flow_work = n / (n - 1.0) * z * R_J_PER_MOL_K * feed.t.value / molar_mass;
    if !flow_work.is_finite() || flow_work <= 0.0 || !n.is_finite() {
        return Err(AzothError::invalid_input(
            "curve",
            format!(
                "the chart's head and the inlet's state give no polytropic ratio: n={n}, \
                 flow work={flow_work}"
            ),
        ));
    }
    let pressure_ratio = ((head_j_per_kg + flow_work) / flow_work).powf(n / (n - 1.0));
    let outlet_pressure = pascals(feed.p.value * pressure_ratio);
    let dh = head_j_per_kg * molar_mass / efficiency;
    let outlet = Stream::from_ph(
        feed.components.clone(),
        feed.z.clone(),
        feed.n,
        outlet_pressure,
        joules_per_mole(feed.h.value + dh),
    )?;

    Ok(ChartOutcome {
        outlet,
        polytropic_head: head_m,
        polytropic_efficiency: efficiency,
        power: watts(feed.n * dh),
    })
}
