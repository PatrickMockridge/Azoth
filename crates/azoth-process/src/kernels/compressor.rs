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
/// The class's anti-surge recycle loop has no cap; this port stops it rather than spinning.
const MAX_SURGE_PASSES: usize = 20;

/// `AntiSurge`'s state that the chart branch reads.
#[derive(Debug, Clone, Copy)]
pub struct SurgeSettings {
    /// `getAntiSurge().setActive(true)`: whether the run-time recycle block runs.
    pub active: bool,
    /// `AntiSurge.surgeControlFactor`, the recycle target's multiple of the surge flow.
    pub control_factor: f64,
    /// `setSurgeControlMargin(f)`, the control-line fraction of the surge flow.
    pub control_margin: f64,
}

impl Default for SurgeSettings {
    fn default() -> Self {
        Self {
            active: false,
            control_factor: 1.05,
            control_margin: 0.0,
        }
    }
}

/// What the compressor's chart branch reports beside the outlet.
#[derive(Debug, Clone)]
pub struct ChartOutcome {
    /// The outlet stream.
    pub outlet: Stream,
    /// The chart head the branch used, m (the class's `polytropicHeadMeter`).
    pub polytropic_head: f64,
    /// The polytropic efficiency it used, as a fraction.
    pub polytropic_efficiency: f64,
    /// The shaft power, W - the boosted flow's power when the recycle block ran.
    pub power: Power,
    /// `getSurgeFlowRate()`, m3/hr, where a surge curve is set.
    pub surge_flow_rate: Option<f64>,
    /// `getSurgeFlowRateMargin()`, m3/hr.
    pub surge_flow_rate_margin: Option<f64>,
    /// `getSurgeFlowRateStd()`, m3/hr at 1.01325 bara and 288.15 K.
    pub surge_flow_rate_std: Option<f64>,
    /// `getDistanceToSurge()`, the inlet flow over the surge flow less one.
    pub distance_to_surge: Option<f64>,
    /// `getControlLineFlow()`, m3/hr.
    pub control_line_flow: Option<f64>,
    /// `getDistanceToControlLine()`.
    pub distance_to_control_line: Option<f64>,
    /// `getRequiredRecycleFractionToControlLine()`.
    pub required_recycle_fraction: Option<f64>,
    /// `getAntiSurgeRecyclePower(required fraction)`, W.
    pub anti_surge_recycle_power: Option<f64>,
    /// `AntiSurge.getCurrentSurgeFraction()`, present only when the block was active.
    pub recycle_fraction: Option<f64>,
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
/// **This is `Compressor.run`'s fixed-speed chart branch**, `Compressor.java:1466-1545`. The
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
/// **The anti-surge block re-runs the machine at a boosted flow.** With `settings.active` and
/// the operating point inside the surge curve, the class sets the system's total flow to
/// `surgeControlFactor * surgeFlow(chart head)`, recomputes the fraction
/// `newMoles/originalMoles - 1`, and repeats until the point is off the curve. The outlet
/// keeps the **fresh** moles (the class restores them after the branch), the reported power is
/// the boosted flow's, and the read-outs evaluate the curve at the *post-run* fluid head
/// `power / (fresh mass flow) * efficiency` - which, after a recycle, is the chart head
/// inflated by `1 + fraction`, because the class divides the recycled power by the restored
/// flow.
///
/// **The class's loop is uncapped and would spin on a curve it cannot leave**; this stops at
/// [`MAX_SURGE_PASSES`] and reports the last pass rather than hanging.
///
/// # Errors
/// `InvalidInput` for a non-positive speed or flow, and whatever the flash, the density or the
/// heat capacity refuses.
pub fn compressor_with_chart(
    feed: &Stream,
    chart: &CompressorChart,
    speed: f64,
    settings: &SurgeSettings,
) -> Result<ChartOutcome> {
    if speed.is_nan() || speed <= 0.0 {
        return Err(AzothError::invalid_input(
            "speed",
            format!("the chart's reduced flow divides by the speed, and {speed} rpm does not"),
        ));
    }

    let fresh_flow_m3hr = system_volumetric_flow(feed)? * 3600.0;
    let fresh_n = feed.n;
    if fresh_n.is_nan() || fresh_n <= 0.0 || fresh_flow_m3hr.is_nan() || fresh_flow_m3hr <= 0.0 {
        return Err(AzothError::invalid_input(
            "inlet_n",
            format!(
                "the chart branch needs a positive flow to scale, and {fresh_n} mol/s gives \
                 {fresh_flow_m3hr} m3/hr"
            ),
        ));
    }
    let molar_mass = feed.molar_mass()?.value;
    // The cubic's own `Z`, which is `getZ()`: `p M/(rho_cubic R T)`.
    let z = feed.p.value * molar_mass / (feed.density()? * R_J_PER_MOL_K * feed.t.value);
    let cp = feed.molar_heat_capacity()?;
    let kappa = cp / (cp - R_J_PER_MOL_K);
    // `getZvolcorr` = `P v/(R T)` at the Peneloux-shifted molar volume.
    let shift = {
        let (mixture, _) = feed.mixture()?;
        mixture.volume_shift(&feed.z)
    };
    let zvolcorr = z - feed.p.value * shift / (R_J_PER_MOL_K * feed.t.value);

    let mut boosted_n = fresh_n;
    let mut passes = 0;
    let (efficiency, head_m, head_j_per_kg, pressure_ratio) = loop {
        let flow_m3hr = fresh_flow_m3hr * (boosted_n / fresh_n);
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

        if settings.active && chart.is_surge(head_m, flow_m3hr) {
            if let Some(target) = chart.surge_flow(head_m) {
                if target.is_finite() && target > 0.0 {
                    boosted_n = fresh_n * target / fresh_flow_m3hr;
                    passes += 1;
                    if passes < MAX_SURGE_PASSES {
                        continue;
                    }
                }
            }
        }
        break (efficiency, head_m, head_j_per_kg, pressure_ratio);
    };

    let outlet_pressure = pascals(feed.p.value * pressure_ratio);
    let dh = head_j_per_kg * molar_mass / efficiency;
    let outlet = Stream::from_ph(
        feed.components.clone(),
        feed.z.clone(),
        fresh_n,
        outlet_pressure,
        joules_per_mole(feed.h.value + dh),
    )?;
    let power = watts(boosted_n * dh);
    let recycle_fraction = boosted_n / fresh_n - 1.0;

    let read = chart
        .surge_flow(power.value / (fresh_n * molar_mass) * efficiency / METRE_TO_JOULE_PER_KG)
        .filter(|surge_flow| *surge_flow > 0.0)
        .map(|surge_flow| {
            let control = surge_flow * (1.0 + settings.control_margin);
            let required = if control > 0.0 {
                ((control - fresh_flow_m3hr) / control).max(0.0)
            } else {
                0.0
            };
            SurgeReadOut {
                surge_flow_rate: surge_flow,
                surge_flow_rate_margin: fresh_flow_m3hr - surge_flow,
                surge_flow_rate_std: surge_flow * (feed.p.value / 1.0e5) / 1.01325 * 288.15
                    / feed.t.value
                    / zvolcorr,
                distance_to_surge: fresh_flow_m3hr / surge_flow - 1.0,
                control_line_flow: control,
                distance_to_control_line: if control > 0.0 {
                    Some(fresh_flow_m3hr / control - 1.0)
                } else {
                    None
                },
                required_recycle_fraction: required,
                // `getAntiSurgeRecyclePower(currentSurgeFraction)`: the power the run-time
                // recycle actually added, which is zero unless the block fired.
                anti_surge_recycle_power: power.value * recycle_fraction.max(0.0),
            }
        });

    Ok(ChartOutcome {
        outlet,
        polytropic_head: head_m,
        polytropic_efficiency: efficiency,
        power,
        surge_flow_rate: read.as_ref().map(|row| row.surge_flow_rate),
        surge_flow_rate_margin: read.as_ref().map(|row| row.surge_flow_rate_margin),
        surge_flow_rate_std: read.as_ref().map(|row| row.surge_flow_rate_std),
        distance_to_surge: read.as_ref().map(|row| row.distance_to_surge),
        control_line_flow: read.as_ref().map(|row| row.control_line_flow),
        distance_to_control_line: read.as_ref().and_then(|row| row.distance_to_control_line),
        required_recycle_fraction: read.as_ref().map(|row| row.required_recycle_fraction),
        anti_surge_recycle_power: read.as_ref().map(|row| row.anti_surge_recycle_power),
        recycle_fraction: settings.active.then_some(recycle_fraction),
    })
}

/// One row of the anti-surge read-outs, computed together with the control line.
struct SurgeReadOut {
    surge_flow_rate: f64,
    surge_flow_rate_margin: f64,
    surge_flow_rate_std: f64,
    distance_to_surge: f64,
    control_line_flow: f64,
    distance_to_control_line: Option<f64>,
    required_recycle_fraction: f64,
    anti_surge_recycle_power: f64,
}
