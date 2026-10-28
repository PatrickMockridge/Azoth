//! `unit_ops.compressor` - the isentropic route, through the kernel the pump landed.

use azoth_core::units::{Power, Pressure, joules_per_mole, pascals, watts};
use azoth_core::{AzothError, Result};

use azoth_eos::{Phase, RootSide, molar_enthalpy_entropy, pt_flash};

use crate::compressor_chart::CompressorChart;
use crate::kernels::pump::pump;
use crate::stream::Stream;

/// `ThermodynamicConstantsInterface.R`, J/(mol*K), as the class's own constant carries it.
const R_J_PER_MOL_K: f64 = 8.314;
/// `Compressor.MIN_EFFICIENCY`: the clamp `setPolytropicEfficiency` applies.
const MIN_EFFICIENCY: f64 = 1.0e-6;
/// `Compressor`'s own metre-to-specific-energy factor (`9.81`, not standard gravity).
const METRE_TO_JOULE_PER_KG: f64 = 9.81;
/// The class's anti-surge recycle loop has no cap; this port stops it rather than spinning.
const MAX_SURGE_PASSES: usize = 20;
/// `Compressor.run`'s speed-solve loop: `while (iteration < 200)`, so 199 passes.
const SPEED_SOLVE_ITERATIONS: usize = 200;
/// The numerical derivative's speed increment, `deltaSpeed`.
const SPEED_SOLVE_DELTA_RPM: f64 = 100.0;
/// `Math.abs(currentPressure - targetPressure) <= 1e-3` bar, in pascals.
const SPEED_SOLVE_TOLERANCE_PA: f64 = 100.0;
/// The derivative floor, `1e-6` bar/rpm, in pascals per rpm.
const SPEED_SOLVE_DERIVATIVE_FLOOR_PA: f64 = 1.0e-4;
/// The fallback slope, `0.01` bar/rpm, in pascals per rpm.
const SPEED_SOLVE_DERIVATIVE_FALLBACK_PA: f64 = 1000.0;
/// `polytropicFluidHead = 0.0001` kJ/kg, in joules per kilogram.
const MIN_POLYTROPIC_HEAD_J_PER_KG: f64 = 0.1;

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
    /// The speed the chart was read at, rpm: the stated one or the solved one.
    pub speed: f64,
    /// Whether the speed-solve mode pegged the speed at a bound.
    pub speed_limited: bool,
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

    let fresh_flow_m3hr = feed_volumetric_flow(feed)? * 3600.0;
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
    // **The feed's own root, at the feed's composition** - `getZ()` and `getCP()` on the one
    // phase the class holds, not the flash's incipient phase.
    let (z, cp) = feed_z_cp(feed)?;
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
            .polytropic_efficiency(flow_m3hr, speed)?
            .clamp(MIN_EFFICIENCY, 1.0);
        // **The chart's head crosses in metres, and the class's metre branch converts it with
        // `9.81`, not standard gravity**: `polytropicFluidHead = head/1000*9.81` kJ/kg, so the
        // specific work is `head * 9.81` J/kg.
        let head_m = chart.polytropic_head(flow_m3hr, speed)?;
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
                    boosted_n = fresh_n * settings.control_factor * target / fresh_flow_m3hr;
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

    let read = surge_read_out(
        chart,
        settings,
        feed,
        zvolcorr,
        fresh_flow_m3hr,
        power.value / (fresh_n * molar_mass) * efficiency / METRE_TO_JOULE_PER_KG,
        power,
        boosted_n,
        fresh_n,
    );

    Ok(ChartOutcome {
        outlet,
        polytropic_head: head_m,
        polytropic_efficiency: efficiency,
        power,
        speed,
        speed_limited: false,
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

/// The surge read-outs at `head_for_curve` and the fresh flow.
///
/// **The head is the caller's**, because the chart branch evaluates the curve at the
/// *post-run* fluid head `power / (fresh mass flow) * efficiency`, which after a recycle is
/// the chart head inflated by `1 + fraction`, while the speed-solve mode reads the head the
/// loop last took from the chart.
#[allow(clippy::too_many_arguments)] // the read-outs are one expression over these locals
fn surge_read_out(
    chart: &CompressorChart,
    settings: &SurgeSettings,
    feed: &Stream,
    zvolcorr: f64,
    fresh_flow_m3hr: f64,
    head_for_curve: f64,
    power: Power,
    boosted_n: f64,
    fresh_n: f64,
) -> Option<SurgeReadOut> {
    let surge_flow = chart
        .surge_flow(head_for_curve)
        .filter(|surge_flow| *surge_flow > 0.0)?;
    let control = surge_flow * (1.0 + settings.control_margin);
    let required = if control > 0.0 {
        ((control - fresh_flow_m3hr) / control).max(0.0)
    } else {
        0.0
    };
    Some(SurgeReadOut {
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
        // `getAntiSurgeRecyclePower(currentSurgeFraction)`: the power the run-time recycle
        // actually added, which is zero unless the block fired.
        anti_surge_recycle_power: power.value * (boosted_n / fresh_n - 1.0).max(0.0),
    })
}

/// `Compressor.run`'s `solveSpeed` branch: the speed that reaches a stated pressure.
#[derive(Debug, Clone, Copy)]
pub struct SpeedSolveSettings {
    /// `getOutletPressure()`: the pressure the chart's ratio is driven to.
    pub target: Pressure,
    /// `getSpeed()` at entry: the Newton's starting point.
    pub initial_speed: f64,
    /// `getMinimumSpeed()`, rpm.
    pub min_speed: f64,
    /// `getMaximumSpeed()`, rpm.
    pub max_speed: f64,
    /// `isLimitSpeed()`: peg at a bound rather than report an out-of-range speed.
    pub limit_speed: bool,
}

/// The speed-solve mode: `Compressor.java:1261-1453`.
///
/// **A Newton on the shaft speed, with a numerical derivative.** Each pass reads the chart's
/// polytropic efficiency and head at the trial speed, forms the polytropic pressure ratio and
/// compares it with the stated pressure; the derivative is the same ratio at `speed + 100` rpm,
/// the step is relaxed by `min(0.8, i/(i+3))` and capped at `max(500, 0.2 N)`, and the loop
/// stops at `1e-3` bar or after 199 passes. **The efficiency clamp is the class's percent
/// form**: a non-positive draw becomes `0.01` and a draw above `100` becomes `100`
/// (a fraction of `100`), which is what `setPolytropicEfficiency` stores.
///
/// **The convergence test reads the pre-update pressure and the update sets the speed.** The
/// class computes the pressure, updates `currentSpeed`, then tests the old pressure - so the
/// reported speed is one Newton step past the head the outlet uses. The port copies that order.
///
/// **`limit_speed` re-runs the fixed-speed branch at the bound**, which is the class's own
/// `setCalcPressureOut(true)` fallback: its intermediate pressure is overwritten by the chart
/// branch, so the answer is that branch at the pegged speed.
///
/// # Errors
/// `InvalidInput` for a non-positive or non-finite starting speed, a chart read that gives no
/// polytropic ratio, or a derivative that is not a number; and whatever the flash refuses.
pub fn compressor_speed_solve(
    feed: &Stream,
    chart: &CompressorChart,
    settings: &SpeedSolveSettings,
    surge: &SurgeSettings,
) -> Result<ChartOutcome> {
    if !settings.initial_speed.is_finite() || settings.initial_speed <= 0.0 {
        return Err(AzothError::invalid_input(
            "speed",
            format!(
                "the speed-solve mode starts from the machine's speed, and {} rpm does not",
                settings.initial_speed
            ),
        ));
    }
    let fresh_flow_m3hr = feed_volumetric_flow(feed)? * 3600.0;
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
    let (z, cp) = feed_z_cp(feed)?;
    let kappa = cp / (cp - R_J_PER_MOL_K);
    let shift = {
        let (mixture, _) = feed.mixture()?;
        mixture.volume_shift(&feed.z)
    };
    let zvolcorr = z - feed.p.value * shift / (R_J_PER_MOL_K * feed.t.value);

    let chart_min = chart.min_curve_speed();
    let chart_max = chart.max_curve_speed();
    if !chart_min.is_finite() || !chart_max.is_finite() || chart_max <= chart_min {
        return Err(AzothError::invalid_input(
            "curve",
            format!(
                "the speed-solve mode needs a measured speed range, and the chart's is \
                 {chart_min} to {chart_max} rpm"
            ),
        ));
    }
    let mut speed = settings.initial_speed;
    if speed < chart_min * 0.5 || speed > chart_max * 2.0 {
        speed = 0.5 * (chart_min + chart_max);
    }
    let mut delta = SPEED_SOLVE_DELTA_RPM;
    let mut efficiency = 0.0;
    let mut head_m = 0.0;
    let mut dh_molar = 0.0;
    let mut converged: Option<f64> = None;
    let mut limited: Option<f64> = None;

    for iteration in 1..SPEED_SOLVE_ITERATIONS {
        let percent = chart.polytropic_efficiency(fresh_flow_m3hr, speed)? * 100.0;
        let (stored, formula) = if percent <= 0.0 {
            (0.01, 0.01)
        } else if percent > 100.0 {
            (100.0, 1.0)
        } else {
            (percent / 100.0, percent / 100.0)
        };
        head_m = chart.polytropic_head(fresh_flow_m3hr, speed)?;
        let mut head_j = head_m * METRE_TO_JOULE_PER_KG;
        if head_j <= 0.0 {
            head_j = MIN_POLYTROPIC_HEAD_J_PER_KG;
        }
        let current = feed.p.value
            * polytropic_pressure_ratio(head_j, z, feed.t.value, molar_mass, kappa, formula)?;
        // The derivative's read is never clamped, as the class's second read is not.
        let percent_delta = chart.polytropic_efficiency(fresh_flow_m3hr, speed + delta)? * 100.0;
        let head_delta =
            chart.polytropic_head(fresh_flow_m3hr, speed + delta)? * METRE_TO_JOULE_PER_KG;
        let ratio_delta = polytropic_pressure_ratio(
            head_delta,
            z,
            feed.t.value,
            molar_mass,
            kappa,
            percent_delta / 100.0,
        )?;
        let next = feed.p.value * ratio_delta;
        let mut derivative = (next - current) / delta;
        if !derivative.is_finite() {
            return Err(AzothError::invalid_input(
                "curve",
                format!(
                    "the chart's derivative at {speed} rpm is not a number: \
                     efficiency {percent_delta}%, head {head_delta} J/kg"
                ),
            ));
        }
        if derivative.abs() < SPEED_SOLVE_DERIVATIVE_FLOOR_PA {
            derivative = SPEED_SOLVE_DERIVATIVE_FALLBACK_PA;
        }
        let relaxation = 0.8_f64.min(iteration as f64 / (iteration as f64 + 3.0));
        let cap = 500.0_f64.max(0.2 * speed);
        let mut update = (settings.target.value - current) / derivative;
        if update.abs() > cap {
            update = update.signum() * cap;
        }
        speed += relaxation * update;
        efficiency = stored;
        dh_molar = head_j * molar_mass / stored;
        if speed < 0.0 {
            speed = if settings.min_speed > 1.0 {
                settings.min_speed
            } else {
                chart_min
            };
        }
        if iteration % 10 == 0 && delta > 10.0 {
            delta /= 2.0;
        }
        if (speed < settings.min_speed || speed > settings.max_speed) && settings.limit_speed {
            limited = Some(if speed > settings.max_speed {
                settings.max_speed
            } else {
                settings.min_speed
            });
            break;
        }
        if (current - settings.target.value).abs() <= SPEED_SOLVE_TOLERANCE_PA {
            converged = Some(speed);
            break;
        }
    }

    if let Some(bound) = limited {
        let mut outcome = compressor_with_chart(feed, chart, bound, surge)?;
        outcome.speed_limited = true;
        return Ok(outcome);
    }
    let speed = converged.unwrap_or(speed);

    // The class's post-convergence surge check: one pass, no re-solve of the speed.
    let mut boosted_n = fresh_n;
    if surge.active && chart.is_surge(head_m, fresh_flow_m3hr) {
        if let Some(target_flow) = chart.surge_flow(head_m) {
            if target_flow.is_finite() && target_flow > 0.0 {
                boosted_n = fresh_n * surge.control_factor * target_flow / fresh_flow_m3hr;
            }
        }
    }
    let outlet = Stream::from_ph(
        feed.components.clone(),
        feed.z.clone(),
        boosted_n,
        settings.target,
        joules_per_mole(feed.h.value + dh_molar),
    )?;
    let power = watts(boosted_n * dh_molar);
    let read = surge_read_out(
        chart,
        surge,
        feed,
        zvolcorr,
        fresh_flow_m3hr,
        head_m,
        power,
        boosted_n,
        fresh_n,
    );

    Ok(ChartOutcome {
        outlet,
        polytropic_head: head_m,
        polytropic_efficiency: efficiency,
        power,
        speed,
        speed_limited: false,
        surge_flow_rate: read.as_ref().map(|row| row.surge_flow_rate),
        surge_flow_rate_margin: read.as_ref().map(|row| row.surge_flow_rate_margin),
        surge_flow_rate_std: read.as_ref().map(|row| row.surge_flow_rate_std),
        distance_to_surge: read.as_ref().map(|row| row.distance_to_surge),
        control_line_flow: read.as_ref().map(|row| row.control_line_flow),
        distance_to_control_line: read.as_ref().and_then(|row| row.distance_to_control_line),
        required_recycle_fraction: read.as_ref().map(|row| row.required_recycle_fraction),
        anti_surge_recycle_power: read.as_ref().map(|row| row.anti_surge_recycle_power),
        recycle_fraction: surge.active.then_some(boosted_n / fresh_n - 1.0),
    })
}

/// The class's polytropic pressure ratio: `((H + W)/W)^(n/(n-1))`.
fn polytropic_pressure_ratio(
    head_j_per_kg: f64,
    z: f64,
    temperature: f64,
    molar_mass: f64,
    kappa: f64,
    efficiency: f64,
) -> Result<f64> {
    let n = 1.0 / (1.0 - (kappa - 1.0) / kappa / efficiency);
    let flow_work = n / (n - 1.0) * z * R_J_PER_MOL_K * temperature / molar_mass;
    if !flow_work.is_finite() || flow_work <= 0.0 || !n.is_finite() {
        return Err(AzothError::invalid_input(
            "curve",
            format!(
                "the chart's head and the inlet's state give no polytropic ratio: n={n}, \
                 flow work={flow_work}"
            ),
        ));
    }
    Ok(((head_j_per_kg + flow_work) / flow_work).powf(n / (n - 1.0)))
}

/// `getFlowRate("m3/sec")` at the feed's own composition.
///
/// **The feed's root, not the flash's incipient one.** `pt_flash` on a superheated stream
/// reports the *incipient* phase's root, whose composition can sit far from the feed's; the
/// class reads the one phase it holds. A two-phase feed is the beta-weighted mean of its
/// phases, which is what the system flow sums.
fn feed_volumetric_flow(feed: &Stream) -> Result<f64> {
    let (mixture, _) = feed.mixture()?;
    let flash = pt_flash(&mixture, feed.t, feed.p, &feed.z)?;
    let reduced = mixture.reduced_parameters(feed.t, feed.p)?;
    let volume = |composition: &[f64], root: f64| -> f64 {
        root * R_J_PER_MOL_K * feed.t.value / feed.p.value - mixture.volume_shift(composition)
    };
    let molar_volume = match flash.phase {
        Phase::AllVapour | Phase::AllLiquid | Phase::Trivial => {
            let side = if matches!(flash.phase, Phase::AllVapour) {
                RootSide::Vapour
            } else {
                RootSide::Liquid
            };
            let root = mixture.phase_state(&reduced, &feed.z, side)?.z;
            volume(&feed.z, root)
        }
        Phase::TwoPhase => {
            let beta = flash.vapour_fraction.unwrap_or(1.0);
            beta * volume(&flash.y, flash.z_vapour)
                + (1.0 - beta) * volume(&flash.x, flash.z_liquid)
        }
    };
    Ok(feed.n * molar_volume)
}

/// `getZ()` and `getCP()/n` at the feed's own composition.
///
/// **The feed's root, not the flash's incipient one**, for the reason
/// [`feed_volumetric_flow`] gives. A two-phase feed is the beta-weighted mean, which is what
/// the class's `getZ()` sums.
fn feed_z_cp(feed: &Stream) -> Result<(f64, f64)> {
    let (mixture, ideal_gas) = feed.mixture()?;
    let flash = pt_flash(&mixture, feed.t, feed.p, &feed.z)?;
    let reduced = mixture.reduced_parameters(feed.t, feed.p)?;
    let single = |composition: &[f64], side: RootSide| -> Result<(f64, f64)> {
        let root = mixture.phase_state(&reduced, composition, side)?.z;
        let cp = molar_enthalpy_entropy(&mixture, &ideal_gas, feed.t, feed.p, composition, root)?
            .cp
            .value;
        Ok((root, cp))
    };
    match flash.phase {
        Phase::AllVapour => single(&feed.z, RootSide::Vapour),
        Phase::AllLiquid | Phase::Trivial => single(&feed.z, RootSide::Liquid),
        Phase::TwoPhase => {
            let beta = flash.vapour_fraction.unwrap_or(1.0);
            let (vapour_z, vapour_cp) = single(&flash.y, RootSide::Vapour)?;
            let (liquid_z, liquid_cp) = single(&flash.x, RootSide::Liquid)?;
            Ok((
                beta * vapour_z + (1.0 - beta) * liquid_z,
                beta * vapour_cp + (1.0 - beta) * liquid_cp,
            ))
        }
    }
}
