//! `TurboExpanderCompressor` - an expander and a compressor on one shaft.
//!
//! **The speed is the unknown, and the balance is the equation.** The expander's isentropic
//! head comes from one `PSflash` to the stated expander outlet pressure; its efficiency is the
//! velocity-ratio parabola (scaled by the design efficiency) or, when a 2-D map is present,
//! the map read at the *raw* `U/C` - the class's own asymmetry, the parabola taking the
//! normalised ratio. The compressor's absorbed power is its design polytropic head scaled by
//! `(N/N_design)²` and a Q/N head factor, over an efficiency that is a Q/N factor times the
//! design value. A damped Newton on `N` drives the expander's power to the compressor's plus
//! the speed-scaled bearing loss.
//!
//! **The replay is two different machines.** The expander outlet is an isentropic expansion at
//! the solved efficiency, which the class allows above one for the inverse mode; the
//! compressor outlet is the polytropic pressure ratio a stated power implies - the class's
//! `setPower` + `setCalcPressureOut` route.

use std::f64::consts::PI;

use azoth_core::units::{
    Pressure, ThermodynamicTemperature, joules_per_mole, joules_per_mole_kelvin, pascals,
};
use azoth_core::{AzothError, Result};
use azoth_eos::{Phase, RootSide, molar_enthalpy_entropy, ph_flash, ps_flash, pt_flash};

use crate::expander_chart::ExpanderChart;
use crate::kernels::expander::expander_replay;
use crate::stream::Stream;

/// `ThermodynamicConstantsInterface.R`, J/(mol*K).
pub const R_J_PER_MOL_K: f64 = 8.314;
/// `run`'s hard-coded speed bounds, rpm.
pub const MIN_SPEED: f64 = 1000.0;
/// See [`MIN_SPEED`].
pub const MAX_SPEED: f64 = 9000.0;
/// `run`'s Newton caps and minimum passes.
const MAX_ITERATIONS: usize = 50;
const MIN_ITERATIONS: usize = 3;
/// `Compressor`'s own metre-to-specific-energy factor.
const METRE_TO_JOULE_PER_KG: f64 = 9.81;
/// The class's no-fit parabola.
const LEGACY_PARABOLA: f64 = -3.56;

/// `TurboExpanderCompressor`'s defaults.
pub const DEFAULT_DESIGN_SPEED: f64 = 6850.0;
pub const DEFAULT_DESIGN_UC: f64 = 0.7;
pub const DEFAULT_DESIGN_QN: f64 = 0.03328;
pub const DEFAULT_IMPELLER_DIAMETER: f64 = 0.424;
pub const DEFAULT_MAXIMUM_IGV_AREA: f64 = 1.637e-2;
pub const DEFAULT_IGV_AREA_INCREASE_FACTOR: f64 = 1.14;
pub const DEFAULT_BEARING_LOSS: f64 = 10.0;
pub const DEFAULT_EXPANDER_DESIGN_EFFICIENCY: f64 = 1.0;
pub const DEFAULT_COMPRESSOR_DESIGN_EFFICIENCY: f64 = 0.81;
pub const DEFAULT_COMPRESSOR_DESIGN_HEAD_M: f64 = 20.47 * 1000.0 / METRE_TO_JOULE_PER_KG;
pub const DEFAULT_IGV_OPENING: f64 = 1.0;

/// The machine's configuration, in the class's units.
#[derive(Debug, Clone)]
pub struct TurboSettings<'a> {
    /// The expander's stated outlet pressure.
    pub expander_out_pressure: Pressure,
    /// The impeller diameter, m.
    pub impeller_diameter: f64,
    /// The design speed, rpm.
    pub design_speed: f64,
    /// The design velocity ratio.
    pub design_uc: f64,
    /// The design flow coefficient, m3/rev.
    pub design_qn: f64,
    /// The maximum IGV area, m2.
    pub maximum_igv_area: f64,
    /// The expanded-area factor.
    pub igv_area_increase_factor: f64,
    /// The stated IGV opening, 0..1.
    pub igv_opening: f64,
    /// Whether the stated opening is held rather than recomputed.
    pub igv_control_mode: bool,
    /// The bearing loss at the design speed, W.
    pub bearing_loss: f64,
    /// The expander's design isentropic efficiency.
    pub expander_design_efficiency: f64,
    /// The compressor's design polytropic efficiency.
    pub compressor_design_efficiency: f64,
    /// The compressor's design polytropic head, m.
    pub compressor_design_head_m: f64,
    /// The legacy velocity-ratio curve, `(uc, efficiency)` pairs.
    pub uc_curve: &'a [(f64, f64)],
    /// The Q/N efficiency correction curve, `(qn, factor)` pairs.
    pub qn_efficiency_curve: &'a [(f64, f64)],
    /// The Q/N head correction curve, `(qn, factor)` pairs.
    pub qn_head_curve: &'a [(f64, f64)],
    /// The 2-D map, where one is supplied.
    pub chart: Option<&'a ExpanderChart>,
    /// The IGV efficiency penalty curve, `(igv, factor)` pairs.
    pub igv_penalty_curve: &'a [(f64, f64)],
    /// A target expander outlet temperature, which selects the inverse.
    pub expander_out_temperature: Option<ThermodynamicTemperature>,
}

/// What the machine solved.
#[derive(Debug, Clone)]
pub struct TurboOutcome {
    /// The expander outlet.
    pub expander_outlet: Stream,
    /// The compressor outlet.
    pub compressor_outlet: Stream,
    /// The shaft speed, rpm.
    pub speed: f64,
    /// The expander's power, W.
    pub power_expander: f64,
    /// The compressor's absorbed power, W.
    pub power_compressor: f64,
    /// The expander's isentropic efficiency.
    pub expander_efficiency: f64,
    /// The expander's design efficiency, back-calculated in the inverse mode.
    pub expander_design_efficiency: f64,
    /// The compressor's polytropic efficiency.
    pub compressor_efficiency: f64,
    /// The compressor's polytropic head, m.
    pub compressor_head_m: f64,
    /// The expander's velocity ratio `U/C/design_uc`.
    pub uc_ratio: f64,
    /// The compressor's Q/N ratio.
    pub qn_ratio: f64,
    /// The IGV opening the machine reports.
    pub igv_opening: f64,
}

/// Run the coupled machine.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] for a non-positive design or geometry value, an
/// efficiency outside `(0, 1]`, a non-positive isentropic drop, or a speed that leaves the
/// class's bounds before the balance converges; and whatever the flash refuses.
pub fn turbo_expander(
    expander_feed: &Stream,
    compressor_feed: &Stream,
    settings: &TurboSettings<'_>,
) -> Result<TurboOutcome> {
    validate(settings)?;

    let (expander_mixture, expander_ideal_gas) = expander_feed.mixture()?;
    let expander_mass = expander_feed.molar_mass()?.value;
    let expander_mass_flow = expander_feed.n * expander_mass;
    let compressor_mass = compressor_feed.molar_mass()?.value;
    let compressor_mass_flow = compressor_feed.n * compressor_mass;
    let expander_state = feed_state(expander_feed)?;
    let compressor_state = feed_state(compressor_feed)?;
    let compressor_flow = compressor_feed.n * compressor_state.molar_volume;
    let expander_density = expander_mass / expander_state.molar_volume;

    // The expander's isentropic head, at the stated outlet pressure.
    let (entropy, _) = ps_flash::entropy_at(
        &expander_mixture,
        &expander_ideal_gas,
        expander_feed.t,
        expander_feed.p,
        &expander_feed.z,
    )?;
    let isentropic = ps_flash::ps_flash(
        &expander_mixture,
        &expander_ideal_gas,
        settings.expander_out_pressure,
        joules_per_mole_kelvin(entropy),
        &expander_feed.z,
    )?;
    let (h_isentropic, _) = ph_flash::enthalpy_at(
        &expander_mixture,
        &expander_ideal_gas,
        isentropic.temperature,
        settings.expander_out_pressure,
        &expander_feed.z,
    )?;
    let isentropic_drop = (expander_feed.h.value - h_isentropic) / expander_mass;
    if isentropic_drop.is_nan() || isentropic_drop <= 0.0 {
        return Err(AzothError::invalid_input(
            "expander_out_pressure",
            format!(
                "the expander's isentropic drop is {isentropic_drop} J/kg at {} Pa, and the \
                 machine needs a positive one",
                settings.expander_out_pressure.value
            ),
        ));
    }
    let spouting = (2.0 * isentropic_drop).sqrt();

    // The inverse's required efficiency, taken from the target temperature.
    let required = match settings.expander_out_temperature {
        Some(target) => {
            let (h_target, _) = ph_flash::enthalpy_at(
                &expander_mixture,
                &expander_ideal_gas,
                target,
                settings.expander_out_pressure,
                &expander_feed.z,
            )?;
            let actual = (expander_feed.h.value - h_target) / expander_mass;
            Some(actual / isentropic_drop)
        }
        None => None,
    };

    let parabola = parabola_coefficient(settings.uc_curve);
    let qn_efficiency = sorted(settings.qn_efficiency_curve);
    let qn_head = sorted(settings.qn_head_curve);
    let penalty = sorted(settings.igv_penalty_curve);

    let evaluate = |speed: f64| -> Pass {
        let tip = PI * settings.impeller_diameter * speed / 60.0;
        let uc_raw = tip / spouting;
        let uc_norm = uc_raw / settings.design_uc;
        let mut efficiency = match settings.chart {
            Some(chart) => chart.efficiency(uc_raw, settings.igv_opening),
            None => settings.expander_design_efficiency * efficiency_from_uc(parabola, uc_norm),
        };
        efficiency *= igv_penalty(&penalty, settings.igv_opening);
        if efficiency < 0.0 {
            efficiency = 0.0;
        }
        let total = if settings.expander_design_efficiency > 1.0e-9 {
            (efficiency / settings.expander_design_efficiency).max(0.0)
        } else {
            0.0
        };
        let efficiency = required.unwrap_or(efficiency);
        let power_expander = expander_mass_flow * isentropic_drop * efficiency;

        let qn_ratio = (compressor_flow * 60.0 / speed) / settings.design_qn;
        let efficiency_factor = spline(&qn_efficiency, qn_ratio).max(1.0e-6);
        let head_factor = spline(&qn_head, qn_ratio).max(1.0e-6);
        let design_head = settings.compressor_design_head_m * METRE_TO_JOULE_PER_KG;
        let head = design_head
            * (speed / settings.design_speed)
            * (speed / settings.design_speed)
            * head_factor;
        let compressor_efficiency = efficiency_factor * settings.compressor_design_efficiency;
        let power_compressor = compressor_mass_flow * head / compressor_efficiency;
        let bearing = settings.bearing_loss
            * (speed / settings.design_speed)
            * (speed / settings.design_speed);

        Pass {
            efficiency,
            total,
            power_expander,
            qn_ratio,
            head,
            compressor_efficiency,
            power_compressor,
            bearing,
            uc_norm,
        }
    };

    let mut speed = settings.design_speed;
    let mut delta = 10.0;
    let mut iteration = 0;
    let mut last: Option<Pass>;
    loop {
        let first = evaluate(speed);
        let second = evaluate(speed + delta);
        let residual = first.residual();
        let derivative = (second.residual() - residual) / delta;
        if derivative.abs() < 1.0e-8 {
            delta += 10.0;
            speed += 10.0;
        } else {
            let damping = (1.0 + iteration as f64) / (iteration as f64 + 5.0);
            speed -= damping * residual / derivative;
        }
        last = Some(first);
        if speed > MAX_SPEED {
            speed = MAX_SPEED;
            break;
        }
        if speed < MIN_SPEED {
            speed = MIN_SPEED;
            break;
        }
        iteration += 1;
        if !((residual.abs() * 100.0 > 1.0e-3 && iteration < MAX_ITERATIONS)
            || iteration < MIN_ITERATIONS)
        {
            break;
        }
    }

    let last = last.expect("the loop runs at least once");
    let design_efficiency = match required {
        Some(target) if last.total > 1.0e-9 => target / last.total,
        _ => settings.expander_design_efficiency,
    };

    let igv_opening = if settings.igv_control_mode {
        settings.igv_opening.clamp(0.0, 1.0)
    } else {
        evaluate_igv(
            isentropic_drop,
            expander_mass_flow,
            expander_density,
            settings.maximum_igv_area,
            settings.igv_area_increase_factor,
        )
    };

    let expander_outlet = expander_replay(
        expander_feed,
        settings.expander_out_pressure,
        last.efficiency,
    )?;
    let compressor_outlet = compressor_replay(
        compressor_feed,
        &compressor_state,
        last.power_compressor,
        last.compressor_efficiency,
    )?;

    Ok(TurboOutcome {
        expander_outlet,
        compressor_outlet,
        speed,
        power_expander: last.power_expander,
        power_compressor: last.power_compressor,
        expander_efficiency: last.efficiency,
        expander_design_efficiency: design_efficiency,
        compressor_efficiency: last.compressor_efficiency,
        compressor_head_m: last.head / METRE_TO_JOULE_PER_KG,
        uc_ratio: last.uc_norm,
        qn_ratio: last.qn_ratio,
        igv_opening,
    })
}

/// One pass's quantities.
struct Pass {
    efficiency: f64,
    total: f64,
    power_expander: f64,
    qn_ratio: f64,
    head: f64,
    compressor_efficiency: f64,
    power_compressor: f64,
    bearing: f64,
    uc_norm: f64,
}

impl Pass {
    fn residual(&self) -> f64 {
        self.power_expander - (self.power_compressor + self.bearing)
    }
}

fn validate(settings: &TurboSettings<'_>) -> Result<()> {
    for (name, value, positive) in [
        ("impeller_diameter", settings.impeller_diameter, true),
        ("design_speed", settings.design_speed, true),
        ("design_uc", settings.design_uc, true),
        ("design_qn", settings.design_qn, true),
        ("maximum_igv_area", settings.maximum_igv_area, true),
        (
            "igv_area_increase_factor",
            settings.igv_area_increase_factor,
            true,
        ),
        ("bearing_loss", settings.bearing_loss, false),
        (
            "compressor_design_head_m",
            settings.compressor_design_head_m,
            true,
        ),
    ] {
        if !value.is_finite() || (positive && value <= 0.0) || (!positive && value < 0.0) {
            return Err(AzothError::invalid_input(
                name,
                format!("the machine needs a valid {name}, and {value} is not"),
            ));
        }
    }
    for (name, value) in [
        (
            "expander_design_isentropic_efficiency",
            settings.expander_design_efficiency,
        ),
        (
            "compressor_design_polytropic_efficiency",
            settings.compressor_design_efficiency,
        ),
    ] {
        if !(value > 0.0 && value <= 1.0) {
            return Err(AzothError::invalid_input(
                name,
                format!("{name} must be in (0, 1], and {value} is not"),
            ));
        }
    }
    Ok(())
}

/// `setUCcurve`'s constrained parabola: `a (uc-1)² + 1`, the vertex fixed at one.
fn parabola_coefficient(curve: &[(f64, f64)]) -> f64 {
    if curve.len() < 2 {
        return 0.0;
    }
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for (uc, efficiency) in curve {
        let dx = uc - 1.0;
        let dx2 = dx * dx;
        numerator += (efficiency - 1.0) * dx2;
        denominator += dx2 * dx2;
    }
    if denominator != 0.0 {
        numerator / denominator
    } else {
        0.0
    }
}

/// `getEfficiencyFromUC`.
fn efficiency_from_uc(coefficient: f64, uc: f64) -> f64 {
    if coefficient != 0.0 {
        coefficient * (uc - 1.0) * (uc - 1.0) + 1.0
    } else {
        LEGACY_PARABOLA * (uc - 1.0) * (uc - 1.0) + 1.0
    }
}

fn sorted(curve: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut rows = curve.to_vec();
    rows.sort_by(|a, b| a.0.partial_cmp(&b.0).expect("finite curve"));
    rows
}

/// `getEfficiencyFromQN`/`getHeadFromQN`: piecewise cubic Hermite, linear extrapolation.
fn spline(curve: &[(f64, f64)], qn: f64) -> f64 {
    if curve.len() < 2 {
        return 1.0;
    }
    let n = curve.len();
    if qn <= curve[0].0 {
        let slope = (curve[1].1 - curve[0].1) / (curve[1].0 - curve[0].0);
        return curve[0].1 + slope * (qn - curve[0].0);
    }
    if qn >= curve[n - 1].0 {
        let slope = (curve[n - 1].1 - curve[n - 2].1) / (curve[n - 1].0 - curve[n - 2].0);
        return curve[n - 1].1 + slope * (qn - curve[n - 1].0);
    }
    let mut index = 0;
    for i in 0..n - 1 {
        if qn >= curve[i].0 && qn <= curve[i + 1].0 {
            index = i;
            break;
        }
    }
    let (x0, y0) = curve[index];
    let (x1, y1) = curve[index + 1];
    let h = x1 - x0;
    let t = (qn - x0) / h;
    let m0 = if index == 0 {
        (y1 - y0) / (x1 - x0)
    } else {
        (y1 - curve[index - 1].1) / (x1 - curve[index - 1].0)
    };
    let m1 = if index == n - 2 {
        (y1 - y0) / (x1 - x0)
    } else {
        (curve[index + 2].1 - y0) / (curve[index + 2].0 - x0)
    };
    let h00 = (1.0 + 2.0 * t) * (1.0 - t) * (1.0 - t);
    let h10 = t * (1.0 - t) * (1.0 - t);
    let h01 = t * t * (3.0 - 2.0 * t);
    let h11 = t * t * (t - 1.0);
    h00 * y0 + h10 * h * m0 + h01 * y1 + h11 * h * m1
}

/// `getIgvEfficiencyPenalty`: edge-clamped linear, unity without a curve.
fn igv_penalty(curve: &[(f64, f64)], igv: f64) -> f64 {
    if curve.len() < 2 {
        return 1.0;
    }
    let n = curve.len();
    if igv <= curve[0].0 {
        return curve[0].1;
    }
    if igv >= curve[n - 1].0 {
        return curve[n - 1].1;
    }
    for i in 0..n - 1 {
        if igv >= curve[i].0 && igv <= curve[i + 1].0 {
            let t = (igv - curve[i].0) / (curve[i + 1].0 - curve[i].0);
            return curve[i].1 + t * (curve[i + 1].1 - curve[i].1);
        }
    }
    curve[n - 1].1
}

/// `evaluateIGV`: the opening the stage drop and flow imply.
fn evaluate_igv(
    stage_drop: f64,
    mass_flow: f64,
    density: f64,
    maximum_area: f64,
    increase_factor: f64,
) -> f64 {
    if mass_flow <= 0.0 || stage_drop <= 0.0 || density <= 0.0 || mass_flow.is_nan() {
        return 0.0;
    }
    let mut half = 0.5 * stage_drop;
    if half <= 0.0 {
        half = (stage_drop * 0.5).max(1.0);
    }
    let velocity = (2.0 * half.max(1.0e-6)).sqrt();
    let required = mass_flow / (density * velocity);
    let expanded_area = maximum_area * increase_factor;
    let available = if required > maximum_area && expanded_area > maximum_area {
        expanded_area
    } else {
        maximum_area
    };
    if available > 0.0 {
        (required / available).clamp(0.0, 1.0)
    } else {
        0.0
    }
}

/// The compressor replay: the polytropic pressure ratio a stated power implies.
fn compressor_replay(
    feed: &Stream,
    state: &FeedState,
    power: f64,
    efficiency: f64,
) -> Result<Stream> {
    if !(efficiency > 0.0 && efficiency <= 1.0) {
        return Err(AzothError::invalid_input(
            "compressor_polytropic_efficiency",
            format!("the replay efficiency is {efficiency}, and it must be in (0, 1]"),
        ));
    }
    let molar_mass = feed.molar_mass()?.value;
    let z = state.z;
    let cp = state.cp;
    let kappa = cp / (cp - R_J_PER_MOL_K);
    let exponent = 1.0 / (1.0 - (kappa - 1.0) / kappa / efficiency);
    let dh = power / feed.n;
    let head = dh * efficiency / molar_mass;
    let flow_work = exponent / (exponent - 1.0) * z * R_J_PER_MOL_K * feed.t.value / molar_mass;
    if !flow_work.is_finite() || flow_work <= 0.0 || !exponent.is_finite() {
        return Err(AzothError::invalid_input(
            "power_compressor",
            format!("the replay has no polytropic ratio: n={exponent}, flow work={flow_work}"),
        ));
    }
    let pressure_ratio = ((head + flow_work) / flow_work).powf(exponent / (exponent - 1.0));
    Stream::from_ph(
        feed.components.clone(),
        feed.z.clone(),
        feed.n,
        pascals(feed.p.value * pressure_ratio),
        joules_per_mole(feed.h.value + dh),
    )
}

/// The feed's own state, as the class reads it: `getZ()`, `getCP()/n`, and the molar volume
/// behind `getFlowRate("m3/sec")` and `getDensity("kg/m3")`.
pub struct FeedState {
    /// The beta-weighted compressibility, `SystemThermo.getZ()`.
    pub z: f64,
    /// The beta-weighted molar heat capacity, `getCP()/n`.
    pub cp: f64,
    /// The beta-weighted molar volume with the Peneloux shift.
    pub molar_volume: f64,
}

/// Read the feed's state at its own composition.
///
/// **A single phase is evaluated at the feed's composition, not the incipient phase's.**
/// `pt_flash` on a superheated stream reports the *incipient* phase's root, and with a large
/// negative vapour fraction that composition sits far from the feed's; the class holds one phase
/// whose composition is the feed's. A genuinely two-phase feed is the beta-weighted mean of its
/// phases, which is what `SystemThermo.getZ()` and `getCP()` sum over.
fn feed_state(feed: &Stream) -> Result<FeedState> {
    let (mixture, ideal_gas) = feed.mixture()?;
    let flash = pt_flash(&mixture, feed.t, feed.p, &feed.z)?;
    let state = |composition: &[f64], root: f64| -> Result<(f64, f64, f64)> {
        let cp = molar_enthalpy_entropy(&mixture, &ideal_gas, feed.t, feed.p, composition, root)?
            .cp
            .value;
        let volume =
            root * R_J_PER_MOL_K * feed.t.value / feed.p.value - mixture.volume_shift(composition);
        Ok((root, cp, volume))
    };
    let (z, cp, molar_volume) = match flash.phase {
        Phase::AllVapour | Phase::AllLiquid | Phase::Trivial => {
            let side = if matches!(flash.phase, Phase::AllVapour) {
                RootSide::Vapour
            } else {
                RootSide::Liquid
            };
            let reduced = mixture.reduced_parameters(feed.t, feed.p)?;
            let root = mixture.phase_state(&reduced, &feed.z, side)?.z;
            state(&feed.z, root)?
        }
        Phase::TwoPhase => {
            let beta = flash.vapour_fraction.unwrap_or(1.0);
            let (_, cp_vapour, volume_vapour) = state(&flash.y, flash.z_vapour)?;
            let (_, cp_liquid, volume_liquid) = state(&flash.x, flash.z_liquid)?;
            (
                beta * flash.z_vapour + (1.0 - beta) * flash.z_liquid,
                beta * cp_vapour + (1.0 - beta) * cp_liquid,
                beta * volume_vapour + (1.0 - beta) * volume_liquid,
            )
        }
    };
    if molar_volume.is_nan() || molar_volume <= 0.0 {
        return Err(AzothError::invalid_input(
            "stream",
            format!(
                "the feed's molar volume is {molar_volume} m3/mol, and a state needs a positive one"
            ),
        ));
    }
    Ok(FeedState {
        z,
        cp,
        molar_volume,
    })
}
