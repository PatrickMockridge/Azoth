//! `process.compressor` - the Compressor's kernel as a registered id.
//!
//! Spec: `specs/models/process/compressor.toml`. The arithmetic is
//! [`crate::kernels::compressor`]; what is here is the boundary a case and a cross-impl
//! test address.
//!
//! **A delegation, because the route is `Pump.run`'s.** `Pump` takes it on `calculateAsCompressor` defaulting to true, and `Compressor.run`'s no-chart branch is the same three lines - so the kernel is `kernels::pump`'s and this adds nothing to it. **A stated curve selects the chart branch**, which is the compressor's own.

use azoth_core::units::{
    Length, Mass, MolarEnergy, Power, Pressure, ThermodynamicTemperature, Velocity, VolumeRate,
    cubic_meter_per_hour, cubic_meters_per_hour, joules_per_mole, kelvins, meter, meters, pascals,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use azoth_eos::Cubic;

use crate::compressor_chart::{ChartKind, CompressorChart, CurvePoint};
use crate::kernels::compressor as kernel;
use crate::kernels::compressor::{
    SpeedSolveSettings, SurgeSettings, compressor_polytropic, compressor_speed_solve,
    compressor_with_chart, compressor_with_power,
};
use crate::model_gen;
use crate::rotating_design::compressor::{
    CompressorDesignInputs, design as design_compressor, head_from_metres, head_from_power,
};
use crate::segment::phase::molar_mass_of;
use crate::stream::Stream;
use crate::surge_curve::SurgeCurve;
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

/// One row of `deposit`: an accumulated mass for a mechanism, in the spec's units.
#[derive(Debug, Clone, PartialEq)]
pub struct DepositMass {
    /// The mechanism's name, as the spec's enum spells it.
    pub mechanism: String,
    /// The accumulated mass.
    pub mass: Mass,
}

/// One row of `surge_curve`: a measured boundary point, in the units the spec declares.
#[derive(Debug, Clone, PartialEq)]
pub struct CompressorSurgePoint {
    /// The measured flow, m3/hr.
    pub flow: VolumeRate,
    /// The measured head, in metres - the same unit as `curve`'s.
    pub head: Length,
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
    /// `getPolytropicExponent`, present in the no-chart polytropic branch.
    pub polytropic_exponent: Option<f64>,
    /// The deposit inventory's summed mass, kg.
    pub total_deposit_mass: Option<f64>,
    /// `sum(mass / density)`, m3.
    pub deposit_volume: Option<f64>,
    /// `volume / wetted_area`, m.
    pub film_thickness: Option<Length>,
    /// `min(film / passage_half_height, 0.95)`.
    pub area_blockage_fraction: Option<f64>,
    /// The polytropic-efficiency multiplier.
    pub efficiency_multiplier: Option<f64>,
    /// The polytropic-head multiplier.
    pub head_multiplier: Option<f64>,
    /// `getDegradationFactor`: the efficiency multiplier the run recorded.
    pub degradation_factor: Option<f64>,
    /// `getFoulingFactor`: `1 - head_multiplier`.
    pub fouling_factor: Option<f64>,
    /// The shaft power, W - the recycled flow's where the anti-surge block ran.
    pub power: Power,
    /// The speed the chart was read at, rpm: the stated one or the solved one.
    pub shaft_speed: f64,
    /// `getSurgeFlowRate()`, m3/hr, where a surge curve is set.
    pub surge_flow_rate: Option<f64>,
    /// `getSurgeFlowRateMargin()`, m3/hr.
    pub surge_flow_rate_margin: Option<f64>,
    /// `getSurgeFlowRateStd()`, m3/hr at 1.01325 bara and 288.15 K.
    pub surge_flow_rate_std: Option<f64>,
    /// `getDistanceToSurge()`: the inlet flow over the surge flow less one.
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
    /// The stage count the impeller screen settled on.
    pub design_number_of_stages: Option<f64>,
    /// The impeller's outer diameter.
    pub design_impeller_diameter: Option<Length>,
    /// The shaft diameter.
    pub design_shaft_diameter: Option<Length>,
    /// The tip speed.
    pub design_tip_speed: Option<Velocity>,
    /// The inlet flow coefficient.
    pub design_flow_coefficient: Option<f64>,
    /// The speed the impeller screen used, rpm.
    pub design_impeller_sizing_speed: Option<f64>,
    /// The runner's inlet volume flow.
    pub design_impeller_sizing_flow: Option<VolumeRate>,
    /// The head the screen divided, kJ/kg.
    pub design_impeller_sizing_head: Option<f64>,
    /// Whether the screen's diameter and flow-coefficient limits are met.
    pub design_impeller_sizing_ok: Option<bool>,
    /// The driver power.
    pub design_driver_power: Option<Power>,
    /// The driver margin.
    pub design_driver_margin: Option<f64>,
    /// The casing design pressure.
    pub design_pressure: Option<Pressure>,
    /// The casing design temperature.
    pub design_temperature: Option<ThermodynamicTemperature>,
    /// The maximum continuous speed, rpm.
    pub design_max_continuous_speed: Option<f64>,
    /// The trip speed, rpm.
    pub design_trip_speed: Option<f64>,
    /// The first critical speed, rpm.
    pub design_first_critical_speed: Option<f64>,
    /// The casing type, by design pressure.
    pub design_casing_type: Option<String>,
    /// The head per stage, kJ/kg.
    pub design_head_per_stage: Option<f64>,
    /// The bearing span.
    pub design_bearing_span: Option<Length>,
    /// The casing weight.
    pub design_casing_weight: Option<Mass>,
    /// The rotor bundle weight.
    pub design_bundle_weight: Option<Mass>,
    /// The rotor weight.
    pub design_rotor_weight: Option<Mass>,
    /// The casing's outer diameter.
    pub design_outer_diameter: Option<Length>,
    /// The casing's inner diameter.
    pub design_inner_diameter: Option<Length>,
    /// The casing's wall thickness.
    pub design_wall_thickness: Option<Length>,
    /// The tangent-to-tangent length.
    pub design_tantan_length: Option<Length>,
    /// The base class's vessel-shell weight.
    pub design_weight_vessel_shell: Option<Mass>,
    /// The base class's internals weight.
    pub design_weight_internals: Option<Mass>,
    /// The base class's nozzle weight, zero by the class's comment.
    pub design_weight_nozzle: Option<Mass>,
    /// The piping weight.
    pub design_weight_piping: Option<Mass>,
    /// The electrical and instrument weight.
    pub design_weight_electro_instrument: Option<Mass>,
    /// The structural steel weight.
    pub design_weight_structural_steel: Option<Mass>,
    /// The total weight.
    pub design_weight_total: Option<Mass>,
    /// The module length.
    pub design_module_length: Option<Length>,
    /// The module width.
    pub design_module_width: Option<Length>,
    /// The module height.
    pub design_module_height: Option<Length>,
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
    use_polytropic_calc: Option<bool>,
    polytropic_efficiency: Option<f64>,
    polytropic_method: Option<&str>,
    compression_ratio: Option<f64>,
    max_outlet_pressure: Option<Pressure>,
    power: Option<Power>,
    deposit: Option<&[DepositMass]>,
    wetted_area: Option<Length>,
    passage_half_height: Option<Length>,
    curve: Option<&[CompressorCurvePoint]>,
    speed: Option<f64>,
    solve_speed: Option<bool>,
    min_speed: Option<f64>,
    max_speed: Option<f64>,
    limit_speed: Option<bool>,
    chart_type: Option<&str>,
    chart_conditions: Option<&[f64]>,
    impeller_diameter: Option<Length>,
    surge_curve: Option<&[CompressorSurgePoint]>,
    surge_control_margin: Option<f64>,
    anti_surge_active: Option<bool>,
    surge_control_factor: Option<f64>,
    mechanical_design: Option<bool>,
    max_stages_per_casing: Option<f64>,
) -> Result<CompressorResult> {
    let spec = &model_gen::COMPRESSOR_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "inlet_t" => Some(inlet_t.value),
            "isentropic_efficiency" => Some(isentropic_efficiency),
            "speed" => speed,
            "surge_control_margin" => surge_control_margin,
            _ => None,
        },
        &mut warnings,
    )?;

    let kind = chart_kind(chart_type)?;

    let feed = Stream::from_pt(
        components.to_vec(),
        inlet_z.to_vec(),
        inlet_n,
        inlet_p,
        inlet_t,
    )?;

    let deposit_state = match deposit {
        Some(rows) if !rows.is_empty() => {
            let pairs: Vec<(String, f64)> = rows
                .iter()
                .map(|row| (row.mechanism.clone(), row.mass.value))
                .collect();
            Some(kernel::deposit_state(
                &pairs,
                wetted_area.map(|value| value.value),
                passage_half_height.map(|value| value.value),
            )?)
        }
        _ => None,
    };

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

    let (surge_flow, surge_head): (Vec<f64>, Vec<f64>) = surge_curve
        .unwrap_or(&[])
        .iter()
        .map(|point| {
            (
                point.flow.get::<cubic_meter_per_hour>(),
                point.head.get::<meter>(),
            )
        })
        .unzip();
    let settings = SurgeSettings {
        active: anti_surge_active.unwrap_or(false),
        control_factor: surge_control_factor.unwrap_or(1.05),
        control_margin: surge_control_margin.unwrap_or(0.0),
    };

    let mut read_outs: [Option<f64>; 9] = [None; 9];
    let field_speed = speed.unwrap_or(3000.0);
    let (out, head, efficiency, power, shaft_speed, polytropic_exponent) = if points.is_empty() {
        if solve_speed.unwrap_or(false) {
            return Err(AzothError::invalid_input(
                "curve",
                "the speed-solve mode reads a chart, and this machine states none",
            ));
        }
        // `setCompressionRatio` flips the mode on and `setOutletPressure` overrides any stated
        // pressure; `setMaximumOutletPressure` caps it.
        let pressure = match compression_ratio {
            Some(ratio) => {
                if !(ratio.is_finite() && ratio > 0.0) {
                    return Err(AzothError::invalid_input(
                        "compression_ratio",
                        format!("a compression ratio is positive, and {ratio} is not"),
                    ));
                }
                let uncapped = feed.p.value * ratio;
                let capped = max_outlet_pressure.map_or(uncapped, |cap| uncapped.min(cap.value));
                pascals(capped)
            }
            None => outlet_pressure.ok_or_else(|| {
                AzothError::invalid_input(
                    "outlet_pressure",
                    "a compressor without a curve reaches only the pressure it is given, and \
                     this one has neither",
                )
            })?,
        };
        let (out, efficiency, exponent) = if let Some(power) = power {
            // `setPower` fixes the shaft work; the polytropic branch keeps its stated
            // efficiency and the isentropic one back-solves it from the isentropic step.
            let stated = if use_polytropic_calc.unwrap_or(false) {
                // `applyDepositDegradation` scales the clean polytropic efficiency before
                // the branch consumes it.
                Some(
                    polytropic_efficiency.unwrap_or(1.0)
                        * deposit_state.map_or(1.0, |state| state.efficiency_multiplier),
                )
            } else {
                None
            };
            let (out, efficiency) = compressor_with_power(&feed, pressure, power, stated)?;
            (out, efficiency, None)
        } else if use_polytropic_calc.unwrap_or(false) {
            let method = polytropic_method.unwrap_or("default");
            match method {
                "default" => {}
                "detailed" => return Err(unported::refuse("polytropic_method=detailed")),
                "schultz" => return Err(unported::refuse("polytropic_method=schultz")),
                other => {
                    return Err(AzothError::invalid_input(
                        "polytropic_method",
                        format!("`{other}` is not one of default, detailed or schultz"),
                    ));
                }
            }
            let clean = polytropic_efficiency.unwrap_or(1.0);
            let degraded = clean * deposit_state.map_or(1.0, |state| state.efficiency_multiplier);
            let (out, exponent) = compressor_polytropic(&feed, pressure, degraded)?;
            (out, degraded, Some(exponent))
        } else {
            (kernel(&feed, pressure, isentropic_efficiency)?, 0.0, None)
        };
        let power = azoth_core::units::watts(feed.n * (out.h.value - feed.h.value));
        (out, 0.0, efficiency, power, field_speed, exponent)
    } else {
        let mut chart = build_chart(
            components,
            inlet_t,
            inlet_p,
            inlet_z,
            &points,
            kind,
            chart_conditions,
            impeller_diameter,
        )?;
        if !surge_flow.is_empty() {
            chart.set_surge_curve(SurgeCurve::new(&surge_flow, &surge_head)?);
        }
        let outcome = if solve_speed.unwrap_or(false) {
            let target = outlet_pressure.ok_or_else(|| {
                AzothError::invalid_input(
                    "outlet_pressure",
                    "the speed-solve mode drives the chart to a stated pressure, and this \
                     machine states none",
                )
            })?;
            compressor_speed_solve(
                &feed,
                &chart,
                &SpeedSolveSettings {
                    target,
                    initial_speed: field_speed,
                    min_speed: min_speed.unwrap_or(0.0),
                    max_speed: max_speed.unwrap_or(30000.0),
                    limit_speed: limit_speed.unwrap_or(false),
                },
                &settings,
            )?
        } else {
            compressor_with_chart(&feed, &chart, field_speed, &settings)?
        };
        // The kernel works in m3/hr, the class's chart unit; the record crosses in SI.
        read_outs = [
            outcome.surge_flow_rate.map(|value| value / 3600.0),
            outcome.surge_flow_rate_margin.map(|value| value / 3600.0),
            outcome.surge_flow_rate_std.map(|value| value / 3600.0),
            outcome.distance_to_surge,
            outcome.control_line_flow.map(|value| value / 3600.0),
            outcome.distance_to_control_line,
            outcome.required_recycle_fraction,
            outcome.anti_surge_recycle_power,
            outcome.recycle_fraction,
        ];
        (
            outcome.outlet,
            outcome.polytropic_head,
            outcome.polytropic_efficiency,
            outcome.power,
            outcome.speed,
            None,
        )
    };

    let design = if mechanical_design.unwrap_or(false) {
        let (mixture, _) = feed.mixture()?;
        let molar_mass = molar_mass_of(&mixture, &feed.z)?;
        let mass_flow_kg_s = feed.n * molar_mass;
        let polytropic_fluid_head = if !points.is_empty() {
            head_from_metres(head)
        } else if use_polytropic_calc.unwrap_or(false) {
            head_from_power(power.value, mass_flow_kg_s, efficiency)
        } else {
            0.0
        };
        if !(polytropic_fluid_head.is_finite() && polytropic_fluid_head > 0.0) {
            return Err(AzothError::invalid_input(
                "mechanical_design",
                "the class's polytropic fluid head is zero on this branch, and `calcDesign` \
                 returns before it sizes anything",
            ));
        }
        Some(design_compressor(
            &feed,
            &out,
            power,
            polytropic_fluid_head,
            shaft_speed,
            &CompressorDesignInputs {
                max_stages_per_casing: max_stages_per_casing.unwrap_or(10.0) as u32,
            },
        )?)
    } else {
        None
    };
    let d = design.as_ref();

    Ok(CompressorResult {
        outlet_n: out.n,
        outlet_z: out.z,
        outlet_p: out.p,
        outlet_t: out.t,
        outlet_h: joules_per_mole(out.h.value),
        polytropic_head: meters(head),
        polytropic_efficiency: efficiency,
        polytropic_exponent,
        total_deposit_mass: deposit_state.map(|state| state.total_mass),
        deposit_volume: deposit_state.map(|state| state.volume),
        film_thickness: deposit_state.map(|state| meters(state.film_thickness)),
        area_blockage_fraction: deposit_state.map(|state| state.blockage),
        efficiency_multiplier: deposit_state.map(|state| state.efficiency_multiplier),
        head_multiplier: deposit_state.map(|state| state.head_multiplier),
        degradation_factor: deposit_state.map(|state| state.efficiency_multiplier),
        fouling_factor: deposit_state.map(|state| 1.0 - state.head_multiplier),
        power,
        shaft_speed,
        surge_flow_rate: read_outs[0],
        surge_flow_rate_margin: read_outs[1],
        surge_flow_rate_std: read_outs[2],
        distance_to_surge: read_outs[3],
        control_line_flow: read_outs[4],
        distance_to_control_line: read_outs[5],
        required_recycle_fraction: read_outs[6],
        anti_surge_recycle_power: read_outs[7],
        recycle_fraction: read_outs[8],
        design_number_of_stages: d.map(|value| f64::from(value.number_of_stages)),
        design_impeller_diameter: d.map(|value| value.impeller_diameter),
        design_shaft_diameter: d.map(|value| value.shaft_diameter),
        design_tip_speed: d.map(|value| value.tip_speed),
        design_flow_coefficient: d.map(|value| value.flow_coefficient),
        design_impeller_sizing_speed: d.map(|value| value.impeller_sizing_speed),
        design_impeller_sizing_flow: d
            .map(|value| cubic_meters_per_hour(value.impeller_sizing_flow)),
        design_impeller_sizing_head: d.map(|value| value.impeller_sizing_head),
        design_impeller_sizing_ok: d.map(|value| value.impeller_sizing_ok),
        design_driver_power: d.map(|value| value.driver_power),
        design_driver_margin: d.map(|value| value.driver_margin),
        design_pressure: d.map(|value| value.design_pressure),
        design_temperature: d.map(|value| value.design_temperature),
        design_max_continuous_speed: d.map(|value| value.max_continuous_speed),
        design_trip_speed: d.map(|value| value.trip_speed),
        design_first_critical_speed: d.map(|value| value.first_critical_speed),
        design_casing_type: d.map(|value| value.casing_type.as_str().to_owned()),
        design_head_per_stage: d.map(|value| value.head_per_stage),
        design_bearing_span: d.map(|value| value.bearing_span),
        design_casing_weight: d.map(|value| value.casing_weight),
        design_bundle_weight: d.map(|value| value.bundle_weight),
        design_rotor_weight: d.map(|value| value.rotor_weight),
        design_outer_diameter: d.map(|value| value.outer_diameter),
        design_inner_diameter: d.map(|value| value.inner_diameter),
        design_wall_thickness: d.map(|value| value.wall_thickness),
        design_tantan_length: d.map(|value| value.tantan_length),
        design_weight_vessel_shell: d.map(|value| value.weight_vessel_shell),
        design_weight_internals: d.map(|value| value.weight_internals),
        design_weight_nozzle: d.map(|value| value.weight_nozzle),
        design_weight_piping: d.map(|value| value.weight_piping),
        design_weight_electro_instrument: d.map(|value| value.weight_electro_instrument),
        design_weight_structural_steel: d.map(|value| value.weight_structural_steel),
        design_weight_total: d.map(|value| value.weight_total),
        design_module_length: d.map(|value| value.module_length),
        design_module_width: d.map(|value| value.module_width),
        design_module_height: d.map(|value| value.module_height),
        warnings,
    })
}

/// `createDefaultFluid`'s phase sound speed: methane/ethane/propane `0.90/0.05/0.05`, rescaled
/// to the chart conditions' molar mass, flashed at their temperature and pressure.
///
/// The class instantiates the *same EOS class* by reflection; azoth's cubic is the one this
/// model runs, so the reference fluid is built the same way and the speed comes from
/// [`azoth_eos::sound_speed`].
pub(crate) fn reference_sound_speed(conditions: &[f64]) -> Result<f64> {
    if conditions.len() != 4 {
        return Err(AzothError::invalid_input(
            "chart_conditions",
            format!(
                "the chart conditions are four numbers (T C, P bara, density, MW g/mol), and {} \
                 were given",
                conditions.len()
            ),
        ));
    }
    let names = ["methane", "ethane", "propane"];
    let (mixture, _) = azoth_eos::databank::mixture_of(&names, Cubic::Pr, None)?;
    let masses: Vec<f64> = mixture
        .components()
        .iter()
        .map(|component| component.molar_mass.unwrap_or(0.0) * 1000.0)
        .collect();
    let (x1, x2, x3) = (0.90_f64, 0.05_f64, 0.05_f64);
    let z = if conditions[3] > 0.0 {
        let a = (conditions[3] - masses[2])
            / (x1 * (masses[0] - masses[2]) + x2 * (masses[1] - masses[2]));
        let first = a * x1;
        let second = a * x2;
        vec![first, second, 1.0 - first - second]
    } else {
        vec![x1, x2, x3]
    };
    let components: Vec<String> = names.iter().map(|name| (*name).to_string()).collect();
    Ok(azoth_eos::sound_speed(
        &components,
        kelvins(conditions[0] + 273.15),
        pascals(conditions[1] * 1.0e5),
        &z,
    )?
    .sound_speed
    .value)
}

/// The `chart_type` input, as the chart's kind.
pub(crate) fn chart_kind(chart_type: Option<&str>) -> Result<ChartKind> {
    match chart_type.unwrap_or("fan_law") {
        "fan_law" => Ok(ChartKind::FanLaw),
        "interpolate" => Ok(ChartKind::MapLookup),
        "interpolate_and_extrapolate" => Ok(ChartKind::MapLookupExtrapolate),
        "khader_2015" => Ok(ChartKind::Khader2015),
        other => Err(AzothError::invalid_input(
            "chart_type",
            format!("'{other}' is not one of the four declared chart types"),
        )),
    }
}

/// The machine's chart, built from its curve and the kind.
///
/// Shared with `process.map_turbo_expander_compressor`, which runs the same map on its
/// compressor side.
#[allow(clippy::too_many_arguments)] // one parameter per declared input
pub(crate) fn build_chart(
    components: &[String],
    inlet_t: ThermodynamicTemperature,
    inlet_p: Pressure,
    inlet_z: &[f64],
    points: &[CurvePoint],
    kind: ChartKind,
    chart_conditions: Option<&[f64]>,
    impeller_diameter: Option<Length>,
) -> Result<CompressorChart> {
    if kind != ChartKind::Khader2015 {
        return CompressorChart::with_kind(points, kind);
    }
    let conditions = chart_conditions.ok_or_else(|| {
        AzothError::invalid_input(
            "chart_conditions",
            "the Khader chart builds its reference fluid at the chart conditions, and this \
             machine states none",
        )
    })?;
    let reference = reference_sound_speed(conditions)?;
    let process = azoth_eos::sound_speed(components, inlet_t, inlet_p, inlet_z)?
        .system_sound_speed
        .value;
    CompressorChart::with_khader(
        points,
        reference,
        process,
        impeller_diameter.map_or(0.3, |diameter| diameter.value),
    )
}
