//! `process.heat_exchanger` - the exchanger's kernel as a registered id.
//!
//! Spec: `specs/models/process/heat_exchanger.toml`. The arithmetic is
//! [`crate::kernels::heat_exchanger`]; what is here is the boundary a case and a
//! cross-impl test address.
//!
//! **The first id whose two ports carry different fluids.** Every other unit operation
//! declares one `components`; this declares `hot_components` and `cold_components`, and
//! the boundary derives the two `Stream`s from them separately.

use azoth_core::units::{
    Area, DynamicViscosity, HeatTransfer, Length, Mass, MolarEnergy, Power, Pressure,
    ThermalConductance, ThermalConductivity, ThermodynamicTemperature, Velocity, joules_per_mole,
    kelvins, kilograms, meters, meters_per_second, pascals, square_meters, watts_per_kelvin,
    watts_per_square_meter_kelvin,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{optional_scalar, scalar, warnings as wire_warnings};
use crate::kernels::heat_exchanger as kernel;
use crate::kernels::heater_mechanical;
use crate::model_gen;
use crate::stream::Stream;

/// One row of the mechanical design's candidate table.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct HeatExchangerMechanicalType {
    /// The configuration's spelling.
    pub mechanical_type: String,
    /// The required area.
    #[serde(serialize_with = "scalar")]
    pub required_area: Area,
    /// The UA the area is computed from.
    #[serde(serialize_with = "scalar")]
    pub required_ua: ThermalConductance,
    /// The type's typical overall coefficient.
    #[serde(serialize_with = "scalar")]
    pub overall_heat_transfer_coefficient: HeatTransfer,
    /// The type's approach.
    #[serde(serialize_with = "scalar")]
    pub approach_temperature: ThermodynamicTemperature,
    /// The estimated length.
    #[serde(serialize_with = "scalar")]
    pub estimated_length: Length,
    /// The equivalent inside diameter.
    #[serde(serialize_with = "scalar")]
    pub inner_diameter: Length,
    /// The equivalent outside diameter.
    #[serde(serialize_with = "scalar")]
    pub outer_diameter: Length,
    /// The wall or sheet thickness.
    #[serde(serialize_with = "scalar")]
    pub wall_thickness: Length,
    /// The tube or plate count.
    pub tube_count: f64,
    /// The tube passes or plate pairs.
    pub tube_passes: f64,
    /// The fin surface area; zero for the unfinned types.
    #[serde(serialize_with = "scalar")]
    pub fin_surface_area: Area,
    /// The class's qualitative pressure-drop proxy, bar.
    pub estimated_pressure_drop: f64,
    /// The estimated weight.
    #[serde(serialize_with = "scalar")]
    pub estimated_weight: Mass,
    /// The module length.
    #[serde(serialize_with = "scalar")]
    pub module_length: Length,
    /// The module width.
    #[serde(serialize_with = "scalar")]
    pub module_width: Length,
    /// The module height.
    #[serde(serialize_with = "scalar")]
    pub module_height: Length,
}

/// Result of `process.heat_exchanger`.
///
/// Two single-multiplicity ports, so the record's five fields per port - the same shape
/// `process.separator` has, under this unit operation's own port names - and the duty and the
/// rating's own numbers beside them.
///
/// **Serialisable, so that a flowsheet's run publishes it.** The duty is on no outlet stream, and
/// neither are the four numbers a rating sized the exchanger by; the wire writes each through
/// [`crate::executor::json::scalar`] and writes `null` for the rating's where there was no rating
/// to size one, which is the same absence a flash's single phase is.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct HeatExchangerResult {
    /// Hot outlet molar flow, mol/s.
    pub hot_out_n: f64,
    /// Hot outlet composition.
    pub hot_out_z: Vec<f64>,
    #[serde(serialize_with = "scalar")]
    /// Hot outlet pressure.
    pub hot_out_p: Pressure,
    #[serde(serialize_with = "scalar")]
    /// Hot outlet temperature.
    pub hot_out_t: ThermodynamicTemperature,
    #[serde(serialize_with = "scalar")]
    /// Hot outlet molar enthalpy.
    pub hot_out_h: MolarEnergy,
    /// Cold outlet molar flow, mol/s.
    pub cold_out_n: f64,
    /// Cold outlet composition.
    pub cold_out_z: Vec<f64>,
    #[serde(serialize_with = "scalar")]
    /// Cold outlet pressure.
    pub cold_out_p: Pressure,
    #[serde(serialize_with = "scalar")]
    /// Cold outlet temperature.
    pub cold_out_t: ThermodynamicTemperature,
    #[serde(serialize_with = "scalar")]
    /// Cold outlet molar enthalpy.
    pub cold_out_h: MolarEnergy,
    /// The heat the hot side released, W: positive when the hot side cools, which is what the
    /// cold side's gain balances. The one figure that means the same thing in both modes.
    #[serde(serialize_with = "scalar")]
    pub duty: Power,
    /// `UA / C_min`, dimensionless, or `null` where one outlet was pinned instead of rated.
    pub ntu: Option<f64>,
    /// The swing the relation scaled, or `null` where one outlet was pinned.
    pub effectiveness: Option<f64>,
    /// The smaller of the two estimated capacities, W/K, or `null` where one outlet was pinned.
    #[serde(serialize_with = "optional_scalar")]
    pub c_min: Option<ThermalConductance>,
    /// The larger of them, W/K, or `null` where one outlet was pinned.
    #[serde(serialize_with = "optional_scalar")]
    pub c_max: Option<ThermalConductance>,
    /// `C_min / C_max`, dimensionless, or `null` where one outlet was pinned.
    pub capacity_ratio: Option<f64>,
    /// `runDeltaT`'s computed overall conductance, W/K; present only in that mode.
    pub ua_value: Option<f64>,
    /// Which side `runDeltaT` set, `"hot"` or `"cold"`; present only in that mode.
    pub delta_t_side: Option<String>,
    /// The selected mechanical configuration.
    pub mechanical_selected_type: Option<String>,
    /// The selected configuration's required area.
    #[serde(serialize_with = "optional_scalar")]
    pub mechanical_selected_area: Option<Area>,
    /// The log-mean temperature difference.
    #[serde(serialize_with = "optional_scalar")]
    pub log_mean_temperature_difference: Option<ThermodynamicTemperature>,
    /// The approach temperature.
    #[serde(serialize_with = "optional_scalar")]
    pub mechanical_approach_temperature: Option<ThermodynamicTemperature>,
    /// The overall coefficient the estimate used.
    #[serde(serialize_with = "optional_scalar")]
    pub used_overall_heat_transfer_coefficient: Option<HeatTransfer>,
    /// The calculated UA.
    #[serde(serialize_with = "optional_scalar")]
    pub calculated_ua: Option<ThermalConductance>,
    /// The rating's correlated overall coefficient.
    #[serde(serialize_with = "optional_scalar")]
    pub rating_u: Option<HeatTransfer>,
    /// `rating_u * rating_area`.
    #[serde(serialize_with = "optional_scalar")]
    pub rating_ua: Option<ThermalConductance>,
    /// The outside area the rating used.
    #[serde(serialize_with = "optional_scalar")]
    pub rating_area: Option<Area>,
    /// The tube side's correlated pressure drop.
    #[serde(serialize_with = "optional_scalar")]
    pub tube_side_pressure_drop: Option<Pressure>,
    /// The shell side's correlated pressure drop.
    #[serde(serialize_with = "optional_scalar")]
    pub shell_side_pressure_drop: Option<Pressure>,
    /// The tube-side film coefficient.
    #[serde(serialize_with = "optional_scalar")]
    pub tube_side_heat_transfer_coefficient: Option<HeatTransfer>,
    /// The shell-side film coefficient.
    #[serde(serialize_with = "optional_scalar")]
    pub shell_side_heat_transfer_coefficient: Option<HeatTransfer>,
    /// The tube-side Reynolds number.
    pub tube_side_reynolds: Option<f64>,
    /// The shell-side Reynolds number.
    pub shell_side_reynolds: Option<f64>,
    /// The tube-side velocity.
    #[serde(serialize_with = "optional_scalar")]
    pub tube_side_velocity: Option<Velocity>,
    /// The shell-side crossflow velocity.
    #[serde(serialize_with = "optional_scalar")]
    pub shell_side_velocity: Option<Velocity>,
    /// The baffle-cut correction factor.
    pub jc: Option<f64>,
    /// The baffle-leakage correction factor.
    pub jl: Option<f64>,
    /// The bundle-bypass correction factor.
    pub jb: Option<f64>,
    /// The unequal-spacing correction factor.
    pub js: Option<f64>,
    /// The adverse-temperature-gradient correction factor.
    pub jr: Option<f64>,
    /// The five candidates, in the class's enum order.
    pub mechanical_types: Vec<HeatExchangerMechanicalType>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

/// Exchange heat between two streams.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] for the mode and arrangement refusals
/// [`crate::kernels::heat_exchanger`] lists, and whatever the databank refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are seventeen
pub fn heat_exchanger(
    hot_components: &[String],
    hot_in_n: f64,
    hot_in_z: &[f64],
    hot_in_p: Pressure,
    hot_in_t: ThermodynamicTemperature,
    cold_components: &[String],
    cold_in_n: f64,
    cold_in_z: &[f64],
    cold_in_p: Pressure,
    cold_in_t: ThermodynamicTemperature,
    ua: Option<ThermalConductance>,
    flow_arrangement: &str,
    hot_outlet_temperature: Option<ThermodynamicTemperature>,
    cold_outlet_temperature: Option<ThermodynamicTemperature>,
    use_delta_t: Option<bool>,
    delta_t: Option<ThermodynamicTemperature>,
    guess_out_temperature: Option<ThermodynamicTemperature>,
    mechanical_design: Option<bool>,
    mechanical_type: Option<&str>,
    mechanical_selection: Option<&str>,
    use_rating_pressure_drop: Option<bool>,
    rating_area: Option<Area>,
    rating_tube_side: Option<&str>,
    rating_tube_outside_diameter: Option<Length>,
    rating_tube_inside_diameter: Option<Length>,
    rating_tube_length: Option<Length>,
    rating_tube_count: Option<f64>,
    rating_tube_passes: Option<f64>,
    rating_tube_pitch: Option<Length>,
    rating_triangular_pitch: Option<bool>,
    rating_shell_inside_diameter: Option<Length>,
    rating_baffle_spacing: Option<Length>,
    rating_baffle_count: Option<f64>,
    rating_baffle_cut: Option<f64>,
    rating_tube_wall_conductivity: Option<ThermalConductivity>,
    rating_fouling_tube: Option<f64>,
    rating_fouling_shell: Option<f64>,
    rating_tube_fouling_layer_thickness: Option<Length>,
    rating_tube_fouling_layer_conductivity: Option<ThermalConductivity>,
    rating_shell_side_method: Option<&str>,
    rating_shell_viscosity_wall: Option<DynamicViscosity>,
    rating_has_sealing: Option<bool>,
    rating_sealing_pairs: Option<f64>,
    rating_tube_to_baffle_clearance: Option<Length>,
    rating_shell_to_baffle_clearance: Option<Length>,
    rating_bypass_area: Option<Area>,
) -> Result<HeatExchangerResult> {
    let spec = &model_gen::HEAT_EXCHANGER_SPEC;
    let mut warnings = Vec::new();
    // One call over every declared input, for the reason `reactions.chemical_equilibrium`
    // gives. `ua` is optional, so an absent one is a *skipped check* rather than a pass -
    // which is the whole reason the bound is declared on it and not on the temperatures.
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "ua" => ua.map(|value| value.value),
            "hot_in_t" => Some(hot_in_t.value),
            "cold_in_t" => Some(cold_in_t.value),
            "rating_area" => rating_area.map(|value| value.value),
            "rating_tube_outside_diameter" => rating_tube_outside_diameter.map(|value| value.value),
            "rating_tube_inside_diameter" => rating_tube_inside_diameter.map(|value| value.value),
            "rating_tube_length" => rating_tube_length.map(|value| value.value),
            "rating_tube_count" => rating_tube_count,
            "rating_tube_passes" => rating_tube_passes,
            "rating_tube_pitch" => rating_tube_pitch.map(|value| value.value),
            "rating_shell_inside_diameter" => rating_shell_inside_diameter.map(|value| value.value),
            "rating_baffle_spacing" => rating_baffle_spacing.map(|value| value.value),
            "rating_baffle_count" => rating_baffle_count,
            "rating_baffle_cut" => rating_baffle_cut,
            "rating_tube_wall_conductivity" => {
                rating_tube_wall_conductivity.map(|value| value.value)
            }
            "rating_fouling_tube" => rating_fouling_tube,
            "rating_fouling_shell" => rating_fouling_shell,
            "rating_tube_fouling_layer_thickness" => {
                rating_tube_fouling_layer_thickness.map(|value| value.value)
            }
            "rating_tube_fouling_layer_conductivity" => {
                rating_tube_fouling_layer_conductivity.map(|value| value.value)
            }
            "rating_sealing_pairs" => rating_sealing_pairs,
            "rating_tube_to_baffle_clearance" => {
                rating_tube_to_baffle_clearance.map(|value| value.value)
            }
            "rating_shell_to_baffle_clearance" => {
                rating_shell_to_baffle_clearance.map(|value| value.value)
            }
            "rating_bypass_area" => rating_bypass_area.map(|value| value.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let mechanical_requested = mechanical_design.unwrap_or(false);
    for (name, stated) in [
        ("mechanical_type", mechanical_type.is_some()),
        ("mechanical_selection", mechanical_selection.is_some()),
    ] {
        if stated && !mechanical_requested {
            return Err(AzothError::invalid_input(
                name,
                "a mechanical-design input is stated without `mechanical_design`",
            ));
        }
    }
    let manual_selection = match mechanical_type {
        Some(spelling) => Some(
            heater_mechanical::MechanicalType::parse(spelling).ok_or_else(|| {
                AzothError::invalid_input(
                    "mechanical_type",
                    format!("{spelling} is not one of the five mechanical configurations"),
                )
            })?,
        ),
        None => None,
    };
    let selection_criterion = match mechanical_selection {
        Some(spelling) => Some(
            heater_mechanical::SelectionCriterion::parse(spelling).ok_or_else(|| {
                AzothError::invalid_input(
                    "mechanical_selection",
                    format!("{spelling} is not one of min_area, min_weight or min_pressure_drop"),
                )
            })?,
        ),
        None => None,
    };
    let mechanical_request = heater_mechanical::MechanicalRequest {
        manual_selection,
        selection_criterion,
    };

    // **A zero flow is not a state**, and every duty-per-flow divides by one.
    for (field, n) in [("hot_in_n", hot_in_n), ("cold_in_n", cold_in_n)] {
        if n <= 0.0 {
            return Err(azoth_core::AzothError::out_of_range(
                field,
                n,
                "an exchanger side needs a positive molar flow",
            ));
        }
    }

    let hot = Stream::from_pt(
        hot_components.to_vec(),
        hot_in_z.to_vec(),
        hot_in_n,
        hot_in_p,
        hot_in_t,
    )?;
    let cold = Stream::from_pt(
        cold_components.to_vec(),
        cold_in_z.to_vec(),
        cold_in_n,
        cold_in_p,
        cold_in_t,
    )?;
    // `runDeltaT`'s `guessOutTemperature` only seeds the pass the `else` branch then
    // overwrites, so it is carried and unused. The declaration keeps the surface whole.
    let _ = guess_out_temperature;
    if use_delta_t.unwrap_or(false)
        && (hot_outlet_temperature.is_some() || cold_outlet_temperature.is_some())
    {
        warnings.push(Warning::for_field(
            azoth_core::WarningCode::FallbackToInput,
            "use_delta_t",
            "a pinned outlet wins in `HeatExchanger.run`: `runDeltaT` delegates to \
             `runSpecifiedStream`, so `use_delta_t` is ignored on this state",
        ));
    }
    // The class's own `UAvalue` basis: the stated one, or its 500 default, which the setter
    // rather than `run` writes - so a pinned run that states none leaves 500 here.
    let class_ua = ua.map_or(500.0, |value| value.value);
    let pinned = hot_outlet_temperature.is_some() || cold_outlet_temperature.is_some();
    let outcome = if use_rating_pressure_drop.unwrap_or(false) {
        if pinned || use_delta_t.unwrap_or(false) {
            return Err(AzothError::invalid_input(
                "use_rating_pressure_drop",
                "the thermal-hydraulic rating requires a UA specification without a pinned \
                 outlet or delta-T",
            ));
        }
        let area = rating_area.ok_or_else(|| {
            AzothError::invalid_input(
                "rating_area",
                "the thermal-hydraulic rating needs a positive outside area",
            )
        })?;
        let tube_is_hot = match rating_tube_side {
            Some("hot") => true,
            Some("cold") => false,
            Some(other) => {
                return Err(AzothError::invalid_input(
                    "rating_tube_side",
                    format!("{other} is not one of hot or cold"),
                ));
            }
            None => {
                return Err(AzothError::invalid_input(
                    "rating_tube_side",
                    "the thermal-hydraulic rating needs the side that is the tube",
                ));
            }
        };
        let shell_side_method = match rating_shell_side_method {
            Some(spelling) => Some(spelling.parse::<azoth_thermal::ShellSideMethod>()?),
            None => None,
        };
        let geometry = kernel::RatingGeometry {
            tube_outside_diameter: rating_tube_outside_diameter.map(|value| value.value),
            tube_inside_diameter: rating_tube_inside_diameter.map(|value| value.value),
            tube_length: rating_tube_length.map(|value| value.value),
            tube_count: rating_tube_count,
            tube_passes: rating_tube_passes,
            tube_pitch: rating_tube_pitch.map(|value| value.value),
            triangular_pitch: rating_triangular_pitch,
            shell_inside_diameter: rating_shell_inside_diameter.map(|value| value.value),
            baffle_spacing: rating_baffle_spacing.map(|value| value.value),
            baffle_count: rating_baffle_count,
            baffle_cut: rating_baffle_cut,
            tube_wall_conductivity: rating_tube_wall_conductivity.map(|value| value.value),
            fouling_tube: rating_fouling_tube,
            fouling_shell: rating_fouling_shell,
            tube_fouling_layer_thickness: rating_tube_fouling_layer_thickness
                .map(|value| value.value),
            tube_fouling_layer_conductivity: rating_tube_fouling_layer_conductivity
                .map(|value| value.value),
            shell_side_method,
            shell_viscosity_wall: rating_shell_viscosity_wall.map(|value| value.value),
            has_sealing: rating_has_sealing,
            sealing_pairs: rating_sealing_pairs,
            tube_to_baffle_clearance: rating_tube_to_baffle_clearance.map(|value| value.value),
            shell_to_baffle_clearance: rating_shell_to_baffle_clearance.map(|value| value.value),
            bypass_area: rating_bypass_area.map(|value| value.value),
        };
        kernel::heat_exchanger_rating(
            &hot,
            &cold,
            tube_is_hot,
            area.value,
            &geometry,
            flow_arrangement,
        )?
    } else {
        kernel(
            &hot,
            &cold,
            ua,
            flow_arrangement,
            hot_outlet_temperature,
            cold_outlet_temperature,
            use_delta_t.unwrap_or(false),
            delta_t.unwrap_or_else(|| azoth_core::units::kelvins(1.0)),
        )?
    };
    let mechanical = if mechanical_requested {
        // **The class's `getDuty()` basis, not this result's `duty`**: the pinned branch
        // writes neither, so `determineDuty` sees zero there.
        let duty_basis = if pinned {
            0.0
        } else {
            outcome.duty.value.abs()
        };
        let effectiveness = outcome
            .rating
            .as_ref()
            .map_or(0.0, |rating| rating.effectiveness);
        Some(heater_mechanical::two_stream_design(
            hot_in_t.value,
            outcome.hot_out.t.value,
            cold_in_t.value,
            outcome.cold_out.t.value,
            duty_basis,
            outcome
                .ua_value
                .or_else(|| {
                    outcome
                        .rating_diagnostics
                        .as_ref()
                        .map(|value| value.rating_ua)
                })
                .unwrap_or(class_ua),
            effectiveness,
            &mechanical_request,
        ))
    } else {
        None
    };
    Ok(HeatExchangerResult::of(
        &outcome,
        mechanical.as_ref(),
        warnings,
    ))
}

impl HeatExchangerResult {
    /// The result of one kernel call, as a flowsheet's dispatcher reads it.
    ///
    /// **The warnings are the caller's, because the two callers have different ones.** A case runs
    /// the spec's input checks through [`apply_checks`]; a flowsheet's equivalent is the
    /// *checker's*, which reports through the envelope's diagnostics rather than through a result.
    /// So an empty list here means "this kernel raised none" and not "nothing checked this call" -
    /// and this kernel raises none of its own, which is why the dispatcher hands it `Vec::new()`.
    #[must_use]
    pub fn of(
        outcome: &kernel::HeatExchangerOutcome,
        mechanical: Option<&heater_mechanical::MechanicalDesign>,
        warnings: Vec<Warning>,
    ) -> Self {
        Self {
            hot_out_n: outcome.hot_out.n,
            hot_out_z: outcome.hot_out.z.clone(),
            hot_out_p: outcome.hot_out.p,
            hot_out_t: outcome.hot_out.t,
            hot_out_h: joules_per_mole(outcome.hot_out.h.value),
            cold_out_n: outcome.cold_out.n,
            cold_out_z: outcome.cold_out.z.clone(),
            cold_out_p: outcome.cold_out.p,
            cold_out_t: outcome.cold_out.t,
            cold_out_h: joules_per_mole(outcome.cold_out.h.value),
            duty: outcome.duty,
            ntu: outcome.rating.as_ref().map(|rating| rating.ntu),
            effectiveness: outcome.rating.as_ref().map(|rating| rating.effectiveness),
            c_min: outcome.rating.as_ref().map(|rating| rating.c_min),
            c_max: outcome.rating.as_ref().map(|rating| rating.c_max),
            capacity_ratio: outcome.rating.as_ref().map(|rating| rating.capacity_ratio),
            ua_value: outcome.ua_value,
            delta_t_side: outcome.delta_t_side.map(str::to_string),
            mechanical_selected_type: mechanical
                .map(|value| value.selected.kind.spelling().to_string()),
            mechanical_selected_area: mechanical
                .map(|value| square_meters(value.selected.required_area_m2)),
            log_mean_temperature_difference: mechanical
                .map(|value| kelvins(value.log_mean_temperature_difference)),
            mechanical_approach_temperature: mechanical
                .map(|value| kelvins(value.approach_temperature)),
            used_overall_heat_transfer_coefficient: mechanical.map(|value| {
                watts_per_square_meter_kelvin(value.used_overall_heat_transfer_coefficient)
            }),
            calculated_ua: mechanical.map(|value| watts_per_kelvin(value.calculated_ua)),
            rating_u: outcome
                .rating_diagnostics
                .as_ref()
                .map(|value| watts_per_square_meter_kelvin(value.rating_u)),
            rating_ua: outcome
                .rating_diagnostics
                .as_ref()
                .map(|value| watts_per_kelvin(value.rating_ua)),
            rating_area: outcome
                .rating_diagnostics
                .as_ref()
                .map(|value| square_meters(value.rating_area)),
            tube_side_pressure_drop: outcome
                .rating_diagnostics
                .as_ref()
                .map(|value| pascals(value.tube_side_pressure_drop)),
            shell_side_pressure_drop: outcome
                .rating_diagnostics
                .as_ref()
                .map(|value| pascals(value.shell_side_pressure_drop)),
            tube_side_heat_transfer_coefficient: outcome.rating_diagnostics.as_ref().map(|value| {
                watts_per_square_meter_kelvin(value.tube_side_heat_transfer_coefficient)
            }),
            shell_side_heat_transfer_coefficient: outcome.rating_diagnostics.as_ref().map(
                |value| watts_per_square_meter_kelvin(value.shell_side_heat_transfer_coefficient),
            ),
            tube_side_reynolds: outcome
                .rating_diagnostics
                .as_ref()
                .map(|value| value.tube_side_reynolds),
            shell_side_reynolds: outcome
                .rating_diagnostics
                .as_ref()
                .map(|value| value.shell_side_reynolds),
            tube_side_velocity: outcome
                .rating_diagnostics
                .as_ref()
                .map(|value| meters_per_second(value.tube_side_velocity)),
            shell_side_velocity: outcome
                .rating_diagnostics
                .as_ref()
                .map(|value| meters_per_second(value.shell_side_velocity)),
            jc: outcome.rating_diagnostics.as_ref().map(|value| value.jc),
            jl: outcome.rating_diagnostics.as_ref().map(|value| value.jl),
            jb: outcome.rating_diagnostics.as_ref().map(|value| value.jb),
            js: outcome.rating_diagnostics.as_ref().map(|value| value.js),
            jr: outcome.rating_diagnostics.as_ref().map(|value| value.jr),
            mechanical_types: mechanical
                .map(|value| {
                    value
                        .types
                        .iter()
                        .map(|row| HeatExchangerMechanicalType {
                            mechanical_type: row.kind.spelling().to_string(),
                            required_area: square_meters(row.required_area_m2),
                            required_ua: watts_per_kelvin(row.required_ua_w_per_k),
                            overall_heat_transfer_coefficient: watts_per_square_meter_kelvin(
                                row.overall_heat_transfer_coefficient,
                            ),
                            approach_temperature: kelvins(row.approach_temperature_k),
                            estimated_length: meters(row.estimated_length_m),
                            inner_diameter: meters(row.inner_diameter_m),
                            outer_diameter: meters(row.outer_diameter_m),
                            wall_thickness: meters(row.wall_thickness_m),
                            tube_count: row.tube_count,
                            tube_passes: row.tube_passes,
                            fin_surface_area: square_meters(row.fin_surface_area_m2),
                            estimated_pressure_drop: row.estimated_pressure_drop_bar,
                            estimated_weight: kilograms(row.estimated_weight_kg),
                            module_length: meters(row.module_length_m),
                            module_width: meters(row.module_width_m),
                            module_height: meters(row.module_height_m),
                        })
                        .collect()
                })
                .unwrap_or_default(),
            warnings,
        }
    }
}
