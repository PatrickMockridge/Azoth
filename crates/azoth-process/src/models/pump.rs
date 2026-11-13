//! `process.pump` - the pump's kernel as a registered id.
//!
//! Spec: `specs/models/process/pump.toml`. The arithmetic is
//! [`crate::kernels::pump`]; what is here is the boundary a case and a cross-impl test
//! address, which is the same split every other namespace makes.

use azoth_core::units::{
    DiffusionCoefficient, Length, Mass, MassDensity, MolarEnergy, Power, Pressure,
    ThermodynamicTemperature, VolumeRate, cubic_meter_per_hour, joules_per_mole, meter, meters,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};

use crate::kernels::pump as kernel;
use crate::kernels::pump::{
    cavitating, npsh_available, npsh_required, pump_with_chart, shaft_power,
};
use crate::model_gen;
use crate::pump_chart::{CurvePoint, PumpChart};
use crate::rotating_design::pump::{
    Api610PumpType, PumpDesignInputs, SealType, design as design_pump,
};
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
    /// The legacy pump family: `overhung`, `between_bearings` or `vertically_suspended`.
    pub design_pump_type: Option<String>,
    /// The API 610 construction type.
    pub design_api610_pump_type: Option<String>,
    /// The seal type.
    pub design_seal_type: Option<String>,
    /// The rated volumetric flow, m3/h.
    pub design_rated_flow: Option<VolumeRate>,
    /// The design's head, m, over the class's own 9.81.
    pub design_rated_head: Option<Length>,
    /// The hydraulic power, W.
    pub design_hydraulic_power: Option<Power>,
    /// The rated speed, rpm.
    pub design_rated_speed: Option<f64>,
    /// The SI specific speed.
    pub design_specific_speed: Option<f64>,
    /// The US suction specific speed.
    pub design_suction_specific_speed: Option<f64>,
    /// The stage count.
    pub design_number_of_stages: Option<f64>,
    /// The vendor best-efficiency flow, m3/h.
    pub design_bep_flow: Option<VolumeRate>,
    /// The vendor best-efficiency head, m.
    pub design_bep_head: Option<Length>,
    /// The BEP's provenance.
    pub design_bep_source: Option<String>,
    /// The NPSHr provenance.
    pub design_npshr_source: Option<String>,
    /// The shutoff head's provenance.
    pub design_shutoff_source: Option<String>,
    /// The effective shutoff head, m.
    pub design_shutoff_head: Option<Length>,
    /// The casing design pressure, Pa.
    pub design_pressure: Option<Pressure>,
    /// The casing design temperature, K.
    pub design_temperature: Option<ThermodynamicTemperature>,
    /// The impeller diameter, m.
    pub design_impeller_diameter: Option<Length>,
    /// The impeller width, m.
    pub design_impeller_width: Option<Length>,
    /// The shaft diameter, m.
    pub design_shaft_diameter: Option<Length>,
    /// The casing wall thickness, m.
    pub design_casing_wall_thickness: Option<Length>,
    /// The maximum allowable working pressure, Pa.
    pub design_max_allowable_working_pressure: Option<Pressure>,
    /// The final driver power, W.
    pub design_driver_power: Option<Power>,
    /// The driver margin.
    pub design_driver_margin: Option<f64>,
    /// The suction nozzle size, m.
    pub design_suction_nozzle_size: Option<Length>,
    /// The discharge nozzle size, m.
    pub design_discharge_nozzle_size: Option<Length>,
    /// The minimum continuous flow, m3/h.
    pub design_minimum_flow: Option<VolumeRate>,
    /// The maximum continuous flow, m3/h.
    pub design_maximum_flow: Option<VolumeRate>,
    /// The casing weight, kg.
    pub design_casing_weight: Option<Mass>,
    /// The motor weight, kg.
    pub design_motor_weight: Option<Mass>,
    /// The baseplate weight, kg.
    pub design_baseplate_weight: Option<Mass>,
    /// The casing inner diameter, m.
    pub design_inner_diameter: Option<Length>,
    /// The casing outer diameter, m.
    pub design_outer_diameter: Option<Length>,
    /// The tangent-to-tangent length, m.
    pub design_tantan_length: Option<Length>,
    /// The base class's wall thickness, m.
    pub design_wall_thickness: Option<Length>,
    /// The piping weight, kg.
    pub design_weight_piping: Option<Mass>,
    /// The electrical and instrument weight, kg.
    pub design_weight_electro_instrument: Option<Mass>,
    /// The structural steel weight, kg.
    pub design_weight_structural_steel: Option<Mass>,
    /// The total weight, kg.
    pub design_weight_total: Option<Mass>,
    /// The module length, m.
    pub design_module_length: Option<Length>,
    /// The module width, m.
    pub design_module_width: Option<Length>,
    /// The module height, m.
    pub design_module_height: Option<Length>,
    /// The NPSH margin, m.
    pub design_npsh_margin: Option<Length>,
    /// The NPSH margin ratio.
    pub design_npsh_margin_ratio: Option<f64>,
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
    mechanical_design: Option<bool>,
    api610_pump_type: Option<String>,
    seal_type: Option<String>,
    maximum_suction_pressure: Option<Pressure>,
    shutoff_head: Option<Length>,
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

    let design_requested = mechanical_design.unwrap_or(false);
    if !design_requested {
        for (name, stated) in [
            ("api610_pump_type", api610_pump_type.is_some()),
            ("seal_type", seal_type.is_some()),
            (
                "maximum_suction_pressure",
                maximum_suction_pressure.is_some(),
            ),
            ("shutoff_head", shutoff_head.is_some()),
        ] {
            if stated {
                return Err(AzothError::invalid_input(
                    name,
                    format!("`{name}` is a `mechanical_design` input, and the design is off"),
                ));
            }
        }
    }
    let api610 = match api610_pump_type.as_deref() {
        None => None,
        Some(value) => Some(Api610PumpType::parse(value).ok_or_else(|| {
            AzothError::invalid_input(
                "api610_pump_type",
                format!("`{value}` is not an API 610 construction type"),
            )
        })?),
    };
    let seal = match seal_type.as_deref() {
        None => SealType::SingleMechanical,
        Some(value) => SealType::parse(value).ok_or_else(|| {
            AzothError::invalid_input("seal_type", format!("`{value}` is not a seal type"))
        })?,
    };

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
    let design = if design_requested {
        let inputs = PumpDesignInputs {
            api610_pump_type: api610,
            seal_type: seal,
            maximum_suction_pressure,
            shutoff_head,
            speed,
            chart: chart.as_ref(),
        };
        Some(design_pump(
            &feed,
            &out,
            power,
            efficiency,
            meters(required),
            meters(available),
            &inputs,
        )?)
    } else {
        None
    };
    let d = design.as_ref();

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
        design_pump_type: d.map(|value| value.pump_type.to_owned()),
        design_api610_pump_type: d.map(|value| value.api610_pump_type.to_owned()),
        design_seal_type: d.map(|value| value.seal_type.to_owned()),
        design_rated_flow: d.map(|value| value.rated_flow),
        design_rated_head: d.map(|value| value.rated_head),
        design_hydraulic_power: d.map(|value| value.hydraulic_power),
        design_rated_speed: d.map(|value| value.rated_speed),
        design_specific_speed: d.map(|value| value.specific_speed),
        design_suction_specific_speed: d.map(|value| value.suction_specific_speed),
        design_number_of_stages: d.map(|value| f64::from(value.number_of_stages)),
        design_bep_flow: d.and_then(|value| value.bep_flow),
        design_bep_head: d.and_then(|value| value.bep_head),
        design_bep_source: d.map(|value| value.bep_source.to_owned()),
        design_npshr_source: d.map(|value| value.npshr_source.to_owned()),
        design_shutoff_source: d.map(|value| value.shutoff_source.to_owned()),
        design_shutoff_head: d.and_then(|value| value.shutoff_head),
        design_pressure: d.map(|value| value.design_pressure),
        design_temperature: d.map(|value| value.design_temperature),
        design_impeller_diameter: d.map(|value| value.impeller_diameter),
        design_impeller_width: d.map(|value| value.impeller_width),
        design_shaft_diameter: d.map(|value| value.shaft_diameter),
        design_casing_wall_thickness: d.map(|value| value.casing_wall_thickness),
        design_max_allowable_working_pressure: d.map(|value| value.max_allowable_working_pressure),
        design_driver_power: d.map(|value| value.driver_power),
        design_driver_margin: d.map(|value| value.driver_margin),
        design_suction_nozzle_size: d.map(|value| value.suction_nozzle_size),
        design_discharge_nozzle_size: d.map(|value| value.discharge_nozzle_size),
        design_minimum_flow: d.map(|value| value.minimum_flow),
        design_maximum_flow: d.map(|value| value.maximum_flow),
        design_casing_weight: d.map(|value| value.casing_weight),
        design_motor_weight: d.map(|value| value.motor_weight),
        design_baseplate_weight: d.map(|value| value.baseplate_weight),
        design_inner_diameter: d.map(|value| value.inner_diameter),
        design_outer_diameter: d.map(|value| value.outer_diameter),
        design_tantan_length: d.map(|value| value.tantan_length),
        design_wall_thickness: d.map(|value| value.wall_thickness),
        design_weight_piping: d.map(|value| value.weight_piping),
        design_weight_electro_instrument: d.map(|value| value.weight_electro_instrument),
        design_weight_structural_steel: d.map(|value| value.weight_structural_steel),
        design_weight_total: d.map(|value| value.weight_total),
        design_module_length: d.map(|value| value.module_length),
        design_module_width: d.map(|value| value.module_width),
        design_module_height: d.map(|value| value.module_height),
        design_npsh_margin: d.map(|value| value.npsh_margin),
        design_npsh_margin_ratio: d.map(|value| value.npsh_margin_ratio),
        warnings,
    })
}
