//! `process.rate_based_packed_column` - the segment model as a registered id.
//!
//! Spec: `specs/models/process/rate_based_packed_column.toml`. **The arithmetic is
//! [`crate::segment`]'s**; what this module adds is the class's own defaults for the eleven
//! palette parameters, the refusal of the four enum values it does not carry, and the flat
//! record a case and a capture can address.

use azoth_core::units::{
    DiffusionCoefficient, DynamicViscosity, Length, MassDensity, MolarEnergy, Power, Pressure,
    ThermodynamicTemperature,
};
use azoth_core::units::{
    kilograms_per_cubic_meter, pascal_seconds, square_meters_per_second, watts,
};
use azoth_core::warning::{Warning, WarningCode};
use azoth_core::{AzothError, CalcResult, Result};

use crate::kernels::rate_based_packed_column::{
    ColumnSolver, FilmModel, HeatTransferModel, MassTransferCorrelation, RateBasedSetup,
    SegmentSolver, rate_based_packed_column as kernel,
};
use crate::segment::{Property, SegmentResult};
use crate::stream::Stream;

/// Result of `process.rate_based_packed_column`.
///
/// **The per-segment profile is parallel vectors over the segments**, the shape
/// `eos.pt_phase_envelope`'s trace points set: one entry per segment for each of the class's
/// carried `SegmentResult` fields. Four of the class's thirty-four are not here - three are
/// constants on the ported path and one is derived - and the module's own doc for
/// [`crate::segment::step::SegmentResult`] says which and why.
#[derive(Debug, Clone, PartialEq)]
pub struct RateBasedPackedColumnResult {
    /// The gas leaving the top segment, mol/s.
    pub gas_out_n: f64,
    /// The gas outlet's composition.
    pub gas_out_z: Vec<f64>,
    /// The gas outlet's pressure.
    pub gas_out_p: Pressure,
    /// The gas outlet's temperature.
    pub gas_out_t: ThermodynamicTemperature,
    /// The gas outlet's molar enthalpy.
    pub gas_out_h: MolarEnergy,
    /// The liquid leaving the bottom segment, mol/s.
    pub liquid_out_n: f64,
    /// The liquid outlet's composition.
    pub liquid_out_z: Vec<f64>,
    /// The liquid outlet's pressure.
    pub liquid_out_p: Pressure,
    /// The liquid outlet's temperature.
    pub liquid_out_t: ThermodynamicTemperature,
    /// The liquid outlet's molar enthalpy.
    pub liquid_out_h: MolarEnergy,
    /// The passes taken.
    pub iterations: u32,
    /// The outlet residual at the last pass, mol/s.
    pub convergence_residual: f64,
    /// Whether the gate was met.
    pub converged: bool,
    /// The sum of every segment's transfers, magnitudes added, mol/s.
    pub total_absolute_molar_transfer: f64,
    /// Each component's net total, mol/s.
    pub component_transfer_totals: Vec<f64>,
    /// Each component the totals are stated over, in the same order.
    pub transfer_components: Vec<String>,
    pub segment_height_from_bottom: Vec<Length>,
    pub segment_gas_temperature: Vec<ThermodynamicTemperature>,
    pub segment_liquid_temperature: Vec<ThermodynamicTemperature>,
    pub segment_gas_pressure: Vec<Pressure>,
    pub segment_liquid_pressure: Vec<Pressure>,
    pub segment_gas_molar_flow: Vec<f64>,
    pub segment_liquid_molar_flow: Vec<f64>,
    pub segment_gas_density: Vec<MassDensity>,
    pub segment_liquid_density: Vec<MassDensity>,
    pub segment_gas_viscosity: Vec<DynamicViscosity>,
    pub segment_liquid_viscosity: Vec<DynamicViscosity>,
    pub segment_gas_diffusivity: Vec<DiffusionCoefficient>,
    pub segment_liquid_diffusivity: Vec<DiffusionCoefficient>,
    pub segment_wetted_area: Vec<f64>,
    pub segment_k_ga: Vec<f64>,
    pub segment_k_la: Vec<f64>,
    pub segment_gas_heat_transfer_coefficient: Vec<f64>,
    pub segment_liquid_heat_transfer_coefficient: Vec<f64>,
    pub segment_overall_heat_transfer_coefficient: Vec<f64>,
    pub segment_interface_temperature: Vec<ThermodynamicTemperature>,
    pub segment_heat_transfer_rate: Vec<Power>,
    pub segment_pressure_drop_per_meter: Vec<Pressure>,
    pub segment_percent_flood: Vec<f64>,
    pub segment_net_molar_transfer: Vec<f64>,
    pub segment_enthalpy_balance_residual: Vec<Power>,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for RateBasedPackedColumnResult {
    const CALC_ID: &'static str = "process.rate_based_packed_column";
    const FIELDS: &'static [&'static str] = &[
        "gas_out_n",
        "gas_out_z",
        "gas_out_p",
        "gas_out_t",
        "gas_out_h",
        "liquid_out_n",
        "liquid_out_z",
        "liquid_out_p",
        "liquid_out_t",
        "liquid_out_h",
        "iterations",
        "convergence_residual",
        "converged",
        "total_absolute_molar_transfer",
        "component_transfer_totals",
        "transfer_components",
        "segment_height_from_bottom",
        "segment_gas_temperature",
        "segment_liquid_temperature",
        "segment_gas_pressure",
        "segment_liquid_pressure",
        "segment_gas_molar_flow",
        "segment_liquid_molar_flow",
        "segment_gas_density",
        "segment_liquid_density",
        "segment_gas_viscosity",
        "segment_liquid_viscosity",
        "segment_gas_diffusivity",
        "segment_liquid_diffusivity",
        "segment_wetted_area",
        "segment_k_ga",
        "segment_k_la",
        "segment_gas_heat_transfer_coefficient",
        "segment_liquid_heat_transfer_coefficient",
        "segment_overall_heat_transfer_coefficient",
        "segment_interface_temperature",
        "segment_heat_transfer_rate",
        "segment_pressure_drop_per_meter",
        "segment_percent_flood",
        "segment_net_molar_transfer",
        "segment_enthalpy_balance_residual",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// The class's own default column diameter, in metres.
const DEFAULT_COLUMN_DIAMETER: f64 = 1.0;
/// Its default packed height, in metres.
const DEFAULT_PACKED_HEIGHT: f64 = 5.0;
/// Its default segment count.
const DEFAULT_SEGMENTS: usize = 10;
/// Its default iteration cap.
const DEFAULT_MAX_ITERATIONS: usize = 30;
/// Its default convergence gate, mol/s.
const DEFAULT_TOLERANCE: f64 = 1.0e-8;

/// Solve a rate-based packed column.
///
/// # Errors
/// [`AzothError::InvalidInput`] for the four enum values the port refuses, each naming the
/// class behind it; for the five conditions `validateSetup` states; and whatever the profile
/// solve raises.
#[allow(clippy::too_many_arguments)] // One per declared input, and the spec declares nineteen.
pub fn rate_based_packed_column(
    gas_components: &[String],
    gas_n: f64,
    gas_z: &[f64],
    gas_p: Pressure,
    gas_t: ThermodynamicTemperature,
    liquid_components: &[String],
    liquid_n: f64,
    liquid_z: &[f64],
    liquid_p: Pressure,
    liquid_t: ThermodynamicTemperature,
    transfer_components: Option<&[String]>,
    column_diameter: Option<f64>,
    packed_height: Option<f64>,
    number_of_segments: Option<usize>,
    packing_type: Option<&str>,
    max_iterations: Option<usize>,
    convergence_tolerance: Option<f64>,
    mass_transfer_correction: Option<f64>,
    mass_transfer_correlation: Option<&str>,
    film_model: Option<&str>,
    heat_transfer_model: Option<&str>,
    segment_solver: Option<&str>,
    column_solver: Option<&str>,
) -> Result<RateBasedPackedColumnResult> {
    // **The four refusals come first, and each names its class** - so a caller who states a
    // solver the port does not carry hears which one would close it rather than reading a
    // number from a machine that silently ran something else.
    MassTransferCorrelation::parse(mass_transfer_correlation.unwrap_or("onda_1968"))?;
    let film = FilmModel::parse(film_model.unwrap_or("maxwell_stefan_matrix"))?;
    let heat = HeatTransferModel::parse(heat_transfer_model.unwrap_or("chilton_colburn_analogy"))?;
    SegmentSolver::parse(segment_solver.unwrap_or("sequential_explicit"))?;
    ColumnSolver::parse(column_solver.unwrap_or("fixed_point_profile"))?;

    if gas_components.len() != gas_z.len() {
        return Err(AzothError::invalid_input(
            "gas_z",
            format!(
                "{} component(s) and {} fraction(s): a composition is one entry per component",
                gas_components.len(),
                gas_z.len()
            ),
        ));
    }
    if liquid_components.len() != liquid_z.len() {
        return Err(AzothError::invalid_input(
            "liquid_z",
            format!(
                "{} component(s) and {} fraction(s): a composition is one entry per component",
                liquid_components.len(),
                liquid_z.len()
            ),
        ));
    }

    let gas = Stream::from_pt(gas_components.to_vec(), gas_z.to_vec(), gas_n, gas_p, gas_t)?;
    let liquid = Stream::from_pt(
        liquid_components.to_vec(),
        liquid_z.to_vec(),
        liquid_n,
        liquid_p,
        liquid_t,
    )?;

    let outcome = kernel(&RateBasedSetup {
        gas,
        liquid,
        column_diameter: column_diameter.unwrap_or(DEFAULT_COLUMN_DIAMETER),
        packed_height: packed_height.unwrap_or(DEFAULT_PACKED_HEIGHT),
        number_of_segments: number_of_segments.unwrap_or(DEFAULT_SEGMENTS),
        packing_type: packing_type.unwrap_or("Pall-Ring-50").to_string(),
        max_iterations: max_iterations.unwrap_or(DEFAULT_MAX_ITERATIONS),
        convergence_tolerance: convergence_tolerance.unwrap_or(DEFAULT_TOLERANCE),
        mass_transfer_correction: mass_transfer_correction.unwrap_or(1.0),
        heat_transfer_correction: 1.0,
        transfer_components: transfer_components.map(<[String]>::to_vec),
        film_model: film,
        heat_transfer_model: heat,
    })?;

    let names: Vec<String> = outcome
        .component_transfer_totals
        .iter()
        .map(|(name, _)| name.clone())
        .collect();
    let totals: Vec<f64> = outcome
        .component_transfer_totals
        .iter()
        .map(|(_, total)| *total)
        .collect();

    let column = |pick: fn(&SegmentResult) -> f64| -> Vec<f64> {
        outcome.segments.iter().map(pick).collect()
    };
    let temperatures = |pick: fn(&SegmentResult) -> f64| -> Vec<ThermodynamicTemperature> {
        outcome
            .segments
            .iter()
            .map(|segment| azoth_core::units::kelvins(pick(segment)))
            .collect()
    };
    let pressures = |pick: fn(&SegmentResult) -> f64| -> Vec<Pressure> {
        outcome
            .segments
            .iter()
            .map(|segment| azoth_core::units::pascals(pick(segment)))
            .collect()
    };
    let lengths = |pick: fn(&SegmentResult) -> f64| -> Vec<Length> {
        outcome
            .segments
            .iter()
            .map(|segment| azoth_core::units::meters(pick(segment)))
            .collect()
    };
    let densities = |pick: fn(&SegmentResult) -> f64| -> Vec<MassDensity> {
        outcome
            .segments
            .iter()
            .map(|segment| kilograms_per_cubic_meter(pick(segment)))
            .collect()
    };
    let viscosities = |pick: fn(&SegmentResult) -> f64| -> Vec<DynamicViscosity> {
        outcome
            .segments
            .iter()
            .map(|segment| pascal_seconds(pick(segment)))
            .collect()
    };
    let diffusivities = |pick: fn(&SegmentResult) -> f64| -> Vec<DiffusionCoefficient> {
        outcome
            .segments
            .iter()
            .map(|segment| square_meters_per_second(pick(segment)))
            .collect()
    };
    let powers = |pick: fn(&SegmentResult) -> f64| -> Vec<Power> {
        outcome
            .segments
            .iter()
            .map(|segment| watts(pick(segment)))
            .collect()
    };

    Ok(RateBasedPackedColumnResult {
        gas_out_n: outcome.gas_out.n,
        gas_out_z: outcome.gas_out.z.clone(),
        gas_out_p: outcome.gas_out.p,
        gas_out_t: outcome.gas_out.t,
        gas_out_h: outcome.gas_out.h,
        liquid_out_n: outcome.liquid_out.n,
        liquid_out_z: outcome.liquid_out.z.clone(),
        liquid_out_p: outcome.liquid_out.p,
        liquid_out_t: outcome.liquid_out.t,
        liquid_out_h: outcome.liquid_out.h,
        iterations: outcome.iterations as u32,
        convergence_residual: outcome.convergence_residual,
        converged: outcome.converged,
        total_absolute_molar_transfer: outcome.total_absolute_molar_transfer,
        component_transfer_totals: totals,
        transfer_components: names,
        segment_height_from_bottom: lengths(|segment| segment.height_from_bottom),
        segment_gas_temperature: temperatures(|segment| segment.gas_temperature),
        segment_liquid_temperature: temperatures(|segment| segment.liquid_temperature),
        segment_gas_pressure: pressures(|segment| segment.gas_pressure),
        segment_liquid_pressure: pressures(|segment| segment.liquid_pressure),
        segment_gas_molar_flow: column(|segment| segment.gas_molar_flow),
        segment_liquid_molar_flow: column(|segment| segment.liquid_molar_flow),
        segment_gas_density: densities(|segment| segment.gas_density),
        segment_liquid_density: densities(|segment| segment.liquid_density),
        segment_gas_viscosity: viscosities(|segment| segment.gas_viscosity),
        segment_liquid_viscosity: viscosities(|segment| segment.liquid_viscosity),
        segment_gas_diffusivity: diffusivities(|segment| segment.gas_diffusivity),
        segment_liquid_diffusivity: diffusivities(|segment| segment.liquid_diffusivity),
        segment_wetted_area: column(|segment| segment.wetted_area),
        segment_k_ga: column(|segment| segment.k_ga),
        segment_k_la: column(|segment| segment.k_la),
        segment_gas_heat_transfer_coefficient: column(|segment| {
            segment.gas_heat_transfer_coefficient
        }),
        segment_liquid_heat_transfer_coefficient: column(|segment| {
            segment.liquid_heat_transfer_coefficient
        }),
        segment_overall_heat_transfer_coefficient: column(|segment| {
            segment.overall_heat_transfer_coefficient
        }),
        segment_interface_temperature: temperatures(|segment| segment.interface_temperature),
        segment_heat_transfer_rate: powers(|segment| segment.heat_transfer_rate),
        segment_pressure_drop_per_meter: pressures(|segment| segment.pressure_drop_per_meter),
        segment_percent_flood: column(|segment| segment.percent_flood),
        segment_net_molar_transfer: column(|segment| segment.net_molar_transfer),
        segment_enthalpy_balance_residual: powers(|segment| segment.enthalpy_balance_residual),
        warnings: warnings_for(outcome.fallbacks, outcome.converged, outcome.iterations),
    })
}

/// **A substitution the class makes silently is a caveat here**, one per property rather than
/// one per segment, and **a profile that missed its gate is published with its residual**
/// rather than refused - which is what the class does and what its own tests expect.
fn warnings_for(
    fallbacks: crate::segment::Fallbacks,
    converged: bool,
    iterations: usize,
) -> Vec<Warning> {
    let mut warnings: Vec<Warning> = fallbacks
        .taken()
        .into_iter()
        .map(|property| Warning {
            code: WarningCode::SolverNotConverged,
            message: format!(
                "the {} is the class's own `DEFAULT_*` constant: the physical-property model \
                 answered nothing, and `RateBasedPackedColumn` substitutes rather than refusing",
                property.as_str()
            ),
            field: Some(property_field(property).to_string()),
        })
        .collect();
    if !converged {
        warnings.push(Warning {
            code: WarningCode::SolverNotConverged,
            message: format!(
                "the profile reached its {iterations}-iteration cap without meeting the gate, \
                 and its last iterate is published - which is what the class does"
            ),
            field: Some("max_iterations".to_string()),
        });
    }
    warnings
}

/// The input a substituted property belongs to, which is the segment vector a caller reads it
/// from.
fn property_field(property: Property) -> &'static str {
    match property {
        Property::Diffusivity => "segment_gas_diffusivity",
        Property::SurfaceTension => "segment_wetted_area",
        Property::ThermalConductivity => "segment_gas_heat_transfer_coefficient",
        Property::HeatCapacity => "segment_gas_heat_transfer_coefficient",
    }
}
