//! `process.gibbs_reactor` - the registry surface over the kernel.
//!
//! The kernel is [`crate::kernels::gibbs_reactor`], which takes a [`Stream`]; this is the same
//! arithmetic as an id the registry addresses, with the declared input list and the checks the
//! spec carries. **Two surfaces, one arithmetic**, for the reason `models`' own module doc gives:
//! the executor calls kernels and a case calls models, and neither wraps the other's types.
//!
//! **The vessel is a second arithmetic beside the solve.** `ReactorMechanicalDesign.calcDesign`
//! reads the *inlet* stream and nothing the solve produces, so it is computed from the feed
//! before the equilibrium and published on the same result.

use azoth_core::units::{
    Length, Mass, MolarEnergy, Pressure, ThermodynamicTemperature, Volume, kelvins, meters, pascals,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{scalar, warnings as wire_warnings};
use crate::kernels::gibbs_reactor::{
    EnergyMode, ReactorNumbers, ReactorSetup, gibbs_reactor as solve_equilibrium,
};
use crate::model_gen;
use crate::reactor::mechanical::{
    CORROSION_ALLOWANCE_MM, DEFAULT_BED_VOID_FRACTION, DEFAULT_CATALYST_BULK_DENSITY,
    DEFAULT_CATALYST_PARTICLE_DIAMETER_MM, DEFAULT_GHSV_PER_HOUR, DEFAULT_LHSV_PER_HOUR,
    DEFAULT_MAX_OPERATION_PRESSURE_BARA, DEFAULT_MAX_OPERATION_TEMPERATURE_K,
    DEFAULT_NUMBER_OF_BEDS, JOINT_EFFICIENCY, ReactorMechanicalGeometry, ReactorMechanicalReport,
    reactor_mechanical_design,
};
use crate::stream::Stream;

/// Result of `process.gibbs_reactor`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct GibbsReactorResult {
    /// Product molar flow, mol/s.
    pub product_n: f64,
    /// Product composition, over the feed's own species order.
    pub product_z: Vec<f64>,
    /// Product pressure.
    #[serde(serialize_with = "scalar")]
    pub product_p: Pressure,
    /// Product temperature.
    #[serde(serialize_with = "scalar")]
    pub product_t: ThermodynamicTemperature,
    /// Product molar enthalpy.
    #[serde(serialize_with = "scalar")]
    pub product_h: MolarEnergy,
    /// `hasConverged()`, which a run that reached the cap reports as false.
    pub converged: bool,
    /// The pass the loop stopped on.
    pub iterations: f64,
    /// The last undamped step norm.
    pub final_error: f64,
    /// The element Lagrange multipliers, J/mol, on the class's own seven element names.
    pub lagrange_multipliers: Vec<f64>,
    /// The outlet element balance less the inlet's, mol/s.
    pub element_balance_difference: Vec<f64>,
    /// The total Gibbs energy at the top of every iteration, W.
    pub gibbs_energy_history: Vec<f64>,
    /// `getVesselDiameter()`, m.
    #[serde(serialize_with = "scalar")]
    pub vessel_diameter: Length,
    /// `getVesselLength()`, m.
    #[serde(serialize_with = "scalar")]
    pub vessel_length: Length,
    /// `getShellThickness()`. SI, shown in mm.
    #[serde(serialize_with = "scalar")]
    pub shell_thickness: Length,
    /// `getHeadThickness()`, which is always the shell's. SI, shown in mm.
    #[serde(serialize_with = "scalar")]
    pub head_thickness: Length,
    /// `getDesignPressureBara()`. SI, shown in bar.
    #[serde(serialize_with = "scalar")]
    pub design_pressure: Pressure,
    /// `getDesignTemperatureC()`, carried in kelvin.
    #[serde(serialize_with = "scalar")]
    pub design_temperature: ThermodynamicTemperature,
    /// `getCatalystVolume()`, m3.
    #[serde(serialize_with = "scalar")]
    pub catalyst_volume: Volume,
    /// `getCatalystMass()`, kg.
    #[serde(serialize_with = "scalar")]
    pub catalyst_mass: Mass,
    /// `getBedPressureDrop()`. SI, shown in bar.
    #[serde(serialize_with = "scalar")]
    pub bed_pressure_drop: Pressure,
    /// `getEmptyVesselWeight()`, kg.
    #[serde(serialize_with = "scalar")]
    pub empty_vessel_weight: Mass,
    /// `getTotalEquippedWeight()`, kg.
    #[serde(serialize_with = "scalar")]
    pub total_equipped_weight: Mass,
    /// `getReactorType()`, in the class's own spelling.
    pub reactor_type: String,
    /// `getNumberOfBeds()`.
    pub number_of_beds: usize,
    /// `getConditionNumberHistory()`, one 2-norm condition number per recorded pass.
    pub condition_number_history: Vec<f64>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

impl GibbsReactorResult {
    /// The result of one kernel call, from the solve's own two answers and the vessel around
    /// them.
    ///
    /// **The warnings are the caller's, because the two callers have different ones.** A case
    /// runs the spec's checks through [`apply_checks`]; a flowsheet's equivalent is the
    /// *checker's*, which reports through the envelope's diagnostics rather than through a result.
    #[must_use]
    pub fn of(
        outlet: &Stream,
        numbers: &ReactorNumbers,
        mechanical: &ReactorMechanicalReport,
        warnings: Vec<Warning>,
    ) -> Self {
        Self {
            product_n: outlet.n,
            product_z: outlet.z.clone(),
            product_p: pascals(outlet.p.value),
            product_t: kelvins(outlet.t.value),
            product_h: outlet.h,
            converged: numbers.converged,
            iterations: f64::from(numbers.iterations),
            final_error: numbers.final_error,
            lagrange_multipliers: numbers
                .lagrange_multipliers
                .iter()
                .map(|v| v * 1000.0)
                .collect(),
            element_balance_difference: numbers.element_balance_difference.to_vec(),
            gibbs_energy_history: numbers
                .gibbs_energy_history
                .iter()
                .map(|v| v * 1000.0)
                .collect(),
            vessel_diameter: mechanical.vessel_diameter,
            vessel_length: mechanical.vessel_length,
            shell_thickness: meters(mechanical.shell_thickness_mm / 1000.0),
            head_thickness: meters(mechanical.head_thickness_mm / 1000.0),
            design_pressure: pascals(mechanical.design_pressure_bara * 1.0e5),
            design_temperature: kelvins(mechanical.design_temperature_c + 273.15),
            catalyst_volume: mechanical.catalyst_volume,
            catalyst_mass: mechanical.catalyst_mass,
            bed_pressure_drop: pascals(mechanical.bed_pressure_drop_bar * 1.0e5),
            empty_vessel_weight: mechanical.empty_vessel_weight,
            total_equipped_weight: mechanical.total_equipped_weight,
            reactor_type: mechanical.reactor_type.clone(),
            number_of_beds: mechanical.number_of_beds,
            condition_number_history: numbers.condition_number_history.clone(),
            warnings,
        }
    }
}

/// The spec's four reactor types, in the class's own spelling.
///
/// **The class stores and echoes whatever string it is given**, and only `CSTR` changes an
/// answer. This id's declaration is the four values, so anything else is refused rather than
/// echoed.
fn reactor_type_name(kind: &str) -> Result<String> {
    match kind {
        "fixed-bed" => Ok("FIXED_BED".to_string()),
        "cstr" => Ok("CSTR".to_string()),
        "packed-bed" => Ok("PACKED_BED".to_string()),
        "fluidized-bed" => Ok("FLUIDIZED_BED".to_string()),
        other => Err(AzothError::invalid_input(
            "reactor_type",
            format!("'{other}' is not one of fixed-bed, cstr, packed-bed, fluidized-bed"),
        )),
    }
}

/// Bring a feed to its Gibbs equilibrium at its own temperature and pressure.
///
/// **The class's own arithmetic is in kJ/mol and two of its outputs cross in other units.** The
/// multipliers are kJ/mol there and J/mol here, and the Gibbs energy history is a sum of `mol/s`
/// against `kJ/mol` - a rate of energy, reported in watts. The spec records both conversions.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] on a feed the spec's ranges refuse, on a fluid that
///   cannot be resolved, or on the singular Newton system the class answers with a pseudo-inverse
///   this port does not carry;
/// * the reactor mechanical design's own refusals, and whatever the vessel's arithmetic raises.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are twenty-six
pub fn gibbs_reactor(
    components: &[String],
    feed_n: f64,
    feed_z: &[f64],
    feed_p: Pressure,
    feed_t: ThermodynamicTemperature,
    energy_mode: &str,
    damping_composition: f64,
    max_iterations: f64,
    convergence_tolerance: f64,
    min_iterations: f64,
    reactor_type: Option<&str>,
    catalyst_bulk_density: Option<f64>,
    catalyst_particle_diameter: Option<f64>,
    bed_void_fraction: Option<f64>,
    ghsv: Option<f64>,
    lhsv: Option<f64>,
    number_of_beds: Option<usize>,
    max_operation_pressure: Option<f64>,
    max_operation_temperature: Option<f64>,
    joint_efficiency: Option<f64>,
    corrosion_allowance: Option<f64>,
    armijo_line_search: Option<bool>,
    armijo_c1: Option<f64>,
    armijo_rho: Option<f64>,
    armijo_max_backtracks: Option<usize>,
    adaptive_step_size: Option<bool>,
    regularization: Option<bool>,
    regularization_threshold: Option<f64>,
    regularization_tau: Option<f64>,
) -> Result<GibbsReactorResult> {
    let spec = &model_gen::GIBBS_REACTOR_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "feed_n" => Some(feed_n),
            "damping_composition" => Some(damping_composition),
            "max_iterations" => Some(max_iterations),
            "convergence_tolerance" => Some(convergence_tolerance),
            _ => None,
        },
        &mut warnings,
    )?;

    // **The four dimensioned values cross in SI and the geometry keeps the class's units.**
    // The spec declares millimetres, bar and pascals' worth of a millimetre, and the boundary
    // hands over the SI magnitude - so the conversions are here, and the defaults are the
    // class's initialisers in the class's own units and are not converted.
    let geometry = ReactorMechanicalGeometry {
        reactor_type: reactor_type_name(reactor_type.unwrap_or("fixed-bed"))?,
        catalyst_bulk_density: catalyst_bulk_density.unwrap_or(DEFAULT_CATALYST_BULK_DENSITY),
        catalyst_particle_diameter_mm: catalyst_particle_diameter
            .map_or(DEFAULT_CATALYST_PARTICLE_DIAMETER_MM, |metres| {
                metres * 1000.0
            }),
        bed_void_fraction: bed_void_fraction.unwrap_or(DEFAULT_BED_VOID_FRACTION),
        ghsv_per_hour: ghsv.unwrap_or(DEFAULT_GHSV_PER_HOUR),
        lhsv_per_hour: lhsv.unwrap_or(DEFAULT_LHSV_PER_HOUR),
        number_of_beds: number_of_beds.unwrap_or(DEFAULT_NUMBER_OF_BEDS),
        max_operation_pressure_bara: max_operation_pressure
            .map_or(DEFAULT_MAX_OPERATION_PRESSURE_BARA, |pascals| {
                pascals / 1.0e5
            }),
        max_operation_temperature_k: max_operation_temperature
            .unwrap_or(DEFAULT_MAX_OPERATION_TEMPERATURE_K),
        joint_efficiency: joint_efficiency.unwrap_or(JOINT_EFFICIENCY),
        corrosion_allowance_mm: corrosion_allowance
            .map_or(CORROSION_ALLOWANCE_MM, |metres| metres * 1000.0),
    };

    let feed = Stream::from_pt(components.to_vec(), feed_z.to_vec(), feed_n, feed_p, feed_t)?;
    // **The vessel is read from the feed, before the equilibrium**, which is the class's own
    // surface: `calcDesign` never looks at the outlet.
    let mechanical = reactor_mechanical_design(&feed, &geometry)?;
    let defaults = ReactorSetup::default();
    let setup = ReactorSetup {
        energy_mode: EnergyMode::named(energy_mode),
        damping_composition,
        max_iterations: max_iterations as u32,
        convergence_tolerance,
        min_iterations: min_iterations as u32,
        armijo_line_search: armijo_line_search.unwrap_or(defaults.armijo_line_search),
        armijo_c1: armijo_c1.unwrap_or(defaults.armijo_c1),
        armijo_rho: armijo_rho.unwrap_or(defaults.armijo_rho),
        armijo_max_backtracks: armijo_max_backtracks
            .unwrap_or(defaults.armijo_max_backtracks as usize)
            as u32,
        adaptive_step_size: adaptive_step_size.unwrap_or(defaults.adaptive_step_size),
        regularization: regularization.unwrap_or(defaults.regularization),
        regularization_threshold: regularization_threshold
            .unwrap_or(defaults.regularization_threshold),
        regularization_tau: regularization_tau.unwrap_or(defaults.regularization_tau),
    };
    let (outlet, numbers) = solve_equilibrium(&feed, &setup)?;

    apply_checks(
        spec.derived_checks(),
        |quantity| match quantity {
            "iterations" => Some(f64::from(numbers.iterations)),
            "final_error" => Some(numbers.final_error),
            _ => None,
        },
        &mut warnings,
    )?;

    Ok(GibbsReactorResult::of(
        &outlet,
        &numbers,
        &mechanical,
        warnings,
    ))
}
