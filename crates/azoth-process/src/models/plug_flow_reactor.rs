//! `process.plug_flow_reactor` - the reactor's kernel as a registered id.
//!
//! Spec: `specs/models/process/plug_flow_reactor.toml`. The arithmetic is
//! [`crate::kernels::plug_flow_reactor`]; what is here is the boundary a case and a
//! cross-impl test address, and the translation from the declaration's names to the kernel's.
//!
//! **The stoichiometry is read and the orders are declared.** `KineticReaction` is a *configured*
//! object in NeqSim - `addReactant(name, coefficient, order)` - and a declaration carries no
//! coefficient map, so the reaction id resolves against `data/reactions/stoichiometry.csv`
//! exactly as `process.stirred_tank_reactor`'s does, and `reaction_orders` supplies the third
//! argument `addReactant` takes. Nothing about the chemistry is invented here.

use azoth_core::units::{
    DiffusionCoefficient, HeatTransfer, MassDensity, MolarEnergy, Power, Pressure,
    ThermodynamicTemperature, joules_per_mole, kelvins, pascals, watts,
};
use azoth_core::{AzothError, CalcResult, Result, Warning, apply_checks};
use azoth_reactions::databank::stoichiometry;
use serde::Serialize;

use crate::executor::json::{scalar, warnings as wire_warnings};
use crate::kernels::plug_flow_reactor::{
    EnergyMode, PlugFlowProfile, ReactorNumbers, ReactorSetup, ThermodynamicCoupling,
    plug_flow_reactor as kernel,
};
use crate::model_gen;
use crate::reactor::catalyst_bed::CatalystBed;
use crate::reactor::kinetic_reaction::{KineticReaction, RateBasis, RateType};
use crate::reactor::stepper::Scheme;
use crate::stream::Stream;

/// Result of `process.plug_flow_reactor`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PlugFlowReactorResult {
    /// Product molar flow, mol/s.
    pub product_n: f64,
    /// Product composition, over the feed's species then the reaction's added ones.
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
    /// The key component's conversion over the whole reactor.
    pub conversion: f64,
    /// The inlet pressure less the outlet's.
    #[serde(serialize_with = "scalar")]
    pub pressure_drop: Pressure,
    /// The outlet temperature the march reports.
    #[serde(serialize_with = "scalar")]
    pub outlet_temperature: ThermodynamicTemperature,
    /// The duty an isothermal reactor supplies, zero on the other branches.
    #[serde(serialize_with = "scalar")]
    pub heat_duty: Power,
    /// Every station's axial position, m.
    pub positions: Vec<f64>,
    /// Every station's temperature, K.
    pub temperature_profile: Vec<f64>,
    /// Every station's pressure, Pa.
    pub pressure_profile: Vec<f64>,
    /// Every station's conversion.
    pub conversion_profile: Vec<f64>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

impl CalcResult for PlugFlowReactorResult {
    const CALC_ID: &'static str = "process.plug_flow_reactor";
    const FIELDS: &'static [&'static str] = &[
        "product_n",
        "product_z",
        "product_p",
        "product_t",
        "product_h",
        "conversion",
        "pressure_drop",
        "outlet_temperature",
        "heat_duty",
        "positions",
        "temperature_profile",
        "pressure_profile",
        "conversion_profile",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// The rate type a declaration's name selects.
fn rate_type(name: &str) -> Result<RateType> {
    match name {
        "power_law" => Ok(RateType::PowerLaw),
        "lhhw" => Ok(RateType::Lhhw),
        "equilibrium" => Ok(RateType::Equilibrium),
        other => Err(AzothError::invalid_input(
            "rate_type",
            format!("`{other}` is not one of power_law, lhhw or equilibrium"),
        )),
    }
}

/// The energy mode a declaration's name selects.
fn energy_mode(name: &str) -> Result<EnergyMode> {
    match name {
        "adiabatic" => Ok(EnergyMode::Adiabatic),
        "isothermal" => Ok(EnergyMode::Isothermal),
        "coolant" => Ok(EnergyMode::Coolant),
        other => Err(AzothError::invalid_input(
            "energy_mode",
            format!("`{other}` is not one of adiabatic, isothermal or coolant"),
        )),
    }
}

/// The coupling a declaration's name selects.
fn coupling(name: &str) -> Result<ThermodynamicCoupling> {
    match name {
        "frozen_properties" => Ok(ThermodynamicCoupling::FrozenProperties),
        "fully_coupled" => Ok(ThermodynamicCoupling::FullyCoupled),
        other => Err(AzothError::invalid_input(
            "thermodynamic_coupling",
            format!("`{other}` is not one of frozen_properties or fully_coupled"),
        )),
    }
}

/// March a feed through a plug-flow reactor and report the whole profile.
///
/// # Errors
/// [`AzothError::InvalidInput`] for an unknown mode, rate type or reaction, a non-positive
/// geometry, an order count the reaction's reactants do not match, and whatever the databank,
/// the flash and the rate law refuse.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are twenty-nine
#[allow(clippy::too_many_lines)] // the declaration is the length
pub fn plug_flow_reactor(
    components: &[String],
    feed_n: f64,
    feed_z: &[f64],
    feed_p: Pressure,
    feed_t: ThermodynamicTemperature,
    length: f64,
    diameter: f64,
    number_of_tubes: f64,
    energy_mode_name: &str,
    coolant_temperature: ThermodynamicTemperature,
    overall_heat_transfer_coefficient: HeatTransfer,
    number_of_steps: f64,
    integration_method: &str,
    property_update_frequency: f64,
    coupling_name: &str,
    reaction: &str,
    reaction_orders: &[f64],
    rate_type_name: &str,
    pre_exponential_factor: f64,
    activation_energy: f64,
    temperature_exponent: f64,
    heat_of_reaction: f64,
    catalyst_bulk_density: Option<MassDensity>,
    catalyst_activity_factor: Option<f64>,
    catalyst_particle_diameter: Option<f64>,
    catalyst_void_fraction: Option<f64>,
    catalyst_molecular_diffusivity: Option<DiffusionCoefficient>,
    catalyst_effectiveness_enabled: Option<bool>,
    key_component: Option<String>,
) -> Result<PlugFlowReactorResult> {
    let spec = &model_gen::PLUG_FLOW_REACTOR_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "length" => Some(length),
            "diameter" => Some(diameter),
            "number_of_steps" => Some(number_of_steps),
            "property_update_frequency" => Some(property_update_frequency),
            "catalyst_activity_factor" => catalyst_activity_factor,
            _ => None,
        },
        &mut warnings,
    )?;

    // **The reaction's stoichiometry is the data's, and only the orders are the caller's.** The
    // rows come back in the table's own order, which is also the order `reaction_orders` is in.
    let rows = stoichiometry(reaction)?;
    if rows.is_empty() {
        return Err(AzothError::invalid_input(
            "reaction",
            format!("`{reaction}` is not a reaction the data carries"),
        ));
    }
    let reactants = rows
        .iter()
        .filter(|(_, coefficient)| *coefficient < 0.0)
        .count();
    if reaction_orders.len() != reactants {
        return Err(AzothError::invalid_input(
            "reaction_orders",
            format!(
                "`{reaction}` names {reactants} reactants, so it needs {reactants} orders, \
                 but {} were given",
                reaction_orders.len()
            ),
        ));
    }

    let mut kinetics = KineticReaction::new(reaction);
    kinetics.rate_type = rate_type(rate_type_name)?;
    kinetics.rate_basis = RateBasis::Volume;
    kinetics.pre_exponential_factor = pre_exponential_factor;
    kinetics.activation_energy = activation_energy;
    kinetics.temperature_exponent = temperature_exponent;
    kinetics.heat_of_reaction = heat_of_reaction;

    let mut next_order = reaction_orders.iter();
    for (component, coefficient) in &rows {
        if *coefficient < 0.0 {
            let order = next_order.next().copied().unwrap_or(0.0);
            kinetics.add_reactant(component, *coefficient, order);
        } else {
            kinetics.add_product(component, *coefficient);
        }
    }

    let bed = catalyst_bulk_density.map(|density| CatalystBed {
        bulk_density: density.value,
        activity_factor: catalyst_activity_factor.unwrap_or(1.0),
        particle_diameter: catalyst_particle_diameter.unwrap_or(0.003),
        void_fraction: catalyst_void_fraction.unwrap_or(0.40),
        ..CatalystBed::default()
    });

    let setup = ReactorSetup {
        length,
        diameter,
        number_of_tubes: number_of_tubes.max(1.0) as usize,
        energy_mode: energy_mode(energy_mode_name)?,
        coolant_temperature: coolant_temperature.value,
        overall_heat_transfer_coefficient: overall_heat_transfer_coefficient.value,
        number_of_steps: number_of_steps.max(1.0) as usize,
        integration_method: Scheme::named(integration_method),
        property_update_frequency: property_update_frequency.max(1.0) as usize,
        thermodynamic_coupling: coupling(coupling_name)?,
        catalyst_bed: bed,
        catalyst_effectiveness_enabled: catalyst_effectiveness_enabled.unwrap_or(false),
        catalyst_molecular_diffusivity: catalyst_molecular_diffusivity.map_or(1.0e-5, |d| d.value),
        reactions: vec![kinetics],
        key_component,
    };

    let feed = Stream::from_pt(components.to_vec(), feed_z.to_vec(), feed_n, feed_p, feed_t)?;
    let (product, numbers, profile) = kernel(&feed, &setup)?;

    Ok(PlugFlowReactorResult::of(
        &product, &numbers, &profile, warnings,
    ))
}

impl PlugFlowReactorResult {
    /// The result of one kernel call, from the march's three pieces.
    ///
    /// **The warnings are the caller's**: a case's are the spec's `apply_checks` and a flowsheet's
    /// are the checker's, which report through the envelope rather than through a result.
    #[must_use]
    pub fn of(
        product: &Stream,
        numbers: &ReactorNumbers,
        profile: &PlugFlowProfile,
        warnings: Vec<Warning>,
    ) -> Self {
        Self {
            product_n: product.n,
            product_z: product.z.clone(),
            product_p: product.p,
            product_t: product.t,
            product_h: joules_per_mole(product.h.value),
            conversion: numbers.conversion,
            pressure_drop: pascals(numbers.pressure_drop_bar * 1.0e5),
            outlet_temperature: kelvins(numbers.outlet_temperature),
            heat_duty: watts(numbers.heat_duty),
            positions: profile.positions.clone(),
            temperature_profile: profile.temperatures.clone(),
            // The march carries bara, as `calculateDerivatives` writes its own row, so the report
            // to the declaration's pascals happens here and nowhere else.
            pressure_profile: profile.pressures.iter().map(|p| p * 1.0e5).collect(),
            conversion_profile: profile.conversions.clone(),
            warnings,
        }
    }
}
