//! The `unit_ops.*` id to kernel table, and the typed coercion a parameter needs.
//!
//! **A palette parameter is a `toml::Value` and a kernel wants a `Pressure`.** That gap is the
//! whole of this module's content, and the middleware's own gap table names it: "the `unit_ops.*`
//! id to kernel mapping, with typed parameter coercion". A parameter arrives as text a user typed
//! on a form, so the unit it is in comes from the *declaration* - `UnitOpSpec::parameters`' own
//! `unit` - and the conversion is the vocabulary's, by name, through
//! [`azoth_core::unit_vocab_gen::conversion`]. Nothing here invents a factor.
//!
//! **A kernel and a model are two call shapes of one arithmetic**, which is why the executor calls
//! kernels: a kernel takes and returns a [`Stream`], which is what a flowsheet's connection
//! carries, while a model takes the record field by field so a case and a capture can address it.
//!
//! **Every entry is either a kernel or a refusal by name**, and the refusals are here rather than
//! absent so that a flowsheet naming one gets a reason and not a missing id.

use std::collections::BTreeMap;

use azoth_core::unit_vocab_gen;
use azoth_core::units::{kelvins, pascals, watts};
use azoth_core::{AzothError, Result};

use crate::kernels;
use crate::model_gen;
use crate::stream::Stream;
use crate::unit_op::UnitOpSpec;

/// The convergence a model's own `[algorithm]` block declares: its tolerance and its cap.
///
/// **A `direct` model has no block, and an absorber's convergence is therefore the column's** -
/// which the absorption spec states itself: `AbsorptionColumn` overrides no `run`, so its stage
/// and its sweeps are `DistillationColumn`'s. That is why both the absorber entries and the
/// column entries below read this one spec.
///
/// **The executor reads the algorithm rather than restating it.** A column's temperature
/// tolerance and iteration cap are the class's own defaults, and
/// `specs/models/process/distillation_column.toml` already carries them in the `[algorithm]`
/// block both implementations run - so reading them here is what keeps a flowsheet's column and
/// a case's column the same machine. A parameter the palette does not declare is not a parameter
/// a form may set, which is why these are not in `Parameters`.
fn algorithm_limits(spec: &'static azoth_core::spec::ModelSpec) -> (f64, usize) {
    let algorithm = spec
        .algorithm
        .expect("a column's spec declares the algorithm it runs");
    (algorithm.tolerance, algorithm.max_iterations as usize)
}

/// A palette entry's arithmetic: its inlets in the order its declaration gives them, its
/// parameters, and its outlets in the same order.
///
/// The order is the declaration's and not the kernel's, because a flowsheet connects by port
/// *name* and this is what turns a name into a position.
pub type Kernel = fn(&[Stream], &Parameters<'_>) -> Result<Vec<Stream>>;

/// The parameters a user typed, read against the declaration that says what they mean.
///
/// A parameter's unit is the declaration's, so a value is converted once, here, and every kernel
/// below is in SI. **An undeclared parameter is refused rather than ignored**, which is the same
/// rule `azoth_process::validate` applies to a whole flowsheet - a form and a checker must not
/// disagree about what a machine takes.
pub struct Parameters<'a> {
    spec: &'a UnitOpSpec,
    values: &'a BTreeMap<String, toml::Value>,
}

impl<'a> Parameters<'a> {
    /// Read a palette entry's parameters against its own declaration.
    #[must_use]
    pub fn new(spec: &'a UnitOpSpec, values: &'a BTreeMap<String, toml::Value>) -> Self {
        Self { spec, values }
    }

    /// The declaration of one parameter, refusing a name the entry does not declare.
    fn declared(&self, name: &str) -> Result<&'a crate::unit_op::Param> {
        self.spec.parameters.get(name).ok_or_else(|| {
            AzothError::invalid_input(
                name,
                format!(
                    "`{}` does not declare a parameter called `{name}`",
                    self.spec.id
                ),
            )
        })
    }

    /// The raw value, or `None` when the instance did not set it.
    fn raw(&self, name: &str) -> Result<Option<&'a toml::Value>> {
        self.declared(name)?;
        Ok(self.values.get(name))
    }

    /// One parameter as an SI magnitude, through the unit its declaration names.
    ///
    /// A declared `dimensionless` parameter is already SI; one with a unit goes through the
    /// vocabulary's own conversion, so a form that says `mm` and a kernel that wants metres
    /// agree because the vocabulary's `conversion` is the one conversion site.
    pub fn si(&self, name: &str) -> Result<f64> {
        let param = self.declared(name)?;
        let value = self.require(name)?;
        let number = as_number(name, value)?;
        match param.unit.as_deref() {
            None | Some("dimensionless") => Ok(number),
            Some(unit) => {
                let convert = unit_vocab_gen::conversion(unit).ok_or_else(|| {
                    AzothError::invalid_input(
                        name,
                        format!("`{unit}` is not a unit the vocabulary carries"),
                    )
                })?;
                Ok(convert(number))
            }
        }
    }

    /// One parameter as an SI magnitude, or `None` when it was not set.
    pub fn optional_si(&self, name: &str) -> Result<Option<f64>> {
        if self.raw(name)?.is_none() {
            return Ok(None);
        }
        self.si(name).map(Some)
    }

    /// One parameter as a bare number. A unit-carrying parameter read this way is refused, so a
    /// caller cannot silently skip the conversion.
    pub fn number(&self, name: &str) -> Result<f64> {
        let param = self.declared(name)?;
        if matches!(param.unit.as_deref(), Some(unit) if unit != "dimensionless") {
            return Err(AzothError::invalid_input(
                name,
                format!(
                    "`{name}` carries the unit `{}`, so it is read as an SI magnitude and not as \
                     a bare number",
                    param.unit.as_deref().unwrap_or_default()
                ),
            ));
        }
        as_number(name, self.require(name)?)
    }

    /// One parameter as a bare number, or `None` when it was not set.
    pub fn optional_number(&self, name: &str) -> Result<Option<f64>> {
        if self.raw(name)?.is_none() {
            return Ok(None);
        }
        self.number(name).map(Some)
    }

    /// One parameter as a boolean.
    pub fn flag(&self, name: &str) -> Result<bool> {
        self.declared(name)?;
        match self.require(name)? {
            toml::Value::Boolean(value) => Ok(*value),
            other => Err(wrong_shape(name, "a boolean", other)),
        }
    }

    /// One parameter as a boolean, or `None` when it was not set.
    pub fn optional_flag(&self, name: &str) -> Result<Option<bool>> {
        if self.raw(name)?.is_none() {
            return Ok(None);
        }
        self.flag(name).map(Some)
    }

    /// One parameter as text.
    pub fn text(&self, name: &str) -> Result<String> {
        self.declared(name)?;
        match self.require(name)? {
            toml::Value::String(value) => Ok(value.clone()),
            other => Err(wrong_shape(name, "a string", other)),
        }
    }

    /// One parameter as text, or `None` when it was not set.
    pub fn optional_text(&self, name: &str) -> Result<Option<String>> {
        if self.raw(name)?.is_none() {
            return Ok(None);
        }
        self.text(name).map(Some)
    }

    /// One parameter as a vector of bare numbers.
    ///
    /// A vector parameter's entries are not converted: every vector a palette declares is
    /// dimensionless - a split factor, a composition, an order - and a dimensional vector would
    /// need a unit per entry rather than one for the whole.
    pub fn vector(&self, name: &str) -> Result<Vec<f64>> {
        self.declared(name)?;
        match self.require(name)? {
            toml::Value::Array(items) => items
                .iter()
                .enumerate()
                .map(|(i, item)| {
                    as_number(name, item).map_err(|_| {
                        AzothError::invalid_input(name, format!("`{name}[{i}]` is not a number"))
                    })
                })
                .collect(),
            other => Err(wrong_shape(name, "an array", other)),
        }
    }

    /// A required parameter, refused by name when the instance left it out.
    fn require(&self, name: &str) -> Result<&'a toml::Value> {
        self.values.get(name).ok_or_else(|| {
            AzothError::invalid_input(
                name,
                format!("`{}` needs a value for `{name}`", self.spec.id),
            )
        })
    }
}

/// A `toml::Value` as a number.
fn as_number(name: &str, value: &toml::Value) -> Result<f64> {
    match value {
        toml::Value::Integer(n) => Ok(*n as f64),
        toml::Value::Float(f) => Ok(*f),
        other => Err(wrong_shape(name, "a number", other)),
    }
}

/// The error for a value of the wrong TOML kind.
fn wrong_shape(name: &str, expected: &str, found: &toml::Value) -> AzothError {
    AzothError::invalid_input(
        name,
        format!("`{name}` is {}, not {expected}", found.type_str()),
    )
}

/// **The table.** One entry per palette id that has a kernel, each returning its outlets in the
/// order its declaration names them.
pub const DISPATCH: &[(&str, Kernel)] = &[
    ("unit_ops.absorption_column", absorption_column),
    ("unit_ops.component_splitter", component_splitter),
    ("unit_ops.compressor", compressor),
    ("unit_ops.cooler", cooler),
    ("unit_ops.distillation_column", distillation_column),
    ("unit_ops.ejector", ejector),
    ("unit_ops.expander", expander),
    ("unit_ops.filter", filter),
    ("unit_ops.flare", flare),
    ("unit_ops.gas_scrubber", gas_scrubber),
    ("unit_ops.heat_exchanger", heat_exchanger),
    ("unit_ops.heater", heater),
    ("unit_ops.manifold", manifold),
    ("unit_ops.mixer", mixer),
    ("unit_ops.pipe", pipe),
    ("unit_ops.plug_flow_reactor", plug_flow_reactor),
    ("unit_ops.pump", pump),
    ("unit_ops.separator", separator),
    (
        "unit_ops.shortcut_distillation_column",
        shortcut_distillation_column,
    ),
    ("unit_ops.splitter", splitter),
    ("unit_ops.stirred_tank_reactor", stirred_tank_reactor),
    ("unit_ops.stripping_column", stripping_column),
    ("unit_ops.tank", tank),
    ("unit_ops.three_phase_separator", three_phase_separator),
    ("unit_ops.throttling_valve", throttling_valve),
];

/// **The entries this executor does not run, and what would close each.** A palette entry with no
/// kernel is refused by name here rather than missing, so a flowsheet naming one is told why.
pub const UNRUNNABLE: &[(&str, &str)] = &[
    (
        "unit_ops.simple_absorber",
        "refused on measured evidence: `SimpleAbsorber` is a fixed-point loop over MDEA/CO2 \
         loading and a faithful port needs the amine electrolyte chemistry P8 declined",
    ),
    (
        "unit_ops.gibbs_reactor",
        "deferred with a measurement: the class carries its own Lagrange-multiplier Newton solve \
         and its own species database, so a port composed from P10 would answer differently",
    ),
    (
        "unit_ops.packed_column",
        "its palette entry declares the packing and not the column: `feed_stage`, \
         `number_of_stages`, the two pressures and the two ends are what `PackedColumn.run` reads \
         through `super.run`, and the entry's own notes say it takes `unit_ops.distillation_column`'s \
         whole declaration - which the file does not carry. A form built from this entry could not \
         configure the machine, so the executor refuses it until the declaration matches the notes",
    ),
    (
        "unit_ops.rate_based_packed_column",
        "a second physics - a segment model with film coefficients and an interphase heat balance \
         - carried by the distillation workstream rather than by this one",
    ),
];

/// The kernel a palette id names, or the reason it has none.
///
/// # Errors
/// [`AzothError::InvalidInput`] for an id the palette does not carry, and for one it carries with
/// no kernel - naming what would close it.
pub fn kernel_for(id: &str) -> Result<Kernel> {
    if let Some((_, kernel)) = DISPATCH.iter().find(|(name, _)| *name == id) {
        return Ok(*kernel);
    }
    if let Some((_, reason)) = UNRUNNABLE.iter().find(|(name, _)| *name == id) {
        return Err(AzothError::invalid_input(
            id,
            format!("`{id}` has no kernel: {reason}"),
        ));
    }
    Err(AzothError::invalid_input(
        id,
        format!("`{id}` is not a unit operation the palette declares"),
    ))
}

/// Call one palette entry's kernel.
///
/// # Errors
/// Whatever the entry's kernel raises, and [`AzothError::InvalidInput`] for an entry that has none.
pub fn dispatch(id: &str, inlets: &[Stream], parameters: &Parameters<'_>) -> Result<Vec<Stream>> {
    kernel_for(id)?(inlets, parameters)
}

// The entries, in the table's order. Each is a translation and nothing else: the arithmetic is
// the kernel's, and what is written here is which declared parameter goes on which argument.

fn component_splitter(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    let (first, second) =
        kernels::component_splitter::component_splitter(&inlets[0], &p.vector("split_factors")?)?;
    Ok(vec![first, second])
}

fn compressor(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    Ok(vec![kernels::compressor::compressor(
        &inlets[0],
        pascals(p.si("outlet_pressure")?),
        p.number("isentropic_efficiency")?,
    )?])
}

fn cooler(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    let out = kernels::cooler::cooler(
        &inlets[0],
        p.optional_si("outlet_temperature")?.map(kelvins),
        p.optional_si("duty")?.map(watts),
        p.optional_si("pressure_drop")?.map(pascals),
    )?;
    Ok(vec![out.outlet])
}

fn expander(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    Ok(vec![kernels::expander::expander(
        &inlets[0],
        pascals(p.si("outlet_pressure")?),
        p.number("isentropic_efficiency")?,
    )?])
}

fn filter(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    Ok(vec![
        kernels::filter::filter(&inlets[0], pascals(p.si("pressure_drop")?))?.outlet,
    ])
}

fn flare(inlets: &[Stream], _p: &Parameters<'_>) -> Result<Vec<Stream>> {
    Ok(vec![kernels::flare::flare(&inlets[0])?.0])
}

fn heater(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    let out = kernels::heater::heater(
        &inlets[0],
        p.optional_si("outlet_temperature")?.map(kelvins),
        p.optional_si("duty")?.map(watts),
        p.optional_si("pressure_drop")?.map(pascals),
    )?;
    Ok(vec![out.outlet])
}

fn manifold(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    kernels::manifold::manifold(inlets, &p.vector("split_factors")?)
}

fn mixer(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    Ok(vec![kernels::mixer::mixer(
        inlets,
        p.optional_si("outlet_pressure")?.map(pascals),
    )?])
}

fn pipe(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    let out = kernels::pipe::pipe(
        &inlets[0],
        azoth_core::units::meters(p.si("length")?),
        azoth_core::units::meters(p.si("diameter")?),
        azoth_core::units::meters(p.si("roughness")?),
    )?;
    Ok(vec![out.outlet])
}

fn pump(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    Ok(vec![kernels::pump::pump(
        &inlets[0],
        pascals(p.si("outlet_pressure")?),
        p.number("isentropic_efficiency")?,
    )?])
}

fn separator(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    let (vapour, liquid) = kernels::separator::separator(
        &inlets[0],
        pascals(p.si("pressure_drop")?),
        p.number("gas_in_liquid")?,
        p.optional_si("heat_input")?.map(watts),
    )?;
    Ok(vec![vapour, liquid])
}

fn gas_scrubber(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    let (gas, liquid) = kernels::gas_scrubber::gas_scrubber(
        &inlets[0],
        pascals(p.si("pressure_drop")?),
        p.number("gas_in_liquid")?,
        p.optional_si("heat_input")?.map(watts),
    )?;
    Ok(vec![gas, liquid])
}

fn splitter(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    kernels::splitter::splitter(&inlets[0], &p.vector("split_factors")?)
}

fn tank(inlets: &[Stream], _p: &Parameters<'_>) -> Result<Vec<Stream>> {
    let (gas, liquid) = kernels::tank::tank(inlets)?;
    Ok(vec![gas, liquid])
}

fn throttling_valve(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    Ok(vec![kernels::throttling_valve::throttling_valve(
        &inlets[0],
        pascals(p.si("outlet_pressure")?),
    )?])
}

fn heat_exchanger(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    let (hot, cold) = kernels::heat_exchanger::heat_exchanger(
        &inlets[0],
        &inlets[1],
        p.optional_si("ua")?
            .map(azoth_core::units::watts_per_kelvin),
        &p.text("flow_arrangement")?,
        p.optional_si("hot_outlet_temperature")?.map(kelvins),
        p.optional_si("cold_outlet_temperature")?.map(kelvins),
    )?;
    Ok(vec![hot, cold])
}

fn ejector(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    Ok(vec![kernels::ejector::ejector(
        &inlets[0],
        &inlets[1],
        kernels::EjectorSetup {
            discharge_pressure: pascals(p.si("discharge_pressure")?),
            motive_nozzle_efficiency: p.number("motive_nozzle_efficiency")?,
            suction_nozzle_efficiency: p.number("suction_nozzle_efficiency")?,
            mixing_efficiency: p.number("mixing_efficiency")?,
            diffuser_efficiency: p.number("diffuser_efficiency")?,
        },
    )?])
}

fn three_phase_separator(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    let (vapour, light, heavy) = kernels::three_phase_separator::three_phase_separator(
        &inlets[0],
        pascals(p.si("pressure_drop")?),
        p.optional_si("heat_input")?.map(watts),
        kernels::Entrainment {
            gas_in_aqueous: p.number("gas_in_aqueous")?,
            gas_in_oil: p.number("gas_in_oil")?,
            oil_in_aqueous: p.number("oil_in_aqueous")?,
            oil_in_gas: p.number("oil_in_gas")?,
            aqueous_in_gas: p.number("aqueous_in_gas")?,
            aqueous_in_oil: p.number("aqueous_in_oil")?,
        },
    )?;
    Ok(vec![vapour, light, heavy])
}

fn stirred_tank_reactor(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    let (product, _duty) = kernels::stirred_tank_reactor::stirred_tank_reactor(
        &inlets[0],
        &kernels::ReactorSetup {
            reaction: p.text("reaction")?,
            limiting_reactant: p.text("limiting_reactant")?,
            conversion: p.number("conversion")?,
            isothermal: p.flag("isothermal")?,
            reactor_temperature: p.optional_si("reactor_temperature")?.map(kelvins),
            reactor_pressure: p.optional_si("reactor_pressure")?.map(pascals),
            pressure_drop: pascals(p.optional_si("pressure_drop")?.unwrap_or(0.0)),
        },
    )?;
    Ok(vec![product])
}

fn plug_flow_reactor(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    // The order vectors are declared over the reactants the reaction data names, and the
    // stoichiometry is read from that data rather than declared - the same route the model takes.
    let orders = p.vector("reaction_orders")?;
    let mut kinetics = crate::reactor::kinetic_reaction::KineticReaction::new(p.text("reaction")?);
    kinetics.rate_type = match p.text("rate_type")?.as_str() {
        "power_law" => crate::reactor::kinetic_reaction::RateType::PowerLaw,
        "lhhw" => crate::reactor::kinetic_reaction::RateType::Lhhw,
        other => {
            return Err(AzothError::invalid_input(
                "rate_type",
                format!("`{other}` is not power_law, lhhw or equilibrium"),
            ));
        }
    };
    kinetics.pre_exponential_factor = p.number("pre_exponential_factor")?;
    kinetics.activation_energy = p.si("activation_energy")?;
    kinetics.temperature_exponent = p.number("temperature_exponent")?;
    kinetics.heat_of_reaction = p.si("heat_of_reaction")?;
    let rows = azoth_reactions::databank::stoichiometry(&kinetics.name)?;
    if rows.is_empty() {
        return Err(AzothError::invalid_input(
            "reaction",
            format!("`{}` is not a reaction the data carries", kinetics.name),
        ));
    }
    let mut next_order = orders.iter();
    for (component, coefficient) in &rows {
        if *coefficient < 0.0 {
            kinetics.add_reactant(
                component,
                *coefficient,
                next_order.next().copied().unwrap_or(0.0),
            );
        } else {
            kinetics.add_product(component, *coefficient);
        }
    }
    let bed = p.optional_si("catalyst_bulk_density")?.map(|density| {
        crate::reactor::catalyst_bed::CatalystBed {
            bulk_density: density,
            activity_factor: p
                .optional_number("catalyst_activity_factor")
                .ok()
                .flatten()
                .unwrap_or(1.0),
            particle_diameter: p
                .optional_si("catalyst_particle_diameter")
                .ok()
                .flatten()
                .unwrap_or(0.003),
            void_fraction: p
                .optional_number("catalyst_void_fraction")
                .ok()
                .flatten()
                .unwrap_or(0.40),
            ..crate::reactor::catalyst_bed::CatalystBed::default()
        }
    });
    let (product, _numbers, _profile) = kernels::plug_flow_reactor::plug_flow_reactor(
        &inlets[0],
        &kernels::plug_flow_reactor::ReactorSetup {
            length: p.si("length")?,
            diameter: p.si("diameter")?,
            number_of_tubes: p.number("number_of_tubes")? as usize,
            energy_mode: match p.text("energy_mode")?.as_str() {
                "adiabatic" => kernels::plug_flow_reactor::EnergyMode::Adiabatic,
                "isothermal" => kernels::plug_flow_reactor::EnergyMode::Isothermal,
                "coolant" => kernels::plug_flow_reactor::EnergyMode::Coolant,
                other => {
                    return Err(AzothError::invalid_input(
                        "energy_mode",
                        format!("`{other}` is not adiabatic, isothermal or coolant"),
                    ));
                }
            },
            coolant_temperature: p.si("coolant_temperature")?,
            overall_heat_transfer_coefficient: p.si("overall_heat_transfer_coefficient")?,
            number_of_steps: p.number("number_of_steps")? as usize,
            integration_method: crate::reactor::stepper::Scheme::named(
                &p.text("integration_method")?,
            ),
            property_update_frequency: p.number("property_update_frequency")? as usize,
            thermodynamic_coupling: match p.text("thermodynamic_coupling")?.as_str() {
                "frozen_properties" => {
                    kernels::plug_flow_reactor::ThermodynamicCoupling::FrozenProperties
                }
                "fully_coupled" => kernels::plug_flow_reactor::ThermodynamicCoupling::FullyCoupled,
                other => {
                    return Err(AzothError::invalid_input(
                        "thermodynamic_coupling",
                        format!("`{other}` is not frozen_properties or fully_coupled"),
                    ));
                }
            },
            catalyst_bed: bed,
            catalyst_effectiveness_enabled: p
                .optional_flag("catalyst_effectiveness_enabled")?
                .unwrap_or(false),
            catalyst_molecular_diffusivity: p
                .optional_si("catalyst_molecular_diffusivity")?
                .unwrap_or(1.0e-5),
            reactions: vec![kinetics],
            key_component: p.optional_text("key_component")?,
        },
    )?;
    Ok(vec![product])
}

fn absorption_column(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    let out = kernels::absorption_column::absorption_column(
        &kernels::absorption_column::AbsorberSetup {
            gas: inlets[0].clone(),
            solvent: inlets[1].clone(),
            number_of_stages: p.number("number_of_stages")? as usize,
            top_pressure: pascals(p.si("top_pressure")?),
            bottom_pressure: pascals(p.si("bottom_pressure")?),
            tray_temperatures: None,
            temperature_tolerance: algorithm_limits(&model_gen::DISTILLATION_COLUMN_SPEC).0,
            max_iterations: algorithm_limits(&model_gen::DISTILLATION_COLUMN_SPEC).1,
            solver_type: kernels::distillation_column::SolverType::DirectSubstitution,
        },
    )?;
    Ok(vec![out.gas_out, out.liquid_out])
}

fn stripping_column(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    // **A stripper is an absorber with its two feeds renamed**, which is the class's own
    // statement: absorption and stripping are the same counter-current stage equations.
    let out = kernels::absorption_column::absorption_column(
        &kernels::absorption_column::AbsorberSetup {
            gas: inlets[0].clone(),
            solvent: inlets[1].clone(),
            number_of_stages: p.number("number_of_stages")? as usize,
            top_pressure: pascals(p.si("top_pressure")?),
            bottom_pressure: pascals(p.si("bottom_pressure")?),
            tray_temperatures: None,
            temperature_tolerance: algorithm_limits(&model_gen::DISTILLATION_COLUMN_SPEC).0,
            max_iterations: algorithm_limits(&model_gen::DISTILLATION_COLUMN_SPEC).1,
            solver_type: kernels::distillation_column::SolverType::DirectSubstitution,
        },
    )?;
    Ok(vec![out.gas_out, out.liquid_out])
}

fn distillation_column(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    let out = kernels::distillation_column::distillation_column(&column_setup(inlets, p)?)?;
    Ok(vec![out.distillate, out.bottoms])
}

/// The column a distillation, packed or stripping entry configures.
fn column_setup(inlets: &[Stream], p: &Parameters<'_>) -> Result<kernels::ColumnSetup> {
    Ok(kernels::ColumnSetup {
        feed: inlets[0].clone(),
        feed_stage: p.number("feed_stage")? as usize,
        number_of_stages: p.number("number_of_stages")? as usize,
        has_reboiler: p.flag("has_reboiler")?,
        has_condenser: p.flag("has_condenser")?,
        top_pressure: pascals(p.si("top_pressure")?),
        bottom_pressure: pascals(p.si("bottom_pressure")?),
        condenser_temperature: p.optional_si("condenser_temperature")?.map(kelvins),
        reboiler_temperature: p.optional_si("reboiler_temperature")?.map(kelvins),
        temperature_tolerance: algorithm_limits(&model_gen::DISTILLATION_COLUMN_SPEC).0,
        max_iterations: algorithm_limits(&model_gen::DISTILLATION_COLUMN_SPEC).1,
        top_specification: None,
        bottom_specification: None,
        top_feed: None,
        tray_temperatures: None,
        solver_type: kernels::distillation_column::SolverType::DirectSubstitution,
    })
}

fn shortcut_distillation_column(inlets: &[Stream], p: &Parameters<'_>) -> Result<Vec<Stream>> {
    let out = kernels::shortcut_distillation_column::shortcut_distillation_column(
        &inlets[0],
        &p.text("light_key")?,
        &p.text("heavy_key")?,
        p.number("light_key_recovery_distillate")?,
        p.number("heavy_key_recovery_bottoms")?,
        p.number("reflux_ratio_multiplier")?,
        p.optional_si("condenser_pressure")?.map(pascals),
        p.optional_si("reboiler_pressure")?.map(pascals),
    )?;
    Ok(vec![out.distillate, out.bottoms])
}
