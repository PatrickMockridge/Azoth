//! `unit_ops.plug_flow_reactor` - a fixed bed marched along its own length.
//!
//! ```text
//! dF_i/dz = A_tot * sum_j nu_ij * r_j      mol/(s*m)
//! dT/dz   = (-Q_gen*A_tot + HT) / sum(F_i Cp_i)   K/m
//! dP/dz   = -Ergun(u, rho, mu) / 1e5       bar/m
//! ```
//!
//! Spec: `specs/models/process/plug_flow_reactor.toml`.
//!
//! **The state vector is the class's, `[F1..Fn, T, P]`, with the pressure in bara.** The
//! derivative's own pressure row converts Pa/m to bar/m rather than the state carrying pascals,
//! so the bara is the marcher's unit and the pascals appear only at the boundary.
//!
//! **Two volumes are read in one derivative, and they are different volumes.** The rate law's
//! concentration is `x_i * getPhase(0).getDensity("mol/m3")`, which is the **Peneloux-corrected**
//! molar density; the superficial velocity is
//! `getVolume("m3") / getNumberOfMoles() * getTotalNumberOfMoles() / A`, and `getVolume` is the
//! **untranslated cubic's**. Measured on the capture's own feed state the two differ by
//! `4.2e-4` relative (`2.782061` against `2.783235` kg/m³), which is small and is not zero, so
//! the port keeps both rather than picking one - the same two-density shape `process.pipe`
//! already carries.
//!
//! **Properties are frozen by default.** `thermodynamicCoupling` is `FROZEN_PROPERTIES`, so the
//! system state is rebuilt only every `property_update_frequency` steps (and never for an RK4
//! sub-state), and the rate law reads a stale temperature and composition in between. That is
//! the class's behaviour, and it is why its two integration schemes agree to thirteen digits on
//! its own default: with the rates frozen for ten steps at a time, the march is near enough to a
//! straight line that RK4 and Euler land together. `FULLY_COUPLED` re-flashes for every residual
//! evaluation instead, and the class pays for it exactly as this does.

use std::cell::RefCell;

use azoth_core::units::pascals;
use azoth_core::{AzothError, Result};
use azoth_eos::viscosity;

use crate::reactor::catalyst_bed::CatalystBed;
use crate::reactor::kinetic_reaction::{KineticReaction, RateBasis};
use crate::reactor::stepper::{Scheme, march};
use crate::stream::Stream;

/// The small mole count `updateSystemState` floors a component at, mol/s.
const MINIMUM_MOLES: f64 = 1.0e-30;

/// The pressure floor the loop holds, bara.
const MINIMUM_PRESSURE: f64 = 0.1;

/// The temperature floor `calculateDerivatives` reads with, K.
const MINIMUM_TEMPERATURE: f64 = 1.0;

/// The empty tube's friction factor, where no bed is configured.
const EMPTY_TUBE_FRICTION_FACTOR: f64 = 0.01;

/// How the reactor's temperature changes along its length.
///
/// The class's `EnergyMode`. **`Adiabatic` is its default**, so a reactor with no stated mode
/// heats or cools from the reaction alone.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnergyMode {
    /// No heat transfer: the temperature moves with the reaction's enthalpy and nothing else.
    Adiabatic,
    /// Held at the inlet's temperature, with the duty reported rather than applied.
    Isothermal,
    /// Heat exchange with a coolant at `Q = U * perimeter * (T_coolant - T)` per unit length.
    Coolant,
}

/// Whether the state is rebuilt for every residual evaluation or reused between updates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThermodynamicCoupling {
    /// Reuse the last flash between property updates, the class's default.
    FrozenProperties,
    /// Rebuild and flash for every residual evaluation.
    FullyCoupled,
}

/// What a plug-flow reactor is configured with.
///
/// Every field is a class setter with the class's own default, so `PlugFlowReactorSetup::default`
/// is the machine a bare `new PlugFlowReactor(name, inlet)` is.
#[derive(Debug, Clone)]
pub struct ReactorSetup {
    /// Tube length, m, the class's default `5.0`.
    pub length: f64,
    /// Tube diameter, m, the class's default `0.10`.
    pub diameter: f64,
    /// How many tubes the flow is split across, the class's default `1`.
    pub number_of_tubes: usize,
    /// How the temperature moves.
    pub energy_mode: EnergyMode,
    /// The coolant's temperature, K, the class's default `298.15`.
    pub coolant_temperature: f64,
    /// The overall heat transfer coefficient, W/(m²·K), the class's default `50.0`.
    pub overall_heat_transfer_coefficient: f64,
    /// How many steps the march takes, the class's default `100`.
    pub number_of_steps: usize,
    /// Which scheme advances it, RK4 by default.
    pub integration_method: Scheme,
    /// How many steps between property updates, the class's default `10`.
    pub property_update_frequency: usize,
    /// Whether the state is rebuilt per residual evaluation.
    pub thermodynamic_coupling: ThermodynamicCoupling,
    /// The bed, or `None` for an empty tube - **the class's default, and the pressure row's
    /// branch**.
    pub catalyst_bed: Option<CatalystBed>,
    /// Whether the pellet's effectiveness factor scales the rate. Off by default.
    pub catalyst_effectiveness_enabled: bool,
    /// The molecular diffusivity the effectiveness factor reads, m²/s, the class's `1.0e-5`.
    pub catalyst_molecular_diffusivity: f64,
    /// The reactions, in the order they are summed.
    pub reactions: Vec<KineticReaction>,
    /// Which component the reported conversion is of. The class defaults to the first reactant
    /// of the first reaction.
    pub key_component: Option<String>,
}

impl Default for ReactorSetup {
    fn default() -> Self {
        Self {
            length: 5.0,
            diameter: 0.10,
            number_of_tubes: 1,
            energy_mode: EnergyMode::Adiabatic,
            coolant_temperature: 298.15,
            overall_heat_transfer_coefficient: 50.0,
            number_of_steps: 100,
            integration_method: Scheme::Rk4,
            property_update_frequency: 10,
            thermodynamic_coupling: ThermodynamicCoupling::FrozenProperties,
            catalyst_bed: None,
            catalyst_effectiveness_enabled: false,
            catalyst_molecular_diffusivity: 1.0e-5,
            reactions: Vec::new(),
            key_component: None,
        }
    }
}

/// Every station's record: **the answer is a curve, not an outlet row.**
#[derive(Debug, Clone, PartialEq)]
pub struct PlugFlowProfile {
    /// The component order the flows are in.
    pub components: Vec<String>,
    /// Axial position, m.
    pub positions: Vec<f64>,
    /// Temperature, K.
    pub temperatures: Vec<f64>,
    /// Pressure, bara.
    pub pressures: Vec<f64>,
    /// The key component's conversion, dimensionless.
    pub conversions: Vec<f64>,
    /// The summed absolute volumetric rate, mol/(m³·s).
    pub rates: Vec<f64>,
    /// Molar flows, mol/s, one row per station and one column per component.
    pub flows: Vec<Vec<f64>>,
}

/// What the reactor reports beside its outlet.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactorNumbers {
    /// The key component's conversion over the whole reactor.
    pub conversion: f64,
    /// The inlet pressure less the outlet's, bar.
    pub pressure_drop_bar: f64,
    /// The outlet temperature, K.
    pub outlet_temperature: f64,
    /// The reactor volume divided by the volumetric flow, s.
    pub residence_time: f64,
    /// The duty an isothermal reactor supplies, W. Zero on the adiabatic branch.
    pub heat_duty: f64,
    /// Which component the conversion is of.
    pub key_component: String,
}

/// The flashed property state the derivatives read.
struct PropertyState {
    stream: Stream,
    molar_mass: f64,
    corrected_density: f64,
    bare_density: f64,
    cp: f64,
    viscosity: f64,
    composition: Vec<f64>,
}

impl PropertyState {
    /// Build the state from a temperature, a pressure in bara and a mole vector.
    fn build(
        components: &[String],
        flows: &[f64],
        temperature: f64,
        pressure_bar: f64,
    ) -> Result<Self> {
        // `updateSystemState` floors each component at `1e-30` mol, so a spent reactant is
        // present but absent from the rate.
        let floored: Vec<f64> = flows.iter().map(|f| f.max(MINIMUM_MOLES)).collect();
        let total: f64 = floored.iter().sum();
        let composition: Vec<f64> = floored.iter().map(|f| f / total).collect();
        let stream = Stream::from_pt(
            components.to_vec(),
            composition.clone(),
            total,
            pascals(pressure_bar * 1.0e5),
            azoth_core::units::kelvins(temperature),
        )?;
        let (mixture, _) = stream.mixture()?;
        Ok(Self {
            molar_mass: stream.molar_mass()?.value,
            corrected_density: stream.corrected_density()?,
            bare_density: stream.density()?,
            cp: stream.molar_heat_capacity()?,
            viscosity: viscosity::viscosity(&mixture, stream.t, stream.p, &composition)?
                .mu
                .value,
            composition,
            stream,
        })
    }

    /// A component's molar concentration, mol/m³: `x_i` times the **corrected** molar density.
    fn concentration(&self, index: usize) -> f64 {
        self.composition[index] * self.corrected_density / self.molar_mass
    }

    /// The superficial velocity through the tubes, m/s.
    ///
    /// `V / n * n_total / A`, where `V` is the system's extensive volume on the **untranslated**
    /// cubic. With the flows as the system's moles the two mole counts cancel and this is
    /// `n_total * v_molar / A`.
    fn superficial_velocity(&self, total_area: f64) -> f64 {
        let n_total = self.stream.n;
        let molar_volume = self.molar_mass / self.bare_density;
        n_total * molar_volume / total_area
    }
}

/// March a feed through a plug-flow reactor.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the geometry or the step count is not positive, if no
///   reaction is configured, or if a state the march reaches cannot be flashed.
pub fn plug_flow_reactor(
    feed: &Stream,
    setup: &ReactorSetup,
) -> Result<(Stream, ReactorNumbers, PlugFlowProfile)> {
    if setup.reactions.is_empty() {
        return Err(AzothError::invalid_input(
            "reactions",
            "a plug-flow reactor with no reaction is a tube, not a reactor",
        ));
    }
    if setup.length <= 0.0 || setup.diameter <= 0.0 || setup.number_of_tubes == 0 {
        return Err(AzothError::invalid_input(
            "geometry",
            "a tube needs a positive length, a positive diameter and at least one tube",
        ));
    }
    let steps = setup.number_of_steps.max(1);
    let update_every = setup.property_update_frequency.max(1);

    // **`ensureProductComponentsExist`.** The species a reaction names but the feed does not
    // carry are added at `1e-20` mol and the fluid is re-initialised, in the reactions' own
    // stoichiometry order.
    let mut components = feed.components.clone();
    let mut flows: Vec<f64> = feed.z.iter().map(|z| z * feed.n).collect();
    for reaction in &setup.reactions {
        for (species, _) in &reaction.stoichiometry {
            if !components.iter().any(|c| c == species) {
                components.push(species.clone());
                flows.push(1.0e-20);
            }
        }
    }
    let n_components = components.len();

    // **The key component defaults to the first reactant of the first reaction**, which is what
    // the class does when none is named.
    let key_index = match &setup.key_component {
        Some(name) => components.iter().position(|c| c == name).ok_or_else(|| {
            AzothError::invalid_input(
                "key_component",
                format!("`{name}` is not in the reactor's species set"),
            )
        })?,
        None => setup.reactions[0]
            .stoichiometry
            .iter()
            .find(|(_, coefficient)| *coefficient < 0.0)
            .and_then(|(species, _)| components.iter().position(|c| c == species))
            .ok_or_else(|| {
                AzothError::invalid_input(
                    "key_component",
                    "no reaction names a reactant, so no conversion can be reported",
                )
            })?,
    };
    let key_component = components[key_index].clone();
    let initial_key_moles = flows[key_index];

    let tube_area = std::f64::consts::PI * setup.diameter * setup.diameter / 4.0;
    let total_area = tube_area * setup.number_of_tubes as f64;
    let perimeter = std::f64::consts::PI * setup.diameter;
    let dz = setup.length / steps as f64;

    let mut temperature = feed.t.value;
    let inlet_temperature = temperature;
    let mut pressure_bar = feed.p.value / 1.0e5;

    let mut state: Vec<f64> = Vec::with_capacity(n_components + 2);
    state.extend_from_slice(&flows);
    state.push(temperature);
    state.push(pressure_bar);

    let properties = RefCell::new(PropertyState::build(
        &components,
        &flows,
        temperature,
        pressure_bar,
    )?);

    let mut profile = PlugFlowProfile {
        components: components.clone(),
        positions: vec![0.0],
        temperatures: vec![temperature],
        pressures: vec![pressure_bar],
        conversions: vec![0.0],
        rates: vec![0.0],
        flows: vec![flows.clone()],
    };

    let mut step_index = 0usize;
    march(
        &mut state,
        steps,
        dz,
        setup.integration_method,
        |trial| {
            derivatives(
                trial,
                &components,
                setup,
                &properties,
                total_area,
                perimeter,
                n_components,
            )
        },
        |settled| {
            // The class's own post-step discipline, in its order.
            for value in settled.iter_mut().take(n_components) {
                *value = value.max(0.0);
            }
            settled[n_components + 1] = settled[n_components + 1].max(MINIMUM_PRESSURE);
            if setup.energy_mode == EnergyMode::Isothermal {
                settled[n_components] = inlet_temperature;
            }
            step_index += 1;
            temperature = settled[n_components];
            pressure_bar = settled[n_components + 1];

            // **The refresh is the class's, including where it does not happen.** Its loop
            // re-flashes at the top of step `s` when `s % frequency == 0 && s > 0`, so a
            // hundred-step march refreshes after steps ten to ninety and **never after the
            // last one** - which is why `calculateResidenceTime`, called before the class's
            // final `updateSystemState`, reads a state at the *ninety*-step temperature rather
            // than the outlet's. Measured: `3.3649` s against the `3.2132` the outlet state
            // gives. The port reproduces that rather than tidying it.
            if setup.thermodynamic_coupling == ThermodynamicCoupling::FrozenProperties
                && step_index % update_every == 0
                && step_index < steps
            {
                if let Ok(refreshed) = PropertyState::build(
                    &components,
                    &settled[..n_components],
                    temperature,
                    pressure_bar,
                ) {
                    properties.replace(refreshed);
                }
            }

            let conversion = if initial_key_moles > 1.0e-30 {
                (1.0 - settled[key_index] / initial_key_moles).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let total_rate = total_reaction_rate(setup, &properties.borrow());
            profile.positions.push(step_index as f64 * dz);
            profile.temperatures.push(temperature);
            profile.pressures.push(pressure_bar);
            profile.conversions.push(conversion);
            profile.rates.push(total_rate);
            profile.flows.push(settled[..n_components].to_vec());
        },
    );

    let final_flows = &state[..n_components];
    let total: f64 = final_flows.iter().sum();
    let final_z: Vec<f64> = final_flows.iter().map(|f| f / total).collect();
    let outlet = Stream::from_pt(
        components.clone(),
        final_z,
        total,
        pascals(pressure_bar * 1.0e5),
        azoth_core::units::kelvins(temperature),
    )?;

    let conversion = if initial_key_moles > 1.0e-30 {
        (1.0 - state[key_index] / initial_key_moles).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let residence_time = {
        let volumetric = properties.borrow().superficial_velocity(total_area) * total_area;
        if volumetric > 0.0 {
            total_area * setup.length / volumetric
        } else {
            0.0
        }
    };

    let numbers = ReactorNumbers {
        conversion,
        pressure_drop_bar: feed.p.value / 1.0e5 - pressure_bar,
        outlet_temperature: temperature,
        residence_time,
        heat_duty: if setup.energy_mode == EnergyMode::Isothermal {
            isothermal_heat_duty(setup, &components, final_flows, feed, &initial_key_moles)
        } else {
            0.0
        },
        key_component,
    };
    Ok((outlet, numbers, profile))
}

/// The state vector's derivative, `[dF1/dz..dFn/dz, dT/dz, dP/dz]`.
fn derivatives(
    trial: &[f64],
    components: &[String],
    setup: &ReactorSetup,
    properties: &RefCell<PropertyState>,
    total_area: f64,
    perimeter: f64,
    n_components: usize,
) -> Vec<f64> {
    let mut derivs = vec![0.0; n_components + 2];
    if setup.thermodynamic_coupling == ThermodynamicCoupling::FullyCoupled {
        let temperature = trial[n_components].max(MINIMUM_TEMPERATURE);
        let pressure = trial[n_components + 1].max(MINIMUM_PRESSURE);
        if let Ok(refreshed) =
            PropertyState::build(components, &trial[..n_components], temperature, pressure)
        {
            properties.replace(refreshed);
        }
    }
    let properties = &properties.borrow();

    let temperature = trial[n_components].max(MINIMUM_TEMPERATURE);
    let total_mol_flow: f64 = trial[..n_components].iter().map(|f| f.max(0.0)).sum();
    if total_mol_flow < MINIMUM_MOLES {
        return derivs;
    }

    let mut heat_generation = 0.0;
    for reaction in &setup.reactions {
        let concentration = |name: &str| {
            components
                .iter()
                .position(|c| c == name)
                .map_or(0.0, |i| properties.concentration(i))
        };
        let Ok(mut rate) = reaction.rate(properties.stream.t.value, concentration) else {
            continue;
        };
        if let Some(bed) = setup.catalyst_bed {
            rate *= bed.activity_factor;
            if setup.catalyst_effectiveness_enabled {
                let diffusivity = bed.effective_diffusivity(setup.catalyst_molecular_diffusivity);
                let thiele = bed.thiele_modulus(reaction.rate_constant(temperature), diffusivity);
                rate *= bed.effectiveness_factor(thiele);
            }
        }
        let volumetric = convert_to_volumetric(rate, reaction, setup.catalyst_bed);
        for (species, coefficient) in &reaction.stoichiometry {
            if let Some(i) = components.iter().position(|c| c == species) {
                derivs[i] += total_area * coefficient * volumetric;
            }
        }
        heat_generation += volumetric * reaction.heat_of_reaction;
    }

    let capacity_flow = properties.cp * total_mol_flow;
    match setup.energy_mode {
        EnergyMode::Adiabatic => {
            if capacity_flow.abs() > 1.0e-20 {
                derivs[n_components] = -heat_generation * total_area / capacity_flow;
            }
        }
        EnergyMode::Coolant => {
            let transfer = setup.overall_heat_transfer_coefficient
                * perimeter
                * setup.number_of_tubes as f64
                * (setup.coolant_temperature - temperature);
            if capacity_flow.abs() > 1.0e-20 {
                derivs[n_components] = (-heat_generation * total_area + transfer) / capacity_flow;
            }
        }
        // Isothermal writes nothing: the temperature is overridden after the step.
        EnergyMode::Isothermal => {}
    }

    let velocity = properties.superficial_velocity(total_area);
    let drop_per_metre = match setup.catalyst_bed {
        Some(bed) => bed.pressure_drop(velocity, properties.bare_density, properties.viscosity),
        // **The empty tube's branch, and the class's default.** A bare `PlugFlowReactor` has no
        // bed, so this is what a pressure row does unless one is set - four orders of magnitude
        // apart from the Ergun drop on the capture's gas.
        None => {
            EMPTY_TUBE_FRICTION_FACTOR * properties.bare_density * velocity * velocity
                / (2.0 * setup.diameter)
        }
    };
    derivs[n_components + 1] = -drop_per_metre / 1.0e5;
    derivs
}

/// A rate in its own basis, as mol/(m³ reactor · s).
fn convert_to_volumetric(rate: f64, reaction: &KineticReaction, bed: Option<CatalystBed>) -> f64 {
    match (reaction.rate_basis, bed) {
        (RateBasis::CatalystMass, Some(bed)) => rate * bed.bulk_density,
        (RateBasis::CatalystArea, Some(bed)) => rate * bed.specific_surface_area * bed.bulk_density,
        _ => rate,
    }
}

/// The summed absolute volumetric rate the profile stores.
fn total_reaction_rate(setup: &ReactorSetup, properties: &PropertyState) -> f64 {
    let mut total = 0.0;
    for reaction in &setup.reactions {
        let components = &properties.stream.components;
        let concentration = |name: &str| {
            components
                .iter()
                .position(|c| c == name)
                .map_or(0.0, |i| properties.concentration(i))
        };
        let Ok(mut rate) = reaction.rate(properties.stream.t.value, concentration) else {
            continue;
        };
        if let Some(bed) = setup.catalyst_bed {
            rate *= bed.activity_factor;
            if setup.catalyst_effectiveness_enabled {
                let diffusivity = bed.effective_diffusivity(setup.catalyst_molecular_diffusivity);
                let thiele = bed.thiele_modulus(
                    reaction.rate_constant(properties.stream.t.value),
                    diffusivity,
                );
                rate *= bed.effectiveness_factor(thiele);
            }
        }
        total += convert_to_volumetric(rate, reaction, setup.catalyst_bed).abs();
    }
    total
}

/// An isothermal reactor's duty, W: `-ΔH` per mole reacted of each reaction's first reactant.
fn isothermal_heat_duty(
    setup: &ReactorSetup,
    components: &[String],
    final_flows: &[f64],
    feed: &Stream,
    initial_key_moles: &f64,
) -> f64 {
    let _ = initial_key_moles;
    let mut duty = 0.0;
    for reaction in &setup.reactions {
        let first_reactant = reaction
            .stoichiometry
            .iter()
            .find(|(_, coefficient)| *coefficient < 0.0);
        if let Some((species, coefficient)) = first_reactant {
            let Some(index) = components.iter().position(|c| c == species) else {
                continue;
            };
            let initial = feed
                .components
                .iter()
                .position(|c| c == species)
                .map_or(0.0, |i| feed.z[i] * feed.n);
            let reacted = initial - final_flows[index];
            duty += -reaction.heat_of_reaction * reacted / coefficient.abs();
        }
    }
    duty
}
