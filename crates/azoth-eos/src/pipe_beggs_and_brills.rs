//! `eos.pipe_beggs_and_brills` - NeqSim's `PipeBeggsAndBrills`, isothermal march.
//!
//! A pipeline is divided into increments; each is flashed at the inlet temperature and its
//! own pressure, its phase properties feed the Beggs-and-Brill segment correlation, and the
//! resulting drop sets the next segment's pressure.

use azoth_core::units::{Length, Pressure, ThermodynamicTemperature, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::beggs_brill_segment::beggs_brill_segment;
use crate::databank;
use crate::mixture::{Component, Mixture};
use crate::model_gen;
use crate::phase_transport::{PhaseKind, viscosity_only};
use crate::results::PipeBeggsAndBrillsResult;

/// The universal gas constant, as the port's volume calculations carry it.
const MOLAR_GAS_CONSTANT: f64 = 8.314_462_618_153_24;
/// SRK's covolume coefficient, for the single-phase phase-type test.
const SRK_OMEGA_B: f64 = 0.086_64;

struct PhaseProperties {
    gas_flow: f64,
    liquid_flow: f64,
    gas_density: f64,
    liquid_density: f64,
    gas_viscosity: f64,
    liquid_viscosity: f64,
}

/// The isothermal pressure profile of a Beggs-and-Brills pipeline.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's hard range.
/// * [`azoth_core::AzothError::InvalidInput`] if a segment's outlet pressure is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, meters, moles_per_second, pascals};
/// use azoth_eos::pipe_beggs_and_brills::pipe_beggs_and_brills;
///
/// let r = pipe_beggs_and_brills(
///     &["methane".to_owned(), "CO2".to_owned(), "water".to_owned()],
///     80.0,
///     &[0.85, 0.05, 0.10],
///     pascals(1.0e7),
///     kelvins(333.15),
///     meters(1000.0),
///     meters(0.1),
///     meters(0.0),
///     10,
///     meters(1.0e-5),
///     Some("srk"),
/// )?;
/// assert_eq!(r.pressure_profile.len(), 11);
/// assert_eq!(r.temperature_profile.len(), 1);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn pipe_beggs_and_brills(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
    length: Length,
    diameter: Length,
    elevation: Length,
    number_of_increments: usize,
    roughness: Length,
    eos: Option<&str>,
) -> Result<PipeBeggsAndBrillsResult> {
    let spec = &model_gen::PIPE_BEGGS_AND_BRILLS_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "inlet_n" => Some(inlet_n),
            "inlet_p" => Some(inlet_p.value),
            "inlet_t" => Some(inlet_t.value),
            "length" => Some(length.value),
            "diameter" => Some(diameter.value),
            "number_of_increments" => Some(number_of_increments as f64),
            _ => None,
        },
        &mut warnings,
    )?;

    let cubic = match eos.unwrap_or("srk") {
        "srk" => crate::Cubic::Srk,
        "pr" => crate::Cubic::Pr,
        other => {
            return Err(AzothError::invalid_input(
                "eos",
                format!("`{other}` is not one of `pr` or `srk`"),
            ));
        }
    };
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, _ideal) = databank::mixture_of(&names, cubic, None)?;

    let segment_length = length.value / number_of_increments as f64;
    let segment_elevation = elevation.value / number_of_increments as f64;

    let mut pressure = inlet_p.value;
    let mut pressure_profile = vec![pressure];
    let mut pressure_drop_profile = vec![0.0];
    let mut holdup_profile = Vec::with_capacity(number_of_increments + 1);
    let mut velocity_profile = Vec::with_capacity(number_of_increments + 1);
    let mut regime_profile = Vec::with_capacity(number_of_increments + 1);

    for _ in 0..number_of_increments {
        let properties =
            phase_properties(&mixture, components, inlet_n, inlet_t, pressure, inlet_z)?;
        let drop = segment_drop(
            &properties,
            inlet_t.value,
            pressure,
            diameter,
            segment_length,
            segment_elevation,
            roughness,
        )?;
        holdup_profile.push(drop.liquid_holdup);
        velocity_profile.push(drop.superficial_mixture_velocity.value);
        regime_profile.push(drop.regime.as_str().to_owned());
        let outlet = pressure - drop.total_drop.value;
        if outlet <= 0.0 {
            return Err(AzothError::invalid_input(
                "inlet_p",
                format!(
                    "a segment's outlet pressure is {outlet} Pa: the line does not deliver at this inlet"
                ),
            ));
        }
        pressure = outlet;
        pressure_profile.push(pressure);
        pressure_drop_profile.push(drop.total_drop.value);
    }

    // The class re-evaluates the correlation once after the loop, so its segment lists are
    // one longer than the pressure profile's segments.
    let outlet_properties =
        phase_properties(&mixture, components, inlet_n, inlet_t, pressure, inlet_z)?;
    let final_drop = segment_drop(
        &outlet_properties,
        inlet_t.value,
        pressure,
        diameter,
        segment_length,
        segment_elevation,
        roughness,
    )?;
    holdup_profile.push(final_drop.liquid_holdup);
    velocity_profile.push(final_drop.superficial_mixture_velocity.value);
    regime_profile.push(final_drop.regime.as_str().to_owned());

    Ok(PipeBeggsAndBrillsResult {
        pressure_profile,
        pressure_drop_profile,
        temperature_profile: vec![inlet_t.value],
        liquid_holdup_profile: holdup_profile,
        mixture_superficial_velocity_profile: velocity_profile,
        flow_regime_profile: regime_profile,
        outlet_pressure: pascals(pressure),
        total_pressure_drop: pascals(inlet_p.value - pressure),
        warnings,
    })
}

#[allow(clippy::too_many_arguments)]
fn segment_drop(
    properties: &PhaseProperties,
    temperature_k: f64,
    pressure_pa: f64,
    diameter: Length,
    segment_length: f64,
    segment_elevation: f64,
    roughness: Length,
) -> Result<crate::results::BeggsBrillSegmentResult> {
    beggs_brill_segment(
        azoth_core::units::cubic_meters_per_second(properties.gas_flow),
        azoth_core::units::cubic_meters_per_second(properties.liquid_flow),
        azoth_core::units::kilograms_per_cubic_meter(properties.gas_density),
        azoth_core::units::kilograms_per_cubic_meter(properties.liquid_density),
        azoth_core::units::pascal_seconds(properties.gas_viscosity),
        azoth_core::units::pascal_seconds(properties.liquid_viscosity),
        azoth_core::units::kelvins(temperature_k),
        pascals(pressure_pa),
        diameter,
        azoth_core::units::meters(segment_length),
        azoth_core::units::meters(segment_elevation),
        roughness,
    )
}

fn molar_mass(composition: &[f64], components: &[Component], names: &[String]) -> Result<f64> {
    let mut total = 0.0;
    for ((fraction, component), name) in composition.iter().zip(components).zip(names) {
        let mass = component.molar_mass.ok_or_else(|| {
            AzothError::invalid_input(
                "components",
                format!("`{name}` carries no molar mass, so a phase density cannot be formed"),
            )
        })?;
        total += fraction * mass;
    }
    Ok(total)
}

fn corrected_molar_volume(
    z_factor: f64,
    temperature_k: f64,
    pressure_pa: f64,
    composition: &[f64],
    components: &[Component],
) -> f64 {
    let ideal_volume = z_factor * MOLAR_GAS_CONSTANT * temperature_k / pressure_pa;
    let shift: f64 = composition
        .iter()
        .zip(components)
        .map(|(fraction, component)| fraction * component.volume_shift_srk)
        .sum();
    ideal_volume - shift
}

fn liquid_phase_kind(composition: &[f64], names: &[String]) -> PhaseKind {
    let water: f64 = composition
        .iter()
        .zip(names)
        .filter(|(_, name)| matches!(name.as_str(), "water" | "water_PC" | "H2O" | "waterpc"))
        .map(|(fraction, _)| fraction)
        .sum();
    if 1.0 - water > water {
        PhaseKind::Oil
    } else {
        PhaseKind::Aqueous
    }
}

fn phase_properties(
    mixture: &Mixture,
    names: &[String],
    inlet_n: f64,
    temperature: ThermodynamicTemperature,
    pressure_pa: f64,
    overall: &[f64],
) -> Result<PhaseProperties> {
    let components = mixture.components();
    let flash =
        crate::tp_multiflash::tp_multiflash(mixture, temperature, pascals(pressure_pa), overall)?;

    let mut gas_flow = 0.0;
    let mut liquid_flow = 0.0;
    let mut liquid_mass = 0.0;
    let mut liquid_weight = 0.0;
    let mut liquid_viscosity_weighted = 0.0;

    for (index, fraction) in flash.phase_fractions.iter().enumerate() {
        if *fraction <= 0.0 {
            continue;
        }
        let composition = &flash.x[index];
        let z_factor = flash.z_factor[index];
        let volume = corrected_molar_volume(
            z_factor,
            temperature.value,
            pressure_pa,
            composition,
            components,
        );
        let mass = molar_mass(composition, components, names)?;
        let density = mass / volume;
        let phase_flow = inlet_n * fraction * volume;
        let kind = phase_kind(
            z_factor,
            temperature.value,
            pressure_pa,
            composition,
            components,
            names,
        );
        let viscosity = viscosity_only(
            mixture,
            kind,
            temperature,
            pascals(pressure_pa),
            composition,
        )?;
        if kind == PhaseKind::Gas {
            gas_flow += phase_flow;
        } else {
            liquid_flow += phase_flow;
            liquid_mass += phase_flow * density;
            liquid_weight += phase_flow;
            liquid_viscosity_weighted += phase_flow * viscosity;
        }
    }

    if liquid_weight > 0.0 {
        let liquid_density = liquid_mass / liquid_weight;
        let liquid_viscosity = liquid_viscosity_weighted / liquid_weight;
        let gas_properties = gas_phase_properties(
            mixture,
            names,
            temperature,
            pressure_pa,
            overall,
            components,
        )?;
        Ok(PhaseProperties {
            gas_flow,
            liquid_flow,
            gas_density: if gas_flow > 0.0 {
                gas_properties.0
            } else {
                liquid_density
            },
            liquid_density,
            gas_viscosity: if gas_flow > 0.0 {
                gas_properties.1
            } else {
                liquid_viscosity
            },
            liquid_viscosity,
        })
    } else {
        let (density, viscosity) = gas_phase_properties(
            mixture,
            names,
            temperature,
            pressure_pa,
            overall,
            components,
        )?;
        Ok(PhaseProperties {
            gas_flow,
            liquid_flow: 0.0,
            gas_density: density,
            liquid_density: density,
            gas_viscosity: viscosity,
            liquid_viscosity: viscosity,
        })
    }
}

/// The single gas state of a feed whose flash left one phase, looked up once.
fn gas_phase_properties(
    mixture: &Mixture,
    names: &[String],
    temperature: ThermodynamicTemperature,
    pressure_pa: f64,
    overall: &[f64],
    components: &[Component],
) -> Result<(f64, f64)> {
    let flash =
        crate::tp_multiflash::tp_multiflash(mixture, temperature, pascals(pressure_pa), overall)?;
    for (index, fraction) in flash.phase_fractions.iter().enumerate() {
        if *fraction <= 0.0 {
            continue;
        }
        let composition = &flash.x[index];
        let z_factor = flash.z_factor[index];
        let volume = corrected_molar_volume(
            z_factor,
            temperature.value,
            pressure_pa,
            composition,
            components,
        );
        let mass = molar_mass(composition, components, names)?;
        let kind = phase_kind(
            z_factor,
            temperature.value,
            pressure_pa,
            composition,
            components,
            names,
        );
        if kind == PhaseKind::Gas {
            let viscosity = viscosity_only(
                mixture,
                kind,
                temperature,
                pascals(pressure_pa),
                composition,
            )?;
            return Ok((mass / volume, viscosity));
        }
    }
    // One liquid phase and no gas: the same liquid fills the absent gas's side.
    let composition = &flash.x[0];
    let z_factor = flash.z_factor[0];
    let volume = corrected_molar_volume(
        z_factor,
        temperature.value,
        pressure_pa,
        composition,
        components,
    );
    let mass = molar_mass(composition, components, names)?;
    let kind = phase_kind(
        z_factor,
        temperature.value,
        pressure_pa,
        composition,
        components,
        names,
    );
    let viscosity = viscosity_only(
        mixture,
        kind,
        temperature,
        pascals(pressure_pa),
        composition,
    )?;
    Ok((mass / volume, viscosity))
}

fn phase_kind(
    z_factor: f64,
    temperature_k: f64,
    pressure_pa: f64,
    composition: &[f64],
    components: &[Component],
    names: &[String],
) -> PhaseKind {
    let volume = z_factor * MOLAR_GAS_CONSTANT * temperature_k / pressure_pa;
    let b: f64 = composition
        .iter()
        .zip(components)
        .map(|(fraction, component)| {
            fraction * SRK_OMEGA_B * MOLAR_GAS_CONSTANT * component.tc.value / component.pc.value
        })
        .sum();
    if b > 0.0 && volume / b > 1.75 {
        PhaseKind::Gas
    } else {
        liquid_phase_kind(composition, names)
    }
}
