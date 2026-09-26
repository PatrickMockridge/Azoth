//! The segment's transport snapshot: the class's `calculateTransportSnapshot`.
//!
//! Every property in it is either read from a phase or substituted for, and the substitutions
//! are the class's own constants - see [`super::fallbacks`].

use azoth_core::Result;
use azoth_core::units::{
    MassRate, kilograms_per_cubic_meter, kilograms_per_second, meters, newtons_per_meter,
    pascal_seconds,
};
use azoth_eos::parachor_mixture_surface_tension::parachor_mixture_surface_tension;
use azoth_hydraulics::packing::packing_or_default;
use azoth_hydraulics::packing_hydraulics::{PackingState, packing_hydraulics};

use super::fallbacks::{
    DEFAULT_GAS_DIFFUSIVITY, DEFAULT_GAS_HEAT_CAPACITY, DEFAULT_GAS_THERMAL_CONDUCTIVITY,
    DEFAULT_LIQUID_DIFFUSIVITY, DEFAULT_LIQUID_HEAT_CAPACITY, DEFAULT_LIQUID_THERMAL_CONDUCTIVITY,
    DEFAULT_SURFACE_TENSION, Fallbacks, MIN_GAS_DIFFUSIVITY, MIN_LIQUID_DIFFUSIVITY,
};
use super::film::finite_positive;
use super::heat::{
    combine_heat_transfer_coefficients, interface_temperature, volumetric_heat_transfer_coefficient,
};
use super::phase::{PhaseView, Pick, phase_view};
use crate::stream::Stream;

/// The geometry and the switches the snapshot reads.
#[derive(Debug, Clone)]
pub struct SnapshotSettings {
    /// The column's internal diameter, in metres.
    pub column_diameter: f64,
    /// The packed section's height, in metres.
    pub packed_height: f64,
    /// The packing's name, which the hydraulics resolves.
    pub packing: String,
    /// Whether the heat-transfer model is `none`, which zeroes both film coefficients.
    pub heat_transfer_none: bool,
    /// `massTransferCorrectionFactor`, which scales both film coefficients.
    pub mass_transfer_correction: f64,
    /// `heatTransferCorrectionFactor`, which scales both heat coefficients.
    pub heat_transfer_correction: f64,
}

/// `calculateTransportSnapshot`: the 17 numbers the segment's transfer and heat steps read.
#[derive(Debug, Clone)]
pub struct TransportSnapshot {
    pub gas_density: f64,
    pub liquid_density: f64,
    pub gas_viscosity: f64,
    pub liquid_viscosity: f64,
    pub gas_diffusivity: f64,
    pub liquid_diffusivity: f64,
    pub wetted_area: f64,
    pub k_ga: f64,
    pub k_la: f64,
    pub gas_heat_capacity: f64,
    pub liquid_heat_capacity: f64,
    pub gas_heat_transfer_coefficient: f64,
    pub liquid_heat_transfer_coefficient: f64,
    pub overall_heat_transfer_coefficient: f64,
    pub interface_temperature: f64,
    pub pressure_drop_per_meter: f64,
    pub percent_flood: f64,
    /// Which of the class's constants stood in for a missing property.
    pub fallbacks: Fallbacks,
}

/// The vector `averageDiffusivity` averages, as NeqSim's own flash path leaves it.
///
/// **It is empty, and that is the measurement.** The class averages
/// `getEffectiveDiffusionCoefficient(i)` over the components whose value is positive, and
/// NeqSim never populates that vector: on the CO2/water absorber's own state every entry of
/// both phases answers `0.0`, so the sum is empty and the class's `DEFAULT_*` constant stands
/// in. What the same state *does* carry is the pair matrix - `getDiffusionCoefficient(0, 1)`
/// is `9.457085848059925e-7` m²/s for the gas - so the film model's ratios are real against a
/// reference that is not.
///
/// Reproduced rather than improved: `Diffusivity.calcEffectiveDiffusionCoefficients` is the
/// class that would close it, and the capture prints both halves.
const NEQSIM_EFFECTIVE_DIFFUSIVITY: &[f64] = &[];

/// `averageDiffusivity`: the mean of the positive-finite effective diffusivities, and the
/// class's constant where there are none.
pub fn average_diffusivity(effective: &[f64], gas_phase: bool) -> (f64, bool) {
    let fallback = if gas_phase {
        DEFAULT_GAS_DIFFUSIVITY
    } else {
        DEFAULT_LIQUID_DIFFUSIVITY
    };
    let mut sum = 0.0;
    let mut count = 0;
    for value in effective {
        if value.is_finite() && *value > 0.0 {
            sum += value;
            count += 1;
        }
    }
    if count > 0 {
        (sum / count as f64, false)
    } else {
        (fallback, true)
    }
}

/// Build the snapshot, the segment's gas phase and its liquid phase.
pub fn calculate_transport_snapshot(
    gas: &Stream,
    liquid: &Stream,
    segment_height: f64,
    settings: &SnapshotSettings,
) -> Result<(TransportSnapshot, PhaseView, PhaseView)> {
    let gas_phase = phase_view(gas, Pick::Gas)?;
    let liquid_phase = phase_view(liquid, Pick::Liquid)?;

    // The class's three sentinels: a non-positive density or viscosity is replaced outright
    // rather than reported, and they are not among the four properties a warning names.
    let gas_density = finite_positive(gas_phase.density, 1.0);
    let liquid_density = finite_positive(liquid_phase.density, 800.0);
    let gas_viscosity = finite_positive(gas_phase.transport.mu.value, 1.0e-5);
    let liquid_viscosity = finite_positive(liquid_phase.transport.mu.value, 1.0e-3);

    let mut fallbacks = Fallbacks::default();
    let (gas_diffusivity, gas_took) = average_diffusivity(NEQSIM_EFFECTIVE_DIFFUSIVITY, true);
    let (liquid_diffusivity, liquid_took) =
        average_diffusivity(NEQSIM_EFFECTIVE_DIFFUSIVITY, false);
    fallbacks.diffusivity = gas_took || liquid_took;

    let (surface_tension, tension_took) = estimate_surface_tension(gas, liquid)?;
    fallbacks.surface_tension = tension_took;

    let gas_heat_capacity = finite_positive(gas_phase.cp_mass, DEFAULT_GAS_HEAT_CAPACITY);
    let liquid_heat_capacity = finite_positive(liquid_phase.cp_mass, DEFAULT_LIQUID_HEAT_CAPACITY);
    fallbacks.heat_capacity = gas_heat_capacity == DEFAULT_GAS_HEAT_CAPACITY
        || liquid_heat_capacity == DEFAULT_LIQUID_HEAT_CAPACITY;

    let gas_conductivity = finite_positive(
        gas_phase.transport.k.value,
        DEFAULT_GAS_THERMAL_CONDUCTIVITY,
    );
    let liquid_conductivity = finite_positive(
        liquid_phase.transport.k.value,
        DEFAULT_LIQUID_THERMAL_CONDUCTIVITY,
    );
    fallbacks.thermal_conductivity = gas_conductivity == DEFAULT_GAS_THERMAL_CONDUCTIVITY
        || liquid_conductivity == DEFAULT_LIQUID_THERMAL_CONDUCTIVITY;

    // ---- The packing hydraulics, on the class's own eleven inputs.
    //
    // **The height it is given is floored at `1e-9`, not passed as it stands** - which is what
    // makes the zero-height state reachable at all, since `setPackedHeight` refuses a bed of
    // no height and the class never asks it to.
    let packing = packing_or_default(&settings.packing);
    let hydraulics = packing_hydraulics(
        &packing.name,
        PackingState {
            column_diameter: meters(settings.column_diameter),
            packed_height: segment_height.max(1.0e-9),
            vapor_mass_flow: mass_rate(gas_phase.n, gas_phase.molar_mass),
            liquid_mass_flow: mass_rate(liquid_phase.n, liquid_phase.molar_mass),
            vapor_density: kilograms_per_cubic_meter(gas_density),
            liquid_density: kilograms_per_cubic_meter(liquid_density),
            vapor_viscosity: pascal_seconds(gas_viscosity),
            liquid_viscosity: pascal_seconds(liquid_viscosity),
            surface_tension: newtons_per_meter(surface_tension),
            vapor_diffusivity: gas_diffusivity,
            liquid_diffusivity,
            hydraulic_capacity_factor: 1.0,
        },
    )?;

    // `BILLET_SCHULTES_1999` is refused by name before a snapshot is built, so both multipliers
    // are one: they exist in the class as a constant scaling of the two coefficients and not as
    // a correlation.
    let k_ga = non_negative(hydraulics.k_ga, 0.0) * settings.mass_transfer_correction;
    let k_la = non_negative(hydraulics.k_la, 0.0) * settings.mass_transfer_correction;

    let gas_heat_transfer_coefficient = volumetric_heat_transfer_coefficient(
        k_ga,
        gas_density,
        gas_heat_capacity,
        gas_viscosity,
        gas_diffusivity,
        gas_conductivity,
        settings.heat_transfer_none,
        settings.heat_transfer_correction,
    );
    let liquid_heat_transfer_coefficient = volumetric_heat_transfer_coefficient(
        k_la,
        liquid_density,
        liquid_heat_capacity,
        liquid_viscosity,
        liquid_diffusivity,
        liquid_conductivity,
        settings.heat_transfer_none,
        settings.heat_transfer_correction,
    );
    let overall_heat_transfer_coefficient = combine_heat_transfer_coefficients(
        gas_heat_transfer_coefficient,
        liquid_heat_transfer_coefficient,
    );
    let interface = interface_temperature(
        gas.t.value,
        liquid.t.value,
        gas_heat_transfer_coefficient,
        liquid_heat_transfer_coefficient,
    );

    Ok((
        TransportSnapshot {
            gas_density,
            liquid_density,
            gas_viscosity,
            liquid_viscosity,
            // **The snapshot lifts both to a floor**, and the floor is reachable only through a
            // negative diffusivity, which the class's own substitutions never produce.
            gas_diffusivity: gas_diffusivity.max(MIN_GAS_DIFFUSIVITY),
            liquid_diffusivity: liquid_diffusivity.max(MIN_LIQUID_DIFFUSIVITY),
            wetted_area: non_negative(hydraulics.wetted_area, 0.0),
            k_ga,
            k_la,
            gas_heat_capacity,
            liquid_heat_capacity,
            gas_heat_transfer_coefficient,
            liquid_heat_transfer_coefficient,
            overall_heat_transfer_coefficient,
            interface_temperature: interface,
            pressure_drop_per_meter: non_negative(hydraulics.pressure_drop_per_meter.value, 0.0),
            percent_flood: non_negative(hydraulics.percent_flood, 0.0),
            fallbacks,
        },
        gas_phase,
        liquid_phase,
    ))
}

/// `estimateSurfaceTension`: a **fourth** mixture, at the segment's own temperature and
/// pressure, asked for its interphase tension.
///
/// The class builds it by adding the liquid's components to a clone of the gas, which is the
/// same construction the interface equilibrium makes - but at the *segment's* state rather
/// than the interface temperature, and it reads the two phases' own densities, molar masses
/// and compositions. A mixture that cannot form two phases has no interface, and the class's
/// answer there is its `DEFAULT_SURFACE_TENSION`.
fn estimate_surface_tension(gas: &Stream, liquid: &Stream) -> Result<(f64, bool)> {
    let names: Vec<&str> = gas
        .components
        .iter()
        .chain(
            liquid
                .components
                .iter()
                .filter(|name| !gas.components.contains(name)),
        )
        .map(String::as_str)
        .collect();
    let total = gas.n + liquid.n;
    let z: Vec<f64> = names
        .iter()
        .map(|name| {
            // **A component the system does not carry contributes nothing, and index zero is
            // not nothing.** `unwrap_or(0)` here read the *first* component's moles for every
            // name the gas does not have - so the solvent's water entered the mixture at
            // methane's flow, and the surface tension came back as the class's fallback.
            let gas_fraction = gas
                .components
                .iter()
                .position(|c| c == name)
                .map_or(0.0, |index| gas.z[index] * gas.n);
            let liquid_fraction = liquid
                .components
                .iter()
                .position(|c| c == name)
                .map_or(0.0, |index| liquid.z[index] * liquid.n);
            (gas_fraction + liquid_fraction) / total
        })
        .collect();
    let mixed = match Stream::from_pt(
        names.iter().map(|name| (*name).to_string()).collect(),
        z,
        total,
        gas.p,
        gas.t,
    ) {
        Ok(mixed) => mixed,
        Err(_) => return Ok((DEFAULT_SURFACE_TENSION, true)),
    };
    // **The class asks for a two-phase mixture with a gas in it, and answers its own default
    // otherwise.** `estimateSurfaceTension` is `if (mixed.hasPhaseType(GAS) &&
    // mixed.getNumberOfPhases() > 1) { ... }` with `return DEFAULT_SURFACE_TENSION` on both the
    // else and the catch - and the gate is *reached*: the class's own absorber state mixes to
    // 68 per cent water, which is one liquid phase, so its surface tension is the `0.025`
    // constant. Measured: feeding that constant to the hydraulics reproduces the capture's
    // wetted area `56.907889197418655` exactly, where the parachor's own `0.05303` gives
    // `32.91880623113832`.
    let phases = match mixed.mixture().and_then(|(mixture, _)| {
        azoth_eos::pt_flash::pt_flash(&mixture, mixed.t, mixed.p, &mixed.z)
    }) {
        Ok(flash) => flash,
        Err(_) => return Ok((DEFAULT_SURFACE_TENSION, true)),
    };
    if !matches!(phases.phase, azoth_eos::Phase::TwoPhase) {
        return Ok((DEFAULT_SURFACE_TENSION, true));
    }
    let (gas_phase, liquid_phase) = match (
        phase_view(&mixed, Pick::Gas),
        phase_view(&mixed, Pick::Liquid),
    ) {
        (Ok(gas_phase), Ok(liquid_phase)) => (gas_phase, liquid_phase),
        _ => return Ok((DEFAULT_SURFACE_TENSION, true)),
    };
    if gas_phase.kind != azoth_eos::phase_transport::PhaseKind::Gas {
        return Ok((DEFAULT_SURFACE_TENSION, true));
    }
    // **A gas-and-aqueous pair answers `0.0` upstream, and the constant is the class's answer.**
    // Measured on this class's own path - the gas cloned, the liquid's positive-mole components
    // added, flashed at the gas's temperature and pressure - `InterfaceProperties.getSurfaceTension`
    // returns `0.0` for the pair at dissolved-CO2 loadings of `0`, `1e-6`, `1e-3`, `1e-2` and
    // `0.04`. `isFinitePositive(0.0)` is false, so `estimateSurfaceTension` falls through to its
    // `DEFAULT_SURFACE_TENSION`, which is what the capture's wetted area is built from.
    //
    // The parachor form is *not* wrong here - it answers `0.05303` on the same mixture, which
    // would move the wetted area from `56.9` to `32.9` - it is a different branch. An oil pair
    // takes it, and the class's own states are all gas-and-aqueous, so this is where its answer
    // comes from.
    if liquid_phase.kind == azoth_eos::phase_transport::PhaseKind::Aqueous {
        return Ok((DEFAULT_SURFACE_TENSION, true));
    }
    let parachors: Vec<f64> = match mixed.mixture() {
        Ok((mixture, _)) => mixture
            .components()
            .iter()
            .map(|component| component.parachor)
            .collect(),
        Err(_) => return Ok((DEFAULT_SURFACE_TENSION, true)),
    };
    if parachors.iter().any(|value| *value <= 0.0) {
        return Ok((DEFAULT_SURFACE_TENSION, true));
    }
    let sigma = parachor_mixture_surface_tension(
        &parachors,
        kilograms_per_cubic_meter(gas_phase.density),
        azoth_core::units::kilograms_per_mole(gas_phase.molar_mass),
        &gas_phase.z,
        kilograms_per_cubic_meter(liquid_phase.density),
        azoth_core::units::kilograms_per_mole(liquid_phase.molar_mass),
        &liquid_phase.z,
    );
    match sigma {
        Ok(result) if result.sigma.value > 0.0 && result.sigma.value.is_finite() => {
            Ok((result.sigma.value, false))
        }
        _ => Ok((DEFAULT_SURFACE_TENSION, true)),
    }
}

/// `finiteNonNegative`.
fn non_negative(value: f64, fallback: f64) -> f64 {
    if value.is_finite() && value >= 0.0 {
        value
    } else {
        fallback
    }
}

/// The phase's mass flow, kg/s: its molar flow times its molar mass.
fn mass_rate(molar_flow: f64, molar_mass: f64) -> MassRate {
    kilograms_per_second(molar_flow * molar_mass)
}
