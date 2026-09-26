//! One segment: the class's `calculateSegment`, and the transfer it applies.

use azoth_core::Result;
use azoth_core::units::kelvins;

use super::equilibrium::{InterfaceEquilibrium, calculate_interface_equilibrium};
use super::fallbacks::Fallbacks;
use super::film::{clamp, combine_film_fluxes, film_coefficient};
use super::heat::{apply_interphase_heat_transfer, heat_capacity_rate};
use super::phase::{PhaseView, index_of};
use super::transport::{SnapshotSettings, TransportSnapshot, calculate_transport_snapshot};
use crate::stream::Stream;

/// `maxTransferFractionPerSegment`: the share of a donor's inventory one segment may move.
pub const MAX_TRANSFER_FRACTION: f64 = 0.35;

/// `maxHeatTransferFractionPerSegment`.
pub const MAX_HEAT_TRANSFER_FRACTION: f64 = 0.50;

/// One segment's answer, as the class's `SegmentResult` carries it.
#[derive(Debug, Clone)]
pub struct SegmentResult {
    /// `segmentNumber`, counted from one at the bottom.
    pub number: usize,
    /// `heightFromBottom`, the segment's mid-point.
    pub height_from_bottom: f64,
    pub gas_temperature: f64,
    pub liquid_temperature: f64,
    pub gas_pressure: f64,
    pub liquid_pressure: f64,
    pub gas_molar_flow: f64,
    pub liquid_molar_flow: f64,
    pub gas_density: f64,
    pub liquid_density: f64,
    pub gas_viscosity: f64,
    pub liquid_viscosity: f64,
    pub gas_diffusivity: f64,
    pub liquid_diffusivity: f64,
    pub wetted_area: f64,
    pub k_ga: f64,
    pub k_la: f64,
    pub gas_heat_transfer_coefficient: f64,
    pub liquid_heat_transfer_coefficient: f64,
    pub overall_heat_transfer_coefficient: f64,
    pub interface_temperature: f64,
    pub heat_transfer_rate: f64,
    pub pressure_drop_per_meter: f64,
    pub percent_flood: f64,
    /// `netMolarTransfer`, the sum of the per-component transfers.
    pub net_molar_transfer: f64,
    /// The per-component transfers, positive from gas to liquid, zero transfers omitted.
    pub component_transfer: Vec<(String, f64)>,
    /// The interface's own compositions and ratios.
    pub interface: InterfaceEquilibrium,
    /// `gas.getEnthalpy() + liquid.getEnthalpy() - inletTotalEnthalpy`, in W.
    pub enthalpy_balance_residual: f64,
    /// Which of the class's constants stood in for a missing property.
    pub fallbacks: Fallbacks,
}

/// A segment's two outlets and its record.
#[derive(Debug, Clone)]
pub struct SegmentComputation {
    pub gas: Stream,
    pub liquid: Stream,
    pub result: SegmentResult,
}

/// `calculateSegment`.
///
/// **The order is load-bearing, and it is the class's.** The transport snapshot is taken
/// *before* the interface equilibrium (which needs the snapshot's interface temperature); the
/// transfers are applied component by component; the heat step follows them; and the record's
/// transport numbers are the **snapshot's**, which are the segment's *inlet* state, while its
/// temperatures and flows are the post-flash outlets'.
#[allow(clippy::too_many_arguments)]
pub fn calculate_segment(
    index: usize,
    gas_in: &Stream,
    liquid_in: &Stream,
    segment_height: f64,
    segment_volume: f64,
    transfer_components: &[String],
    settings: &SnapshotSettings,
    matrix_model: bool,
) -> Result<SegmentComputation> {
    let mut gas = gas_in.clone();
    let mut liquid = liquid_in.clone();
    let inlet_enthalpy = total_enthalpy(&gas) + total_enthalpy(&liquid);

    let (snapshot, gas_phase, liquid_phase) =
        calculate_transport_snapshot(&gas, &liquid, segment_height, settings)?;

    let components: Vec<String> = if transfer_components.is_empty() {
        union(&gas.components, &liquid.components)
    } else {
        transfer_components.to_vec()
    };
    let interface = calculate_interface_equilibrium(
        &gas,
        &liquid,
        &components,
        snapshot.interface_temperature,
        0.5 * (gas.p.value + liquid.p.value),
    )?;

    let mut transfers: Vec<(String, f64)> = Vec::new();
    let mut heat_transfer_rate = 0.0;
    // **A bed of no height transfers nothing and heats nothing**, and the whole block is
    // behind the class's own `if (segmentHeight > 0.0)`.
    if segment_height > 0.0 {
        for component in &components {
            let transfer = component_transfer(
                component,
                &gas,
                &liquid,
                &gas_phase,
                &liquid_phase,
                &snapshot,
                &interface,
                segment_volume,
                matrix_model,
            )?;
            if transfer != 0.0 {
                apply_component_transfer(&mut gas, component, -transfer)?;
                apply_component_transfer(&mut liquid, component, transfer)?;
                transfers.push((component.clone(), transfer));
            }
        }
        let gas_capacity = heat_capacity_rate(
            gas_phase.n * gas_phase.molar_mass,
            snapshot.gas_heat_capacity,
        );
        let liquid_capacity = heat_capacity_rate(
            liquid_phase.n * liquid_phase.molar_mass,
            snapshot.liquid_heat_capacity,
        );
        let step = apply_interphase_heat_transfer(
            gas.t.value,
            liquid.t.value,
            gas_capacity,
            liquid_capacity,
            snapshot.overall_heat_transfer_coefficient,
            segment_volume,
            settings.heat_transfer_none,
            MAX_HEAT_TRANSFER_FRACTION,
        );
        heat_transfer_rate = step.rate;
        gas = at_temperature(&gas, step.gas_temperature)?;
        liquid = at_temperature(&liquid, step.liquid_temperature)?;
    }

    let net_molar_transfer: f64 = transfers.iter().map(|(_, value)| *value).sum();
    let result = SegmentResult {
        number: index + 1,
        height_from_bottom: (index as f64 + 0.5) * segment_height,
        gas_temperature: gas.t.value,
        liquid_temperature: liquid.t.value,
        gas_pressure: gas.p.value,
        liquid_pressure: liquid.p.value,
        gas_molar_flow: gas.n,
        liquid_molar_flow: liquid.n,
        gas_density: snapshot.gas_density,
        liquid_density: snapshot.liquid_density,
        gas_viscosity: snapshot.gas_viscosity,
        liquid_viscosity: snapshot.liquid_viscosity,
        gas_diffusivity: snapshot.gas_diffusivity,
        liquid_diffusivity: snapshot.liquid_diffusivity,
        wetted_area: snapshot.wetted_area,
        k_ga: snapshot.k_ga,
        k_la: snapshot.k_la,
        gas_heat_transfer_coefficient: snapshot.gas_heat_transfer_coefficient,
        liquid_heat_transfer_coefficient: snapshot.liquid_heat_transfer_coefficient,
        overall_heat_transfer_coefficient: snapshot.overall_heat_transfer_coefficient,
        interface_temperature: interface.temperature,
        heat_transfer_rate,
        pressure_drop_per_meter: snapshot.pressure_drop_per_meter,
        percent_flood: snapshot.percent_flood,
        net_molar_transfer,
        component_transfer: transfers,
        interface,
        enthalpy_balance_residual: total_enthalpy(&gas) + total_enthalpy(&liquid) - inlet_enthalpy,
        fallbacks: snapshot.fallbacks,
    };
    Ok(SegmentComputation {
        gas,
        liquid,
        result,
    })
}

/// `calculateComponentTransfer`: the unbounded film transfer, then the class's inventory cap.
#[allow(clippy::too_many_arguments)]
fn component_transfer(
    component: &str,
    gas: &Stream,
    liquid: &Stream,
    gas_phase: &PhaseView,
    liquid_phase: &PhaseView,
    snapshot: &TransportSnapshot,
    interface: &InterfaceEquilibrium,
    segment_volume: f64,
    matrix_model: bool,
) -> Result<f64> {
    let k_value = interface.ratio(component);
    let gas_fraction = gas_phase.mole_fraction(&gas.components, component);
    let liquid_fraction = liquid_phase.mole_fraction(&liquid.components, component);
    let gas_interface =
        interface.gas_fraction(component, clamp(k_value * liquid_fraction, 0.0, 0.999999));
    let liquid_interface = interface.liquid_fraction(
        component,
        if k_value > 1.0e-12 {
            clamp(gas_interface / k_value, 0.0, 0.999999)
        } else {
            liquid_fraction
        },
    );
    let gas_driving_force = gas_fraction - gas_interface;
    // A vanishing driving force, or a film that is not there, is no transfer - and the guard
    // is the class's, so a state with one zero coefficient still reports the other.
    if gas_driving_force.abs() < 1.0e-12 || snapshot.k_ga <= 0.0 || snapshot.k_la <= 0.0 {
        return Ok(0.0);
    }

    let gas_index = index_of(&gas.components, component).unwrap_or(0);
    let gas_film = film_coefficient(
        gas_phase,
        gas_index,
        snapshot.k_ga,
        snapshot.gas_diffusivity,
        true,
        matrix_model,
    );
    let liquid_film = film_coefficient(
        liquid_phase,
        gas_index,
        snapshot.k_la,
        snapshot.liquid_diffusivity,
        false,
        matrix_model,
    );
    let gas_flux = gas_film * concentration(gas_phase) * gas_driving_force;
    let liquid_flux =
        liquid_film * concentration(liquid_phase) * (liquid_interface - liquid_fraction);

    let mut transfer_density = combine_film_fluxes(gas_flux, liquid_flux);
    // **The harmonic mean's collapse is not a zero transfer**: where the two films disagree in
    // sign the class falls to the overall two-resistance form, which is a different expression
    // on the same state.
    if transfer_density == 0.0 {
        let y_star = clamp(k_value * liquid_fraction, 0.0, 0.999999);
        let driving_force = gas_fraction - y_star;
        let overall = 1.0 / (1.0 / gas_film + k_value.max(1.0e-12) / liquid_film);
        transfer_density = overall * concentration(gas_phase) * driving_force;
    }
    let proposed = transfer_density * segment_volume;
    Ok(limit_transfer(
        component,
        proposed,
        component_moles(gas),
        &gas.components,
        component_moles(liquid),
        &liquid.components,
    ))
}

/// `limitTransfer`: the cap on a donor's inventory.
fn limit_transfer(
    component: &str,
    proposed: f64,
    gas_moles: Vec<f64>,
    gas_components: &[String],
    liquid_moles: Vec<f64>,
    liquid_components: &[String],
) -> f64 {
    let inventory = |moles: &Vec<f64>, components: &[String]| {
        index_of(components, component).map_or(0.0, |index| moles[index].max(0.0))
    };
    if proposed > 0.0 {
        let available = inventory(&gas_moles, gas_components) * MAX_TRANSFER_FRACTION;
        proposed.min(available.max(0.0))
    } else if proposed < 0.0 {
        let available = inventory(&liquid_moles, liquid_components) * MAX_TRANSFER_FRACTION;
        -(-proposed).min(available.max(0.0))
    } else {
        0.0
    }
}

/// `molarConcentration`: the phase's molar density, mol/m³.
fn concentration(view: &PhaseView) -> f64 {
    view.density / super::film::finite_positive(view.molar_mass, 0.020)
}

/// The system's per-component molar flow, in its own component order.
pub fn component_moles(stream: &Stream) -> Vec<f64> {
    stream
        .z
        .iter()
        .map(|fraction| fraction * stream.n)
        .collect()
}

/// The union of two systems' components, in the order the class builds it.
fn union(gas: &[String], liquid: &[String]) -> Vec<String> {
    let mut names: Vec<String> = gas.to_vec();
    for name in liquid {
        if !names.contains(name) {
            names.push(name.clone());
        }
    }
    names
}

/// `addComponent`: move moles into or out of a system, holding its temperature and pressure.
fn apply_component_transfer(stream: &mut Stream, component: &str, delta: f64) -> Result<()> {
    let index = match index_of(&stream.components, component) {
        Some(index) => index,
        None => {
            stream.components.push(component.to_string());
            stream.z.push(0.0);
            stream.components.len() - 1
        }
    };
    let moles = stream.z[index] * stream.n + delta;
    if moles < 0.0 {
        return Err(azoth_core::AzothError::invalid_input(
            component,
            format!(
                "a transfer of {delta} mol/s would leave {} with {moles} mol of {component}",
                stream.n
            ),
        ));
    }
    let total = stream.n + delta;
    let z: Vec<f64> = stream
        .z
        .iter()
        .enumerate()
        .map(|(component_index, fraction)| {
            let amount = if component_index == index {
                moles
            } else {
                fraction * stream.n
            };
            if total > 0.0 { amount / total } else { 0.0 }
        })
        .collect();
    let rebuilt = Stream::from_pt(stream.components.clone(), z, total, stream.p, stream.t)?;
    *stream = rebuilt;
    Ok(())
}

/// The same system at a different temperature, with its enthalpy re-derived.
fn at_temperature(stream: &Stream, temperature: f64) -> Result<Stream> {
    Stream::from_pt(
        stream.components.clone(),
        stream.z.clone(),
        stream.n,
        stream.p,
        kelvins(temperature),
    )
}

/// A system's total enthalpy, in W: its molar enthalpy times its molar flow.
pub fn total_enthalpy(stream: &Stream) -> f64 {
    stream.h.value * stream.n
}
