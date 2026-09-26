//! The interface equilibrium: a third mixture, flashed at the interface temperature.

use azoth_core::Result;
use azoth_core::units::{kelvins, pascals};

use super::phase::{Pick, index_of, phase_view};
use crate::stream::Stream;

/// The interface's equilibrium: the compositions the two films drive against, and the ratios
/// their driving forces are built from.
#[derive(Debug, Clone)]
pub struct InterfaceEquilibrium {
    /// The temperature the mixture was flashed at, K.
    pub temperature: f64,
    /// The components the interface is *reported* over, which is the active transfer list.
    pub components: Vec<String>,
    /// The interface gas phase's mole fractions, one per reported component.
    pub gas: Vec<f64>,
    /// The interface liquid phase's.
    pub liquid: Vec<f64>,
    /// `y/x`, floored at `1e-12`, and `1.0` where there is no liquid to divide by.
    pub ratios: Vec<f64>,
}

impl InterfaceEquilibrium {
    /// A component's equilibrium ratio, or `1.0` where the interface does not carry it.
    pub fn ratio(&self, component: &str) -> f64 {
        index_of(&self.components, component).map_or(1.0, |index| self.ratios[index])
    }

    /// The interface's gas mole fraction, and the caller's fallback where it is absent.
    ///
    /// **The fallback is nearly unreachable and it is still the class's structure**: the
    /// interface is reported over the same list the transfer walks, so the lookup answers -
    /// and the clamped `K x` the caller offers as a default is what would stand in for a
    /// component the interface does not carry.
    pub fn gas_fraction(&self, component: &str, fallback: f64) -> f64 {
        index_of(&self.components, component).map_or(fallback, |index| self.gas[index])
    }

    /// The interface's liquid mole fraction, with the caller's fallback.
    pub fn liquid_fraction(&self, component: &str, fallback: f64) -> f64 {
        index_of(&self.components, component).map_or(fallback, |index| self.liquid[index])
    }
}

/// `calculateInterfaceEquilibrium`: the segment's third mixture.
///
/// **Two component lists, and they are not the same one.** The mixture is built by cloning the
/// gas system and adding every component of the liquid that carries positive moles, so it is
/// stated over the **union**; the maps it answers are then read over the **active transfer
/// list**, which is the caller's `setTransferComponents` when there is one. Collapsing the two
/// is a real defect and not a tidy-up: a transfer list of `["CO2"]` read against the union's
/// indices picks methane's moles for CO2, which moves the transfer by three orders of
/// magnitude.
///
/// **On a failure the class falls back to the two unmixed phases**, which is the branch that
/// keeps a segment alive when the combined mixture will not flash. It is reproduced: the
/// fallback is the two systems the segment already holds, and it is not an error.
pub fn calculate_interface_equilibrium(
    gas: &Stream,
    liquid: &Stream,
    components: &[String],
    temperature: f64,
    pressure: f64,
) -> Result<InterfaceEquilibrium> {
    // ---- The union, in the class's own order: the gas's components, then whatever the liquid
    // adds. Only a component the liquid carries in positive moles is added.
    let mut union: Vec<String> = gas.components.clone();
    for (index, name) in liquid.components.iter().enumerate() {
        let moles = liquid.z[index] * liquid.n;
        if moles > 0.0 && !union.contains(name) {
            union.push(name.clone());
        }
    }

    let amounts: Vec<f64> = union
        .iter()
        .map(|name| moles_of(gas, name) + moles_of(liquid, name))
        .collect();
    let total: f64 = amounts.iter().sum();

    let mixed = if total > 0.0 {
        let z: Vec<f64> = amounts.iter().map(|amount| amount / total).collect();
        Stream::from_pt(
            union.clone(),
            z,
            total,
            pascals(pressure),
            kelvins(temperature),
        )
        .ok()
        .and_then(|stream| {
            let gas_phase = phase_view(&stream, Pick::Gas).ok()?;
            let liquid_phase = phase_view(&stream, Pick::Liquid).ok()?;
            Some((union.clone(), gas_phase, liquid_phase))
        })
    } else {
        None
    };

    let (names, gas_phase, liquid_phase) = match mixed {
        Some(triple) => triple,
        None => (
            gas.components.clone(),
            phase_view(gas, Pick::Gas)?,
            phase_view(liquid, Pick::Liquid)?,
        ),
    };

    let mut gas_fractions = Vec::with_capacity(components.len());
    let mut liquid_fractions = Vec::with_capacity(components.len());
    let mut ratios = Vec::with_capacity(components.len());
    for component in components {
        let x = liquid_phase.mole_fraction(&names, component);
        let y = gas_phase.mole_fraction(&names, component);
        gas_fractions.push(y);
        liquid_fractions.push(x);
        // `y/x` is a partial function at a zero denominator, and the class's policy is a
        // **substitution rather than a clamp**: where there is no liquid the ratio is one.
        if x > 1.0e-12 && y >= 0.0 {
            ratios.push((y / x).max(1.0e-12));
        } else {
            ratios.push(1.0);
        }
    }
    Ok(InterfaceEquilibrium {
        temperature,
        components: components.to_vec(),
        gas: gas_fractions,
        liquid: liquid_fractions,
        ratios,
    })
}

/// `componentMoles(system, name)`: a component's molar flow in the system, by name.
pub fn moles_of(stream: &Stream, name: &str) -> f64 {
    index_of(&stream.components, name).map_or(0.0, |index| (stream.z[index] * stream.n).max(0.0))
}
