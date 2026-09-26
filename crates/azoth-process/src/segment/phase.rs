//! One phase of a flashed system: the class's `PhaseInterface`, as a value.
//!
//! `RateBasedPackedColumn` holds two `SystemInterface`s - a gas one and a liquid one - and
//! reads a *phase* out of each. This module is that read: flash the system, pick the phase the
//! class would pick, and answer the state and the transport the segment model asks for.

use azoth_core::units::{Pressure, ThermodynamicTemperature, kilograms_per_mole, pascals};
use azoth_core::{AzothError, Result};
use azoth_eos::hydrate_inhibitor_wt::{PhaseLabel, label};
use azoth_eos::phase_transport::{PhaseKind, phase_transport};
use azoth_eos::results::PhaseTransportResult;
use azoth_eos::{
    Cubic, Mixture, Phase, databank, molar_enthalpy_entropy, pr_mass_density, pr_molar_volume,
    pt_flash,
};

use crate::stream::Stream;

/// Which of a flashed system's phases the class is asking for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Pick {
    /// `getGasPhase`: the phase labelled `GAS`, else phase zero.
    Gas,
    /// `getLiquidPhase`: `AQUEOUS`, then `OIL`, then the first phase that is not the gas.
    Liquid,
}

/// A phase, at the state the system's flash put it in.
#[derive(Debug, Clone)]
pub struct PhaseView {
    /// The phase's composition.
    pub z: Vec<f64>,
    /// The phase's molar flow, mol/s.
    pub n: f64,
    /// The phase's temperature, which is its system's.
    pub t: ThermodynamicTemperature,
    /// The phase's pressure, which is its system's.
    pub p: Pressure,
    /// The phase's molar mass, kg/mol.
    pub molar_mass: f64,
    /// The phase's mass density, kg/m³ - the cubic's volume with the Peneloux shift applied,
    /// which is `getDensity("kg/m3")`.
    pub density: f64,
    /// The phase's heat capacity on a mass basis, J/(kg·K) - `getCp("J/kgK")`.
    pub cp_mass: f64,
    /// Which of NeqSim's three phase kinds it carries.
    pub kind: PhaseKind,
    /// Its viscosity, conductivity and diffusivities.
    pub transport: PhaseTransportResult,
}

impl PhaseView {
    /// A component's mole fraction in the phase, or zero where it is absent.
    pub fn mole_fraction(&self, components: &[String], name: &str) -> f64 {
        index_of(components, name)
            .and_then(|index| self.z.get(index).copied())
            .map_or(0.0, |fraction| fraction.max(0.0))
    }
}

/// The index of a component in a system's own order.
pub fn index_of(components: &[String], name: &str) -> Option<usize> {
    components.iter().position(|component| component == name)
}

/// The phase the class picks out of `stream`'s system.
pub fn phase_view(stream: &Stream, pick: Pick) -> Result<PhaseView> {
    let (mixture, ideal_gas) = stream.mixture()?;
    let flash = pt_flash(&mixture, stream.t, stream.p, &stream.z)?;
    let reduced = mixture.reduced_parameters(stream.t, stream.p)?;

    // **The phase, and it is the label rule that names it, not `beta`.** `getGasPhase` and
    // `getLiquidPhase` look for a `PhaseType`, and `PhaseEos.init` is what assigns one - the
    // volume ratio past `1.75` is the gas, and of the two liquids the one whose hydrocarbons
    // outweigh its aqueous components is the oil. It is the same rule `three_phase_separator`
    // and `pipe` already read, and it is public in `azoth-eos` for exactly this reason.
    let labelled: Vec<(PhaseLabel, Vec<f64>, f64)> = match flash.phase {
        // **A single phase has the system's composition, and the flash's own `x` or `y` is not
        // it.** The flash answers a trial phase there - `column/tray.rs` records the
        // measurement - so the one phase *is* the system.
        Phase::AllVapour => vec![(PhaseLabel::Gas, stream.z.clone(), flash.z_vapour)],
        Phase::AllLiquid | Phase::Trivial => vec![(
            label(&mixture, &reduced, &stream.z, flash.z_liquid)?,
            stream.z.clone(),
            flash.z_liquid,
        )],
        Phase::TwoPhase => {
            if flash.beta.is_none() {
                return Err(AzothError::invalid_input(
                    "flash",
                    "a two-phase answer with no vapour fraction, so neither phase's share is known",
                ));
            }
            vec![
                (PhaseLabel::Gas, flash.y.clone(), flash.z_vapour),
                (
                    label(&mixture, &reduced, &flash.x, flash.z_liquid)?,
                    flash.x.clone(),
                    flash.z_liquid,
                ),
            ]
        }
    };

    let chosen = match pick {
        // `getGasPhase` asks for the gas phase and falls back to phase zero, which on a system
        // with no vapour is the liquid's own state.
        Pick::Gas => labelled
            .iter()
            .find(|(kind, _, _)| *kind == PhaseLabel::Gas)
            .or_else(|| labelled.first()),
        // `getLiquidPhase` is a ladder: AQUEOUS, then OIL, then LIQUID, then the first phase
        // that is not the gas, then phase zero. The first three are the label rule's two
        // liquid branches in that order.
        Pick::Liquid => labelled
            .iter()
            .find(|(kind, _, _)| *kind == PhaseLabel::Aqueous)
            .or_else(|| {
                labelled
                    .iter()
                    .find(|(kind, _, _)| *kind == PhaseLabel::Oil)
            })
            .or_else(|| {
                labelled
                    .iter()
                    .find(|(kind, _, _)| *kind != PhaseLabel::Gas)
            })
            .or_else(|| labelled.first()),
    }
    .ok_or_else(|| {
        AzothError::invalid_input("flash", "the flash answered no phase at all".to_string())
    })?;

    let (kind, z, root) = chosen.clone();
    let share = if labelled.len() == 1 {
        1.0
    } else if kind == PhaseLabel::Gas {
        flash.beta.unwrap_or(1.0)
    } else {
        1.0 - flash.beta.unwrap_or(0.0)
    };

    let count = mixture.components().len();
    if z.len() != count {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "a mixture of {count} components needs {count} fractions, got {}",
                z.len()
            ),
        ));
    }

    let molar_mass = molar_mass_of(&mixture, &z)?;
    let transport = phase_transport(
        &mixture,
        &ideal_gas,
        phase_kind(kind),
        stream.t,
        stream.p,
        &z,
    )?;
    let density = mass_density(&mixture, stream.t, stream.p, &z, root, molar_mass)?;
    let cp_mass = molar_enthalpy_entropy(&mixture, &ideal_gas, stream.t, stream.p, &z, root)?
        .cp
        .value
        / molar_mass;

    Ok(PhaseView {
        z,
        n: stream.n * share,
        t: stream.t,
        p: stream.p,
        molar_mass,
        density,
        cp_mass,
        kind: phase_kind(kind),
        transport,
    })
}

/// The three kinds `eos.phase_transport` dispatches on, from the label rule's three.
pub const fn phase_kind(label: PhaseLabel) -> PhaseKind {
    match label {
        PhaseLabel::Gas => PhaseKind::Gas,
        PhaseLabel::Oil => PhaseKind::Oil,
        PhaseLabel::Aqueous => PhaseKind::Aqueous,
    }
}

/// The mixture's molar mass at a stated composition, kg/mol.
pub fn molar_mass_of(mixture: &Mixture, z: &[f64]) -> Result<f64> {
    let mut total = 0.0;
    for (index, (component, fraction)) in mixture.components().iter().zip(z).enumerate() {
        let mass = component.molar_mass.ok_or_else(|| {
            // `Component` carries no name by design - it is a caller-supplied record - so the
            // refusal names the position and the caller's own list names the substance.
            AzothError::property_unavailable(
                format!("component {index}"),
                "molar_mass",
                "a component of this mixture carries no molar mass, so the mixture has none",
            )
        })?;
        total += fraction * mass;
    }
    Ok(total)
}

/// `getDensity("kg/m3")`: the cubic's volume at the phase's own root, with the Peneloux shift.
///
/// **The root is the flash's own, and it is not re-derived.** A phase of a split sits on the
/// root the split put it on; asking the cubic again from the phase's composition alone can
/// answer a different one, which is the divergence [`Stream::from_side`] exists for.
pub fn mass_density(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    root: f64,
    molar_mass: f64,
) -> Result<f64> {
    let volume = pr_molar_volume(root, t, p)?.v.value;
    let shift = mixture.volume_shift(z);
    Ok(pr_mass_density(
        kilograms_per_mole(molar_mass),
        azoth_core::units::cubic_meters_per_mole(volume - shift),
    )?
    .rho
    .value)
}

/// A fluid resolves, which every caller needs before it can build a system.
pub fn resolve(names: &[&str]) -> Result<(Mixture, azoth_eos::IdealGasModel)> {
    databank::mixture_of(names, Cubic::Pr, None)
}

/// A pressure in pascals, which the flash and the correlations take.
pub fn pascals_of(value: f64) -> Pressure {
    pascals(value)
}
