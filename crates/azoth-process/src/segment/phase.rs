//! One phase of a flashed system: the class's `PhaseInterface`, as a value.
//!
//! `RateBasedPackedColumn` holds two `SystemInterface`s - a gas one and a liquid one - and
//! reads a *phase* out of each. This module is that read: flash the system, pick the phase the
//! class would pick, and answer the state and the transport the segment model asks for.

use azoth_core::units::{Pressure, ThermodynamicTemperature, kilograms_per_mole, pascals};
use azoth_core::{AzothError, Result};
use azoth_eos::hydrate_inhibitor_wt::{PhaseLabel, label};
use azoth_eos::phase_transport::{PhaseKind, phase_transport};
use azoth_eos::results::{PhaseTransportResult, PtFlashResult};
use azoth_eos::{
    Cubic, IdealGasModel, Mixture, Phase, databank, molar_enthalpy_entropy, pr_mass_density,
    pr_molar_volume, pt_flash,
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

    // **The phase, and it is the label rule that names it, not `beta`.** `getGasPhase` and
    // `getLiquidPhase` look for a `PhaseType`, and `PhaseEos.init` is what assigns one - the
    // volume ratio past `1.75` is the gas, and of the two liquids the one whose hydrocarbons
    // outweigh its aqueous components is the oil. It is the same rule `three_phase_separator`
    // and `pipe` already read, and it is public in `azoth-eos` for exactly this reason.
    let labelled = labelled_phases(stream, &mixture, &flash)?;

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
        flash.vapour_fraction.unwrap_or(1.0)
    } else {
        1.0 - flash.vapour_fraction.unwrap_or(0.0)
    };

    view(&mixture, &ideal_gas, stream, kind, z, root, share)
}

/// One phase, built from its label, its composition and its cubic root.
fn view(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    stream: &Stream,
    kind: PhaseLabel,
    z: Vec<f64>,
    root: f64,
    share: f64,
) -> Result<PhaseView> {
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

    let molar_mass = molar_mass_of(mixture, &z)?;
    let transport = phase_transport(mixture, ideal_gas, phase_kind(kind), stream.t, stream.p, &z)?;
    let density = mass_density(mixture, stream.t, stream.p, &z, root, molar_mass)?;
    let cp_mass = molar_enthalpy_entropy(mixture, ideal_gas, stream.t, stream.p, &z, root)?
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

/// The **vapour** of a stream that carries one, where a re-flash of a saturated composition can
/// answer a single *liquid* phase.
///
/// **The class reads a tray's own system and asks `hasPhaseType("gas")`, and a tray with a vapour
/// flow has a gas phase whatever a re-flash of its composition answers.** A pinned end is where
/// that matters: the reboiler boils at the temperature it is pinned to, so its vapour's
/// composition sits on the dew point *by construction*, and a re-flash of it can land on the
/// other side of the knife edge. [`phase_view`]'s `Pick::Gas` falls back to phase zero there and
/// would hand back the liquid; this builds the vapour from the flash's own **vapour root**.
///
/// # Errors
/// Whatever the flash, the label rule or the phase's own properties raise.
pub fn vapour_view(stream: &Stream) -> Result<PhaseView> {
    let (mixture, ideal_gas) = stream.mixture()?;
    let flash = pt_flash(&mixture, stream.t, stream.p, &stream.z)?;
    let labelled = labelled_phases(stream, &mixture, &flash)?;
    if let Some((_, z, root)) = labelled
        .iter()
        .find(|(kind, _, _)| *kind == PhaseLabel::Gas)
    {
        let share = if labelled.len() == 1 {
            1.0
        } else {
            flash.vapour_fraction.unwrap_or(1.0)
        };
        return view(
            &mixture,
            &ideal_gas,
            stream,
            PhaseLabel::Gas,
            z.clone(),
            *root,
            share,
        );
    }
    // A single phase the label rule did not call a gas, on a tray that has vapour all the same.
    view(
        &mixture,
        &ideal_gas,
        stream,
        PhaseLabel::Gas,
        stream.z.clone(),
        flash.z_vapour,
        1.0,
    )
}

/// The **liquid** of a stream that carries one, on [`vapour_view`]'s reasoning and measured on
/// the same tray.
///
/// **The capture's reboiler is the case.** Its liquid composition flashes all-vapour in one
/// implementation and all-liquid in the other, so `Pick::Liquid`'s ladder falls through to the
/// gas and answers the *vapour's* density - `44.17832800507808` against the oil's
/// `435.74535196919646` on `spec_top_flow_rate`'s state - which is a factor of `208` on that
/// tray's flood and moves the report's maximum from `11.67` to `2427.5`.
///
/// # Errors
/// Whatever the flash, the label rule or the phase's own properties raise.
pub fn liquid_view(stream: &Stream) -> Result<PhaseView> {
    let (mixture, ideal_gas) = stream.mixture()?;
    let flash = pt_flash(&mixture, stream.t, stream.p, &stream.z)?;
    let labelled = labelled_phases(stream, &mixture, &flash)?;
    if let Some((kind, z, root)) = labelled
        .iter()
        .find(|(kind, _, _)| *kind != PhaseLabel::Gas)
    {
        let share = if labelled.len() == 1 {
            1.0
        } else {
            1.0 - flash.vapour_fraction.unwrap_or(0.0)
        };
        return view(&mixture, &ideal_gas, stream, *kind, z.clone(), *root, share);
    }
    // A single *vapour* phase where the caller knows the tray carries liquid: the class's own
    // system would have its oil phase here, and the flash's liquid root is still its own.
    view(
        &mixture,
        &ideal_gas,
        stream,
        PhaseLabel::Oil,
        stream.z.clone(),
        flash.z_liquid,
        1.0,
    )
}

/// The stream's phases, labelled by the class's own rule and in the order its system carries
/// them: **the gas first, then the liquids**. That order is what makes [`Pick::Gas`] the spelling
/// of `getPhase(0)` - the class indexes its phase array, and the array is gas-first.
fn labelled_phases(
    stream: &Stream,
    mixture: &Mixture,
    flash: &PtFlashResult,
) -> Result<Vec<(PhaseLabel, Vec<f64>, f64)>> {
    labelled_phases_at(mixture, stream.t, stream.p, &stream.z, flash)
}

/// [`labelled_phases`] at a state rather than at a stream.
///
/// **The pair a capacity read needs is the state's and not the stream's.** `process.gas_scrubber`
/// asks about the system `run` left behind - the outlet's pressure, the outlet's temperature and
/// the *feed's* composition - which no one of its three streams carries together.
pub(crate) fn labelled_phases_at(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    flash: &PtFlashResult,
) -> Result<Vec<(PhaseLabel, Vec<f64>, f64)>> {
    let reduced = mixture.reduced_parameters(t, p)?;
    Ok(match flash.phase {
        // **A single phase has the system's composition, and the flash's own `x` or `y` is not
        // it.** The flash answers a trial phase there - `column/tray.rs` records the
        // measurement - so the one phase *is* the system.
        Phase::AllVapour => vec![(PhaseLabel::Gas, z.to_vec(), flash.z_vapour)],
        Phase::AllLiquid | Phase::Trivial => {
            vec![(
                label(mixture, &reduced, z, flash.z_liquid)?,
                z.to_vec(),
                flash.z_liquid,
            )]
        }
        Phase::TwoPhase => {
            if flash.vapour_fraction.is_none() {
                return Err(AzothError::invalid_input(
                    "flash",
                    "a two-phase answer with no vapour fraction, so neither phase's share is known",
                ));
            }
            vec![
                (PhaseLabel::Gas, flash.y.clone(), flash.z_vapour),
                (
                    label(mixture, &reduced, &flash.x, flash.z_liquid)?,
                    flash.x.clone(),
                    flash.z_liquid,
                ),
            ]
        }
    })
}

/// The whole system's mass and volume, per mole of the system: `(kg/mol, m³/mol)`.
///
/// **The sum over the phases, and every capacity-limit reading needs it.** `getDensity("kg/m3")`
/// and `getFlowRate("m3/sec")` are *system* quantities in NeqSim, not phase ones, and on a
/// two-phase outlet neither is either phase's: the capture's over-cooled overhead row reads a
/// system density of `15.34351809893286` beside a gas phase of `14.4273646532561`.
fn system_per_mole(stream: &Stream) -> Result<(f64, f64)> {
    let (mixture, _) = stream.mixture()?;
    let flash = pt_flash(&mixture, stream.t, stream.p, &stream.z)?;
    let phases = labelled_phases(stream, &mixture, &flash)?;
    let mut mass = 0.0;
    let mut volume = 0.0;
    for (kind, z, root) in &phases {
        let share = if phases.len() == 1 {
            1.0
        } else if *kind == PhaseLabel::Gas {
            flash.vapour_fraction.unwrap_or(1.0)
        } else {
            1.0 - flash.vapour_fraction.unwrap_or(0.0)
        };
        let molar_mass = molar_mass_of(&mixture, z)?;
        let density = mass_density(&mixture, stream.t, stream.p, z, *root, molar_mass)?;
        mass += share * molar_mass;
        volume += share * molar_mass / density;
    }
    Ok((mass, volume))
}

/// The stream's mass density, kg/m³ - `SystemInterface.getDensity("kg/m3")`.
///
/// # Errors
/// Whatever the flash, the label rule or a phase's own density raises.
pub fn system_mass_density(stream: &Stream) -> Result<f64> {
    let (mass, volume) = system_per_mole(stream)?;
    Ok(mass / volume)
}

/// The whole stream's volumetric flow, m³/s - `SystemInterface.getFlowRate("m3/sec")`.
///
/// **The identity the capacity capture holds**: `n*M/rho` equals NeqSim's own `m3/sec` on every
/// one of its rows, and this is that identity summed over the phases.
///
/// # Errors
/// Whatever the flash, the label rule or a phase's own density raises.
pub fn system_volumetric_flow(stream: &Stream) -> Result<f64> {
    let (_, volume) = system_per_mole(stream)?;
    Ok(stream.n * volume)
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
