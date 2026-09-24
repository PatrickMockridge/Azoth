//! `unit_ops.three_phase_separator` - a multiphase flash, split three ways.

use azoth_core::units::{
    Power, Pressure, ThermodynamicTemperature, joules_per_mole, kelvins, pascals,
};
use azoth_core::{AzothError, Result};
use azoth_eos::hydrate_inhibitor_wt::{PhaseLabel, label};
use azoth_eos::{
    IdealGasModel, Mixture, RootSide, molar_enthalpy_entropy, ph_flash, tp_multiflash,
};

use crate::stream::Stream;

/// The six entrainment directions, **in the order `run` applies them**.
///
/// The order is load-bearing rather than tidy: each transfer moves a fraction of the *current*
/// from-phase's moles, so a later pair reads what an earlier one left behind. NeqSim runs
/// `gas -> aqueous`, `gas -> oil`, `oil -> aqueous`, `oil -> gas`, `aqueous -> gas`,
/// `aqueous -> oil`, and this does the same.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Entrainment {
    /// A fraction of the vapour's moles carried into the aqueous phase.
    pub gas_in_aqueous: f64,
    /// A fraction of the vapour's moles carried into the oil.
    pub gas_in_oil: f64,
    /// A fraction of the oil's moles carried into the aqueous phase.
    pub oil_in_aqueous: f64,
    /// A fraction of the oil's moles carried into the vapour.
    pub oil_in_gas: f64,
    /// A fraction of the aqueous phase's moles carried into the vapour.
    pub aqueous_in_gas: f64,
    /// A fraction of the aqueous phase's moles carried into the oil.
    pub aqueous_in_oil: f64,
}

/// One of the three outlets, as `run` names it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Slot {
    Vapour,
    Oil,
    Aqueous,
}

impl Slot {
    /// Which side of the cubic a phase with this label sits on.
    fn side(self) -> RootSide {
        match self {
            Slot::Vapour => RootSide::Vapour,
            Slot::Oil | Slot::Aqueous => RootSide::Liquid,
        }
    }
}

impl From<PhaseLabel> for Slot {
    fn from(label: PhaseLabel) -> Self {
        match label {
            PhaseLabel::Gas => Slot::Vapour,
            PhaseLabel::Oil => Slot::Oil,
            PhaseLabel::Aqueous => Slot::Aqueous,
        }
    }
}

impl Entrainment {
    /// The six fractions with the pair they act on, in `run`'s order.
    fn directions(&self) -> [(Slot, Slot, f64); 6] {
        [
            (Slot::Vapour, Slot::Aqueous, self.gas_in_aqueous),
            (Slot::Vapour, Slot::Oil, self.gas_in_oil),
            (Slot::Oil, Slot::Aqueous, self.oil_in_aqueous),
            (Slot::Oil, Slot::Vapour, self.oil_in_gas),
            (Slot::Aqueous, Slot::Vapour, self.aqueous_in_gas),
            (Slot::Aqueous, Slot::Oil, self.aqueous_in_oil),
        ]
    }

    /// The fractions beside their names, so a refusal can name the one that is not one.
    fn named(&self) -> [(&'static str, f64); 6] {
        [
            ("gas_in_aqueous", self.gas_in_aqueous),
            ("gas_in_oil", self.gas_in_oil),
            ("oil_in_aqueous", self.oil_in_aqueous),
            ("oil_in_gas", self.oil_in_gas),
            ("aqueous_in_gas", self.aqueous_in_gas),
            ("aqueous_in_oil", self.aqueous_in_oil),
        ]
    }
}

/// One phase of the split.
#[derive(Debug, Clone, PartialEq)]
struct Split {
    /// Moles of each component **per mole of feed**, which is the basis the transfers move.
    amounts: Vec<f64>,
    /// The compressibility the cubic reported for it.
    z_factor: f64,
    /// Which of the three outlets it is, from `PhaseEos.init`'s label.
    slot: Slot,
}

impl Split {
    /// The share of the feed's moles in this phase.
    fn beta(&self) -> f64 {
        self.amounts.iter().sum()
    }

    /// The phase's composition, or all zero for a phase the transfers emptied.
    fn composition(&self) -> Vec<f64> {
        let beta = self.beta();
        if beta <= 0.0 {
            return vec![0.0; self.amounts.len()];
        }
        self.amounts.iter().map(|amount| amount / beta).collect()
    }
}

/// The phases a split reports, named the way NeqSim's outlet getters are.
#[derive(Debug, Clone, Default, PartialEq)]
struct Phases {
    vapour: Option<Split>,
    oil: Option<Split>,
    aqueous: Option<Split>,
}

impl Phases {
    fn get(&self, slot: Slot) -> Option<&Split> {
        match slot {
            Slot::Vapour => self.vapour.as_ref(),
            Slot::Oil => self.oil.as_ref(),
            Slot::Aqueous => self.aqueous.as_ref(),
        }
    }

    fn get_mut(&mut self, slot: Slot) -> Option<&mut Split> {
        match slot {
            Slot::Vapour => self.vapour.as_mut(),
            Slot::Oil => self.oil.as_mut(),
            Slot::Aqueous => self.aqueous.as_mut(),
        }
    }

    /// **Move each entrainment fraction, component by component, in `run`'s order.**
    ///
    /// `addPhaseFractionToPhase` computes `change_i = moles_i(from) * fraction` for every
    /// component and adds it to one phase and subtracts it from the other, so the transfer is
    /// proportional to the from-phase's *own* composition - it needs no equilibrium and moves
    /// no energy. It returns unchanged unless both phases are present and the fraction is
    /// above `1e-30`, which is the guard here.
    fn carry(&mut self, entrainment: &Entrainment) {
        for (from, to, fraction) in entrainment.directions() {
            if fraction < 1.0e-30 || self.get(from).is_none() || self.get(to).is_none() {
                continue;
            }
            let moved: Vec<f64> = {
                let source = self.get(from).expect("checked above");
                source
                    .amounts
                    .iter()
                    .map(|amount| amount * fraction)
                    .collect()
            };
            for (amount, change) in self
                .get_mut(from)
                .expect("checked above")
                .amounts
                .iter_mut()
                .zip(&moved)
            {
                *amount -= change;
            }
            for (amount, change) in self
                .get_mut(to)
                .expect("checked above")
                .amounts
                .iter_mut()
                .zip(&moved)
            {
                *amount += change;
            }
        }
    }
}

/// Flash a feed into vapour, oil and aqueous outlets.
///
/// `run` reduces the pressure, **enables `multiPhaseCheck` for the flash and disables it
/// again**, and splits: a `TPflash` at the feed's temperature, or a `PHflash` at the enthalpy
/// a `heatInput` implies. The three phases are the flash's own, named by `PhaseEos.init`'s
/// label rule - a volume ratio past `1.75` is the gas, and of the two liquids the one whose
/// hydrocarbons outweigh its aqueous components is the oil.
///
/// `pressure_drop` is not a throttling: the vessel holds the feed's temperature.
///
/// # Errors
/// [`AzothError::InvalidInput`] if the pressure drop would take the outlet to a non-positive
/// pressure, or if a fraction is outside `[0, 1]`; [`AzothError::OutOfRange`] if no
/// temperature carries the enthalpy a `heat_input` asks for. The flash's own refusals pass
/// through.
pub fn three_phase_separator(
    feed: &Stream,
    pressure_drop: Pressure,
    heat_input: Option<Power>,
    entrainment: Entrainment,
) -> Result<(Stream, Stream, Stream)> {
    for (name, fraction) in entrainment.named() {
        if !(0.0..=1.0).contains(&fraction) {
            return Err(AzothError::invalid_input(
                name,
                format!("an entrainment fraction is a fraction, and {fraction} is not in [0, 1]"),
            ));
        }
    }
    let p_out = pascals(feed.p.value - pressure_drop.value);
    if p_out.value <= 0.0 {
        return Err(AzothError::invalid_input(
            "pressure_drop",
            format!(
                "a pressure drop of {} Pa takes a {} Pa inlet to {} Pa, which is not a pressure",
                pressure_drop.value, feed.p.value, p_out.value
            ),
        ));
    }

    let (mixture, ideal_gas) = feed.mixture()?;

    let (temperature, mut phases) = match heat_input.filter(|power| power.value != 0.0) {
        Some(power) => {
            // One second of flow, for the reason `Heater.run` gives: the duty is W and the
            // enthalpy it moves is per mole, so it divides by the same flow a case states in
            // mol/s.
            let target = feed.h.value + power.value / feed.n;
            let t = solve_temperature(&mixture, &ideal_gas, p_out, target, &feed.z)?;
            (t, split(&mixture, t, p_out, &feed.z)?)
        }
        None => (feed.t.value, split(&mixture, feed.t.value, p_out, &feed.z)?),
    };

    phases.carry(&entrainment);
    let temperature = kelvins(temperature);

    // **The outlet is the phase, except where the class re-runs it - and here that is most of
    // them.** `run` leaves the vapour as the flash's gas phase unless `oilInGas` or
    // `aqueousInGas` is set, and re-runs the *oil and aqueous* outlets without a condition
    // wherever those phases exist (`liquidOutStream.run(id)` and `waterOutStream.run(id)`).
    // A stream `run` is a `TPflash` of that phase's composition at its own temperature and
    // pressure, which is a different state from the phase's own root wherever the composition
    // separates once its parent's other phases are gone - see `Stream::from_side`.
    let reflashed_vapour = entrainment.oil_in_gas != 0.0 || entrainment.aqueous_in_gas != 0.0;
    Ok((
        outlet(
            feed,
            phases.vapour.as_ref(),
            p_out,
            temperature,
            reflashed_vapour,
        )?,
        outlet(feed, phases.oil.as_ref(), p_out, temperature, true)?,
        outlet(feed, phases.aqueous.as_ref(), p_out, temperature, true)?,
    ))
}

/// One outlet: the phase, or the state its composition settles on where the class re-runs it.
///
/// **An absent phase is a zero-flow stream** at the outlet state. NeqSim answers with a
/// `1e-20 kg/hr` clone of the whole system instead, whose composition need not sum to one and
/// whose enthalpy is that clone's own; the capture carries those rows and the case declares
/// them rather than pinning them.
fn outlet(
    feed: &Stream,
    phase: Option<&Split>,
    p_out: Pressure,
    t: ThermodynamicTemperature,
    reflashed: bool,
) -> Result<Stream> {
    let Some(phase) = phase else {
        return Stream::from_pt(feed.components.clone(), feed.z.clone(), 0.0, p_out, t);
    };
    let n = feed.n * phase.beta();
    let composition = phase.composition();
    if reflashed {
        settled(feed, composition, n, p_out, t, phase.slot.side())
    } else {
        Stream::from_side(
            feed.components.clone(),
            composition,
            n,
            p_out,
            t,
            phase.slot.side(),
        )
    }
}

/// **The state a composition settles on at `(p, t)`, which is what a stream `run` reports.**
///
/// `Stream::from_pt` is that state for a composition whose own flash resolves, and it is what
/// this reaches for whenever two phases or more come back. Where **one** phase comes back the
/// state *is* that phase, and its root is the answer - the same arithmetic through
/// [`Stream::from_side`], which is the other branch of the same rule rather than a second one.
///
/// **That branch is not decoration: the two-phase flash cannot always resolve the question.**
/// Measured on this model's own aqueous outlet - `0.999999986` water with `1e-8` methane and
/// `1e-15` n-butane at 300 K and 20 bar - `eos.pt_flash` reports a **two-phase split** with a
/// vapour fraction of `2.6e-15` and a vapour composition summing to `0.9602200424212375`, which
/// is not a composition, and `molar_enthalpy_entropy` refuses it: `Stream::from_pt` fails on a
/// state the class asks for. `eos.tp_multiflash` answers the same question with one liquid
/// phase, which is what NeqSim's own flash reports - its outlet's enthalpy is `-44702.486`
/// J/mol against the liquid root's `-44702.552`, the two libraries' ideal-gas offset. The
/// defect is `azoth-eos`'s and is recorded rather than hidden: this asks the general flash the
/// question and takes the state it names.
fn settled(
    feed: &Stream,
    composition: Vec<f64>,
    n: f64,
    p_out: Pressure,
    t: ThermodynamicTemperature,
    side: RootSide,
) -> Result<Stream> {
    let (mixture, _) = feed.mixture()?;
    let flash = tp_multiflash(&mixture, t, p_out, &composition)?;
    if flash.phase_count > 1 {
        return Stream::from_pt(feed.components.clone(), composition, n, p_out, t);
    }
    Stream::from_side(feed.components.clone(), composition, n, p_out, t, side)
}

/// The phases a multiphase flash at `(t, p_out)` reports, named the way `run` names them.
fn split(mixture: &Mixture, t: f64, p_out: Pressure, z: &[f64]) -> Result<Phases> {
    let flash = tp_multiflash(mixture, kelvins(t), p_out, z)?;
    let reduced = mixture.reduced_parameters(kelvins(t), p_out)?;
    let mut phases = Phases::default();
    for index in 0..flash.phase_count as usize {
        let slot = Slot::from(label(
            mixture,
            &reduced,
            &flash.x[index],
            flash.z_factor[index],
        )?);
        let phase = Split {
            amounts: flash.x[index]
                .iter()
                .map(|fraction| fraction * flash.beta[index])
                .collect(),
            z_factor: flash.z_factor[index],
            slot,
        };
        *match slot {
            Slot::Vapour => &mut phases.vapour,
            Slot::Oil => &mut phases.oil,
            Slot::Aqueous => &mut phases.aqueous,
        } = Some(phase);
    }
    Ok(phases)
}

/// The temperature at which a multiphase flash carries the target molar enthalpy.
///
/// NeqSim's `PHflash` iterates on its own flash; this solves the same equation on
/// [`tp_multiflash`], which is the flash `run` enabled `multiPhaseCheck` for. The iteration
/// starts at the two-phase `PH` flash's own temperature and bisects a bracket that widens
/// from there: a mixture's enthalpy is monotone in the temperature at fixed pressure and
/// composition, so a bracket that straddles the target holds the root.
fn solve_temperature(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    p_out: Pressure,
    target: f64,
    z: &[f64],
) -> Result<f64> {
    let at = |t: f64| -> Result<f64> {
        let flash = tp_multiflash(mixture, kelvins(t), p_out, z)?;
        let mut total = 0.0;
        for index in 0..flash.phase_count as usize {
            let state = molar_enthalpy_entropy(
                mixture,
                ideal_gas,
                kelvins(t),
                p_out,
                &flash.x[index],
                flash.z_factor[index],
            )?;
            total += flash.beta[index] * state.h.value;
        }
        Ok(total)
    };

    let seed = ph_flash(mixture, ideal_gas, p_out, joules_per_mole(target), z)?
        .temperature
        .value;
    let (mut low, mut high) = (seed - 25.0, seed + 25.0);
    let (mut h_low, mut h_high) = (at(low)?, at(high)?);
    for _ in 0..40 {
        if h_low <= target && target <= h_high {
            break;
        }
        let step = 0.5 * (high - low);
        if h_low > target {
            if low - step <= 0.0 {
                break;
            }
            low -= step;
            h_low = at(low)?;
        } else {
            high += step;
            h_high = at(high)?;
        }
    }
    if !(h_low <= target && target <= h_high) {
        return Err(AzothError::out_of_range(
            "heat_input",
            target,
            "no temperature in the bracket carries this enthalpy",
        ));
    }

    for _ in 0..200 {
        let middle = 0.5 * (low + high);
        let h_middle = at(middle)?;
        if (h_middle - target).abs() < 1.0e-9 {
            return Ok(middle);
        }
        if h_middle < target {
            low = middle;
            h_low = h_middle;
        } else {
            high = middle;
            h_high = h_middle;
        }
    }
    let _ = (h_low, h_high);
    Ok(0.5 * (low + high))
}
