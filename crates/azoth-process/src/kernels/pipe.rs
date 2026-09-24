//! `unit_ops.pipe` - a line's pressure drop, with its outlet pressure solved from the
//! geometry rather than stated.

use azoth_core::units::{Length, Pressure, pascals};
use azoth_core::{AzothError, Result};
use azoth_eos::hydrate_inhibitor_wt::PhaseLabel;
use azoth_eos::{aqueous_viscosity, hydrate_inhibitor_wt, pt_flash, viscosity};

use crate::stream::Stream;

/// The outlet's pressure tolerance and pass cap, from `AdiabaticPipe.run`.
///
/// `1e-2` bar - the class's own figure, and `1e-3` in the Pa this crate works in.
const PRESSURE_TOLERANCE: f64 = 1.0e3;
/// `iter < 25`, the class's own cap.
const MAX_PASSES: usize = 25;
/// The Rackett-free gas constant, as `eos` uses it everywhere.
const R: f64 = 8.31446261815324;

/// What a line did: the outlet record, and the interior a hydraulic calculation reads.
pub struct PipeOutcome {
    /// The outlet record.
    pub outlet: Stream,
    /// The velocity the pressure drop was computed at, m/s.
    pub velocity: f64,
    /// The Reynolds number it was computed at.
    pub reynolds: f64,
    /// The Darcy friction factor.
    pub friction_factor: f64,
    /// Which of the three regimes the friction factor came from.
    pub regime: FlowRegime,
    /// The passes the outlet pressure took to settle.
    pub passes: usize,
}

/// The friction factor's regime, as `AdiabaticPipe` names them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlowRegime {
    /// Below `Re = 2300`: `64 / Re`.
    Laminar,
    /// The `2300`-`4000` band, where the class **interpolates** between the laminar and
    /// turbulent values rather than refusing.
    Transition,
    /// Above `4000`: the Haaland form.
    Turbulent,
}

impl FlowRegime {
    /// The label `AdiabaticPipe.getFlowRegime` reports.
    pub fn as_str(self) -> &'static str {
        match self {
            FlowRegime::Laminar => "laminar",
            FlowRegime::Transition => "transition",
            FlowRegime::Turbulent => "turbulent",
        }
    }
}

/// Drop a stream's pressure along a line, solving the outlet pressure it implies.
///
/// **The outlet pressure is a function of the line, not an input** - which is what makes a
/// pipe a different shape of unit operation from the two-port machines. `AdiabaticPipe.run`
/// iterates `calcPressureOut()` against the state it is evaluating at until the two agree to
/// `1e-2` bar or twenty-five passes run out, and the *outlet* is a flash at the pressure it
/// settled on, at the inlet's temperature.
///
/// **Two equations, one per phase**, and the branch is the phase's label - NeqSim's own
/// `PhaseEos.init` rule, which `eos` ports for the hydrate dose pair. A **gas** takes the
/// compressible `P1^2 - P2^2` form; anything else takes Darcy-Weisbach. The two differ in
/// more than algebra: the gas form reads the cubic's `Z` and molar mass, while the liquid
/// form deliberately recomputes the velocity from the *physical-properties* density, because
/// "cubic EOS liquid volumes are inaccurate for polar fluids (e.g. water)".
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the geometry is not positive, and whatever
/// the flash, the density, the viscosity or the label refuses.
pub fn pipe(
    feed: &Stream,
    length: Length,
    diameter: Length,
    roughness: Length,
) -> Result<PipeOutcome> {
    for (name, value) in [
        ("length", length.value),
        ("diameter", diameter.value),
        ("roughness", roughness.value),
    ] {
        // `NaN` is not a length either, and this says so without clippy reading the test as
        // a negated comparison on a partial order.
        if value.is_nan() || value <= 0.0 {
            return Err(AzothError::invalid_input(
                name,
                format!("a pipe's {name} is a positive length, and {value} is not"),
            ));
        }
    }

    // The ideal-gas model is the Pitzer departure's, which this kernel does not read: the
    // outlet's enthalpy is a `TP` state, and both branches work from the flash's own `Z`.
    let (mixture, _ideal_gas) = feed.mixture()?;
    let relative_roughness = roughness.value / diameter.value;
    let area = std::f64::consts::PI / 4.0 * diameter.value * diameter.value;

    // The inlet's own pressure, which the compressible form's `P1^2` term holds fixed while
    // the state it is evaluated at walks down to the outlet.
    let inlet_pressure = feed.p.value;
    let mut pressure = feed.p.value;
    let mut passes = 0;
    let mut state = State::at(
        &mixture,
        feed,
        pressure,
        area,
        diameter.value,
        relative_roughness,
    )?;

    loop {
        passes += 1;
        let next = state.outlet_pressure(inlet_pressure, length.value, diameter.value);
        let settled = (next - pressure).abs() <= PRESSURE_TOLERANCE;
        pressure = next;
        if settled || passes >= MAX_PASSES {
            break;
        }
        state = State::at(
            &mixture,
            feed,
            pressure,
            area,
            diameter.value,
            relative_roughness,
        )?;
    }

    // The outlet is a flash at the pressure the loop settled on, at the inlet's temperature:
    // `run` ends with a `TPflash` over the system it has been moving, and never sets a
    // temperature.
    Ok(PipeOutcome {
        outlet: Stream::from_pt(
            feed.components.clone(),
            feed.z.clone(),
            feed.n,
            pascals(pressure),
            feed.t,
        )?,
        velocity: state.velocity,
        reynolds: state.reynolds,
        friction_factor: state.friction,
        regime: state.regime,
        passes,
    })
}

/// One pass's state, and the outlet pressure it implies.
struct State {
    /// The velocity the drop was computed at, m/s.
    velocity: f64,
    /// The Reynolds number.
    reynolds: f64,
    /// The Darcy friction factor.
    friction: f64,
    /// Which regime produced it.
    regime: FlowRegime,
    /// Whether the phase took the compressible branch.
    is_gas: bool,
    /// The compressible form's `Z`, molar mass and temperature.
    z_factor: f64,
    molar_mass: f64,
    temperature: f64,
    /// The molar flow the compressible form's leading term carries.
    moles: f64,
    /// The liquid form's density, from the physical properties.
    density: f64,
}

impl State {
    /// Evaluate the line at one pressure.
    #[allow(clippy::too_many_arguments)] // one per quantity the class reads, and there are six
    fn at(
        mixture: &azoth_eos::Mixture,
        feed: &Stream,
        pressure: f64,
        area: f64,
        diameter: f64,
        relative_roughness: f64,
    ) -> Result<State> {
        let p = pascals(pressure);
        let flash = pt_flash(mixture, feed.t, p, &feed.z)?;
        // Phase 0 of the flashed system, which is the liquid when there are two and the
        // single phase when there is one. `pt_flash` reports `x` for the liquid root and `y`
        // for the vapour one, so a gas state's composition is `y`.
        let (z_factor, composition) = match flash.phase {
            azoth_eos::Phase::TwoPhase => (flash.z_liquid, flash.x.clone()),
            azoth_eos::Phase::AllVapour => (flash.z_vapour, flash.y.clone()),
            azoth_eos::Phase::AllLiquid | azoth_eos::Phase::Trivial => {
                (flash.z_liquid, flash.x.clone())
            }
        };

        let reduced = mixture.reduced_parameters(feed.t, p)?;
        let label = hydrate_inhibitor_wt::label(mixture, &reduced, &composition, z_factor)?;
        let is_gas = label == PhaseLabel::Gas;

        // **Two densities, and the class uses them in different places.** The liquid
        // branch's *velocity* takes the physical-properties density, which applies the
        // Peneloux translation - and its *Reynolds number* takes `getKinematicViscosity()`,
        // which divides the same viscosity by the **untranslated cubic** `M/(Z R T / P)`.
        // Measured on the captured butane row, the two differ by 6%: `567.33` against
        // `604.0`, which is the whole of the Reynolds gap between this port and NeqSim
        // before the distinction is made.
        let density = density_at(mixture, feed, p, z_factor)?;
        let cubic_density =
            molar_mass_of(mixture, &feed.z)? / (z_factor * R * feed.t.value / pressure);
        // **The viscosity is the phase's, and NeqSim dispatches on the phase type.** A gas or
        // a hydrocarbon liquid takes `PFCTViscosityMethodHeavyOil` - `eos.viscosity` - and an
        // **aqueous** phase takes `WaterPhysicalProperties`, whose correlation is the liquid
        // `Viscosity` class: `eos.aqueous_viscosity`. Water at 300 K is `8.5510e-4` Pa s
        // through the second and `5.3097e-4` through the first.
        let mu = match label {
            PhaseLabel::Aqueous => {
                aqueous_viscosity::aqueous_viscosity(mixture, feed.t, p, &composition)?
                    .viscosity
                    .value
            }
            PhaseLabel::Gas | PhaseLabel::Oil => {
                viscosity::viscosity(mixture, feed.t, p, &composition)?
                    .mu
                    .value
            }
        };
        let kinematic = mu / cubic_density;

        // The gas branch takes the velocity from the cubic's total volume; the liquid branch
        // recomputes it from the mass flow and the physical-properties density. They are the
        // same number in the limit and not otherwise, which is why both are written out.
        let molar_mass = molar_mass_of(mixture, &feed.z)?;
        let velocity = if is_gas {
            let volume_rate = feed.n * z_factor * R * feed.t.value / pressure;
            volume_rate / area
        } else {
            feed.n * molar_mass / density / area
        };

        let reynolds = velocity * diameter / kinematic;
        let (friction, regime) = friction_factor(reynolds, relative_roughness);

        Ok(State {
            velocity,
            reynolds,
            friction,
            regime,
            is_gas,
            z_factor,
            molar_mass,
            temperature: feed.t.value,
            moles: feed.n,
            density,
        })
    }

    /// The outlet pressure this state implies, in Pa.
    fn outlet_pressure(&self, inlet_pressure: f64, length: f64, diameter: f64) -> f64 {
        if self.is_gas {
            // The compressible general flow equation, `P1^2 - P2^2`, exactly as `run` writes
            // it - including that its leading term is the *mass* flow, `4 n M / pi`.
            let leading = (4.0 * self.moles * self.molar_mass / std::f64::consts::PI).powi(2);
            let dp = leading * self.friction * length * self.z_factor * R / self.molar_mass
                * self.temperature
                / diameter.powi(5);
            // No gravity term: the palette declares no elevations, and `run`'s
            // `dp_gravity` is `rho g (z_in - z_out)` - zero for a level line. A line with
            // a rise is a column the entry cannot describe, which is stated in its spec
            // rather than defaulted.
            (inlet_pressure.powi(2) - dp).max(0.0).sqrt()
        } else {
            let dp =
                self.friction * length / diameter * self.density * self.velocity * self.velocity
                    / 2.0;
            inlet_pressure - dp
        }
    }
}

/// The Darcy friction factor and the regime it came from, as `calcWallFrictionFactor` states
/// them.
///
/// **The transition band is interpolated, not refused.** Between `2300` and `4000` the class
/// linearly blends the laminar value at `2300` with the Haaland value at `4000`, and this
/// reproduces it: where the class draws a straight line between two points, a port that
/// refused the band would be answering a smaller question than the class does.
fn friction_factor(reynolds: f64, relative_roughness: f64) -> (f64, FlowRegime) {
    let re = reynolds.abs();
    if re < 1e-10 {
        return (0.0, FlowRegime::Laminar);
    }
    if re < 2300.0 {
        return (64.0 / re, FlowRegime::Laminar);
    }
    if re < 4000.0 {
        let laminar = 64.0 / 2300.0;
        let turbulent = haaland(4000.0, relative_roughness);
        return (
            laminar + (turbulent - laminar) * (re - 2300.0) / 1700.0,
            FlowRegime::Transition,
        );
    }
    (haaland(re, relative_roughness), FlowRegime::Turbulent)
}

/// Haaland's explicit approximation to Colebrook, as the class writes it.
fn haaland(reynolds: f64, relative_roughness: f64) -> f64 {
    let inner = 6.9 / reynolds + (relative_roughness / 3.7).powf(1.11);
    (1.0 / (-1.8 * inner.log10())).powi(2)
}

/// The physical-properties density at a state, which is the cubic's volume with the volume
/// translation applied.
fn density_at(
    mixture: &azoth_eos::Mixture,
    feed: &Stream,
    p: Pressure,
    z_factor: f64,
) -> Result<f64> {
    let volume =
        azoth_eos::pr_molar_volume(z_factor, feed.t, p)?.v.value - mixture.volume_shift(&feed.z);
    Ok(azoth_eos::pr_mass_density(
        azoth_core::units::kilograms_per_mole(molar_mass_of(mixture, &feed.z)?),
        azoth_core::units::cubic_meters_per_mole(volume),
    )?
    .rho
    .value)
}

/// The mixture's molar mass, kg/mol.
fn molar_mass_of(mixture: &azoth_eos::Mixture, z: &[f64]) -> Result<f64> {
    let mut total = 0.0;
    for (component, zi) in mixture.components().iter().zip(z) {
        let mass = component.molar_mass.ok_or_else(|| {
            AzothError::property_unavailable(
                "a pipe's fluid",
                "molar_mass",
                "a component of this stream carries no molar mass, so the mixture has none",
            )
        })?;
        total += zi * mass;
    }
    Ok(total)
}
