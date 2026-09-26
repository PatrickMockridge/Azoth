//! A stream: the value a channel carries.
//!
//! The palette's shared record is `n, z, P, T, h` — molar flow, composition,
//! pressure, temperature and molar enthalpy. This module is that value: the thing
//! a unit operation's kernel reads on its inlets and writes on its outlets.

use azoth_core::units::{
    MolarEnergy, MolarMass, Pressure, ThermodynamicTemperature, cubic_meters_per_mole,
    joules_per_mole, kilograms_per_mole,
};
use azoth_core::{AzothError, Result};
use azoth_eos::{
    Cubic, IdealGasModel, Mixture, Phase, RootSide, databank, molar_enthalpy_entropy, ph_flash,
    pr_mass_density, pr_molar_volume, ps_flash, pt_flash,
};

/// A material stream on the shared record.
#[derive(Debug, Clone)]
pub struct Stream {
    /// Component names, resolved against the databank with no keycard overlay.
    pub components: Vec<String>,
    /// Mole fractions, one per component, summing to one.
    pub z: Vec<f64>,
    /// Molar flow in mol/s.
    pub n: f64,
    pub p: Pressure,
    pub t: ThermodynamicTemperature,
    /// Molar enthalpy in J/mol, on the databank's reference.
    pub h: MolarEnergy,
    /// The vapour fraction of the state, when the construction knew one.
    ///
    /// **Carried, not derived, and this is the one field where that is the honest choice.** `s`,
    /// `cp` and the density are functions of `(T, P, z)` a caller can cheaply recompute - see
    /// [`Stream::entropy`] - so they are not fields. The vapour fraction is a function of the same
    /// three, but the accessor that would compute it *refuses* a two-phase stream rather than
    /// answering ([`Stream::single_phase_root`]), and a two-phase stream is exactly the one whose
    /// vapour fraction a reader wants. So whichever constructor already ran the flash keeps its
    /// answer: [`Stream::from_pt`] and [`Stream::from_ph`] take [`vapour_fraction`] of the flash and
    /// [`Stream::from_side`] reports the side it was given. `None` means the state is not known to
    /// be one phase and no split was solved, and it is not a `0.5` by default.
    pub vapour_fraction: Option<f64>,
}

/// The vapour fraction a flash established, from its phase and its `beta`.
///
/// **`beta` outside `[0, 1]` is not a vapour fraction and is not published.** A flash at a state
/// that is not two-phase converges a Rachford-Rice root that does not lie in the interval, and
/// `PtFlashResult` carries it with a warning rather than withholding it - its own doc says `None` is
/// the honest answer, "the flash extrapolates a split that does not exist, and reporting it would
/// invite a caller to use it". This is that policy at the boundary where the number becomes a
/// stream's field, and the shipped demo is where it was found: a subcooled methane feed flashed to
/// `1.9847` and a heater's outlet to `1.2878`, both of which a front end would have drawn as a
/// vapour fraction.
///
/// **Where the state *is* one phase the endpoint is a fact and is reported**: a separator's vapour
/// outlet is all vapour whatever its composition would settle on by itself, and a reader comparing
/// two streams wants that `1` rather than a blank.
#[must_use]
pub fn vapour_fraction(phase: Phase, beta: Option<f64>) -> Option<f64> {
    match phase {
        Phase::AllVapour => Some(1.0),
        Phase::AllLiquid => Some(0.0),
        Phase::TwoPhase => beta.filter(|fraction| (0.0..=1.0).contains(fraction)),
        // The solution is trivial: the same composition on both sides, so there is no split to
        // report and no phase to point at.
        Phase::Trivial => None,
    }
}

impl Stream {
    /// Resolve the fluid this stream names.
    pub fn mixture(&self) -> azoth_core::Result<(Mixture, IdealGasModel)> {
        let names: Vec<&str> = self.components.iter().map(String::as_str).collect();
        databank::mixture_of(&names, Cubic::Pr, None)
    }

    /// A stream at a known temperature, with its enthalpy computed at that state.
    pub fn from_pt(
        components: Vec<String>,
        z: Vec<f64>,
        n: f64,
        p: Pressure,
        t: ThermodynamicTemperature,
    ) -> azoth_core::Result<Self> {
        let names: Vec<&str> = components.iter().map(String::as_str).collect();
        let (mixture, ideal_gas) = databank::mixture_of(&names, Cubic::Pr, None)?;
        // The flash was already run for the enthalpy, so its vapour fraction is the one this
        // constructor does not have to pay for - and a two-phase feed, which is the case a reader
        // most wants it for, is one this constructor reaches every day.
        let (h, flash) = ph_flash::enthalpy_at(&mixture, &ideal_gas, t, p, &z)?;
        Ok(Stream {
            components,
            z,
            n,
            p,
            t,
            h: joules_per_mole(h),
            vapour_fraction: vapour_fraction(flash.phase, flash.beta),
        })
    }

    /// A stream at a known temperature and pressure, **on a stated side of the cubic**.
    ///
    /// **This is how an outlet that is *a phase* is built, and it is not the same question
    /// as [`Stream::from_pt`].** `from_pt` re-flashes the composition, which answers with
    /// the state that composition settles on *on its own*; a phase of a split is on the
    /// root the split put it on. Where the two differ, the difference is not small. On a
    /// water-bearing gas that leaves a tank, the outlet's composition re-flashes to two
    /// phases at a vapour fraction of `0.767` and an enthalpy of `-10043.91` J/mol, where
    /// the phase's own root gives `441.77` - and NeqSim's `setThermoSystemFromPhase`
    /// reports `441.83`.
    ///
    /// Where the composition *is* stable on its own the two agree, which is why this was
    /// invisible until a feed whose phase separates when its parent does not.
    ///
    /// # Errors
    /// Whatever the databank and the equation of state raise.
    pub fn from_side(
        components: Vec<String>,
        z: Vec<f64>,
        n: f64,
        p: Pressure,
        t: ThermodynamicTemperature,
        side: RootSide,
    ) -> Result<Self> {
        let names: Vec<&str> = components.iter().map(String::as_str).collect();
        let (mixture, ideal_gas) = databank::mixture_of(&names, Cubic::Pr, None)?;
        let reduced = mixture.reduced_parameters(t, p)?;
        let root = mixture.phase_state(&reduced, &z, side)?.z;
        let state = molar_enthalpy_entropy(&mixture, &ideal_gas, t, p, &z, root)?;
        Ok(Stream {
            components,
            z,
            n,
            p,
            t,
            h: joules_per_mole(state.h.value),
            // The side *is* the phase, so this is not a flash's estimate of one: a vapour outlet
            // is all vapour whatever its composition would settle on by itself.
            vapour_fraction: Some(match side {
                RootSide::Vapour => 1.0,
                RootSide::Liquid => 0.0,
            }),
        })
    }

    /// The molar entropy at the stream's own state, J/(mol·K).
    ///
    /// **Derived, not carried.** The record has five fields and this is not one of them,
    /// for the reason `s` is a function of `(T, P, z)` exactly as `h` is: a sixth field
    /// would have to be kept in step with the other five at every port of every unit
    /// operation, and could then disagree with them.
    /// The molar heat capacity at the stream's own state, J/(mol·K).
    ///
    /// **Derived, not carried**, for the reason [`Stream::entropy`] gives, and it is the same
    /// call: `eos.molar_enthalpy_entropy` returns `cp` beside `h` and `s`, so this is a read of
    /// the model the stream's enthalpy already comes from rather than a second surface.
    ///
    /// **It is not the ideal-gas capacity.** The departure is carried, and at the state
    /// `PlugFlowReactor`'s own capture uses it is `0.0535` of `31.5162` J/(mol·K). Comparing it
    /// against NeqSim's `getCp("J/molK")` gives `6.85e-4` relative, which is the library
    /// divergence `crates/azoth-process/tests/reactor.rs` pins and not a difference in kind.
    pub fn molar_heat_capacity(&self) -> Result<f64> {
        let (mixture, ideal_gas) = self.mixture()?;
        let root = self.single_phase_root(&mixture)?;
        let state = molar_enthalpy_entropy(&mixture, &ideal_gas, self.t, self.p, &self.z, root)?;
        Ok(state.cp.value)
    }

    pub fn entropy(&self) -> Result<f64> {
        let (mixture, ideal_gas) = self.mixture()?;
        let (s, _) = ps_flash::entropy_at(&mixture, &ideal_gas, self.t, self.p, &self.z)?;
        Ok(s)
    }

    /// The mixture's molar mass, kg/mol.
    ///
    /// # Errors
    /// [`azoth_core::AzothError::PropertyUnavailable`] if any component carries no molar
    /// mass, which is true of a component built from critical constants alone. A
    /// correlation that needs the mass refuses rather than defaulting it to zero.
    pub fn molar_mass(&self) -> Result<MolarMass> {
        let (mixture, _) = self.mixture()?;
        let mut total = 0.0;
        for (component, zi) in mixture.components().iter().zip(&self.z) {
            let mass = component.molar_mass.ok_or_else(|| {
                AzothError::property_unavailable(
                    self.components.join(" / "),
                    "molar_mass",
                    "a component of this stream carries no molar mass, so the mixture has none",
                )
            })?;
            total += zi * mass;
        }
        Ok(kilograms_per_mole(total))
    }

    /// The mass flow, kg/s: the molar flow times the molar mass.
    pub fn mass_flow(&self) -> Result<f64> {
        Ok(self.n * self.molar_mass()?.value)
    }

    /// The mass density at the stream's own state, kg/m³.
    ///
    /// # Errors
    /// [`azoth_core::AzothError::InvalidInput`] **if the stream is two-phase**. A stream
    /// that is two phases has no one density, and the two candidate conventions - the
    /// vapour's, the liquid's, or a homogeneous average - are three different answers
    /// rather than one. NeqSim's `Stream.getDensity()` reports phase 0's without saying
    /// so; this refuses, and a caller with a two-phase line has to say which density they
    /// mean.
    pub fn density(&self) -> Result<f64> {
        let (mixture, _) = self.mixture()?;
        let root = self.single_phase_root(&mixture)?;
        let v = pr_molar_volume(root, self.t, self.p)?.v;
        Ok(pr_mass_density(self.molar_mass()?, v)?.rho.value)
    }

    /// The mass density a *physical-properties* model reports, kg/m³: the cubic's volume
    /// with the Peneloux volume translation applied.
    ///
    /// **The second of NeqSim's three densities, and the one a hydraulic calculation
    /// reads.** [`Stream::density`] is the cubic's own `M/(Z R T / P)`; NeqSim's
    /// `getPhase(0).getDensity("kg/m3")` adds the translation, and it is *that* one
    /// `getPhysicalProperties().getDensity()` returns - the density
    /// `AdiabaticPipe.calcPressureOut` reads for its liquid branch and, through the
    /// kinematic viscosity, for both.
    ///
    /// Measured on the captured states, the translation moves the density `6.4%` for
    /// liquid n-butane (`601.26` untranslated against `565.04`) and `14%` for water
    /// (`848.23` against `983.93`), because water's databank row carries a Rackett
    /// compressibility the fallback correlation does not reproduce.
    ///
    /// # Errors
    /// [`azoth_core::AzothError::InvalidInput`] **if the stream is two-phase**, for the
    /// reason [`Stream::density`] gives: two phases have no one density.
    pub fn corrected_density(&self) -> Result<f64> {
        let (mixture, _) = self.mixture()?;
        let root = self.single_phase_root(&mixture)?;
        let v = pr_molar_volume(root, self.t, self.p)?.v.value;

        // The mixture's own rule, `v_corr = v - sum_i x_i c_i`, which lives on `Mixture`
        // because the components' shifts do - `databank::Entry::component` is where each
        // one is filled in.
        let shift = mixture.volume_shift(&self.z);
        Ok(
            pr_mass_density(self.molar_mass()?, cubic_meters_per_mole(v - shift))?
                .rho
                .value,
        )
    }

    /// The cubic root a single-phase stream sits on.
    ///
    /// Shared by the two accessors that need one, so the refusal is stated once: a
    /// two-phase stream has two roots and nothing here may choose between them.
    fn single_phase_root(&self, mixture: &Mixture) -> Result<f64> {
        let flash = pt_flash(mixture, self.t, self.p, &self.z)?;
        match flash.phase {
            Phase::AllVapour => Ok(flash.z_vapour),
            Phase::AllLiquid | Phase::Trivial => Ok(flash.z_liquid),
            Phase::TwoPhase => Err(AzothError::invalid_input(
                "stream",
                format!(
                    "this stream is two phases at {} K and {} Pa, so it has no one density                      or viscosity; a caller has to say which phase's they mean",
                    self.t.value, self.p.value
                ),
            )),
        }
    }

    /// A stream at a known molar enthalpy, with its temperature solved for.
    pub fn from_ph(
        components: Vec<String>,
        z: Vec<f64>,
        n: f64,
        p: Pressure,
        h: MolarEnergy,
    ) -> azoth_core::Result<Self> {
        let names: Vec<&str> = components.iter().map(String::as_str).collect();
        let (mixture, ideal_gas) = databank::mixture_of(&names, Cubic::Pr, None)?;
        let r = ph_flash::ph_flash(&mixture, &ideal_gas, p, h, &z)?;
        Ok(Stream {
            components,
            z,
            n,
            p,
            t: r.temperature,
            h,
            vapour_fraction: vapour_fraction(r.phase, r.beta),
        })
    }
}
