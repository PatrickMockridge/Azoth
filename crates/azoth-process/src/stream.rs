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
    Cubic, IdealGasModel, Mixture, Phase, databank, ph_flash, pr_mass_density, pr_molar_volume,
    ps_flash, pt_flash,
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
        let (h, _) = ph_flash::enthalpy_at(&mixture, &ideal_gas, t, p, &z)?;
        Ok(Stream {
            components,
            z,
            n,
            p,
            t,
            h: joules_per_mole(h),
        })
    }

    /// The molar entropy at the stream's own state, J/(mol·K).
    ///
    /// **Derived, not carried.** The record has five fields and this is not one of them,
    /// for the reason `s` is a function of `(T, P, z)` exactly as `h` is: a sixth field
    /// would have to be kept in step with the other five at every port of every unit
    /// operation, and could then disagree with them.
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
        })
    }
}
