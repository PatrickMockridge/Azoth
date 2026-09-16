//! A stream: the value a channel carries.
//!
//! The palette's shared record is `n, z, P, T, h` — molar flow, composition,
//! pressure, temperature and molar enthalpy. This module is that value: the thing
//! a unit operation's kernel reads on its inlets and writes on its outlets.

use azoth_core::units::{MolarEnergy, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_eos::{IdealGasModel, Mixture, databank, ph_flash};

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
        databank::mixture_of(&names, None)
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
        let (mixture, ideal_gas) = databank::mixture_of(&names, None)?;
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

    /// A stream at a known molar enthalpy, with its temperature solved for.
    pub fn from_ph(
        components: Vec<String>,
        z: Vec<f64>,
        n: f64,
        p: Pressure,
        h: MolarEnergy,
    ) -> azoth_core::Result<Self> {
        let names: Vec<&str> = components.iter().map(String::as_str).collect();
        let (mixture, ideal_gas) = databank::mixture_of(&names, None)?;
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
