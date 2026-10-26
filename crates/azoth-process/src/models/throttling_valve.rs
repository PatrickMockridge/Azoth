//! `process.throttling_valve` - the valve's kernel as a registered id.
//!
//! Spec: `specs/models/process/throttling_valve.toml`. The arithmetic is
//! [`crate::kernels::throttling_valve`]; what is here is the boundary a case and a
//! cross-impl test address.
//!
//! **The kernel needed no change to become a port.** NeqSim's `ThrottlingValve.run`
//! flashes `PHflash` at the inlet's enthalpy, which is what `Stream::from_ph` does, and
//! `acceptNegativeDP` defaults to `true` - so a stated outlet pressure above the inlet is
//! taken at face value, exactly as this does. The entry was wrong about its *parameters*,
//! not about its arithmetic.

use azoth_core::units::{MolarEnergy, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{AzothError, Result, Warning, apply_checks};

use crate::kernels::throttling_valve as kernel;
use crate::model_gen;
use crate::stream::Stream;
use crate::unported;
use crate::valve::{
    DEFAULT_OPENING_PERCENT, DEFAULT_XT, ValveCharacteristic, ValveSizingGeometry, sizing_report,
};

/// Result of `process.throttling_valve`.
#[derive(Debug, Clone, PartialEq)]
pub struct ThrottlingValveResult {
    /// Outlet molar flow, mol/s: the inlet's, since a valve adds no moles.
    pub outlet_n: f64,
    /// Outlet composition, which is the inlet's.
    pub outlet_z: Vec<f64>,
    /// Outlet pressure, which is the parameter the valve drops the stream to.
    pub outlet_p: Pressure,
    /// Outlet temperature, solved from the unchanged enthalpy at the outlet pressure.
    pub outlet_t: ThermodynamicTemperature,
    /// Outlet molar enthalpy: the inlet's on the isenthalpic branch, and the flashed value on
    /// the isothermal one.
    pub outlet_h: MolarEnergy,
    /// `getKv()`: the rated flow coefficient at 100% opening, stated or sized.
    pub kv: f64,
    /// `getCv()`: `Kv * 1.156`.
    pub cv: f64,
    /// `getCg()`: `Cv * 1360`.
    pub cg: f64,
    /// `calcValveSize`'s choked flag.
    pub choked: bool,
    /// `isGasValve()`, which selects the sizing's branch.
    pub gas_valve: bool,
    /// `calculateMolarFlow()`: the flow the adjusted `Kv` and the drop give, mol/s.
    pub calculated_molar_flow: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

/// Drop a stream to a lower pressure without heat or work.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the inlet pressure is not positive and
/// whatever the databank refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are fifteen
pub fn throttling_valve(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
    outlet_pressure: Pressure,
    valve_opening: Option<f64>,
    kv: Option<f64>,
    fouling_fraction: Option<f64>,
    valve_characterization: Option<&str>,
    xt: Option<f64>,
    allow_choked: Option<bool>,
    is_calc_pressure: Option<bool>,
    is_isothermal: Option<bool>,
    sizing_standard: Option<&str>,
) -> Result<ThrottlingValveResult> {
    let spec = &model_gen::THROTTLING_VALVE_SPEC;
    let mut warnings = Vec::new();
    let opening = valve_opening.unwrap_or(DEFAULT_OPENING_PERCENT);
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "outlet_pressure" => Some(outlet_pressure.value),
            "inlet_t" => Some(inlet_t.value),
            "valve_opening" => Some(opening),
            _ => None,
        },
        &mut warnings,
    )?;

    // **The eight other standards are second machines**, each refused by name from the
    // declaration's `[[unported]]` rows.
    match sizing_standard.unwrap_or("default") {
        "default" => {}
        "iec_60534" => return Err(unported::refuse("sizing_standard=iec_60534")),
        "iec_60534_full" => return Err(unported::refuse("sizing_standard=iec_60534_full")),
        "prod_choke" => return Err(unported::refuse("sizing_standard=prod_choke")),
        "sachdeva" => return Err(unported::refuse("sizing_standard=sachdeva")),
        "gilbert" => return Err(unported::refuse("sizing_standard=gilbert")),
        "baxendell" => return Err(unported::refuse("sizing_standard=baxendell")),
        "ros" => return Err(unported::refuse("sizing_standard=ros")),
        "achong" => return Err(unported::refuse("sizing_standard=achong")),
        other => {
            return Err(AzothError::invalid_input(
                "sizing_standard",
                format!("'{other}' is not one of the nine declared standards"),
            ));
        }
    }

    let geometry = ValveSizingGeometry {
        kv,
        opening_percent: opening,
        fouling_fraction: fouling_fraction.unwrap_or(0.0),
        characteristic: ValveCharacteristic::named(valve_characterization.unwrap_or("linear"))?,
        xt: xt.unwrap_or(DEFAULT_XT),
        allow_choked: allow_choked.unwrap_or(false),
    };

    let feed = Stream::from_pt(
        components.to_vec(),
        inlet_z.to_vec(),
        inlet_n,
        inlet_p,
        inlet_t,
    )?;
    // **The sizing runs on the inlet and the stated drop**, which is what the class's `calcKv`
    // reads; the pressure mode then replaces the drop with the one the trim answers.
    let report = sizing_report(&feed, outlet_pressure.value / 1.0e5, &geometry)?;
    let out_pressure = if is_calc_pressure.unwrap_or(false) {
        crate::valve::outlet_pressure(&report)
    } else {
        outlet_pressure
    };
    let out = if is_isothermal.unwrap_or(false) {
        // `isIsoThermal` flashes at the inlet temperature rather than isenthalpically.
        Stream::from_pt(
            components.to_vec(),
            inlet_z.to_vec(),
            inlet_n,
            out_pressure,
            inlet_t,
        )?
    } else {
        kernel(&feed, out_pressure)?
    };

    Ok(ThrottlingValveResult {
        outlet_n: out.n,
        outlet_z: out.z,
        outlet_p: out.p,
        outlet_t: out.t,
        outlet_h: joules_per_mole(out.h.value),
        kv: report.kv,
        cv: report.cv,
        cg: report.cg,
        choked: report.choked,
        gas_valve: report.gas_valve,
        calculated_molar_flow: report.calculated_molar_flow,
        warnings,
    })
}
