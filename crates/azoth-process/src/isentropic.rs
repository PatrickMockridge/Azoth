//! The procedure a compressor, a pump and an expander have in common.
//!
//! The three machines are the same eleven statements with one difference: which way the
//! efficiency scales the ideal enthalpy change, which is what [`Direction`] names.
//! Dividing puts the outlet further from the inlet than ideal - a machine consuming work;
//! multiplying puts it closer - a machine producing it.
//!
//! Each machine's spec carries its own provenance and its own list of what is not
//! ported: `specs/models/process/compressor.yaml`, `pump.yaml` and `expander.yaml`.

use azoth_core::units::{
    Pressure, ThermodynamicTemperature, joules_per_mole, joules_per_mole_kelvin,
};
use azoth_core::{AzothError, Result, Warning};
use azoth_eos::Phase;
use azoth_eos::mixture::Mixture;
use azoth_eos::molar_enthalpy_entropy::IdealGasModel;
use azoth_eos::ph_flash::{enthalpy_at, ph_flash};
use azoth_eos::ps_flash::{entropy_at, ps_flash};

/// Which way a machine moves the pressure, and so which way the efficiency scales.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    /// Work in. The outlet is at a higher pressure and a real machine needs *more* than
    /// the ideal enthalpy rise, so the efficiency **divides**.
    Consuming,
    /// Work out. The outlet is at a lower pressure and a real machine delivers *less*
    /// than the ideal drop, so the efficiency **multiplies**.
    Producing,
}

/// Everything the three machines compute, before it is put into their own result types.
#[derive(Debug, Clone)]
pub struct Solved {
    /// The real outlet temperature.
    pub temperature: ThermodynamicTemperature,
    /// The outlet pressure, which was an input.
    pub pressure: Pressure,
    /// The vapour fraction at the outlet.
    pub beta: Option<f64>,
    /// Which phase the stream is in at the outlet.
    pub phase: Phase,
    /// The shaft power in watts, signed: positive into the fluid.
    pub power: f64,
    /// The temperature at unit efficiency.
    pub isentropic_temperature: ThermodynamicTemperature,
    /// Iterations, summed over the flashes this ran.
    pub iterations: u32,
    /// Caveats from the flashes.
    pub warnings: Vec<Warning>,
}

/// Run one machine: a stated pressure change at a stated isentropic efficiency.
///
/// `efficiency` is the **isentropic** efficiency, in `(0, 1]`. The spec's range check
/// refuses a value outside it rather than clamping, because a clamp turns a caller's
/// mistake into a slightly different answer.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if the pressure moves the wrong way for
///   the direction claimed - a compressor whose outlet is not above its inlet, or an
///   expander whose outlet is not below it. The schema's range checks are per quantity
///   and cannot express a relation *between* two inputs, so this is the only place the
///   check can live.
/// * Whatever the flashes raise: a state with no root, or a temperature outside
///   `eos.ps_flash`'s or `eos.ph_flash`'s bracket.
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn run(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t_in: ThermodynamicTemperature,
    p_in: Pressure,
    p_out: Pressure,
    n: f64,
    z: &[f64],
    efficiency: f64,
    direction: Direction,
) -> Result<Solved> {
    let consistent = match direction {
        Direction::Consuming => p_out.value > p_in.value,
        Direction::Producing => p_out.value < p_in.value,
    };
    if !consistent {
        let (machine, sense) = match direction {
            Direction::Consuming => ("a compressor or a pump", "above"),
            Direction::Producing => ("an expander", "below"),
        };
        return Err(AzothError::InvalidInput {
            field: "outlet_pressure".to_string(),
            reason: format!(
                "{} must leave at a pressure {sense} its inlet; this one goes from {} Pa \
                 to {} Pa. A machine that moves the pressure the other way is a different \
                 unit operation, and admitting it here would return an answer for a \
                 process nobody can build",
                machine, p_in.value, p_out.value
            ),
        });
    }

    // The state the fluid arrives in. Both properties come from the same model at the
    // same state, so they describe one inlet and not two.
    let (h_in, _) = enthalpy_at(mixture, ideal_gas, t_in, p_in, z)?;
    let (s_in, _) = entropy_at(mixture, ideal_gas, t_in, p_in, z)?;

    // Where it would get to if the machine were perfect.
    let ideal = ps_flash(mixture, ideal_gas, p_out, joules_per_mole_kelvin(s_in), z)?;
    let (h_ideal, _) = enthalpy_at(mixture, ideal_gas, ideal.temperature, p_out, z)?;

    // The one expression that differs between the three machines.
    let d_h = match direction {
        Direction::Consuming => (h_ideal - h_in) / efficiency,
        Direction::Producing => (h_ideal - h_in) * efficiency,
    };

    let flash = ph_flash(mixture, ideal_gas, p_out, joules_per_mole(h_in + d_h), z)?;

    let mut warnings = ideal.warnings;
    warnings.extend(flash.warnings);

    Ok(Solved {
        temperature: flash.temperature,
        pressure: p_out,
        beta: flash.beta,
        phase: flash.phase,
        power: n * d_h,
        isentropic_temperature: ideal.temperature,
        iterations: ideal.iterations.saturating_add(flash.iterations),
        warnings,
    })
}
