//! `process.mixer` - several feeds blended into one.
//!
//! Spec: `specs/models/process/mixer.yaml`
//!
//! # The port
//!
//! NeqSim's `Mixer.run` is lines 642-729 of a 1,235-line `Mixer.java`, and its physics
//! is four statements:
//!
//! ```text
//! n_out = sum_s n_s                              (mixStream, :192-320)
//! z_out = sum_s (n_s * z_s) / n_out              (:292, addComponent accumulating moles)
//! P_out = min_s(P_s)                             (:681)
//! T_out = PHflash(P_out, sum_s H_s)              (:702-711)
//! ```
//!
//! The rest of the file is the component-matching machinery - NeqSim accumulates moles
//! component by component and has to handle an inlet whose pseudo-fraction has a
//! different molar mass from the destination's (`:285-291`), and an inlet carrying a
//! component the others do not. None of that exists here: azoth takes one `Mixture` for
//! every stream, so all the inlets already share a component set and the accumulation is
//! a vector sum.
//!
//! # One expression that had to change, and why
//!
//! NeqSim's `calcMixStreamEnthalpy` (`:537-544`) sums **total** enthalpies in joules,
//! and its `PHflash` takes joules, so the sum is the whole story. `eos.ph_flash` takes a
//! **molar** enthalpy, so the translation is the flow-weighted *mean*:
//! `H = sum(n_s H_s) / sum(n_s)`. Same physics - an energy balance on an open system
//! with no heat and no work - and a different expression.
//!
//! **This is the one place in the port where a formula looks different from the source
//! it came from**, which is why it is written down rather than left to look like a
//! transcription error. It is also worth knowing that the two are equal when the flows
//! are: the mean and the sum differ only by the factor `sum(n_s)` that divides out
//! again when the flash is given a molar enthalpy.
//!
//! # The outlet pressure is the lowest inlet pressure
//!
//! Not the highest and not an average. A mixer is a vessel: nothing in it can be above
//! the lowest pressure any feed arrives at, and a feed arriving at 20 bar into a vessel
//! held at 5 bar flashes as it enters. Taking the maximum, or the arithmetic mean, would
//! give an answer that is arithmetically fine and describes a pump.

use azoth_core::units::{Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{AzothError, Result, apply_checks};
use azoth_eos::mixture::Mixture;
use azoth_eos::molar_enthalpy_entropy::IdealGasModel;
use azoth_eos::ph_flash::{enthalpy_at, ph_flash};

use crate::model_gen;
use crate::results::MixerResult;

/// Several feeds blended into one, at the lowest inlet pressure.
///
/// `t_in` and `p_in` are the inlets' states and `n_in` their molar flows in mol/s. `z`
/// is every inlet's composition, **flattened row-major** - inlet `s` occupies
/// `z[s * N .. (s + 1) * N]` with `N = mixture.len()`. The flattened shape is the same
/// one `kij` crosses in, for the same reason: the boundary carries one list rather than
/// a list of lists, and the component order is the one thing the two sides have to agree
/// about.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if the inlet vectors disagree about how many
///   streams there are, if `z` is not `S * N` long, or if any inlet's composition is not
///   a composition.
/// * [`azoth_core::AzothError::OutOfRange`] if an inlet state is non-positive.
/// * Whatever the flash raises.
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn mixer(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t_in: &[ThermodynamicTemperature],
    p_in: &[Pressure],
    n_in: &[f64],
    z: &[f64],
) -> Result<MixerResult> {
    let spec = &model_gen::MIXER_SPEC;
    let mut warnings = Vec::new();

    let streams = t_in.len();
    let components = mixture.len();
    if p_in.len() != streams || n_in.len() != streams {
        return Err(AzothError::InvalidInput {
            field: "P".to_string(),
            reason: format!(
                "a mixer's inlets must be the same number of streams; got {} temperature(s), \
                 {} pressure(s) and {} flow(s)",
                streams,
                p_in.len(),
                n_in.len()
            ),
        });
    }
    if z.len() != streams * components {
        return Err(AzothError::InvalidInput {
            field: "z".to_string(),
            reason: format!(
                "`z` must carry one composition per inlet - {streams} x {components} = {} \
                 number(s) - and carries {}",
                streams * components,
                z.len()
            ),
        });
    }

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "n" => Some(n_in.iter().copied().fold(f64::INFINITY, f64::min)),
            "P" => Some(p_in.iter().map(|p| p.value).fold(f64::INFINITY, f64::min)),
            "T" => Some(t_in.iter().map(|t| t.value).fold(f64::INFINITY, f64::min)),
            _ => None,
        },
        &mut warnings,
    )?;

    let total_flow: f64 = n_in.iter().sum();
    let outlet_pressure = p_in.iter().map(|p| p.value).fold(f64::INFINITY, f64::min);

    // The blend, and the enthalpy balance, in one pass over the inlets.
    let mut z_out = vec![0.0; components];
    let mut total_enthalpy = 0.0;
    for stream in 0..streams {
        let composition = &z[stream * components..(stream + 1) * components];
        let (h, _) = enthalpy_at(mixture, ideal_gas, t_in[stream], p_in[stream], composition)?;
        total_enthalpy += n_in[stream] * h;
        for (out, inlet) in z_out.iter_mut().zip(composition) {
            *out += n_in[stream] * inlet;
        }
    }
    let molar_enthalpy = total_enthalpy / total_flow;
    for fraction in &mut z_out {
        *fraction /= total_flow;
    }

    let flash = ph_flash(
        mixture,
        ideal_gas,
        azoth_core::units::pascals(outlet_pressure),
        joules_per_mole(molar_enthalpy),
        &z_out,
    )?;
    warnings.extend(flash.warnings);

    Ok(MixerResult {
        temperature: flash.temperature,
        pressure: azoth_core::units::pascals(outlet_pressure),
        flow: total_flow,
        z_out,
        beta: flash.beta,
        phase: flash.phase,
        iterations: flash.iterations,
        warnings,
    })
}
