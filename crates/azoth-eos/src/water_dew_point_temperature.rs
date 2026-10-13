//! `eos.water_dew_point_temperature` - the temperature a gas condenses water at.
//!
//! NeqSim's `WaterDewPointTemperatureFlash`, reached through
//! `ThermodynamicOperations.waterDewPointTemperatureFlash()`. Spec:
//! `specs/models/eos/water_dew_point_temperature.toml`. Oracle:
//! `validation/neqsim/captures/water_dew_point_probe.tsv`.
//!
//! # The second phase is *stated*, not found
//!
//! Every other saturation operation this crate ports asks a flash where a phase appears.
//! This one does not: it calls `setNumberOfPhases(2)`, sets phase 1 to `1e-30` on every
//! component but one `water`-or-`MEG` at `1.0`, and walks the temperature until that
//! phase's water fugacity matches the vapour's. So it is the one water-family operation
//! that does **not** rest on the aqueous seeds `crate::tp_multiflash` declines
//! (`seedHydrocarbonLiquidFromFeed`, `seedAdditionalPhaseFromFeed`); the phase is
//! prescribed, and nothing has to find it.
//!
//! # It is a fixed point with a fixed gain, and that is not a detail
//!
//! `funk = z_k - beta_0 * phi_1k / phi_0k`, then `T += 100.0 * funk` - **a hundred
//! kelvin per unit of residual**, not a Newton step. The root `funk = 0` is the
//! incipient-dew condition `z_k phi_vapour = 1 * phi_liquid` with the aqueous phase's
//! water at unit mole fraction, and the fixed point walks to it at a rate this gain
//! fixes. A port that solved the same condition by Newton would return the same
//! temperature and not the same path, so the seed - the component's own
//! `melting_point_temperature`, `273.15` K for water against a triple point of
//! `273.16` - and the gain are both carried.
//!
//! # The maximum, and the zero
//!
//! The class loops over *every* component named `water` or `MEG`, whatever its mole
//! fraction, and returns the **largest** temperature it reached. When a fluid has
//! neither, `maxTemperature` stays at its initialised `0.0` and the class writes a
//! zero temperature - which this port refuses rather than returns, because a zero
//! kelvin is not a state [`ThermodynamicTemperature`] represents. That refusal is
//! stated in the spec rather than patched over.

use azoth_core::units::{Pressure, kelvins};
use azoth_core::{AzothError, Result, apply_checks};

use crate::Cubic;
use crate::mixture::{Mixture, RootSide};
use crate::model_gen;
use crate::results::WaterDewPointTemperatureResult;

/// The residual the fixed point stops on, in mole fraction. NeqSim's own `1e-7`.
const TOLERANCE: f64 = 1.0e-7;

/// NeqSim's own `maxNumberOfIterations`.
const MAXIMUM_ITERATIONS: u32 = 10_000;

/// The temperature step per unit of residual, in kelvin: NeqSim's own `100.0 * funk`.
const GAIN: f64 = 100.0;

/// The mole fraction every non-target component of the incipient phase is given:
/// NeqSim's own `1e-30`.
const TRACE: f64 = 1.0e-30;

/// The two components the class seeds an incipient phase from.
///
/// NeqSim matches them with `equals("water") || equals("MEG")`, and its own table spells
/// the second one in capitals. This library's spells it `meg`, so the match is made
/// case-insensitively: the *substances* are the same two, and a caller cannot express
/// NeqSim's spelling of the glycol here.
const TARGETS: [&str; 2] = ["water", "meg"];

/// One component's fixed point, and what it cost.
struct Walk {
    temperature: f64,
    residual: f64,
    iterations: u32,
    converged: bool,
}

/// The incipient phase's composition: every component at [`TRACE`], the target at one.
fn incipient(index: usize, count: usize) -> Vec<f64> {
    let mut x = vec![TRACE; count];
    x[index] = 1.0;
    x
}

/// The fixed point on one target component.
///
/// `feed` is the overall composition, which is phase 0's - the class never changes it.
fn walk(mixture: &Mixture, p: Pressure, feed: &[f64], index: usize, seed: f64) -> Result<Walk> {
    let x = incipient(index, mixture.len());
    let mut temperature = seed;
    let mut residual = 0.0;
    let mut iterations = 0;
    let mut converged = false;
    while iterations < MAXIMUM_ITERATIONS {
        iterations += 1;
        let reduced = mixture.reduced_parameters(kelvins(temperature), p)?;
        let vapour = mixture.phase_state(&reduced, feed, RootSide::Vapour)?;
        let aqueous = mixture.phase_state(&reduced, &x, RootSide::Liquid)?;
        // `beta_0` is phase 0's own fraction, which `setNumberOfPhases(2)` leaves at one:
        // the second phase is incipient, so its amount is zero and the first phase is the
        // whole of the fluid. The ratio is taken in logs because the two coefficients span
        // orders of magnitude at a water mole fraction of `3.5e-5`.
        residual = feed[index] - (aqueous.ln_phi[index] - vapour.ln_phi[index]).exp();
        temperature += GAIN * residual;
        if residual.abs() < TOLERANCE {
            converged = true;
            break;
        }
    }
    Ok(Walk {
        temperature,
        residual,
        iterations,
        converged,
    })
}

/// The temperature a gas condenses water at, by NeqSim's `WaterDewPointTemperatureFlash`.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `components` and `z` are not the same length, if `z`
///   does not name a `water` or `MEG` component at a positive mole fraction, or if that
///   component's `melting_point_temperature` is absent - which is a card's substance.
pub fn water_dew_point_temperature(
    components: &[String],
    z: &[f64],
    p: Pressure,
    cubic: Cubic,
) -> Result<WaterDewPointTemperatureResult> {
    let spec = &model_gen::WATER_DEW_POINT_TEMPERATURE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "P" => Some(p.value),
            _ => None,
        },
        &mut warnings,
    )?;

    if components.len() != z.len() {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "{} components and {} mole fractions",
                components.len(),
                z.len()
            ),
        ));
    }

    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, _) = crate::databank::mixture_of(&names, cubic, None)?;

    let mut best: Option<Walk> = None;
    for (index, name) in components.iter().enumerate() {
        let lowered = name.trim().to_lowercase();
        if !TARGETS.contains(&lowered.as_str()) {
            continue;
        }
        // The class seeds from every `water` or `MEG` component whatever its mole
        // fraction, so a trace one walks too; a zero would make the residual `-ratio` and
        // is refused here rather than walked, which the spec states.
        if !z[index].is_finite() || z[index] <= 0.0 {
            return Err(AzothError::invalid_input(
                "z",
                format!(
                    "`{name}` is the phase this model seeds from and its mole fraction is \
                     {}, so there is no dew point to find - NeqSim walks it anyway and \
                     answers a number that is not one",
                    z[index]
                ),
            ));
        }
        let seed = mixture.components()[index].melting_point_temperature;
        if !seed.is_finite() || seed <= 0.0 {
            return Err(AzothError::invalid_input(
                "components",
                format!(
                    "`{name}` carries no melting-point temperature, which is the \
                     temperature the fixed point starts from; a card's substance has none, \
                     and seeding at zero kelvin is not a state this library represents"
                ),
            ));
        }
        let result = walk(&mixture, p, z, index, seed)?;
        if best
            .as_ref()
            .is_none_or(|held| result.temperature > held.temperature)
        {
            best = Some(result);
        }
    }

    let walk = best.ok_or_else(|| {
        AzothError::invalid_input(
            "components",
            "the fluid names neither `water` nor `MEG`, which are the two components this \
             model seeds an incipient phase from. NeqSim writes a zero temperature there; \
             this refuses."
                .to_string(),
        )
    })?;

    Ok(WaterDewPointTemperatureResult {
        temperature: kelvins(walk.temperature),
        residual: walk.residual,
        iterations: walk.iterations,
        converged: walk.converged,
        warnings,
    })
}
