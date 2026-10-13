//! `eos.water_dew_point_multiphase` - the water dew point a multiphase flash brackets.
//!
//! NeqSim's `WaterDewPointTemperatureMultiphaseFlash`, reached through
//! `ThermodynamicOperations.waterDewPointTemperatureMultiphaseFlash()`. Spec:
//! `specs/models/eos/water_dew_point_multiphase.toml`. Oracle:
//! `validation/neqsim/captures/water_dew_point_probe.tsv`.
//!
//! # The other route, and it took a seeding to become possible
//!
//! [`crate::water_dew_point_temperature`] **prescribes** its incipient phase and walks the
//! temperature to the fugacity equality. This one does the opposite: it asks a *multiphase
//! flash* where an aqueous phase appears and brackets on the answer, so it is only as good as
//! the flash's seeding underneath. When this class was first read, this tree's flash answered
//! a single gas phase at the very temperature NeqSim's dew point names - because
//! [`crate::tp_multiflash`] was seeding from [`crate::stability_test`], the *base* method's
//! two-trial test, where `TPmultiflash` overrides it with a pure-component fallback.
//! `eos.multiflash_stability_test` is that override, and with it the two classes agree.
//!
//! # The answer is a bracket's midpoint
//!
//! Upstream starts at `500` K and steps **down** by `20` K until an aqueous phase appears,
//! expanding upward instead if `500` K already has one; then bisects to `1e-4` K; then asks
//! for `0.5 * (Thigh + Tlow)` once more. So the answer is a dyadic midpoint on a grid the
//! bracket's own history fixes, and the last digit of it is a fact about that history rather
//! than about the fluid - which is why `bracket_width` is reported beside it.
//!
//! **The two caps are upstream's own and they are not the same counter.** The coarse search
//! runs while `i < 50` and the bisection while `i < 100`, over a counter neither resets - so
//! a feed that spends many steps in the first has fewer left for the second, and a port that
//! counted them separately would bisect to a finer bracket and answer a different number.

use azoth_core::units::{Pressure, kelvins};
use azoth_core::{AzothError, Result, apply_checks};

use crate::Cubic;
use crate::hydrate_inhibitor_wt::{PhaseLabel, label};
use crate::model_gen;
use crate::results::WaterDewPointMultiphaseResult;
use crate::tp_multiflash::tp_multiflash;

/// The temperature the coarse search starts from, upstream's `Thigh`.
const START_TEMPERATURE: f64 = 500.0;

/// The floor the coarse search steps down to, upstream's `Tlow`.
const FLOOR_TEMPERATURE: f64 = 200.0;

/// The temperature the search expands to when `500` K already has an aqueous phase.
const CEILING_TEMPERATURE: f64 = 600.0;

/// The coarse step, upstream's `step`.
const STEP: f64 = 20.0;

/// The bracket width the bisection stops at, upstream's `tolerance`.
const TOLERANCE: f64 = 1.0e-4;

/// The coarse search's step cap, upstream's `i < 50`.
const COARSE_CAP: u32 = 50;

/// The bisection's step cap, upstream's `i < 100` - **the same counter** as the coarse cap.
const TOTAL_CAP: u32 = 100;

/// A trial state, and whether it has an aqueous phase.
///
/// The label is NeqSim's own `PhaseEos.init` rule, the one [`crate::hydrate_inhibitor_wt`]
/// already carries: a volume ratio past `1.75` is a gas, and past that it is aqueous exactly
/// where its hydrocarbon content does not outweigh the rest.
struct Trial {
    aqueous: bool,
    phases: u32,
}

fn trial(mixture: &crate::mixture::Mixture, t: f64, p: Pressure, z: &[f64]) -> Result<Trial> {
    let result = tp_multiflash(mixture, kelvins(t), p, z)?;
    let reduced = mixture.reduced_parameters(kelvins(t), p)?;
    let mut aqueous = false;
    for (composition, z_factor) in result.x.iter().zip(&result.z_factor) {
        if label(mixture, &reduced, composition, *z_factor)? == PhaseLabel::Aqueous {
            aqueous = true;
        }
    }
    Ok(Trial {
        aqueous,
        phases: result.phase_count,
    })
}

/// The temperature an aqueous phase first appears at, by NeqSim's own bracket.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `components` and `z` are not the same length, or if `z`
///   does not name a `water` component at a positive mole fraction - a fluid with no water has
///   no aqueous phase to find and the bracket would step to its floor and answer that.
/// * Propagates the flash's and the kernels' range checks.
pub fn water_dew_point_multiphase(
    components: &[String],
    z: &[f64],
    p: Pressure,
    cubic: Cubic,
) -> Result<WaterDewPointMultiphaseResult> {
    let spec = &model_gen::WATER_DEW_POINT_MULTIPHASE_SPEC;
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
    let water = components
        .iter()
        .position(|name| name.trim().eq_ignore_ascii_case("water"))
        .ok_or_else(|| {
            AzothError::invalid_input(
                "components",
                "the fluid names no `water`, and a water dew point is the temperature an \
                 aqueous phase appears at - without water there is none to find"
                    .to_string(),
            )
        })?;
    if z[water] <= 0.0 {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "`water` is present and its mole fraction is {}, so there is no aqueous phase \
                 for the bracket to find",
                z[water]
            ),
        ));
    }

    let mut iterations = 0;
    let mut high = START_TEMPERATURE;
    let mut low = FLOOR_TEMPERATURE;

    let at_start = trial(&mixture, high, p, z)?;
    iterations += 1;
    if at_start.aqueous {
        low = high;
        high = CEILING_TEMPERATURE;
        trial(&mixture, high, p, z)?;
        iterations += 1;
    } else {
        let mut t = high - STEP;
        let mut found = false;
        while t >= low && iterations < COARSE_CAP {
            let step = trial(&mixture, t, p, z)?;
            iterations += 1;
            if step.aqueous {
                low = t;
                high = t + STEP;
                found = true;
                break;
            }
            t -= STEP;
        }
        if !found {
            // Upstream answers the floor here. This refuses instead: a dew point that is not
            // the temperature of an appearance is not one, and the caller cannot tell the two
            // apart from the number alone.
            return Err(AzothError::invalid_input(
                "T",
                format!(
                    "no aqueous phase appears anywhere in {FLOOR_TEMPERATURE} to \
                     {START_TEMPERATURE} K at this pressure, so there is no water dew point to \
                     bracket. Upstream returns the floor, which reads as an answer"
                ),
            ));
        }
    }

    while (high - low) > TOLERANCE && iterations < TOTAL_CAP {
        let mid = 0.5 * (high + low);
        if trial(&mixture, mid, p, z)?.aqueous {
            low = mid;
        } else {
            high = mid;
        }
        iterations += 1;
    }

    let answer = 0.5 * (high + low);
    let width = high - low;
    let at = trial(&mixture, answer, p, z)?;

    Ok(WaterDewPointMultiphaseResult {
        temperature: kelvins(answer),
        bracket_width: kelvins(width).value,
        phase_count: at.phases,
        iterations,
        warnings,
    })
}
