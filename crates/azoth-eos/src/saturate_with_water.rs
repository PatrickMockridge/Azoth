//! `eos.saturate_with_water`: the water a fluid holds when it is water-saturated.
//!
//! NeqSim's `SaturateWithWater` floods the feed with water, flashes, and reads the water
//! dissolved in the non-aqueous phases - a fluid coexisting with a free aqueous phase is, by
//! definition, saturated - then resets the overall water to that amount and re-flashes.
//!
//! **The class mutates a system; this is a function.** The three paths that restore the caller's
//! composition instead of returning a flooded fluid - the flood flash diverged, no aqueous phase
//! formed, a residual incipient aqueous phase survived the correction - become "return the feed",
//! with the reason in [`SaturateWithWaterResult::warnings`].

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result, Warning, WarningCode, apply_checks};

use crate::cubic::Cubic;
use crate::databank::mixture_of;
use crate::hydrate_inhibitor_wt::{PhaseLabel, label};
use crate::mixture::Mixture;
use crate::model_gen;
use crate::results::{SaturateWithWaterResult, TpMultiflashResult};
use crate::tp_multiflash::tp_multiflash;

/// The water a fluid holds at saturation, with the phase state at that composition.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `moles` does not match `components`, `eos` is not a cubic
///   this library carries, or the fluid cannot be built.
/// * Whatever the final flash raises. The *flood* flash's error is recovered from, as the
///   class's `runFlashSafely` does.
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn saturate_with_water(
    components: &[String],
    moles: &[f64],
    T: ThermodynamicTemperature,
    P: Pressure,
    eos: &str,
    associating: bool,
) -> Result<SaturateWithWaterResult> {
    let spec = &model_gen::SATURATE_WITH_WATER_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T.value),
            "P" => Some(P.value),
            _ => None,
        },
        &mut warnings,
    )?;

    if moles.len() != components.len() {
        return Err(AzothError::invalid_input(
            "moles",
            format!(
                "{} component(s) and {} mole number(s)",
                components.len(),
                moles.len()
            ),
        ));
    }
    // **The associating fluid is refused, and the reason is measured.** The operation is a
    // *multiphase* flash, and `eos.tp_multiflash` has no associating case anywhere in this
    // tree - so an associating fluid's answer here would rest on the one combination no case
    // covers. Measured against the source's own two CPA states, it is 19% high on one and
    // 15% low on the other, in opposite directions: a plausible number and a wrong one.
    if associating {
        return Err(AzothError::invalid_input(
            "associating",
            "this operation floods and flashes multiphase, and the multiphase flash has no \
             associating case in this tree; the two source states it was measured on disagree \
             by 15-19% in opposite directions, so an associating fluid is refused rather than \
             answered"
                .to_string(),
        ));
    }
    let cubic: Cubic = eos.parse().map_err(|_| {
        AzothError::invalid_input(
            "eos",
            format!("`{eos}` is not a cubic this library carries"),
        )
    })?;

    let mut names: Vec<String> = components.to_vec();
    let mut feed: Vec<f64> = moles.to_vec();

    // **The trace of water, and it is not an optimisation to skip.** `addComponent("water",
    // totalMoles / 100)` runs before the flood, so a feed without water arrives at the flood
    // with 1% of a feed more moles in it - and the saturation amount is reported against that
    // total. The measured check values are the ones this reproduces.
    let total: f64 = feed.iter().sum();
    if total <= 0.0 || !total.is_finite() {
        return Err(AzothError::out_of_range(
            "moles",
            total,
            "the feed's mole numbers do not describe a feed",
        ));
    }
    if !names.iter().any(|name| name.eq_ignore_ascii_case("water")) {
        names.push("water".to_string());
        feed.push(total / 100.0);
    }
    let total: f64 = feed.iter().sum();
    let input_moles = feed.clone();
    let input_total = total;

    let refs: Vec<&str> = names.iter().map(String::as_str).collect();
    // A classical fluid only: `associating` is refused above, and the branch that built a CPA
    // mixture comes back when the multiphase flash has an associating case behind it.
    let (mixture, _) = mixture_of(&refs, cubic, None)?;
    let water = mixture.index_of("water").ok_or_else(|| {
        AzothError::invalid_input(
            "components",
            "no water: this reports the water a fluid holds at saturation, so a feed without one \
             has nothing to saturate"
                .to_string(),
        )
    })?;

    // The flood: one feed's worth of water.
    let flooded_total = total * 2.0;
    let mut flooded = feed.clone();
    flooded[water] += total;
    let flooded_fractions: Vec<f64> = flooded.iter().map(|value| value / flooded_total).collect();

    let flood = match tp_multiflash(&mixture, T, P, &flooded_fractions) {
        Ok(flash) => flash,
        Err(error) => {
            warnings.push(Warning::new(
                WarningCode::FallbackToInput,
                format!("the flood flash failed ({error}); the feed is returned unchanged"),
            ));
            return at_the_feed(&mixture, &input_moles, input_total, water, T, P, warnings);
        }
    };
    let flood_labels = labels_of(&mixture, &flood, T, P)?;

    if !flood_labels.contains(&PhaseLabel::Aqueous) {
        warnings.push(Warning::new(
            WarningCode::FallbackToInput,
            "flooding produced no free aqueous phase, so saturation is undefined here and the \
             feed is returned unchanged",
        ));
        return at_the_feed(&mixture, &input_moles, input_total, water, T, P, warnings);
    }
    let aqueous = flood_labels
        .iter()
        .position(|phase| *phase == PhaseLabel::Aqueous)
        .expect("an aqueous phase was just found");

    // The water dissolved in every non-aqueous phase, which is the saturation amount.
    let saturated: f64 = (0..flood.phase_count as usize)
        .filter(|phase| *phase != aqueous)
        .map(|phase| flood.phase_fractions[phase] * flood.x[phase][water] * flooded_total)
        .sum();

    let mut saturated_moles = input_moles.clone();
    saturated_moles[water] = saturated;
    let saturated_total: f64 = saturated_moles.iter().sum();
    if saturated_total <= 0.0 || !saturated_total.is_finite() {
        return Err(AzothError::out_of_range(
            "moles",
            saturated_total,
            "the saturated composition has no moles, which the correction produced",
        ));
    }
    let saturated_fractions: Vec<f64> = saturated_moles
        .iter()
        .map(|value| value / saturated_total)
        .collect();

    let saturated_flash = tp_multiflash(&mixture, T, P, &saturated_fractions)?;
    let saturated_labels = labels_of(&mixture, &saturated_flash, T, P)?;

    // **A residual incipient aqueous phase is dropped, not flashed away.** The class calls
    // `removePhase` and re-flashes; a flash here is a function of the composition alone, so the
    // re-flash is the same answer and the removal is the part that shows.
    let kept: Vec<usize> = (0..saturated_flash.phase_count as usize)
        .filter(|phase| {
            !(saturated_labels[*phase] == PhaseLabel::Aqueous && saturated_flash.phase_count > 1)
        })
        .collect();
    let kept_fraction: f64 = kept
        .iter()
        .map(|phase| saturated_flash.phase_fractions[*phase])
        .sum();

    Ok(SaturateWithWaterResult {
        water_mole_fraction: saturated / saturated_total,
        phase_count: kept.len() as u32,
        phase_fractions: kept
            .iter()
            .map(|phase| saturated_flash.phase_fractions[*phase] / kept_fraction)
            .collect(),
        x: kept
            .iter()
            .map(|phase| saturated_flash.x[*phase].clone())
            .collect(),
        z_factor: kept
            .iter()
            .map(|phase| saturated_flash.z_factor[*phase])
            .collect(),
        aqueous_present: true,
        warnings,
    })
}

/// The feed's own state, which is what each recovery path returns.
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
fn at_the_feed(
    mixture: &Mixture,
    moles: &[f64],
    total: f64,
    water: usize,
    T: ThermodynamicTemperature,
    P: Pressure,
    warnings: Vec<Warning>,
) -> Result<SaturateWithWaterResult> {
    let fractions: Vec<f64> = moles.iter().map(|value| value / total).collect();
    let flash = tp_multiflash(mixture, T, P, &fractions)?;
    Ok(SaturateWithWaterResult {
        water_mole_fraction: moles[water] / total,
        phase_count: flash.phase_count,
        phase_fractions: flash.phase_fractions.clone(),
        x: flash.x.clone(),
        z_factor: flash.z_factor.clone(),
        aqueous_present: false,
        warnings,
    })
}

/// The label of every phase in a flash, in the flash's own order.
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
fn labels_of(
    mixture: &Mixture,
    flash: &TpMultiflashResult,
    T: ThermodynamicTemperature,
    P: Pressure,
) -> Result<Vec<PhaseLabel>> {
    let reduced = mixture.reduced_parameters(T, P)?;
    (0..flash.phase_count as usize)
        .map(|phase| label(mixture, &reduced, &flash.x[phase], flash.z_factor[phase]))
        .collect()
}
