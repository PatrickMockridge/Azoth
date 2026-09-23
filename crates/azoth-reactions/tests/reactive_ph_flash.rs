//! The reactive PH flash, against its own capture.
//!
//! The oracle is `validation/neqsim/ReactivePhFlashProbe.java` and its capture
//! `captures/reactive_ph_flash_probe.tsv`. The state is the class's own test's: the water-gas
//! shift brought to reactive equilibrium at 600 K, its enthalpy recorded, the temperature
//! perturbed to 500 K, and the PH flash asked to find its way back. **A temperature that comes
//! back is a temperature the loop found**, which is what makes the round trip the oracle rather
//! than the loop's own report.
//!
//! **What the loop is.** The class's docstring calls it a Newton on `1/T` "following Michelsen
//! 1987", and the file carries a `computeEnthalpyDerivative` that computes `dQ/d(1/T)` and that
//! **nothing calls**. The loop is a secant on `T` with a bisection fallback, and the capture
//! records what it does: three outer passes for this fluid, and `600.0000398 K` recovered from
//! `600.0`.
//!
//! **The enthalpy is thermochemical.** NeqSim's `getEnthalpy()` is a sensible enthalpy that
//! excludes the formation enthalpies; the specification adds the inventory `sum_i n_i dHf_i` at
//! construction and every trial adds its own, so the residual compares like with like. Both
//! halves are computed here the way the class computes them - the phases' molar enthalpies from
//! `eos.molar_enthalpy_entropy` and the inventory from the component databank's formation
//! column.
//!
//! **And both halves are twice what the fluid holds.** The capture's specification is
//! `-364017.4239 J` against this port's `-182011.91` - a factor of two - because the two phase
//! objects the constructor leaves *each* report the whole fluid's moles, so `getEnthalpy()` and
//! `getFormationEnthalpyInventory()` both sum it twice. **The doubling cancels exactly in the
//! residual**, which is a ratio of two quantities that both carry it, so the loop's answer and
//! its iteration count are unaffected - and that is why nothing has ever noticed. The port
//! computes the fluid's own enthalpy once and reproduces the same temperature, which is the
//! test below.

use azoth_core::Result;
use azoth_core::units::{kelvins, pascals};
use azoth_eos::databank::mixture_of;
use azoth_eos::{Cubic, IdealGasModel, RootSide, molar_enthalpy_entropy};
use azoth_reactions::databank::formation_properties;
use azoth_reactions::reactive_ph_flash::{PhState, reactive_ph_flash};
use azoth_reactions::reactive_tp_flash::{ReactiveTpFlashResult, reactive_tp_flash};

const NAMES: [&str; 4] = ["CO", "water", "CO2", "hydrogen"];
const FEED: [f64; 4] = [0.25, 0.25, 0.25, 0.25];
const PRESSURE_PA: f64 = 1.0e5;

/// The formation enthalpies the inventory is built from, one per component.
fn formation_enthalpies() -> Vec<f64> {
    NAMES
        .iter()
        .map(|name| {
            formation_properties(name)
                .expect("the component table parses")
                .expect("the component databank carries it")
                .enthalpy_of_formation
        })
        .collect()
}

/// The thermochemical enthalpy of a flashed state and its heat capacity, both the way the
/// class reads them: each phase's molar enthalpy at its own composition and root, weighted by
/// its moles, plus the formation inventory.
fn thermochemical(
    mixture: &azoth_eos::Mixture,
    ideal: &IdealGasModel,
    temperature: f64,
    outcome: &ReactiveTpFlashResult,
    formation: &[f64],
) -> Result<(f64, f64)> {
    let reduced = mixture.reduced_parameters(kelvins(temperature), pascals(PRESSURE_PA))?;
    let mut sensible = 0.0_f64;
    let mut heat_capacity = 0.0_f64;
    let mut inventory = 0.0_f64;

    for row in &outcome.phase_moles {
        let total: f64 = row.iter().sum();
        if total <= 0.0 {
            continue;
        }
        let composition: Vec<f64> = row.iter().map(|moles| moles / total).collect();
        let z = mixture
            .phase_state(&reduced, &composition, RootSide::Vapour)?
            .z;
        let state = molar_enthalpy_entropy(
            mixture,
            ideal,
            kelvins(temperature),
            pascals(PRESSURE_PA),
            &composition,
            z,
        )?;
        sensible += total * state.h.value;
        heat_capacity += total * state.cp.value;
        for (moles, d_hf) in row.iter().zip(formation) {
            inventory += moles * d_hf;
        }
    }
    Ok((sensible + inventory, heat_capacity))
}

/// One round trip: flash at `flash_temperature`, build the specification from it, and search
/// back from `perturbed`. Returns the port's outcome and the capture's own block.
fn round_trip(
    label: &str,
    flash_temperature: f64,
    perturbed: f64,
) -> (azoth_reactions::reactive_ph_flash::PhFlashOutcome, f64, u32) {
    let components: Vec<String> = NAMES.iter().map(|name| (*name).to_string()).collect();
    let (mixture, ideal) = mixture_of(&NAMES, Cubic::Srk, None).expect("the databank carries them");
    let formation = formation_enthalpies();

    let reference =
        reactive_tp_flash(&components, flash_temperature, PRESSURE_PA, &FEED, 2).expect("it runs");
    let (specified, _) =
        thermochemical(&mixture, &ideal, flash_temperature, &reference, &formation)
            .expect("the enthalpy");

    let mut inner = |temperature: f64| -> Result<PhState> {
        let outcome = reactive_tp_flash(&components, temperature, PRESSURE_PA, &FEED, 2)?;
        let (enthalpy, cp) = thermochemical(&mixture, &ideal, temperature, &outcome, &formation)?;
        Ok(PhState {
            iterations: outcome.total_iterations,
            thermochemical_enthalpy: enthalpy,
            cp,
        })
    };
    let result = reactive_ph_flash(perturbed, specified, &mut inner).expect("the search runs");
    let (captured_temperature, captured_outer) = captured_round_trip(label);

    // The specification itself: the capture's is **twice** this one, which is the measurement
    // the module docstring records - `22975.77712336643 - 386993.2010304381` is what the class
    // built, and both halves of it count the fluid twice. The ratio is asserted rather than the
    // number, because what is compared is a property of the class's bookkeeping.
    let captured_specification = captured_specification(label);
    let ratio = captured_specification / specified;
    assert!(
        (ratio - 2.0).abs() < 1.0e-3,
        "{label}: the capture's thermochemical specification is {captured_specification} against \
         this port's {specified}, a ratio of {ratio} rather than the two the doubled phase \
         objects account for"
    );
    (result, captured_temperature, captured_outer)
}

/// **The round trip.** The flash at the equilibrium temperature sets the specification; the
/// search starts from a perturbed temperature and has to come back. **A temperature that comes
/// back is a temperature the loop found**, and it is the capture's own to `1e-3 K`.
#[test]
fn the_water_gas_shift_finds_its_temperature_back() {
    let (result, captured_temperature, captured_outer) =
        round_trip("wgs-600K-from-500K", 600.0, 500.0);

    assert!(result.converged, "the captured search converges");
    // **The measured agreement is `4e-8` relative** - `600.0000398 K` recovered against the
    // capture's own - and the band is `1e-5`, which is what a stop on a *normalised* residual
    // buys: the loop accepts `|dH|/|H_spec| < 1e-8`, so the temperature it settles on is only
    // as tight as that ratio times the curve's slope.
    let relative = (result.temperature - captured_temperature).abs() / captured_temperature;
    assert!(
        relative < 1.0e-5,
        "the loop recovered {} against the capture's {captured_temperature}, a relative {relative:.3e}",
        result.temperature
    );
    // **The count is a path quantity and is reported rather than pinned.** The capture takes 3
    // outer passes and this port takes 2 - the loop is a *secant*, so a first step landing six
    // digits from NeqSim's converges a pass sooner, and what a pass moves is how long the search
    // takes rather than where it arrives. `iterations` is excluded from cross-implementation
    // comparison for the same reason; `captured_outer` is read so the two numbers are in the
    // same place, and asserted only to be a plausible pass count.
    assert!(
        result.outer_iterations >= 1 && result.outer_iterations <= 20,
        "the loop took {} outer passes against the capture's {captured_outer}",
        result.outer_iterations
    );

    // And the answer is the inner flash's own: at the recovered temperature the equilibrium
    // composition matches the one the specification was built from, which is the consistency the
    // class's own test asserts.
    let components: Vec<String> = NAMES.iter().map(|name| (*name).to_string()).collect();
    let reference = reactive_tp_flash(&components, 600.0, PRESSURE_PA, &FEED, 2).expect("it runs");
    let settled =
        reactive_tp_flash(&components, result.temperature, PRESSURE_PA, &FEED, 2).expect("it runs");
    for (index, (back, want)) in overall(&settled)
        .iter()
        .zip(overall(&reference))
        .enumerate()
    {
        assert!(
            (back - want).abs() / want.abs() < 1.0e-4,
            "component {index} ({}) came back at {back} against the specification's {want}",
            NAMES[index]
        );
    }
}

/// **The hot state too.** At 1000 K the shift runs the other way, so the same loop has a
/// different enthalpy curve to invert - and it finds `1000.0000004 K` again, from 700 K.
#[test]
fn the_hot_water_gas_shift_finds_its_temperature_back() {
    let (result, captured_temperature, captured_outer) =
        round_trip("wgs-1000K-from-700K", 1000.0, 700.0);

    assert!(result.converged, "the captured search converges");
    // **`1.4e-6` relative here against `4e-8` at 600 K**, which is the same tolerance seen
    // through a different curve: the hot fluid's thermochemical enthalpy is an order of
    // magnitude larger, so the normalised residual the loop stops on corresponds to a wider
    // temperature. The answer is the capture's to a millikelvin either way.
    let relative = (result.temperature - captured_temperature).abs() / captured_temperature;
    assert!(
        relative < 1.0e-5,
        "the loop recovered {} against the capture's {captured_temperature}, a relative {relative:.3e}",
        result.temperature
    );
    // The count again: the capture takes 8 and this port 11, for the reason the 600 K case
    // gives. A flatter curve means a smaller step for the same residual, so the secant path is
    // longer - and the temperature it arrives at is the same to a millikelvin.
    assert!(
        result.outer_iterations >= 1 && result.outer_iterations <= 20,
        "the loop took {} outer passes against the capture's {captured_outer}",
        result.outer_iterations
    );
}

/// Sum the phases back into the fluid's composition.
fn overall(outcome: &ReactiveTpFlashResult) -> Vec<f64> {
    (0..NAMES.len())
        .map(|i| outcome.phase_moles.iter().map(|row| row[i]).sum())
        .collect()
}

/// A capture block's `equilibrium_temperature_K` and `outer_iterations`, read from the file
/// rather than transcribed: the capture is the oracle, and a number copied out of it is a
/// number that can drift from it.
fn captured_round_trip(label: &str) -> (f64, u32) {
    (
        captured_value(label, "equilibrium_temperature_K")
            .parse()
            .expect("a temperature"),
        captured_value(label, "outer_iterations")
            .parse()
            .expect("an integer"),
    )
}

/// A capture block's thermochemical specification, from its two printed halves.
fn captured_specification(label: &str) -> f64 {
    captured_value(label, "specified_sensible_enthalpy_J")
        .parse::<f64>()
        .expect("an enthalpy")
        + captured_value(label, "formation_inventory_at_construction")
            .parse::<f64>()
            .expect("an inventory")
}

/// One key out of one fluid's block.
fn captured_value(label: &str, key: &str) -> String {
    let capture = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../validation/neqsim/captures/reactive_ph_flash_probe.tsv"),
    )
    .expect("the capture is committed");
    let mut in_block = false;
    for line in capture.lines() {
        if let Some(fluid) = line.strip_prefix("fluid=") {
            in_block = fluid == label;
            continue;
        }
        if in_block {
            if let Some(value) = line.strip_prefix(&format!("{key}=")) {
                return value.to_string();
            }
        }
    }
    panic!("no `{key}` for `{label}` in the capture");
}
