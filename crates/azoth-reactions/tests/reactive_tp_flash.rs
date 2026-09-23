//! The reactive flash as a model, against the same capture the stack was built on.
//!
//! The oracle is `validation/neqsim/ReactiveFlashProbe.java` and its capture
//! `captures/reactive_flash_probe.tsv`. What is asserted here is the *model's* boundary rather
//! than the driver's: the components, the state and the overall moles go in, and the phases
//! come back with their mole numbers and the driver's own numbers beside them.
//!
//! **The composition is the pin and the split is not.** Both captured fluids end with two
//! phases of the *same* composition, so the split is a flat direction of the Gibbs energy and
//! is decided by the path; the moles summed over the phases are what the element balance and
//! the equilibrium fix. That is why the overall is asserted here and the per-phase numbers are
//! not - the driver's own tests do that, where the two codes are shown to land on different
//! points of the same line.

use azoth_reactions::reactive_tp_flash::reactive_tp_flash;

const BAR: f64 = 1.0e5;

/// Sum the phases back into the fluid's composition.
fn overall(result: &azoth_reactions::reactive_tp_flash::ReactiveTpFlashResult) -> Vec<f64> {
    let components = result.phase_moles[0].len();
    (0..components)
        .map(|i| result.phase_moles.iter().map(|phase| phase[i]).sum())
        .collect()
}

/// **The class's own fluid, at 600 K.** `CO + H2O = CO2 + H2` has one independent reaction
/// over three elements and its products are strongly favoured. The capture's
/// `equilibrium_moles[0]` is the whole fluid's composition, because both its phases hold it.
#[test]
fn the_water_gas_shift_is_the_captured_composition() {
    let components: Vec<String> = ["CO", "water", "CO2", "hydrogen"]
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    let result = reactive_tp_flash(&components, 600.0, BAR, &[0.25; 4], 2).expect("the flash runs");

    assert!(result.converged, "the captured state converges");
    assert_eq!(result.phase_count, 2, "the constructor's pair survives");
    let captured = [
        0.079_140_502_548_023_4,
        0.079_140_502_560_693_7,
        0.420_859_061_906_977_14,
        0.420_859_061_880_687_67,
    ];
    for (index, (got, want)) in overall(&result).iter().zip(captured).enumerate() {
        let relative = (got - want).abs() / want.abs();
        assert!(
            relative < 1.0e-5,
            "component {index}: {got} against the capture's {want}, a relative {relative:.3e}"
        );
    }

    // `final_gibbs_energy`, **twice** one phase's worth: the driver weighs each phase by the
    // fraction its list carries, and the constructor leaves both at one.
    let captured_gibbs = -2.259_535_542_715_054_7_f64;
    assert!(
        (result.gibbs_energy - captured_gibbs).abs() / captured_gibbs.abs() < 1.0e-5,
        "the Gibbs measure is {} against the capture's {captured_gibbs}",
        result.gibbs_energy
    );
    // The total moves away from the feed's one mole only by the element residual, and it is
    // `getEquilibriumTotalMoles` either way.
    assert!((result.equilibrium_total_moles - 1.0).abs() < 1.0e-4);
}

/// **The four-component fluid at 1000 K**, where more than one reaction runs - and where the
/// comparison has to say what it is comparing.
///
/// **Both codes stop on the relaxed multiphase tolerance, and at different points inside it**:
/// this port at `5.5e-5` and NeqSim at `8.0e-6`, against a tolerance of `1e-4`. The species
/// that carry most of the fluid agree to `3.6e-5`, and the small species are where a stopping
/// point shows - **water at `3.6e-2` off by `7.2e-4`, hydrogen at `5.3e-1` off by `1.1e-4`**.
/// So the composition band here is `1e-3` and the *residual* is the pin: `5.5e-5` against the
/// `1e-4` the class accepts. A caller comparing compositions on a relaxed stop is comparing
/// stopping points.
#[test]
fn the_hot_four_component_fluid_is_the_captured_composition() {
    let components: Vec<String> = ["methane", "water", "CO2", "hydrogen"]
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    let result =
        reactive_tp_flash(&components, 1000.0, BAR, &[0.4, 0.2, 0.2, 0.2], 2).expect("it runs");

    assert!(result.converged);
    // The capture's residual, and this port's: both inside the `1e-4` the class accepts, and
    // that is the band the composition comparison lives in.
    assert!(
        result.residual < 1.0e-4,
        "the solve stopped at {} against the class's relaxed tolerance of 1e-4",
        result.residual
    );
    let captured = [
        0.318_153_488_941_055_16,
        0.036_306_383_641_311_11,
        0.281_847_077_069_566_24,
        0.527_387_691_620_036_2,
    ];
    for (index, (got, want)) in overall(&result).iter().zip(captured).enumerate() {
        let relative = (got - want).abs() / want.abs();
        assert!(
            relative < 1.0e-3,
            "component {index}: {got} against the capture's {want}, a relative {relative:.3e}"
        );
    }
    let captured_gibbs = -2.329_146_095_859_745_6_f64;
    assert!(
        (result.gibbs_energy - captured_gibbs).abs() / captured_gibbs.abs() < 1.0e-4,
        "the Gibbs measure is {} against the capture's {captured_gibbs}",
        result.gibbs_energy
    );
}

/// **A fluid with no reaction to run** takes the conventional fallback, and the model reports
/// it the same way: the phases, their fractions, and a driver that counts no iterations.
#[test]
fn the_non_reactive_fluid_takes_the_fallback() {
    let components: Vec<String> = ["methane", "water"]
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    let result =
        reactive_tp_flash(&components, 300.0, 50.0 * BAR, &[0.5, 0.5], 2).expect("the flash runs");

    assert!(result.converged);
    assert_eq!(result.total_iterations, 0, "the fallback counts none");
    assert!(
        (result.phase_fraction[0] - 0.500_334_895_161_083_2).abs() < 1.0e-12,
        "the vapour fraction is {} against the capture's 0.5003348951610832",
        result.phase_fraction[0]
    );
    // The phase's *composition*, from its own moles: the gas is methane-rich, and the
    // capture's `phase_x[0]` is that number.
    let held: f64 = result.phase_moles[0].iter().sum();
    let methane = result.phase_moles[0][0] / held;
    assert!(
        (methane - 0.999_330_366_296_782).abs() < 1.0e-4,
        "the gas is {methane} methane against the capture's 0.999330366296782"
    );
}

/// **The ionic branch is refused, not approximated.** NeqSim pins a gas-phase ion to `EPS` and
/// corrects an electrolyte phase's reference state through `getLogInfiniteDiluteFugacity`;
/// neither is ported, so a charged component is an error rather than a neutral fluid's answer.
#[test]
fn a_charged_component_is_refused() {
    let components: Vec<String> = ["water", "OH-"]
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    let refused = reactive_tp_flash(&components, 298.15, BAR, &[1.0, 1.0e-7], 2);
    assert!(refused.is_err(), "a charged component is refused");
}
