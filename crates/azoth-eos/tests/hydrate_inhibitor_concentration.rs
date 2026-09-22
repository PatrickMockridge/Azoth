//! Spec-driven tests for the `eos.hydrate_inhibitor_concentration` model.
//!
//! The oracle is `validation/neqsim/captures/hydrate_inhibitor_probe.tsv`: NeqSim's own
//! `HydrateInhibitorConcentrationFlash` on the composition its `main` uses, run on two fluids.
//! This file pins the **plain `SystemSrkEos`** column, which is the fluid this library can
//! build, and records the CPA one - see the case's `source` for why the two are `1.663` and
//! `0.326` mol of MEG apart.

// The oracle values are NeqSim's ground truth, so their full digits are the point.
#![allow(clippy::excessive_precision)]

use azoth_core::AzothError;
use azoth_core::units::{kelvins, pascals};
use azoth_eos::hydrate::{self, HydrateModel};
use azoth_eos::{Cubic, hydrate_inhibitor_concentration, model_gen};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.hydrate_inhibitor_concentration";

/// NeqSim's own `main` for this class: methane, ethane, propane, i-butane, MEG and water.
const NAMES: [&str; 6] = ["methane", "ethane", "propane", "i-butane", "MEG", "water"];

const MOLES: [f64; 6] = [1.0, 0.10, 0.050, 0.0050, 0.1, 1.0];

/// `(target in K, the capture's plain-SRK inhibitor moles, its mass fraction)`.
///
/// **The mole number agrees to `2e-8` and the mass fraction to `4e-9`, and the two bars
/// differ for a reason.** The secant stops when its residual is inside `1e-3` K, so the mole
/// number it lands on is only fixed to that band divided by `dT/dC` - which is why the looser
/// bar is the mole number's. Measured across the three: `2.4e-9`, `1.1e-8` and `2.1e-8` for
/// the moles, `3.6e-10`, `1.3e-9` and `3.8e-9` for the fraction.
const STATES: [(f64, f64, f64); 3] = [
    (270.9, 1.663_215_479_419_76, 0.851_421_604_021_146),
    (265.0, 2.163_911_968_852_98, 0.881_734_574_278_325),
    (275.0, 1.347_410_172_703_67, 0.822_769_695_596_269),
];

fn mixture() -> azoth_eos::mixture::Mixture {
    hydrate::hydrate_mixture_of(&NAMES, Cubic::Srk, None, HydrateModel::Pvtsim)
        .expect("resolves")
        .0
}

#[test]
fn the_plain_cubic_column_reproduces_the_capture() {
    let fluid = mixture();
    for (target, want_moles, want_fraction) in STATES {
        let result =
            hydrate_inhibitor_concentration(&fluid, "MEG", &MOLES, kelvins(target), pascals(1.0e7))
                .unwrap_or_else(|e| panic!("{target} K: {e:?}"));
        assert!(
            (result.inhibitor_moles / want_moles - 1.0).abs() < 1.0e-6,
            "{target} K: {} mol against NeqSim's {want_moles}",
            result.inhibitor_moles
        );
        common::assert_close(
            result.weight_fraction,
            want_fraction,
            1.0e-8,
            &format!("{target} K: the mass fraction"),
        );
        // The residual is the secant's own stopping rule, and NeqSim's tolerance is `1e-3`.
        assert!(
            result.residual.abs() <= 1.0e-3,
            "{target} K: the residual is {}",
            result.residual
        );
        assert!(
            (result.hydrate_temperature.value - target).abs() <= 1.0e-3,
            "{target} K: reached {}",
            result.hydrate_temperature.value
        );
    }
}

/// **NeqSim's secant always runs three steps**, whatever the residual is, so the count is a
/// floor rather than a convergence history.
#[test]
fn the_secant_takes_at_least_three_steps() {
    let fluid = mixture();
    let result =
        hydrate_inhibitor_concentration(&fluid, "MEG", &MOLES, kelvins(290.0), pascals(1.0e7))
            .expect("a target close to the uninhibited temperature still lands");
    assert!(result.iterations >= 3, "{} steps", result.iterations);
}

/// A target the feed already meets needs no inhibitor, and the secant finds that.
#[test]
fn a_target_near_the_uninhibited_temperature_needs_little_inhibitor() {
    let fluid = mixture();
    let warm =
        hydrate_inhibitor_concentration(&fluid, "MEG", &MOLES, kelvins(290.0), pascals(1.0e7))
            .expect("solves");
    let cold =
        hydrate_inhibitor_concentration(&fluid, "MEG", &MOLES, kelvins(265.0), pascals(1.0e7))
            .expect("solves");
    assert!(
        cold.inhibitor_moles > warm.inhibitor_moles,
        "a colder target needed {} mol and a warmer one {}",
        cold.inhibitor_moles,
        warm.inhibitor_moles
    );
}

/// An inhibitor the feed does not carry is refused: there is nothing to add to it.
#[test]
fn an_inhibitor_the_feed_does_not_have_is_refused() {
    let error = hydrate_inhibitor_concentration(
        &mixture(),
        "methanol",
        &MOLES,
        kelvins(270.9),
        pascals(1.0e7),
    )
    .expect_err("the feed has no methanol");
    assert!(
        matches!(error, AzothError::InvalidInput { .. }),
        "{error:?}"
    );
}

/// A feed with no water has nothing to inhibit.
#[test]
fn a_feed_with_no_water_is_refused() {
    // A hydrate mixture refuses a water-free feed before this model is reached, so the only
    // way to hold one is through `Mixture` directly - which is what a caller with a dry feed
    // would have to do, and what this model's own check is for.
    let (dry, _) = azoth_eos::databank::mixture_of(&["methane", "ethane"], Cubic::Srk, None)
        .expect("the pair resolves");
    let error = hydrate_inhibitor_concentration(
        &dry,
        "methane",
        &[0.5, 0.5],
        kelvins(270.9),
        pascals(1.0e7),
    )
    .expect_err("a mixture with no hydrate tables has no inhibitor to walk");
    assert!(
        matches!(error, AzothError::InvalidInput { .. }),
        "{error:?}"
    );
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty());
    for case in spec.cases {
        let names = case
            .list("components")
            .expect("the case declares components");
        let fluid = hydrate::hydrate_mixture_of(names, Cubic::Srk, None, HydrateModel::Pvtsim)
            .expect("the case's components resolve for a hydrate")
            .0;
        let result = hydrate_inhibitor_concentration(
            &fluid,
            common::input_str(case, "inhibitor"),
            case.vector("moles").expect("moles"),
            kelvins(common::input(case, "T_target")),
            pascals(common::input(case, "P")),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));
        common::assert_close(
            result.inhibitor_moles,
            common::expected(case, "inhibitor_moles"),
            case.tolerance,
            &format!("{}::{} (inhibitor moles)", spec.id, case.id),
        );
        common::assert_close(
            result.weight_fraction,
            common::expected(case, "weight_fraction"),
            case.tolerance,
            &format!("{}::{} (mass fraction)", spec.id, case.id),
        );
    }
}
