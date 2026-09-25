//! The reactive stage, against NeqSim's own claims about it.

use azoth_core::units::{kelvins, pascals, watts};
use azoth_process::Stream;
use azoth_process::column::{SideDraws, tray};

fn relative(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let scale = 1.0 + actual.abs().max(expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{what}: {actual} vs {expected}, {} relative",
        (actual - expected).abs() / scale
    );
}

/// A fluid with no independent reaction - methane and ethane - which is NeqSim's own
/// `testReactiveColumnMassBalanceNR0` fluid.
fn unreactive(t: f64, p_bar: f64) -> Stream {
    Stream::from_pt(
        vec!["methane".to_string(), "ethane".to_string()],
        vec![0.7, 0.3],
        1.0,
        pascals(p_bar * 1.0e5),
        kelvins(t),
    )
    .expect("the fluid builds")
}

/// The water-gas shift, which has one independent reaction.
fn wgs(t: f64, p_bar: f64) -> Stream {
    Stream::from_pt(
        ["CO", "water", "CO2", "hydrogen"]
            .iter()
            .map(|name| (*name).to_string())
            .collect(),
        vec![0.25; 4],
        1.0,
        pascals(p_bar * 1.0e5),
        kelvins(t),
    )
    .expect("the fluid builds")
}

/// **The class's own claim about the reactive route: where the fluid has no reaction, it is the
/// plain flash.** `testReactiveColumnMassBalanceNR0` asserts a reactive column matching a
/// standard one to `0.01` kg/hr on exactly this fluid, and the flash probe measures the `NR = 0`
/// delegation returning in zero iterations. So the two stages must land on the same state, and
/// they do - through two different flashes, which is why the gate is loose.
#[test]
fn a_stage_with_no_reaction_to_run_is_the_plain_stage() {
    let inlets = [unreactive(300.0, 15.0)];
    let plain = tray(&inlets, None, None, watts(0.0), SideDraws::NONE, false)
        .expect("the plain tray solves");
    let reactive = tray(&inlets, None, None, watts(0.0), SideDraws::NONE, true)
        .expect("the reactive tray solves");

    relative(
        reactive.temperature.value,
        plain.temperature.value,
        1.0e-4,
        "the two stages' temperatures",
    );
    let (a, b) = (
        reactive.gas.as_ref().map_or(0.0, |g| g.n),
        plain.gas.as_ref().map_or(0.0, |g| g.n),
    );
    assert!(a > 0.0 && b > 0.0, "both stages hold a vapour: {a}, {b}");
    relative(a, b, 1.0e-3, "the two vapours");
    let (c, d) = (
        reactive.liquid.as_ref().map_or(0.0, |l| l.n),
        plain.liquid.as_ref().map_or(0.0, |l| l.n),
    );
    relative(c, d, 1.0e-3, "the two liquids");
    assert!(
        reactive.warnings.is_empty(),
        "nothing fell back on this fluid: {:?}",
        reactive.warnings
    );
}

/// **Where the reactive answer has two phases of the same composition, the stage refuses.** The
/// water-gas shift is the class's own reactive example and its two phases converge to one
/// composition there - measured, on both cubics - so there is nothing to say which row leaves up
/// the column and which leaves down. A tray that picked one would be inventing an outlet, and
/// this port reports the absence as an absence.
#[test]
fn a_two_phase_reactive_answer_without_a_label_is_refused() {
    let inlets = [wgs(600.0, 1.0)];
    let error = tray(&inlets, None, None, watts(0.0), SideDraws::NONE, true)
        .expect_err("the two phases are the same state, so no outlet pair is determined");
    let message = format!("{error}");
    assert!(
        message.contains("same composition") || message.contains("label"),
        "{message}"
    );
    // The plain stage on the same fluid is an ordinary tray and is not refused.
    let plain = tray(&inlets, None, None, watts(0.0), SideDraws::NONE, false)
        .expect("the plain tray solves");
    assert!(plain.gas.is_some() || plain.liquid.is_some());
}
