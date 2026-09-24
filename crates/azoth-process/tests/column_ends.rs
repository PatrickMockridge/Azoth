//! The column's two ends, against `validation/neqsim/captures/process_column_condenser.tsv`
//! and `..._reboiler.tsv`.
//!
//! Each end is a stage with one specification added, so these rows also exercise the tray: a
//! reboiler with no ratio is `Reboiler.run`'s `super.run`, and a condenser with no reflux is
//! `Condenser.run`'s. What the ends add is the *ratio* - and the ratio is the whole point, in
//! the parking note's words: it is a specification at the end and not a duty parameter.
//!
//! # What the capture cannot be an oracle for, and why it is the duty
//!
//! **NeqSim's own `duty` field and its own published phase streams disagree with each other
//! on every ratio row.** Its duty is `mixedStream.getFluid().getEnthalpy() -
//! calcMixStreamEnthalpy0()` - the flashed *system* against the inlets - while the outlets a
//! caller reads are `phaseToSystem` extracts of that system, scaled. Measured on these four
//! rows, the two differ by **1801, 930, 759 and 4340 J/mol** respectively, and azoth agrees
//! with whichever of the two is the property-consistent one, to `0.02-0.11` J/mol:
//!
//! | row | azoth's `h` at the captured state | the duty field's `h` | the streams' `h` |
//! |---|---|---|---|
//! | reboiler `V/B=2` | **-2715.907** | -4516.863 | -2715.757 |
//! | reboiler `V/B=0.5` | **-7924.534** | -6993.877 | -7924.435 |
//! | condenser `R=1.5` | **-7857.970** | -8616.551 | -7857.891 |
//! | condenser total `R=1.5` | **-16305.110** | -16305.129 | -11964.862 |
//!
//! So each row below asserts the duty against the figure that is a property evaluation of the
//! state the flash reports - the published streams' balance on the ratio rows, the duty field
//! on the total one - and names the other. Its own columns in the capture are all there, so a
//! reader can redo the arithmetic.

use azoth_core::units::{kelvins, pascals, watts};
use azoth_process::Stream;
use azoth_process::column::{
    CondenserMode, CondenserOutcome, ReboilerMode, ReboilerOutcome, condenser, reboiler,
};

fn feed(t: f64, p_bar: f64) -> Stream {
    Stream::from_pt(
        vec![
            "methane".into(),
            "ethane".into(),
            "propane".into(),
            "n-butane".into(),
        ],
        vec![0.1, 0.3, 0.4, 0.2],
        1.0,
        pascals(p_bar * 1.0e5),
        kelvins(t),
    )
    .expect("the fluid resolves")
}

fn relative(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let scale = 1.0 + actual.abs().max(expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{what}: {actual} vs {expected}, {} relative",
        (actual - expected).abs() / scale
    );
}

fn absolute(actual: f64, expected: f64, tolerance: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{what}: {actual} vs {expected}, {} absolute",
        (actual - expected).abs()
    );
}

fn vapour_n(out: &ReboilerOutcome) -> f64 {
    out.vapour.as_ref().expect("a boilup").n
}

fn liquid_n(out: &ReboilerOutcome) -> f64 {
    out.liquid.as_ref().expect("a bottoms stream").n
}

fn distillate_n(out: &CondenserOutcome) -> f64 {
    out.distillate.as_ref().expect("a distillate").n
}

fn reflux_n(out: &CondenserOutcome) -> f64 {
    out.reflux.as_ref().expect("a reflux").n
}

fn temperature_of(out: &CondenserOutcome) -> f64 {
    out.temperature.value
}

/// **The equilibrium reboiler, which is the stage and nothing else.** No ratio set reaches
/// `super.run`, so this is the tray's own flash at the mixture's enthalpy - and its duty is
/// zero by construction. Measured on the capture: `1.4e-11` W, the flash's residual rather
/// than a duty.
#[test]
fn an_equilibrium_reboiler_is_the_stage() {
    let out = reboiler(
        &[feed(320.0, 25.0)],
        None,
        None,
        watts(0.0),
        ReboilerMode::Equilibrium,
    )
    .expect("the reboiler solves");

    absolute(
        out.temperature.value,
        320.0000000000001,
        1e-8,
        "temperature",
    );
    relative(vapour_n(&out), 0.48390324619327096, 1e-9, "boilup");
    relative(liquid_n(&out), 0.516096753806729, 1e-9, "bottoms");
    // **Zero, and as a bound rather than an equality.** NeqSim's residual is `1.4e-11` W and
    // azoth's is `1.1e-6` W: both are the flash's own round-off, and the second is the phase
    // round trip the tray measured at `7e-6` J/mol. Writing either as the expected value
    // would be asserting a library's rounding.
    assert!(
        out.duty.value.abs() < 1.0e-4,
        "an equilibrium reboiler is adiabatic; its duty is {} W",
        out.duty.value
    );
}

/// **The boilup ratio, and the fraction it lands on exactly.** `PVrefluxflash(ratio, 1)`
/// searches for the temperature at which the *liquid's* fraction satisfies
/// `ratio = 1/beta_L - 1`; at `ratio = 2` that is `beta_L = 1/3`, and the capture's boilup is
/// `0.6666666666662084`. So the ratio is not approximated by the flash - it is what the flash
/// solves for, to twelve digits, which is why it can be the column's second specification.
#[test]
fn a_boilup_ratio_is_solved_for_and_not_approximated() {
    let two = reboiler(
        &[feed(320.0, 25.0)],
        None,
        None,
        watts(0.0),
        ReboilerMode::VaporBoilupRatio(2.0),
    )
    .expect("the reboiler solves");
    relative(
        vapour_n(&two),
        0.6666666666662084,
        1e-6,
        "boilup at V/B = 2",
    );
    relative(
        liquid_n(&two),
        0.33333333333379156,
        1e-6,
        "bottoms at V/B = 2",
    );
    relative(
        vapour_n(&two) / liquid_n(&two),
        2.0,
        1e-9,
        "the ratio itself",
    );
    relative(
        two.temperature.value,
        327.72921512704346,
        1e-6,
        "temperature",
    );
    // The duty against the published phase streams' own balance:
    // `0.6666666666662084 * 1337.483 + 0.33333333333379156 * (-10822.236) - (-5511.230)`.
    // **NeqSim's `duty_W` field says `994.37`**, which is its flashed system's enthalpy less
    // the inlets' and not the sum of the phases it publishes - a 1801 J/mol disagreement with
    // its own capture. See this file's header.
    absolute(two.duty.value, 2795.4727, 0.5, "duty");

    let half = reboiler(
        &[feed(320.0, 25.0)],
        None,
        None,
        watts(0.0),
        ReboilerMode::VaporBoilupRatio(0.5),
    )
    .expect("the reboiler solves");
    relative(
        vapour_n(&half) / liquid_n(&half),
        0.5,
        1e-6,
        "the ratio itself",
    );
    relative(
        half.temperature.value,
        312.20960118690476,
        1e-6,
        "temperature",
    );
    // The published balance, as above; the field says `-1482.65`, 930 J/mol away.
    absolute(half.duty.value, -2413.2053, 0.5, "duty");
    // A negative duty on a reboiler is a *cooler*: the ratio asked for less vapour than the
    // feed already carried, so the end has to take heat out to reach it. The class reports it
    // with the sign the balance gives, and so does this.
    assert!(half.duty.value < 0.0);
}

/// The equilibrium partial condenser: the stage again, so its duty is zero.
#[test]
fn an_equilibrium_partial_condenser_is_the_stage() {
    let out = condenser(&[feed(300.0, 20.0)], None, CondenserMode::Equilibrium)
        .expect("the condenser solves");

    absolute(
        out.temperature.value,
        300.0000000000223,
        1e-9,
        "temperature",
    );
    relative(distillate_n(&out), 0.3106080504317247, 1e-9, "distillate");
    relative(reflux_n(&out), 0.6893919495682753, 1e-9, "reflux");
    absolute(out.duty.value, 2.457454684190452e-9, 1e-6, "duty");
    assert!(out.liquid_product.is_none(), "no mode makes one here");
}

/// **A partial condenser at a reflux ratio**: the same search as the reboiler's, on the other
/// phase. `R = 1.5` puts the vapour fraction at `1/(1 + 1.5) = 0.4`, and the capture's
/// distillate is `0.4000000000009474`.
#[test]
fn a_partial_condensers_reflux_ratio_is_the_vapour_fraction() {
    let out = condenser(&[feed(300.0, 20.0)], None, CondenserMode::RefluxRatio(1.5))
        .expect("the condenser solves");

    relative(distillate_n(&out), 0.4000000000009474, 1e-9, "distillate");
    relative(reflux_n(&out), 0.5999999999990526, 1e-9, "reflux");
    relative(temperature_of(&out), 305.4051147269737, 1e-6, "temperature");
    // The published balance; the field says `815.26`, 759 J/mol away.
    absolute(out.duty.value, 1573.9208, 0.5, "duty");
}

/// **A total condenser, where the reflux is a split of the condensate and both products are
/// liquid.** The flash is a bubble-point one at 263.24 K, and the split is `R/(1 + R) = 0.6`
/// to the reflux. The capture reaches the *distillate* through `getGasOutStream()`, which for
/// this mode returns the splitter's second branch - a liquid under a gas getter's name.
#[test]
fn a_total_condenser_splits_its_condensate() {
    let out = condenser(&[feed(300.0, 20.0)], None, CondenserMode::Total(1.5))
        .expect("the condenser solves");

    relative(
        temperature_of(&out),
        263.2395720086368,
        1e-6,
        "bubble temperature",
    );
    relative(reflux_n(&out), 0.5999999999999999, 1e-9, "reflux");
    relative(distillate_n(&out), 0.4, 1e-9, "distillate");
    relative(
        reflux_n(&out) / distillate_n(&out),
        1.5,
        1e-9,
        "the ratio itself",
    );
    // **This is the row where the field is the consistent one.** The class's duty is the
    // whole condensate's enthalpy less the inlets', `1.0 * (-16305.13) - (-9431.81)`, and it
    // agrees with azoth's `h` at the bubble point to `0.02` J/mol. Its *published* streams
    // say `-11964.86` J/mol for that same condensate - `4340` J/mol above a liquid at 263 K
    // that was `-13696.96` J/mol at 300 K - which is the stale figure, not this one.
    absolute(out.duty.value, -6873.316636201762, 0.5, "duty");

    // A ratio of zero returns no reflux at all: the class publishes `1e-50` mol/s, which is
    // its `Splitter`'s zero-flow branch and not a state.
    let none = condenser(&[feed(300.0, 20.0)], None, CondenserMode::Total(0.0))
        .expect("the condenser solves");
    assert!(none.reflux.is_none(), "a zero reflux is absent, not empty");
    relative(distillate_n(&none), 1.0, 1e-9, "all of it is distillate");
    absolute(none.duty.value, -6873.316636201762, 0.5, "the same duty");
}

/// **A fixed liquid reflux, which is a flow and not a ratio.** The class takes it from the
/// liquid outlet and leaves the remainder as a *separate liquid product* - the one mode with
/// three products - and it keeps a shortfall residual rather than refusing a reflux larger
/// than the condensate.
#[test]
fn a_fixed_liquid_reflux_leaves_a_liquid_product() {
    let out = condenser(
        &[feed(300.0, 20.0)],
        None,
        CondenserMode::LiquidRefluxSplit(0.5),
    )
    .expect("the condenser solves");

    relative(
        reflux_n(&out),
        0.49999999999999994,
        1e-9,
        "the reflux is the flow asked for",
    );
    relative(
        distillate_n(&out),
        0.3106080504317247,
        1e-9,
        "the vapour is still the vapour",
    );
    let product = out.liquid_product.as_ref().expect("a liquid product").n;
    relative(product, 0.18939194956827537, 1e-9, "the remainder");
    // The condensate is conserved across the two withdrawals.
    absolute(
        reflux_n(&out) + product,
        0.6893919495682753,
        1e-9,
        "the liquid balance",
    );
    // The reflux is a *branch* of the liquid outlet, so it carries the liquid's own molar
    // enthalpy and the duty stays zero. **NeqSim's published reflux says `-17762.60` J/mol
    // against the `-13696.96` of the liquid it was taken from** - its
    // `scalePhaseSystemToNormalizedMoles` path, which a split of a stream should not move at
    // all.
    absolute(out.duty.value, 2.457454684190452e-9, 1e-6, "duty");
    // The liquid's own molar enthalpy, to the libraries' ideal-gas offset (`0.064` J/mol,
    // `tests/stream.rs`'s) - so the branch is the same state and not a re-flash of it.
    absolute(
        out.reflux.as_ref().expect("a reflux").h.value,
        -13696.958208647447,
        0.1,
        "the reflux carries the liquid's own molar enthalpy",
    );
}

/// A ratio is a ratio and a reflux flow is a flow.
#[test]
fn the_ends_refuse_what_is_not_a_ratio() {
    assert!(
        condenser(
            &[feed(300.0, 20.0)],
            None,
            CondenserMode::RefluxRatio(f64::NAN)
        )
        .is_err()
    );
    assert!(condenser(&[feed(300.0, 20.0)], None, CondenserMode::Total(-1.0)).is_err());
    assert!(
        condenser(
            &[feed(300.0, 20.0)],
            None,
            CondenserMode::LiquidRefluxSplit(-0.5)
        )
        .is_err()
    );
    assert!(
        reboiler(
            &[feed(320.0, 25.0)],
            None,
            None,
            watts(0.0),
            ReboilerMode::VaporBoilupRatio(-2.0)
        )
        .is_err()
    );
}
