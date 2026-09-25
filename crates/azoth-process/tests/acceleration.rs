//! The tear's two accelerators, held to NeqSim step by step.
//!
//! `validation/neqsim/captures/process_acceleration.tsv` is the oracle, and it exists because a
//! flowsheet cannot be one for this: on `specs/flowsheets/` the only loop accumulates in its
//! *flow*, and `Recycle.applyStreamValues` writes back the **composition and nothing else** - so
//! an acceleration there is invisible however wrong it is. What the probe drives instead is the
//! arithmetic, on a scripted sequence, where a step one call early is a different number.
//!
//! **The two methods are checkable to different depths, and the difference is NeqSim's.**
//! `BroydenAccelerator.accelerate` is public and self-contained, so the port is held to its
//! *answer* at every step. Wegstein's accelerator is `private`, and what the probe can read back
//! through `Recycle.run` is the `q` factors - so the port is held to those, and **not** to the
//! outlet the class published: `applyStreamValues` writes normalised fractions through
//! `Component.setx` on both phases, and the two-phase system reads back *different* numbers from
//! the ones it was handed (pass 4 writes `[1.0, 0.0]` and reads `[1.0, 0.4]`). That is a `setx`
//! bookkeeping property of a two-phase fluid, and a port whose `Stream` carries one composition
//! has nothing to model it with.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use azoth_process::recycle::{BroydenAccelerator, WEGSTEIN_Q_MAX, WEGSTEIN_Q_MIN, wegstein};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// The capture, as `key -> the numbers in its value`.
///
/// A `null` value - Wegstein's `q` before the branch is taken - is absent from the map rather
/// than present as a zero, because "no `q` was computed" and "a `q` of zero" are the two states
/// the delay boundary is drawn between.
fn capture() -> BTreeMap<String, Vec<f64>> {
    let text =
        std::fs::read_to_string(root().join("validation/neqsim/captures/process_acceleration.tsv"))
            .expect("the acceleration capture is there");
    text.lines()
        .filter_map(|line| line.split_once('='))
        .filter(|(_, value)| *value != "null")
        .map(|(key, value)| {
            let numbers = value
                .split_whitespace()
                .map(|entry| entry.parse().expect("a number"))
                .collect();
            (key.to_string(), numbers)
        })
        .collect()
}

fn vector(row: &BTreeMap<String, Vec<f64>>, key: &str) -> Vec<f64> {
    row.get(key)
        .unwrap_or_else(|| panic!("the capture has no `{key}`"))
        .clone()
}

#[test]
fn broyden_reproduces_the_captured_step_at_every_call() {
    let row = capture();

    // **The scalar sequence first**, because its fixed point is checkable by hand: `g(x) = 0.5x +
    // 1` has `x* = 2`, and the whole sequence is five numbers a reader can follow.
    let mut scalar = BroydenAccelerator::new();
    for step in 1..=5 {
        let x = vector(&row, &format!("broyden_scalar_x_{step}"));
        let g = vector(&row, &format!("broyden_scalar_g_{step}"));
        let accelerated = scalar.accelerate(&x, &g);
        assert_eq!(
            accelerated,
            vector(&row, &format!("broyden_scalar_acc_{step}")),
            "the scalar case, call {step}"
        );
    }

    // And the three-variable one, where the Jacobian update actually has a matrix to correct.
    let dimension: usize = vector(&row, "broyden_dim")[0] as usize;
    assert_eq!(dimension, 3);
    let mut accelerator = BroydenAccelerator::new();
    for step in 1..=8 {
        let x = vector(&row, &format!("broyden_x_{step}"));
        let g = vector(&row, &format!("broyden_g_{step}"));
        assert_eq!(x.len(), dimension);
        let accelerated = accelerator.accelerate(&x, &g);
        assert_eq!(
            accelerated,
            vector(&row, &format!("broyden_acc_{step}")),
            "the coupled case, call {step}"
        );
        assert_eq!(
            f64::from(accelerator.iteration_count()),
            vector(&row, &format!("broyden_count_{step}"))[0],
            "the counter the delay is measured in"
        );
    }
}

#[test]
fn wegstein_reproduces_the_captured_q_factors() {
    let row = capture();
    let passes = vector(&row, "wegstein_passes")[0] as usize;

    // `previous_input`/`previous_output` are the pair of the *last call the branch was taken on*,
    // which is what `applyWegsteinToStream` stores after each one - and which is `None` before the
    // first, since the class only reaches the store through that method.
    let mut previous: Option<(Vec<f64>, Vec<f64>)> = None;
    for pass in 1..=passes {
        let output = vector(&row, &format!("wegstein_in_{pass}"));

        // The branch is not taken on the first two passes: `iterations > wegsteinDelayIterations`
        // excludes them, and the capture records `null` rather than a `q` vector for both.
        if pass <= 2 {
            assert!(
                !row.contains_key(&format!("wegstein_q_{pass}")),
                "pass {pass} took the branch, and the class's delay says it cannot"
            );
            // What the recycle publishes on those passes is the flash's own answer, which is what
            // the next pass reads as its input.
            previous = None;
            continue;
        }

        // The input a pass sees is the previous pass's *published* vector.
        let input = vector(&row, &format!("wegstein_out_{}", pass - 1));
        let step = wegstein(
            &input,
            &output,
            previous.as_ref().map(|(input, _)| input.as_slice()),
            previous.as_ref().map(|(_, output)| output.as_slice()),
            WEGSTEIN_Q_MIN,
            WEGSTEIN_Q_MAX,
        );

        assert_eq!(
            step.q_factors,
            vector(&row, &format!("wegstein_q_{pass}")),
            "pass {pass}: the q factors"
        );
        assert_eq!(
            f64::from(passes as u32),
            vector(&row, "wegstein_passes")[0],
            "a run's pass count is the class's own field"
        );
        previous = Some((input, output));
    }

    // **The third pass is the seed and the fourth is the first step**, and the capture shows both:
    // an all-zero `q` is the identity that stores the secant, and pass 4 is the first one with a
    // real slope in it.
    assert_eq!(vector(&row, "wegstein_q_3"), vec![0.0; 5]);
    assert!(
        vector(&row, "wegstein_q_4").iter().any(|q| *q != 0.0),
        "the fourth pass is where the secant first moves something"
    );
}

/// **The class's Broyden step has the wrong sign, and the port reproduces it.**
///
/// On `g(x) = 0.5x + 1` the fixed point is `x* = 2`. The residual is `f(x) = g(x) - x = 1 - 0.5x`,
/// so `f'(x) = -0.5` and **one** Newton step from anywhere lands on `2`. The class's own comment
/// says the step it takes is `dx = -B^{-1} f(x)`; the code computes `B^{-1} f` and **adds** it.
///
/// The capture is what makes that a measurement rather than a reading: from `x = 1.5` the class
/// answers `1.0`, then `0.0`, then `-2.0` - walking away from a fixed point it could have reached
/// exactly, in one step, in the other direction. The port must reproduce that sequence, because
/// the rule is that the class wins where the two disagree - so this test asserts **both** halves:
/// the divergence, and that the port diverges with it.
///
/// **It is an upstream defect and it is worth reporting**: a `Recycle` declaring
/// `accelerationMethod = BROYDEN` is pushed away from its solution by the accelerator, where
/// `DIRECT_SUBSTITUTION` would converge and `WEGSTEIN` would accelerate.
#[test]
fn the_class_s_broyden_step_walks_away_from_the_fixed_point() {
    let row = capture();
    let mut accelerator = BroydenAccelerator::new();

    // What the capture says, replayed through the port.
    let mut value = vector(&row, "broyden_scalar_x_1")[0];
    for step in 1..=5 {
        let g = 0.5 * value + 1.0;
        assert_eq!(
            g,
            vector(&row, &format!("broyden_scalar_g_{step}"))[0],
            "the map the probe scripted"
        );
        value = accelerator.accelerate(&[value], &[g])[0];
        assert_eq!(
            value,
            vector(&row, &format!("broyden_scalar_acc_{step}"))[0],
            "call {step}"
        );
    }

    // The fixed point is 2, the capture's fourth call answers 0.0, and its fifth answers -2.0.
    // `DIRECT_SUBSTITUTION` from the same three-call history would be walking *towards* it: the
    // sequence `0, 1, 1.5` is `g` applied three times, and the next direct step is `1.75`.
    assert_eq!(value, -2.0, "the fifth call");
    let direct = 0.5 * vector(&row, "broyden_scalar_x_5")[0] + 1.0;
    assert!(
        (direct - 2.0).abs() < (value - 2.0).abs(),
        "direct substitution at {direct} is closer to the fixed point 2 than the accelerator's \
         {value}, which is the whole of the finding"
    );
    // And the one Newton step the class's own comment describes, for the record.
    let newton: f64 = 1.5 - (1.0 - 0.5 * 1.5) / -0.5;
    assert!(
        (newton - 2.0).abs() < 1e-12,
        "Newton reaches it exactly: {newton}"
    );
}

/// **`applyStreamValues`' own two guards, which nothing else exercises.**
///
/// The class clamps each accelerated fraction at zero and normalises the sum to one - and
/// **skips** the write entirely when that sum is at or below `1e-15`, rather than dividing by it.
/// Both are reachable from an accelerated step that overshoots past a physical composition, which
/// is what a secant extrapolation does, and both are the class's own behaviour rather than a
/// tidying.
#[test]
fn the_composition_write_back_clamps_normalises_and_declines_to_divide() {
    use azoth_core::units::{kelvins, pascals};
    use azoth_process::Stream;
    use azoth_process::recycle::{apply_composition, extract};

    let stream = || {
        Stream::from_pt(
            vec!["methane".to_string(), "n-butane".to_string()],
            vec![0.6, 0.4],
            1.0,
            pascals(5.0e5),
            kelvins(300.0),
        )
        .expect("the two resolve")
    };

    // The vector's first three slots are T, P and the flow, and they are read but not written.
    let values = extract(&stream());
    assert_eq!(values.len(), 5);
    assert_eq!(values[0], 300.0);
    assert_eq!(values[1], 5.0e5);
    assert_eq!(values[2], 1.0);

    // An overshoot past zero: component 1 is clamped, and the sum renormalised to one.
    let written = apply_composition(&stream(), &[300.0, 5.0e5, 1.0, 1.2, -0.2]);
    assert_eq!(written.z, vec![1.0, 0.0]);
    assert_eq!(written.n, 1.0, "the flow is left alone");
    assert_eq!(written.t.value, 300.0, "and so is the temperature");

    // **A vector summing to at most `1e-15` leaves the stream as it was** rather than dividing by
    // a zero - the class's guard, and the reason the two states are distinguishable at all.
    let untouched = apply_composition(&stream(), &[300.0, 5.0e5, 1.0, 0.0, 0.0]);
    assert_eq!(untouched.z, vec![0.6, 0.4], "it declined to normalise");

    // A vector too short for the component count is left alone too, which the class reaches by
    // its own `values.length >= 3 + numComponents` test.
    let short = apply_composition(&stream(), &[300.0, 5.0e5, 1.0, 1.0]);
    assert_eq!(short.z, vec![0.6, 0.4]);
}
