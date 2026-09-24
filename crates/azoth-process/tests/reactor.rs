//! The reactor tier's internals: the axial stepper, on its own.
//!
//! **A stepper's oracle is an integral with a closed form.** Every other kernel in this crate
//! is pinned against a NeqSim capture, but a capture of a *stepper* would only be another
//! stepper's answer - so this holds `reactor::stepper::march` to arithmetic instead: the same
//! equations integrated exactly, at two step sizes, and the two schemes held to their own
//! orders. It is the one test in the tier that no NeqSim row can replace.

use azoth_process::reactor::stepper::{Scheme, march};

/// The exact solution of `dy/dz = -k y`, `y(0) = y0`.
fn decaying(z: f64, k: f64, y0: f64) -> f64 {
    y0 * (-k * z).exp()
}

/// March `dy/dz = -k y` over `length` in `steps` steps and return the final `y`.
fn integrate(k: f64, y0: f64, length: f64, steps: usize, scheme: Scheme) -> f64 {
    let mut state = vec![y0];
    let dz = length / steps as f64;
    march(
        &mut state,
        steps,
        dz,
        scheme,
        |state| vec![-k * state[0]],
        |_| {},
    );
    state[0]
}

/// **RK4 is fourth order and Euler is first, and the port's weights are the classical ones.**
///
/// If the `1:2:2:1` weights or the `1/6` were wrong, RK4 would fall to first or second order
/// and the ratio below would not be `16`.
#[test]
fn rk4_is_fourth_order_and_euler_is_first() {
    let (k, y0, length) = (0.7, 2.0, 1.0);
    let exact = decaying(length, k, y0);

    let coarse = (integrate(k, y0, length, 8, Scheme::Rk4) - exact).abs();
    let fine = (integrate(k, y0, length, 16, Scheme::Rk4) - exact).abs();
    let ratio = coarse / fine;
    assert!(
        (ratio - 16.0).abs() < 1.0,
        "halving the step should cut RK4's error sixteenfold, not {ratio} (coarse {coarse:e}, fine {fine:e})"
    );

    let euler_coarse = (integrate(k, y0, length, 8, Scheme::Euler) - exact).abs();
    let euler_fine = (integrate(k, y0, length, 16, Scheme::Euler) - exact).abs();
    let euler_ratio = euler_coarse / euler_fine;
    assert!(
        (euler_ratio - 2.0).abs() < 0.2,
        "halving the step should halve Euler's error, not scale it by {euler_ratio}"
    );

    // The two schemes are not the same arithmetic, which is the point of carrying both.
    assert!(
        euler_coarse > 100.0 * coarse,
        "Euler ({euler_coarse:e}) should be far behind RK4 ({coarse:e}) on the same step"
    );
}

/// **A march terminates on its own step count, and the length is the steps times the step.**
///
/// The class computes `dz = length / numberOfSteps` and loops exactly `numberOfSteps` times, so
/// the profile ends at `length` and nowhere else.
#[test]
fn the_profile_ends_at_the_length_it_was_given() {
    let mut state = vec![0.0];
    let mut stations = vec![0.0];
    let (length, steps) = (2.5, 5);
    let dz = length / steps as f64;

    march(
        &mut state,
        steps,
        dz,
        Scheme::Rk4,
        |_| vec![1.0],
        |state| stations.push(state[0]),
    );

    assert_eq!(
        stations.len(),
        steps + 1,
        "one station per step, plus the start"
    );
    // `dy/dz = 1` from `y(0) = 0` is `y = z`, so the state is the position.
    assert!(
        (state[0] - length).abs() < 1e-12,
        "the march should end at {length}, not {}",
        state[0]
    );
    assert!(
        (stations[steps] - length).abs() < 1e-12,
        "and the last station should be the last state"
    );
}

/// **`settle` runs after every advance, and what it writes is what the next step reads.**
///
/// The reactor clamps its molar flows at zero and its pressure at `0.1` there, and overrides an
/// isothermal reactor's temperature - so a `settle` that did not take effect before the next
/// derivative would silently change the answer.
#[test]
fn settle_rewrites_the_state_before_the_next_step() {
    let mut state = vec![0.0];
    let mut saw = Vec::new();

    // `dy/dz = 1` makes the state the step index; `settle` then doubles it, so each step starts
    // from twice the last. Holding 1, 3, 7, 15 is that rule and not the bare march's 1, 2, 3.
    march(
        &mut state,
        4,
        1.0,
        Scheme::Euler,
        |state| {
            saw.push(state[0]);
            vec![1.0]
        },
        |state| state[0] *= 2.0,
    );

    assert_eq!(
        saw,
        vec![0.0, 2.0, 6.0, 14.0],
        "the next derivative must see the settled state"
    );
    assert_eq!(state[0], 30.0, "and the last settle is the answer");
}

/// **The scheme name follows the class's setter, fallback and all.**
///
/// `PlugFlowReactor.setIntegrationMethod` matches `"EULER"` and answers RK4 otherwise - an
/// unknown name is the default rather than a refusal, which is a behaviour and not an
/// oversight.
#[test]
fn an_unknown_scheme_name_is_the_default_and_not_a_refusal() {
    assert_eq!(Scheme::named("EULER"), Scheme::Euler);
    assert_eq!(Scheme::named("RK4"), Scheme::Rk4);
    assert_eq!(Scheme::named("runge"), Scheme::Rk4);
    assert_eq!(
        Scheme::named("euler"),
        Scheme::Rk4,
        "the class's match is case-sensitive"
    );
}
