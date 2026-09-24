//! The two schemes `PlugFlowReactor.run` marches with, and the loop that applies them.
//!
//! **A fixed step and no step-size control, because that is what the class has.** `dz` is
//! `length / numberOfSteps`, computed once (`PlugFlowReactor.java:326`), and the loop runs
//! exactly `numberOfSteps` times. Nothing here adapts, retries or estimates an error - a port
//! that added a controller would answer a different question than the class does, which is the
//! one thing a port may not do.
//!
//! **The state vector is the class's.** `[F₁ … Fₙ, T, P]` - component molar flows in mol/s,
//! temperature in K and pressure in **bara** (`PlugFlowReactor.java:752-780`'s own javadoc).
//! The pressure row's derivative converts Pa/m to bar/m rather than the step converting the
//! state, so the bara is carried through the whole march and only the caller decides what to
//! do with it.
//!
//! **What is not here is the fixups, and they are the caller's.** After each step the class
//! clamps the molar flows at zero and the pressure at `0.1`, overrides an isothermal reactor's
//! temperature back to its inlet, and records the station. Those are the reactor's state
//! discipline rather than the integrator's, so [`march`] takes them as a `settle` callback
//! that runs in the loop's own order - after the advance, before the next derivative.

/// Which scheme advances the profile.
///
/// The class's `IntegrationMethod`. RK4 is its default (`PlugFlowReactor.java:177`), and Euler
/// is reachable only by asking for it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scheme {
    /// First-order forward Euler: one derivative evaluation per step.
    Euler,
    /// Classical fourth-order Runge-Kutta: four evaluations, weights `1:2:2:1` over six.
    Rk4,
}

impl Scheme {
    /// The scheme a class's `setIntegrationMethod` name selects.
    ///
    /// NeqSim's setter matches `"EULER"` case-sensitively and falls back to RK4 for anything
    /// else (`PlugFlowReactor.java:1028`), which is the rule here: an unknown name is not a
    /// refusal but the default, and the caller that cares compares against [`Scheme::Rk4`].
    #[must_use]
    pub fn named(name: &str) -> Self {
        if name == "EULER" {
            Self::Euler
        } else {
            Self::Rk4
        }
    }
}

/// Advance `state` in place by `steps` steps of `dz`, evaluating `derivatives` as the class
/// does and calling `settle` after every advance.
///
/// `derivatives` is `FnMut` and not `Fn` because the class's own is: under
/// `FULLY_COUPLED` it writes the trial state back into the thermodynamic system before
/// reading properties off it, so the four RK4 evaluations are not independent and each one
/// mutates the system the next reads. Flattening that into a pure function would change the
/// arithmetic.
///
/// `settle` sees the state after each step in the class's own order, which is advance, then
/// clamp, then override, so a caller that needs the class's profile records it there. It runs
/// *before* the step count is tested, so a one-step march settles once.
pub fn march(
    state: &mut [f64],
    steps: usize,
    dz: f64,
    scheme: Scheme,
    mut derivatives: impl FnMut(&[f64]) -> Vec<f64>,
    mut settle: impl FnMut(&mut [f64]),
) {
    let mut trial = state.to_vec();
    for _ in 0..steps {
        match scheme {
            Scheme::Euler => {
                let k = derivatives(state);
                for (value, derivative) in state.iter_mut().zip(&k) {
                    *value += derivative * dz;
                }
            }
            Scheme::Rk4 => {
                let k1 = derivatives(state);
                scaled_into(state, &k1, 0.5 * dz, &mut trial);
                let k2 = derivatives(&trial);
                scaled_into(state, &k2, 0.5 * dz, &mut trial);
                let k3 = derivatives(&trial);
                scaled_into(state, &k3, dz, &mut trial);
                let k4 = derivatives(&trial);
                for (i, value) in state.iter_mut().enumerate() {
                    *value += dz / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
                }
            }
        }
        settle(state);
    }
}

/// `out = base + scale * delta`, the class's `addScaled`.
///
/// It writes a *fresh* trial state from the pre-step state every time, which is what makes the
/// four RK4 stages the class's: each is taken from the step's own origin and not from the last
/// stage's endpoint.
fn scaled_into(base: &[f64], delta: &[f64], scale: f64, out: &mut [f64]) {
    for (slot, (start, step)) in out.iter_mut().zip(base.iter().zip(delta)) {
        *slot = start + scale * step;
    }
}
