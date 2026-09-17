//! The solid argon Helmholtz equation, checked against NeqSim.

// The oracle values are NeqSim's ground truth, so their full digits are the point, and
// the per-state tuple is clearer than a factored type.
#![allow(clippy::excessive_precision, clippy::type_complexity)]

use azoth_eos::argon_solid;
use azoth_eos::dual;

fn assert_close(actual: f64, expected: f64, rel_tol: f64) {
    let bound = (expected.abs() * rel_tol).max(1e-12);
    assert!(
        (actual - expected).abs() <= bound,
        "actual {actual} != expected {expected} (diff {})",
        (actual - expected).abs()
    );
}

/// The molar state reproduces NeqSim's `ArgonSolidHelmholtzEquation.evaluate` from the
/// triple point to the 16 GPa limit.
#[test]
fn properties_reproduce_neqsim() {
    // (T in K, p in Pa, v, a, u, s, h, g, cv, cp, ln_phi).
    let cases: [(f64, f64, f64, f64, f64, f64, f64, f64, f64, f64, f64); 4] = [
        (
            70.0,
            1.0e6,
            2.395_456_902_490_360e-5,
            -9.007_458_032_090_453e3,
            -6.703_430_469_514_010e3,
            3.291_467_946_537_777e1,
            -6.679_475_900_489_107e3,
            -8.983_503_463_065_550e3,
            2.291_388_515_385_048e1,
            3.028_612_487_583_490e1,
            -1.773_783_059_830_211e1,
        ),
        (
            20.0,
            1.0e3,
            2.263_457_708_484_870e-5,
            -7.980_305_903_625_450e3,
            -7.856_979_648_119_364e3,
            6.166_312_775_304_262e0,
            -7.856_957_013_542_279e3,
            -7.980_283_269_048_365e3,
            1.199_711_265_904_549e1,
            1.234_705_268_500_559e1,
            -4.338_520_593_857_422e1,
        ),
        (
            300.0,
            1.6e10,
            1.220_095_641_718_160e-5,
            2.411_509_891_692_012e4,
            3.306_667_782_235_993e4,
            2.983_859_635_146_606e1,
            2.282_819_804_972_656e5,
            2.193_304_015_918_258e5,
            2.414_247_552_045_469e1,
            2.485_501_023_110_843e1,
            7.594_835_615_867_119e1,
        ),
        (
            83.8058,
            6.8891e4,
            2.460_993_959_848_729e-5,
            -9.503_072_959_250_616e3,
            -6.250_543_211_454_313e3,
            3.881_031_799_465_315e1,
            -6.248_847_808_105_434e3,
            -9.501_377_555_901_736e3,
            2.340_016_227_517_340e1,
            3.532_737_724_838_895e1,
            -1.326_308_374_591_701e1,
        ),
    ];
    for (t, p, v, a, u, s, h, g, cv, cp, ln_phi) in cases {
        let state = argon_solid::properties(t, p);
        assert_close(state.v, v, 1e-8);
        assert_close(state.a, a, 1e-8);
        assert_close(state.u, u, 1e-8);
        assert_close(state.s, s, 1e-8);
        assert_close(state.h, h, 1e-8);
        assert_close(state.g, g, 1e-8);
        assert_close(state.cv, cv, 1e-8);
        assert_close(state.cp, cp, 1e-8);
        assert_close(state.ln_phi, ln_phi, 1e-8);
    }
}

/// The normalized third-order Debye function reproduces NeqSim's values.
#[test]
fn debye_function_reproduces_neqsim() {
    assert_close(dual::debye_function3(1.0), 0.224_805_188_025_938_2, 1e-12);
    assert_close(dual::debye_function3(1.0e-8), 1.0 / 3.0, 1e-8);
}

/// The state satisfies the Gibbs, internal-energy and enthalpy identities, and the
/// `cp - cv` relation holds against finite differences of the pressure.
#[test]
fn properties_satisfy_the_thermodynamic_identities() {
    for (t, p) in [(70.0, 1.0e6), (20.0, 1.0e3), (300.0, 1.6e10)] {
        let state = argon_solid::properties(t, p);
        assert_close(state.g, state.h - t * state.s, 1e-9);
        assert_close(state.u, state.h - p * state.v, 1e-9);
        assert_close(state.u, state.a + t * state.s, 1e-9);
        assert_close(state.g, state.a + p * state.v, 1e-9);

        // cp - cv = T (dP/dT)_V^2 / -(dP/dV)_T, by central differences of P(T, V).
        let dt = t * 1.0e-5;
        let dv = state.v * 1.0e-5;
        let dp_dt = (argon_solid::pressure(t + dt, state.v)
            - argon_solid::pressure(t - dt, state.v))
            / (2.0 * dt);
        let dp_dv = (argon_solid::pressure(t, state.v + dv)
            - argon_solid::pressure(t, state.v - dv))
            / (2.0 * dv);
        assert_close(state.cp - state.cv, t * dp_dt * dp_dt / -dp_dv, 1e-5);
    }
}
