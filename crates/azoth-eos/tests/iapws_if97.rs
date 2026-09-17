//! The IAPWS-IF97 (water) steam formulation, checked against NeqSim.

// The oracle values are NeqSim's ground truth, so their full digits are the point, and
// the per-state tuple is clearer than a factored type.
#![allow(clippy::excessive_precision, clippy::type_complexity)]

use azoth_eos::iapws_if97::{self, R};

fn assert_close(actual: f64, expected: f64, tol: f64) {
    let bound = (expected.abs() * tol).max(1e-12);
    assert!(
        (actual - expected).abs() <= bound,
        "actual {actual} != expected {expected} (diff {})",
        (actual - expected).abs()
    );
}

/// The specific properties reproduce NeqSim's `Iapws_if97` for two liquid and two vapour
/// states. NeqSim exposes no `cv`, so it is checked separately below.
#[test]
fn properties_reproduce_neqsim() {
    // (p in MPa, T in K, v, h, s, cp, w).
    let cases: [(f64, f64, f64, f64, f64, f64, f64); 4] = [
        (
            1.0,
            300.0,
            1.003_048_947_481_398e-3,
            1.134_923_020_764_579e2,
            3.928_488_884_700_410e-1,
            4.178_572_324_566_614,
            1.504_557_549_675_154e3,
        ),
        (
            5.0,
            500.0,
            1.199_733_158_655_553e-3,
            9.759_873_996_062_025e2,
            2.576_505_158_601_740,
            4.638_420_445_058_126,
            1.249_713_885_375_161e3,
        ),
        (
            0.1,
            400.0,
            1.826_205_557_277_711,
            2.730_397_845_967_862e3,
            7.502_400_892_087_548,
            2.008_059_686_205_170,
            4.903_066_732_262_722e2,
        ),
        (
            10.0,
            600.0,
            2.009_298_821_831_874e-2,
            2.819_826_634_144_003e3,
            5.775_376_080_918_124,
            5.140_697_425_375_575,
            5.033_474_443_380_647e2,
        ),
    ];
    for (p, t, v, h, s, cp, w) in cases {
        let props = iapws_if97::properties(p, t);
        assert_close(props.v, v, 1e-9);
        assert_close(props.h, h, 1e-9);
        assert_close(props.s, s, 1e-9);
        assert_close(props.cp, cp, 1e-9);
        assert_close(props.w, w, 1e-9);
    }
}

/// The Region 4 saturation equations reproduce NeqSim's `psat_t` and `tsat_p`.
#[test]
fn saturation_equations_reproduce_neqsim() {
    assert_close(
        iapws_if97::saturation_pressure(300.0),
        0.003_536_589_413,
        1e-9,
    );
    assert_close(
        iapws_if97::saturation_temperature(1.0),
        453.035_632_391_467,
        1e-9,
    );
}

/// The property set satisfies the exact Gibbs, internal-energy and compressibility
/// identities, and `cv` matches the finite-difference relation `cv = cp - T v a_p^2 / k_T`.
#[test]
fn properties_satisfy_the_thermodynamic_identities() {
    for (p, t) in [(1.0, 300.0), (5.0, 500.0), (0.1, 400.0), (10.0, 600.0)] {
        let props = iapws_if97::properties(p, t);

        assert_close(props.g, props.h - t * props.s, 1e-12);
        assert_close(props.u, props.h - 1000.0 * p * props.v, 1e-9);
        assert_close(props.z, 1000.0 * p * props.v / (R * t), 1e-9);

        // cv = cp - T v (dv/dT)_p^2 / (dv/dp)_T, by central differences of v(p, T) in SI.
        let p_pa = p * 1e6;
        let dt = 0.05;
        let dp = 1e-4 * p_pa;
        let v = props.v;
        let dv_dt = (iapws_if97::properties(p, t + dt).v - iapws_if97::properties(p, t - dt).v)
            / (2.0 * dt);
        let alpha = dv_dt / v;
        let dv_dp = (iapws_if97::properties((p_pa - dp) / 1e6, t).v
            - iapws_if97::properties((p_pa + dp) / 1e6, t).v)
            / (-2.0 * dp);
        let kappa = -dv_dp / v;
        let cv_fd = (props.cp * 1000.0 - t * v * alpha * alpha / kappa) / 1000.0;
        assert_close(props.cv, cv_fd, 1e-4);
    }
}
