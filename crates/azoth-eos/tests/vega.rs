//! The Vega (NIST helium) reference equation of state, checked against NeqSim.

// The oracle values are NeqSim's ground truth, so their full digits are the point, and
// the per-state tuple is clearer than a factored type for four rows.
#![allow(clippy::excessive_precision, clippy::type_complexity)]

use azoth_eos::vega::{self, R, RHO_CRIT, T_CRIT};

fn assert_close(actual: f64, expected: f64) {
    let tol = (expected.abs() * 1e-9).max(1e-12);
    assert!(
        (actual - expected).abs() <= tol,
        "actual {actual} != expected {expected} (diff {})",
        (actual - expected).abs()
    );
}

/// The molar properties reproduce NeqSim's `propertiesVega`.
#[test]
fn properties_reproduce_neqsim() {
    // (T, P in Pa, rho in mol/m³, Z, u, h, s, cv, cp, g).
    let cases: [(f64, f64, f64, f64, f64, f64, f64, f64, f64, f64); 4] = [
        (
            273.15,
            1e7,
            4185.601151595783,
            1.051977081289438,
            3438.131816891616,
            5827.274890446882,
            71.90291916965856,
            12.57753705703601,
            20.81172148386873,
            -13813.00748074535,
        ),
        (
            300.0,
            1e5,
            40.07172609531133,
            1.000474491201717,
            3761.829372397934,
            6257.354515541211,
            112.1030803176974,
            12.47268631699359,
            20.78632637100741,
            -27373.56957976801,
        ),
        (
            150.0,
            5e5,
            398.9611042993242,
            1.004878404053713,
            1890.709478499191,
            3143.964481585582,
            84.31155939301351,
            12.48167815020072,
            20.79925708136838,
            -9502.769427366446,
        ),
        (
            400.0,
            2e7,
            5634.294345182381,
            1.067322826621424,
            5042.388106074379,
            8592.078408836254,
            74.10811629278300,
            12.60878442140509,
            20.76592409159910,
            -21051.16810827695,
        ),
    ];
    for (t, p, rho_e, z, u, h, s, cv, cp, g) in cases {
        // Density solve: Newton from ideal gas finds the same root NeqSim's does.
        assert_close(vega::solve_density(t, p), rho_e);
        // Properties at the oracle density: a tight check of the Helmholtz and property set.
        let props = vega::properties(t, rho_e);
        assert_close(props.z, z);
        assert_close(props.u, u);
        assert_close(props.h, h);
        assert_close(props.s, s);
        // NeqSim's derivative code has two defects: its Gaussian tau-derivative drops the
        // `(de/dtau)^2` term (~5e-9 in `Cv`) and its exponential delta-derivative flips the
        // sign of `d^2 e/ddelta^2` (~2e-5 in `Cp`). The closed form here is checked against
        // a finite difference in the derivative test below.
        assert!(
            (props.cv - cv).abs() < 1e-4 * cv.abs(),
            "cv {} far from NeqSim {}",
            props.cv,
            cv
        );
        assert!(
            (props.cp - cp).abs() < 1e-4 * cp.abs(),
            "cp {} far from NeqSim {}",
            props.cp,
            cp
        );
        assert_close(props.g, g);
    }
}

/// The residual derivatives match a finite difference of the Helmholtz energy, the
/// independent check that NeqSim's dropped `(de/dtau)^2` term fails.
#[test]
fn residual_derivatives_match_finite_difference() {
    let t = 273.15;
    let rho = 4185.601151595783;
    let delta = rho / RHO_CRIT;
    let tau = T_CRIT / t;
    let res = vega::residual(delta, tau);

    let h = 1e-5;
    let a = |d: f64, tt: f64| vega::residual(d, tt).alpha;
    // Second tau derivative by central difference. Helium's Tc is 5.2 K, so `tau` at room
    // temperature is ~0.019 and the fourth derivative of the `tau^t` terms is large; the
    // tolerance is sized to the finite difference's truncation, not to a defect.
    let alpha_tt_fd = (a(delta, tau + h) - 2.0 * a(delta, tau) + a(delta, tau - h)) / (h * h);
    assert!(
        (res.alpha_tau_tau - alpha_tt_fd).abs() < 1e-4,
        "alpha_tau_tau"
    );
    // Second delta derivative by central difference.
    let alpha_dd_fd = (a(delta + h, tau) - 2.0 * a(delta, tau) + a(delta - h, tau)) / (h * h);
    assert!(
        (res.alpha_delta_delta - alpha_dd_fd).abs() < 1e-4,
        "alpha_delta_delta"
    );
    // Cross derivative by a mixed difference.
    let cross_fd = (a(delta + h, tau + h) - a(delta + h, tau - h) - a(delta - h, tau + h)
        + a(delta - h, tau - h))
        / (4.0 * h * h);
    assert!(
        (res.alpha_delta_tau - cross_fd).abs() < 1e-4,
        "alpha_delta_tau"
    );
}

/// The model assembles the density solve and the property set into the phase state the
/// spec's worked examples name.
#[test]
fn helium_phase_model_matches_the_worked_example() {
    let r = azoth_eos::helium_phase(
        azoth_core::units::kelvins(273.15),
        azoth_core::units::pascals(10_000_000.0),
    )
    .unwrap();
    assert_close(r.z_factor, 1.051977081289438);
    assert_close(r.u.value, 3438.1318168916159);
    assert_close(r.h.value, 5827.2748904468817);
    assert_close(r.s.value, 71.90291916965856);
    assert_close(r.cv.value, 12.577536993905438);
    assert_close(r.cp.value, 20.811314291462409);
    assert_close(r.g.value, -13813.007480745353);
}

/// The property set satisfies the Gibbs identity and the `cp - cv` relation.
#[test]
fn properties_satisfy_the_thermodynamic_identities() {
    for t in [10.0, 100.0, 300.0, 500.0] {
        let rho = 0.5 * RHO_CRIT;
        let p = vega::properties(t, rho);
        assert_close(p.g, p.h - t * p.s);

        let delta = rho / RHO_CRIT;
        let tau = T_CRIT / t;
        let id = vega::ideal(delta, tau);
        let res = vega::residual(delta, tau);
        let numer = 1.0 + delta * res.alpha_delta - delta * tau * res.alpha_delta_tau;
        let denom = 1.0 + 2.0 * delta * res.alpha_delta + delta * delta * res.alpha_delta_delta;
        assert_close(p.cp - p.cv, R * numer * numer / denom);
        let _ = id;
    }
}
