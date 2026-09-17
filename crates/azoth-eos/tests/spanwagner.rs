//! The Span-Wagner CO2 reference equation of state, checked against NeqSim.

// The oracle values are NeqSim's ground truth, so their full digits are the point, and
// the per-state tuple is clearer than a factored type for four rows.
#![allow(clippy::excessive_precision, clippy::type_complexity)]

use azoth_eos::spanwagner::{self, R, RHO_CRIT, T_CRIT};

fn assert_close(actual: f64, expected: f64) {
    let tol = (expected.abs() * 1e-9).max(1e-12);
    assert!(
        (actual - expected).abs() <= tol,
        "actual {actual} != expected {expected} (diff {})",
        (actual - expected).abs()
    );
}

/// The molar properties reproduce NeqSim's `getProperties`.
#[test]
fn properties_reproduce_neqsim() {
    // (T, P, rho, Z, h, s, cp, cv, u, g, w, phi) from NeqSim's NeqSimSpanWagner.
    let cases: [(f64, f64, f64, f64, f64, f64, f64, f64, f64, f64, f64, f64); 4] = [
        (
            300.0,
            1e6,
            422.16459308,
            0.94964279824,
            21953.758222,
            100.75454630,
            40.527617019,
            30.021598422,
            19585.013859,
            -8272.6056698,
            262.43079659,
            0.95173812800,
        ),
        (
            250.0,
            1e5,
            48.541636158,
            0.99108047578,
            20508.579556,
            114.14630660,
            35.428373006,
            26.766071175,
            18448.492425,
            -8027.9970925,
            247.79417196,
            0.99115342043,
        ),
        (
            350.0,
            5e6,
            2036.3488287,
            0.84374871586,
            22656.859248,
            90.848281624,
            52.264553620,
            34.213615841,
            20201.484250,
            -9140.0393205,
            266.44681196,
            0.86040884202,
        ),
        (
            300.0,
            7.3e6,
            16462.678094,
            0.17777245059,
            12103.768623,
            54.803064734,
            219.29282355,
            45.119130069,
            11660.341378,
            -4337.150797212196,
            302.72588156,
            0.63154003367,
        ),
    ];
    for (t, p, rho_e, z, h, s, cp, cv, u, g, w, phi) in cases {
        let rho = spanwagner::solve_density(t, p, false);
        assert_close(rho, rho_e);
        let props = spanwagner::properties(t, rho);
        assert_close(props.z, z);
        assert_close(props.h, h);
        assert_close(props.s, s);
        assert_close(props.cp, cp);
        assert_close(props.cv, cv);
        assert_close(props.u, u);
        assert_close(props.g, g);
        assert_close(props.sound, w);
        assert_close(props.phi, phi);
    }
}

/// The property set satisfies the Gibbs identity and the `cp - cv` relation.
#[test]
fn properties_satisfy_the_thermodynamic_identities() {
    for t in [230.0, 280.0, 320.0, 400.0] {
        let rho = 0.5 * RHO_CRIT;
        let p = spanwagner::properties(t, rho);
        assert_close(p.g, p.h - t * p.s);

        let delta = rho / RHO_CRIT;
        let tau = T_CRIT / t;
        let id = spanwagner::ideal(delta, tau);
        let res = spanwagner::residual(delta, tau);
        let numer = 1.0 + delta * res.alpha_delta - delta * tau * res.alpha_delta_tau;
        let denom = 1.0 + 2.0 * delta * res.alpha_delta + delta * delta * res.alpha_delta_delta;
        assert_close(p.cp - p.cv, R * numer * numer / denom);
        let _ = id;
    }
}
