//! The ammonia reference equation of state, checked against NeqSim's `Ammonia2023`.

// The oracle values are NeqSim's ground truth, so their full digits are the point, and
// the per-state tuple is clearer than a factored type for four rows.
#![allow(clippy::excessive_precision, clippy::type_complexity)]

use azoth_eos::ammonia::{self, MOLAR_MASS, R, RHO_CRIT, T_CRIT};

fn assert_close(actual: f64, expected: f64) {
    let tol = (expected.abs() * 1e-10).max(1e-12);
    assert!(
        (actual - expected).abs() <= tol,
        "actual {actual} != expected {expected} (diff {})",
        (actual - expected).abs()
    );
}

/// The reduced variables from NeqSim's reported density (kg/m³) and temperature.
fn reduced(t: f64, rho_kg: f64) -> (f64, f64) {
    let rho = rho_kg / MOLAR_MASS;
    (rho / RHO_CRIT, T_CRIT / t)
}

/// The Helmholtz energy and its derivatives reproduce NeqSim's `getAlpha0`/`getAlphares`.
#[test]
fn helmholtz_reproduces_neqsim() {
    // (T, rho_kg/m³, a0[3], ar[6]) from NeqSim's SystemAmmoniaEos.
    let cases: [(f64, f64, [f64; 3], [f64; 6]); 4] = [
        (
            298.15,
            7.776761102217432,
            [-1.462137932877399, 10.66822021087482, -3.274456863694800],
            [
                -0.1165083493057211,
                -0.1165943917180288,
                -0.0001331828250870903,
                -0.3588857203672393,
                -0.3699886564477126,
                -1.117044420257865,
            ],
        ),
        (
            300.0,
            0.6897588419899237,
            [-3.950541034688663, 10.62264423975695, -3.280594753637151],
            [
                -0.01013636813623184,
                -0.01013676891622481,
                -7.800689433490305e-7,
                -0.03019448186838622,
                -0.03027677318620893,
                -0.08641772177530725,
            ],
        ),
        (
            250.0,
            13.29852827223677,
            [1.076013454732514, 12.10646030495096, -3.133737757052752],
            [
                -0.3685814966865799,
                -0.3839015247188007,
                -0.02995739987445334,
                -1.415807628721409,
                -1.581965260610990,
                -5.417288869795858,
            ],
        ),
        (
            350.0,
            525.5704288290348,
            [1.130151467998181, 9.586286601405450, -3.459162818976216],
            [
                -2.428630312781604,
                -0.8886488369508564,
                3.829219815162751,
                -5.465363413825743,
                -3.916476559512957,
                -2.142567372321456,
            ],
        ),
    ];
    for (t, rho_kg, a0, ar) in cases {
        let (delta, tau) = reduced(t, rho_kg);
        let id = ammonia::ideal(delta, tau);
        let res = ammonia::residual(delta, tau);
        assert_close(id.alpha, a0[0]);
        assert_close(tau * id.alpha_tau, a0[1]);
        assert_close(tau * tau * id.alpha_tau_tau, a0[2]);
        assert_close(res.alpha, ar[0]);
        assert_close(delta * res.alpha_delta, ar[1]);
        assert_close(delta * delta * res.alpha_delta_delta, ar[2]);
        assert_close(tau * res.alpha_tau, ar[3]);
        assert_close(tau * delta * res.alpha_delta_tau, ar[4]);
        assert_close(tau * tau * res.alpha_tau_tau, ar[5]);
    }
}

/// The molar properties reproduce NeqSim's `properties()`.
///
/// `g` is absent: NeqSim's Ammonia2023 subtracts an extra `u` from its Gibbs energy,
/// so the physical `g = h - T s` (checked by the identity test) does not match it.
#[test]
fn properties_reproduce_neqsim() {
    // (T, rho_kg/m³, u, h, s, cv, cp, W, kappa).
    let cases: [(f64, f64, f64, f64, f64, f64, f64, f64, f64); 4] = [
        (
            298.15,
            7.776761102217432,
            25556.39720576612,
            27746.32174840689,
            98.84217173986131,
            36.51297326299638,
            53.55010850788054,
            404.5606452042579,
            1.152251100731124e-6,
        ),
        (
            300.0,
            0.6897588419899237,
            26421.15826398924,
            28890.21251356837,
            121.0014317549826,
            27.99489936167249,
            36.82669984618634,
            434.4681629539527,
            1.010347334126684e-5,
        ),
        (
            250.0,
            13.29852827223677,
            22221.75801004262,
            23502.38994549931,
            83.00511547011064,
            71.09719123576569,
            269.7289352577949,
            306.0157191402957,
            3.046379766892541e-6,
        ),
        (
            350.0,
            525.5704288290348,
            11992.14162839433,
            12316.18040732534,
            45.05941561082694,
            46.57537627252494,
            90.77336517601113,
            1008.149449481734,
            3.648558445961311e-9,
        ),
    ];
    for (t, rho_kg, u, h, s, cv, cp, w, kappa) in cases {
        let rho = rho_kg / MOLAR_MASS;
        let p = ammonia::properties(t, rho);
        assert_close(p.u, u);
        assert_close(p.h, h);
        assert_close(p.s, s);
        assert_close(p.cv, cv);
        assert_close(p.cp, cp);
        assert_close(p.sound, w);
        assert_close(p.kappa, kappa);
    }
}

/// The property set satisfies the Gibbs identity `g = h - T s` and the isobaric relation.
#[test]
fn properties_satisfy_the_thermodynamic_identities() {
    for t in [250.0, 300.0, 350.0, 405.0] {
        let rho = 0.3 * RHO_CRIT;
        let p = ammonia::properties(t, rho);
        assert_close(p.g, p.h - t * p.s);

        // cp - cv = R numer^2 / denom, re-derived from the Helmholtz directly.
        let delta = rho / RHO_CRIT;
        let tau = T_CRIT / t;
        let id = ammonia::ideal(delta, tau);
        let res = ammonia::residual(delta, tau);
        let numer = 1.0 + delta * res.alpha_delta - delta * tau * res.alpha_delta_tau;
        let denom = 1.0 + 2.0 * delta * res.alpha_delta + delta * delta * res.alpha_delta_delta;
        assert_close(p.cp - p.cv, R * numer * numer / denom);
        let _ = id;
    }
}
