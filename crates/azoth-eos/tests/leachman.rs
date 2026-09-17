//! The Leachman (hydrogen) reference equation of state, checked against NeqSim.

// The oracle values are NeqSim's ground truth, so their full digits are the point, and
// the per-state tuple is clearer than a factored type.
#![allow(clippy::excessive_precision, clippy::type_complexity)]

use azoth_eos::leachman::{self, HydrogenType, R};

fn assert_close(actual: f64, expected: f64) {
    let tol = (expected.abs() * 1e-9).max(1e-12);
    assert!(
        (actual - expected).abs() <= tol,
        "actual {actual} != expected {expected} (diff {})",
        (actual - expected).abs()
    );
}

/// The molar properties reproduce NeqSim's `propertiesLeachman` for all three hydrogen types.
#[test]
fn properties_reproduce_neqsim() {
    // (type, T, P in Pa, rho in mol/m³, Z, u, h, s, cv, cp, g).
    let cases: [(
        HydrogenType,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
    ); 6] = [
        (
            HydrogenType::Normal,
            300.0,
            1e5,
            40.06713592000857,
            1.000584534646334,
            5483.632550699886,
            7979.443586448574,
            107.8880264448619,
            20.53458433576209,
            28.85299070530393,
            -24386.96434700999,
        ),
        (
            HydrogenType::Normal,
            100.0,
            1e6,
            1204.673127522504,
            0.9983759652275699,
            1773.315295196198,
            2603.415989860626,
            59.74983629167724,
            14.32576936590632,
            23.20854934562025,
            -3371.567639307098,
        ),
        (
            HydrogenType::Para,
            300.0,
            1e5,
            40.06714615284991,
            1.000584279104730,
            6486.525588313916,
            8982.335986651638,
            114.5962420804385,
            21.60828133155501,
            29.92679427636361,
            -25396.53663747990,
        ),
        (
            HydrogenType::Para,
            100.0,
            1e6,
            1205.075901321336,
            0.9980422769679870,
            1860.386500056668,
            2690.209749283979,
            60.86600528152406,
            18.76284030968661,
            27.65274576878174,
            -3396.390778868428,
        ),
        (
            HydrogenType::Ortho,
            300.0,
            1e5,
            40.06725460425209,
            1.000581570789196,
            5147.766152414441,
            7643.569795257185,
            105.6755296171790,
            20.17836590088738,
            28.49650418378744,
            -24059.08908989651,
        ),
        (
            HydrogenType::Ortho,
            100.0,
            1e6,
            1203.653857425213,
            0.9992214032751717,
            1744.259968773596,
            2575.063603748141,
            59.41295206256161,
            12.84476109421034,
            21.70740983796728,
            -3366.231602508020,
        ),
    ];
    for (ht, t, p, rho_e, z, u, h, s, cv, cp, g) in cases {
        assert_close(leachman::solve_density(t, p, ht), rho_e);
        let props = leachman::properties(t, rho_e, ht);
        assert_close(props.z, z);
        assert_close(props.u, u);
        assert_close(props.h, h);
        assert_close(props.s, s);
        assert_close(props.cv, cv);
        assert_close(props.cp, cp);
        assert_close(props.g, g);
    }
}

/// The model assembles the density solve and the property set into the phase state the
/// spec's worked examples name.
#[test]
fn hydrogen_phase_model_matches_the_worked_example() {
    let r = azoth_eos::hydrogen_phase(
        azoth_core::units::kelvins(300.0),
        azoth_core::units::pascals(100_000.0),
        "normal",
    )
    .unwrap();
    assert_close(r.z_factor, 1.000584534646334);
    assert_close(r.u.value, 5483.632550699886);
    assert_close(r.h.value, 7979.443586448574);
    assert_close(r.s.value, 107.8880264448619);
    assert_close(r.cv.value, 20.53458433576209);
    assert_close(r.cp.value, 28.85299070530393);
    assert_close(r.g.value, -24386.96434700999);
}

/// The property set satisfies the Gibbs identity and the `cp - cv` relation for all types.
#[test]
fn properties_satisfy_the_thermodynamic_identities() {
    for ht in [
        HydrogenType::Normal,
        HydrogenType::Para,
        HydrogenType::Ortho,
    ] {
        for t in [25.0, 100.0, 300.0] {
            let rho = 0.5 * ht.rhoc();
            let p = leachman::properties(t, rho, ht);
            assert_close(p.g, p.h - t * p.s);

            let delta = rho / ht.rhoc();
            let tau = ht.tc() / t;
            let id = leachman::ideal(delta, tau, ht);
            let res = leachman::residual(delta, tau, ht);
            let numer = 1.0 + delta * res.alpha_delta - delta * tau * res.alpha_delta_tau;
            let denom = 1.0 + 2.0 * delta * res.alpha_delta + delta * delta * res.alpha_delta_delta;
            assert_close(p.cp - p.cv, R * numer * numer / denom);
            let _ = id;
        }
    }
}
