//! The GERG-2008 multi-fluid Helmholtz equation of state, checked against NeqSim.

// The oracle values are NeqSim's ground truth, so their full digits are the point.
#![allow(
    clippy::excessive_precision,
    clippy::type_complexity,
    clippy::inconsistent_digit_grouping
)]

use azoth_eos::gerg2008::{self, R};

fn assert_close(actual: f64, expected: f64, rel_tol: f64) {
    let bound = (expected.abs() * rel_tol).max(1e-12);
    assert!(
        (actual - expected).abs() <= bound,
        "actual {actual} != expected {expected} (diff {})",
        (actual - expected).abs()
    );
}

/// Build a 1-indexed composition from `(index, fraction)` pairs.
fn composition(pairs: &[(usize, f64)]) -> [f64; 22] {
    let mut x = [0.0; 22];
    for &(i, f) in pairs {
        x[i] = f;
    }
    x
}

/// The molar properties reproduce NeqSim's `propertiesGERG` at fixed densities, for two
/// pure fluids and a natural-gas mixture. Fixed density rather than solved: the
/// comparison is then bit-for-bit rather than bounded by the density-root precision.
#[test]
fn properties_reproduce_neqsim() {
    // (T, D in mol/L, composition, P in kPa, Z, u, h, s, cv, cp, g).
    let cases: [(
        f64,
        f64,
        &[(usize, f64)],
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
    ); 4] = [
        (
            298.15,
            0.4,
            &[(3, 1.0)],
            943.5208942839346,
            0.9515290285098041,
            -2759.642604325134,
            -400.8403686152976,
            -19.50018178895629,
            29.89772681600762,
            40.31437707616752,
            5413.138831762020,
        ),
        (
            298.15,
            4.0,
            &[(3, 1.0)],
            5810.674211394409,
            0.5859992312678176,
            -5160.139401297455,
            -3707.470848448853,
            -43.33281878056451,
            42.40948861026042,
            128.8144674528486,
            9212.209070976454,
        ),
        (
            298.15,
            0.4,
            &[(1, 1.0)],
            974.9852068220410,
            0.9832603944217748,
            -2591.051469809317,
            -153.5884527542141,
            -19.20024144051945,
            27.55524414639745,
            36.64679167352956,
            5570.963532736659,
        ),
        (
            250.0,
            3.0,
            &[(1, 0.8), (4, 0.1), (3, 0.05), (2, 0.05)],
            4953.584314999493,
            0.7943714389399581,
            -4863.716278014667,
            -3212.521506348171,
            -37.05825083337178,
            29.58437176674896,
            51.25654419196839,
            6052.041201994774,
        ),
    ];
    for (t, d, pairs, p, z, u, h, s, cv, cp, g) in cases {
        let x = composition(pairs);
        let props = gerg2008::properties(t, d, &x);
        assert_close(props.z, z, 1e-12);
        assert_close(props.u, u, 1e-12);
        assert_close(props.h, h, 1e-12);
        assert_close(props.s, s, 1e-12);
        assert_close(props.cv, cv, 1e-12);
        assert_close(props.cp, cp, 1e-12);
        assert_close(props.g, g, 1e-12);
        assert_close(props.pressure_kpa, p, 1e-12);
    }
}

/// The density solve inverts the pressure to the molar density NeqSim reports.
#[test]
fn density_solve_reproduces_neqsim() {
    for (t, p, pairs, z) in [
        (298.15, 1000.0, &[(3, 1.0)][..], 0.9485173228128433),
        (298.15, 100.0, &[(3, 1.0)][..], 0.9950151247616322),
        (
            250.0,
            5000.0,
            &[(1, 0.8), (4, 0.1), (3, 0.05), (2, 0.05)][..],
            0.7923095385298691,
        ),
    ] {
        let x = composition(pairs);
        let d_oracle = p / (R * t) / z;
        assert_close(gerg2008::solve_density(t, p, &x), d_oracle, 1e-6);
    }
}

/// The property set satisfies the Gibbs identity and the molar-mass sum.
#[test]
fn properties_satisfy_the_thermodynamic_identities() {
    let x = composition(&[(1, 0.8), (4, 0.1), (3, 0.05), (2, 0.05)]);
    let t = 250.0;
    let p = 5000.0;
    let d = gerg2008::solve_density(t, p, &x);
    let props = gerg2008::properties(t, d, &x);
    assert_close(props.g, props.h - t * props.s, 1e-9);
    assert_close(props.h, props.u + p / d, 1e-9);
}
