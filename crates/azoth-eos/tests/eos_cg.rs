//! The EOS-CG multi-fluid Helmholtz equation of state, checked against NeqSim.

// The oracle values are NeqSim's ground truth, so their full digits are the point.
#![allow(
    clippy::excessive_precision,
    clippy::type_complexity,
    clippy::inconsistent_digit_grouping
)]

use azoth_eos::eos_cg::{self, R};

fn assert_close(actual: f64, expected: f64, rel_tol: f64) {
    let bound = (expected.abs() * rel_tol).max(1e-12);
    assert!(
        (actual - expected).abs() <= bound,
        "actual {actual} != expected {expected} (diff {})",
        (actual - expected).abs()
    );
}

/// Build a 1-indexed composition from `(index, fraction)` pairs.
fn composition(pairs: &[(usize, f64)]) -> [f64; 29] {
    let mut x = [0.0; 29];
    for &(i, f) in pairs {
        x[i] = f;
    }
    x
}

/// The molar properties reproduce NeqSim's `propertiesEOSCG` at fixed densities, for two
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
            943.519_829_620_654_6,
            0.951_529_028_509_804_1,
            -2759.639_490_361_096,
            -400.839_916_309_459_2,
            -19.500_150_403_074_36,
            29.897_561_338_576_97,
            40.314_199_844_651_63,
            5413.129_926_367_160,
        ),
        (
            298.15,
            4.0,
            &[(3, 1.0)],
            5810.667_654_664_715,
            0.585_999_231_267_817_6,
            -5160.133_578_627_326,
            -3707.466_664_961_147,
            -43.332_760_502_078_91,
            42.409_309_014_633_32,
            128.814_190_358_365_1,
            9212.195_878_733_677,
        ),
        (
            298.15,
            0.4,
            &[(1, 1.0)],
            974.984_106_654_620_9,
            0.983_260_394_421_774_8,
            -2591.048_546_082_660,
            -153.588_279_446_107_5,
            -19.200_210_393_087_98,
            27.555_087_857_303_41,
            36.646_625_125_588_19,
            5570.954_449_253_073,
        ),
        (
            250.0,
            3.0,
            &[(1, 0.8), (4, 0.1), (3, 0.05), (2, 0.05)],
            4953.578_704_696_020,
            0.794_371_435_619_022_9,
            -4863.704_686_057_494,
            -3212.511_784_492_154,
            -37.058_177_300_818_80,
            29.584_215_023_987_57,
            51.256_363_406_401_12,
            6052.032_540_712_545,
        ),
    ];
    for (t, d, pairs, p, z, u, h, s, cv, cp, g) in cases {
        let x = composition(pairs);
        let props = eos_cg::properties(t, d, &x);
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
        (298.15, 1000.0, &[(3, 1.0)][..], 0.948_517_262_506_362_8),
        (298.15, 100.0, &[(3, 1.0)][..], 0.995_015_119_152_256_5),
        (
            250.0,
            5000.0,
            &[(1, 0.8), (4, 0.1), (3, 0.05), (2, 0.05)][..],
            0.792_309_282_412_182_0,
        ),
    ] {
        let x = composition(pairs);
        let d_oracle = p / (R * t) / z;
        assert_close(eos_cg::solve_density(t, p, &x), d_oracle, 1e-6);
    }
}

/// The property set satisfies the Gibbs identity and the molar-mass sum.
#[test]
fn properties_satisfy_the_thermodynamic_identities() {
    let x = composition(&[(1, 0.8), (4, 0.1), (3, 0.05), (2, 0.05)]);
    let t = 250.0;
    let p = 5000.0;
    let d = eos_cg::solve_density(t, p, &x);
    let props = eos_cg::properties(t, d, &x);
    assert_close(props.g, props.h - t * props.s, 1e-9);
    assert_close(props.h, props.u + p / d, 1e-9);
}
