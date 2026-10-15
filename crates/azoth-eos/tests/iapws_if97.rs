//! The IAPWS-IF97 (water) steam formulation, checked against NeqSim and against the standard.
//!
//! **Two oracles, and they answer different questions.** `properties_reproduce_neqsim` reads
//! `validation/neqsim/captures/iapws_if97_probe.tsv`, written by
//! `validation/neqsim/IapwsIf97Probe.java` against the pinned jar, and says *this port is what
//! NeqSim computes*. `the_release_check_values_hold` reads the IAPWS release's own verification
//! tables and says *what NeqSim computes is the standard* - so a pin bump that moved NeqSim away
//! from R7-97 would fail here rather than pass everywhere.
//!
//! The capture became the first oracle on 2026-10-15; before that the digits were NeqSim's own
//! shipped `Iapws_if97Test` literals copied by hand, which an earlier reader could not re-run and
//! so could not say whether a pin move had changed.

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

/// The specific properties reproduce NeqSim's `Iapws_if97` for three liquid and three vapour
/// states. NeqSim exposes no `cv`, so it is checked separately below.
#[test]
fn properties_reproduce_neqsim() {
    // (p in MPa, T in K, v, h, s, cp, w), from `captures/iapws_if97_probe.tsv`. The last row is
    // the state NeqSim's own shipped test states, kept so a reader can tie this to `h = 3479.00`.
    let cases: [(f64, f64, f64, f64, f64, f64, f64); 6] = [
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
        (
            3.0,
            500.0,
            1.202_418_003_378_336e-3,
            9.755_422_390_972_245e2,
            2.580_419_120_051_806,
            4.655_806_822_111_205,
            1.240_713_373_101_727e3,
        ),
        (
            1.0,
            773.15,
            3.541_102_656_247_171e-1,
            3.479_003_674_758_910e3,
            7.763_965_038_557_445,
            2.168_203_745_568_266,
            6.723_448_068_553_923e2,
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

/// **The release's own verification values, which are not NeqSim's.**
///
/// IAPWS R7-97 states check values for Region 1 and Region 2 to nine significant digits, and they
/// are what any implementation of the formulation is judged by. They are written here rather than
/// read from a capture on purpose: the test above says *this port agrees with NeqSim*, and this
/// one says *NeqSim agrees with the standard* - so a pin bump that moved the class away from
/// R7-97 fails here instead of passing everywhere.
///
/// Measured against the pinned jar, NeqSim reproduces all thirty of these numbers to `3.1e-9`
/// relative, which is the release's own nine-digit rounding; the bound below is a decade above it.
#[test]
fn the_release_check_values_hold() {
    // (p in MPa, T in K, v, h, s, cp, w), R7-97 Table 5 and Table 15.
    let region_1: [(f64, f64, f64, f64, f64, f64, f64); 3] = [
        (
            3.0,
            300.0,
            0.100_215_168e-2,
            0.115_331_273e3,
            0.392_294_792,
            0.417_301_218e1,
            0.150_773_921e4,
        ),
        (
            80.0,
            300.0,
            0.971_180_894e-3,
            0.184_142_828e3,
            0.368_563_852,
            0.401_008_987e1,
            0.163_469_054e4,
        ),
        (
            3.0,
            500.0,
            0.120_241_800e-2,
            0.975_542_239e3,
            0.258_041_912e1,
            0.465_580_682e1,
            0.124_071_337e4,
        ),
    ];
    let region_2: [(f64, f64, f64, f64, f64, f64, f64); 3] = [
        (
            0.0035,
            300.0,
            0.394_913_866e2,
            0.254_991_145e4,
            0.852_238_967e1,
            0.191_300_162e1,
            0.427_920_172e3,
        ),
        (
            0.0035,
            700.0,
            0.923_015_898e2,
            0.333_568_375e4,
            0.101_749_996e2,
            0.208_141_274e1,
            0.644_289_068e3,
        ),
        (
            30.0,
            700.0,
            0.542_946_619e-2,
            0.263_149_474e4,
            0.517_540_298e1,
            0.103_505_092e2,
            0.480_386_523e3,
        ),
    ];

    for (name, states) in [("Region 1", region_1), ("Region 2", region_2)] {
        for (p, t, v, h, s, cp, w) in states {
            let context = &format!("{name} at p = {p} MPa, T = {t} K");
            // The region the release states this state is in is the one the class has to select,
            // which is a separate reading of the formulation from the properties themselves.
            assert_eq!(
                t <= iapws_if97::saturation_temperature(p),
                name == "Region 1",
                "{context}: the class selected the other region"
            );
            let props = iapws_if97::properties(p, t);
            for (field, actual, expected) in [
                ("v", props.v, v),
                ("h", props.h, h),
                ("s", props.s, s),
                ("cp", props.cp, cp),
                ("w", props.w, w),
            ] {
                let relative = (actual / expected - 1.0).abs();
                assert!(
                    relative <= 1.0e-8,
                    "{context}: {field} = {actual} against the release's {expected} \
                     ({relative:e} from it)"
                );
            }
        }
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

/// The model converts the mass-based property set to the molar phase state the spec's
/// worked examples name.
#[test]
fn water_phase_model_matches_the_worked_example() {
    let r = azoth_eos::water_phase(
        azoth_core::units::kelvins(300.0),
        azoth_core::units::pascals(1_000_000.0),
    )
    .unwrap();
    assert_close(r.z_factor, 0.007_244_437_998_302_719, 1e-9);
    assert_close(r.u.value, 2_026.525_392_109_388, 1e-9);
    assert_close(r.h.value, 2_044.595_599_751_971, 1e-9);
    assert_close(r.s.value, 7.077_282_723_476_562, 1e-9);
    assert_close(r.cv.value, 74.362_176_806_013_20, 1e-9);
    assert_close(r.cp.value, 75.278_150_427_318_45, 1e-9);
    assert_close(r.g.value, -78.589_217_290_997_03, 1e-9);
}
