//! The solid para-hydrogen Helmholtz equation, checked against NeqSim.

// The oracle values are NeqSim's ground truth, so their full digits are the point, and
// the per-state tuple is clearer than a factored type.
#![allow(clippy::excessive_precision, clippy::type_complexity)]

use azoth_eos::parahydrogen_solid;

fn assert_close(actual: f64, expected: f64, rel_tol: f64) {
    let bound = (expected.abs() * rel_tol).max(1e-12);
    assert!(
        (actual - expected).abs() <= bound,
        "actual {actual} != expected {expected} (diff {})",
        (actual - expected).abs()
    );
}

/// The molar state reproduces NeqSim's `ParaHydrogenSolidHelmholtzEquation.evaluate` from
/// the triple point to the 1 GPa limit.
#[test]
fn properties_reproduce_neqsim() {
    // (T in K, p in Pa, v, a, u, s, h, g, cv, cp, ln_phi).
    let cases: [(f64, f64, f64, f64, f64, f64, f64, f64, f64, f64, f64); 4] = [
        (
            4.2,
            7042.0,
            2.313_809_319_665_469e-5,
            -1.452_957_411_925_972e1,
            -1.440_767_785_183_070e1,
            2.902_292_081_643_306e-2,
            -1.424_473_939_953_986e1,
            -1.436_663_566_696_888e1,
            8.911_452_118_340_994e-2,
            8.917_908_176_775_872e-2,
            2.241_871_009_686_271e0,
        ),
        (
            20.0,
            1.0e8,
            1.821_137_965_908_958e-5,
            1.606_947_184_910_250e2,
            1.781_288_404_282_875e2,
            8.717_060_968_631_262e-1,
            1.999_266_806_337_245e3,
            1.981_832_684_399_983e3,
            2.760_332_287_046_142e0,
            2.802_738_495_196_803e0,
            5.010_229_676_412_527e0,
        ),
        (
            80.0,
            1.0e9,
            1.212_191_771_436_695e-5,
            2.343_350_943_535_304e3,
            2.898_023_914_771_470e3,
            6.933_412_140_452_078e0,
            1.501_994_162_913_842e4,
            1.446_526_865_790_226e4,
            1.526_010_893_661_340e1,
            1.571_316_527_176_994e1,
            1.253_681_008_094_105e1,
        ),
        (
            13.8033,
            7042.0,
            2.339_884_191_704_222e-5,
            -1.956_387_135_714_953e1,
            2.571_932_404_200_801e0,
            1.603_660_266_845_634e0,
            2.736_707_048_980_613e0,
            -1.939_909_671_236_972e1,
            4.662_105_736_423_913e0,
            5.293_679_708_064_397e0,
            2.484_247_719_393_687e0,
        ),
    ];
    for (t, p, v, a, u, s, h, g, cv, cp, ln_phi) in cases {
        let state = parahydrogen_solid::properties(t, p);
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

/// The state satisfies the Gibbs, internal-energy and enthalpy identities, and the
/// `cp - cv` relation holds against finite differences of the pressure.
#[test]
fn properties_satisfy_the_thermodynamic_identities() {
    for (t, p) in [(20.0, 1.0e8), (80.0, 1.0e9)] {
        let state = parahydrogen_solid::properties(t, p);
        assert_close(state.g, state.h - t * state.s, 1e-9);
        assert_close(state.u, state.h - p * state.v, 1e-9);
        assert_close(state.u, state.a + t * state.s, 1e-9);
        assert_close(state.g, state.a + p * state.v, 1e-9);

        let dt = t * 1.0e-5;
        let dv = state.v * 1.0e-5;
        let dp_dt = (parahydrogen_solid::pressure(t + dt, state.v)
            - parahydrogen_solid::pressure(t - dt, state.v))
            / (2.0 * dt);
        let dp_dv = (parahydrogen_solid::pressure(t, state.v + dv)
            - parahydrogen_solid::pressure(t, state.v - dv))
            / (2.0 * dv);
        assert_close(state.cp - state.cv, t * dp_dt * dp_dt / -dp_dv, 1e-5);
    }
}
