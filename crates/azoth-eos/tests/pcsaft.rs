//! PC-SAFT's layers, against NeqSim's own.
//!
//! Every expected number comes from `validation/neqsim/PcsaftProbe.java`, which prints
//! `PhasePCSAFT`'s intermediates with the method each came from. A layer-by-layer
//! comparison is the point: a single Helmholtz energy that agrees can be two errors
//! cancelling, and this tranche has already spent a session on a quantity read at the
//! wrong scale.

use azoth_eos::pcsaft::{PcsaftComponent, d_pressure_over_rt_dv, pressure_over_rt, state};

/// Methane at 300 K and 50 bara, and methane/n-butane at 350 K and 30 bara - both from
/// the probe, whose molar volumes are its `volumeSAFT`.
fn methane() -> PcsaftComponent {
    PcsaftComponent {
        m: 1.0,
        sigma: 3.703900119e-10,
        epsik: 150.029998779,
    }
}

fn n_butane() -> PcsaftComponent {
    PcsaftComponent {
        m: 2.3316,
        sigma: 3.7086e-10,
        epsik: 222.86,
    }
}

/// The layers that are a pure function of the parameters and the state, which reproduce
/// the probe to its printed digits.
///
/// **The packing fraction and everything built on it agree to only `1e-9`**, and the
/// residual is inside NeqSim's own bookkeeping rather than in the formula: `d`, `md3` and
/// `m_bar` reproduce exactly, while the `eta` its `volInit` reports needs an `md3` `1.1e-9`
/// *smaller* than the one its `getDSAFT()` returns on the same phase. So a layer that
/// depends only on the parameters is held tightly and a layer that depends on `eta` is held
/// to the agreement the oracle can offer - which is the same shape as this tranche's other
/// two accessor-versus-value findings, `getMolarVolume` against `getZ`.
fn exact(actual: f64, expected: f64, context: &str) {
    if expected == 0.0 {
        assert!(actual.abs() < 1.0e-15, "{context}: {actual}, expected zero");
        return;
    }
    let relative = (actual / expected - 1.0).abs();
    assert!(
        relative < 1.0e-12,
        "{context}: {actual} against the probe's {expected}, a relative {relative:e}"
    );
}

fn via_eta(actual: f64, expected: f64, context: &str) {
    let relative = (actual / expected - 1.0).abs();
    assert!(
        relative < 1.0e-8,
        "{context}: {actual} against the probe's {expected}, a relative {relative:e}"
    );
}

#[test]
fn pure_methane_reproduces_every_layer() {
    let s = state(&[methane()], &[0.0], &[1.0], 300.0, 4.55032070500990e-4).expect("a state");

    exact(s.d[0], 3.60475564638566e-10, "d");
    exact(s.m_bar, 1.0, "m_bar");
    exact(s.m_minus_1, 0.0, "m_minus_1");
    exact(s.md3, 4.68411438936926e-29, "md3");
    via_eta(s.eta, 0.0324636218320417, "eta");
    via_eta(s.a_hs, 0.135337273046932, "a_hs");
    via_eta(s.g_hs, 1.08615265067751, "g_hs");
    exact(s.s1, 2.54117545218473e-29, "s1");
    exact(s.s2, 1.27084183329500e-29, "s2");
    via_eta(s.i1, 0.933239716943836, "I1");
    via_eta(s.i2, 0.791857336220105, "I2");
    via_eta(s.c1, 0.772825283837257, "C1");
    via_eta(s.f_hc, 0.135337273046932, "F_hc");
    via_eta(s.f_disp1, -0.197232549697721, "F_disp1");
    via_eta(s.f_disp2, -0.0323400359099549, "F_disp2");
    via_eta(s.f(), -0.0942353125607442, "F");
}

#[test]
fn a_binary_reproduces_every_layer() {
    // The `k_ij` is the `KIJPCSAFT` column, 0.022 for this pair - which is what NeqSim's
    // classic mixing rule reads, and the only rule its PC-SAFT can be driven with.
    let s = state(
        &[methane(), n_butane()],
        &[0.0, 0.022, 0.022, 0.0],
        &[0.6, 0.4],
        350.0,
        8.14109734626316e-4,
    )
    .expect("a state");

    exact(s.d[0], 3.58105717328466e-10, "d methane");
    via_eta(s.d[1], 3.64271455064515e-10, "d n-butane");
    exact(s.m_bar, 1.53264, "m_bar");
    exact(s.m_minus_1, 0.53264, "m_minus_1");
    exact(s.md3, 7.26345992024210e-29, "md3");
    via_eta(s.eta, 0.0281366301313126, "eta");
    via_eta(s.a_hs, 0.116643052723215, "a_hs");
    via_eta(s.g_hs, 1.07406652907799, "g_hs");
    exact(s.s1, 6.52038038734557e-29, "s1");
    exact(s.s2, 3.62241907207195e-29, "s2");
    via_eta(s.i1, 0.832366049443170, "I1");
    via_eta(s.i2, 0.582385191477821, "I2");
    via_eta(s.c1, 0.766029681691526, "C1");
    via_eta(s.f_hc, 0.140713647373789, "F_hc");
    via_eta(s.f_disp1, -0.252288278624795, "F_disp1");
    via_eta(s.f_disp2, -0.0575671271038317, "F_disp2");
    via_eta(s.f(), -0.169141758354838, "F");
}

/// A component the table has no PC-SAFT set for is refused, not solved for.
///
/// The compiled table spells absence as zeros in all three columns - 47 of the 286 rows -
/// so a model that took them would answer for a fluid with no segments.
#[test]
fn a_component_without_a_set_is_refused() {
    let absent = PcsaftComponent {
        m: 0.0,
        sigma: 0.0,
        epsik: 0.0,
    };
    let error = state(&[absent], &[0.0], &[1.0], 300.0, 4.55e-4).unwrap_err();
    let message = error.to_string();
    assert!(
        message.contains("no PC-SAFT set"),
        "the failure should say what is missing: {message}"
    );
}

/// The pressure the kernel predicts at NeqSim's own converged volume.
///
/// `P/(RT) = 1/v - F_V` and `F_V = -eta F_eta / v`, so `P v/(RT)` must be the state's
/// compressibility - and the probe reports both `v` and `Z` on the same phase. That makes
/// this the check on `F_eta`: it ties the derivative to the *state* the oracle converged
/// to, so a sign, a term or a power wrong in the chain shows up as a pressure that is not
/// the one the volume came from.
///
/// It also cannot be satisfied by a wrong `F` at a wrong `v`: the volume and the
/// compressibility are both NeqSim's.
#[test]
fn the_pressure_at_neqsims_volume_is_neqsims_compressibility() {
    let pure = pressure_over_rt(&[methane()], &[0.0], &[1.0], 300.0, 4.55032070500990e-4)
        .expect("a pressure");
    via_eta(pure * 4.55032070500990e-4, 0.912129702491155, "methane Z");

    let binary = pressure_over_rt(
        &[methane(), n_butane()],
        &[0.0, 0.022, 0.022, 0.0],
        &[0.6, 0.4],
        350.0,
        8.14109734626316e-4,
    )
    .expect("a pressure");
    via_eta(binary * 8.14109734626316e-4, 0.839270581279057, "binary Z");
}

/// `F_eta` is what it says, by finite difference of the energy in the packing fraction.
///
/// Self-contained rather than oracle-pinned: the energy is a function of `eta` and the
/// volume enters only through it, so a central difference in `v` divided by `-v/eta`
/// gives the derivative the pressure is built from.
#[test]
fn the_eta_derivative_is_the_energy_s_own() {
    let components = [methane(), n_butane()];
    let kij = [0.0, 0.022, 0.022, 0.0];
    let x = [0.6, 0.4];
    let (t, v) = (350.0, 8.14109734626316e-4);
    let h = 1.0e-8 * v;

    let up = state(&components, &kij, &x, t, v + h).expect("a state");
    let down = state(&components, &kij, &x, t, v - h).expect("a state");
    let base = state(&components, &kij, &x, t, v).expect("a state");
    // `dF/dv = -eta F_eta / v`, so `F_eta = -v/eta dF/dv`.
    let f_eta = -v / base.eta * (up.f() - down.f()) / (2.0 * h);

    let pressure = pressure_over_rt(&components, &kij, &x, t, v).expect("a pressure");
    let wanted = (pressure * v - 1.0) / base.eta;
    assert!(
        (f_eta / wanted - 1.0).abs() < 1.0e-6,
        "the finite difference gives F_eta = {f_eta} and the pressure's chain rule needs \
         {wanted}"
    );
}

/// The second derivative, tied to the energy rather than to the first derivative's own
/// code: a central difference of `A^R/(RT)` in `v` is `f''(v)`, and the pressure's
/// curvature is `-(f'' + 1/v^2)`.
#[test]
fn the_second_derivative_is_the_energy_s_own() {
    let components = [methane(), n_butane()];
    let kij = [0.0, 0.022, 0.022, 0.0];
    let x = [0.6, 0.4];
    let (t, v) = (350.0, 8.14109734626316e-4);
    // **`1e-4` of the volume, not of the energy's scale.** The second difference is
    // `f'' h^2` against an energy of order `0.1`, so it loses twelve digits to
    // cancellation before the truncation error is worth anything; measured over six
    // decades of `h`, this is where the two meet and the difference is good to `1e-8`.
    let h = 1.0e-4 * v;

    let up = state(&components, &kij, &x, t, v + h).expect("a state");
    let down = state(&components, &kij, &x, t, v - h).expect("a state");
    let base = state(&components, &kij, &x, t, v).expect("a state");
    let f_d2 = (up.f() - 2.0 * base.f() + down.f()) / (h * h);

    let d = d_pressure_over_rt_dv(&components, &kij, &x, t, v).expect("a derivative");
    let wanted = -f_d2 - 1.0 / (v * v);
    assert!(
        (d / wanted - 1.0).abs() < 1.0e-6,
        "the finite difference gives d2F/dv2 = {f_d2} and the pressure's chain rule needs \
         {d} against {wanted}"
    );
}

/// And the two public functions are consistent: the curvature of `P/(RT)` is the central
/// difference of `pressure_over_rt` itself. The first test ties the second derivative to
/// the energy and this one ties it to the pressure, so an error would have to be in both
/// the energy and the pressure to pass.
#[test]
fn the_second_derivative_is_the_pressure_s_own() {
    let components = [methane(), n_butane()];
    let kij = [0.0, 0.022, 0.022, 0.0];
    let x = [0.6, 0.4];
    let (t, v) = (350.0, 8.14109734626316e-4);
    let h = 1.0e-5 * v;

    let up = pressure_over_rt(&components, &kij, &x, t, v + h).expect("a pressure");
    let down = pressure_over_rt(&components, &kij, &x, t, v - h).expect("a pressure");
    let numeric = (up - down) / (2.0 * h);

    let d = d_pressure_over_rt_dv(&components, &kij, &x, t, v).expect("a derivative");
    assert!(
        (d / numeric - 1.0).abs() < 1.0e-6,
        "the finite difference gives {numeric} and the closed form gives {d}"
    );
}
