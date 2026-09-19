//! PC-SAFT's layers, against NeqSim's own.
//!
//! Every expected number comes from `validation/neqsim/PcsaftProbe.java`, which prints
//! `PhasePCSAFT`'s intermediates with the method each came from. A layer-by-layer
//! comparison is the point: a single Helmholtz energy that agrees can be two errors
//! cancelling, and this tranche has already spent a session on a quantity read at the
//! wrong scale.

use azoth_eos::pcsaft::{PcsaftComponent, state};

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
