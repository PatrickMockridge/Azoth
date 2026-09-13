//! The model layer's own arithmetic, checked by reduction.
//!
//! The mixture fugacity coefficient and the mixture departure functions have **no
//! registered spec**: the registry's inputs are scalars, and there is nowhere in it to
//! put a composition vector. That makes them the one piece of arithmetic in this crate
//! that no kernel checks directly, and the discipline that replaces a spec is
//! *reduction* - each must reproduce a registered calculation where the registered
//! calculation applies.
//!
//! This file is where those reductions live. The alternative, asserting the mixture
//! form against values it produced itself, would check nothing.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::mixture::{Component, Mixture, RootSide};

fn mix(tc: &[f64], pc: &[f64], omega: &[f64], kij: Vec<f64>) -> Mixture {
    let components = (0..tc.len())
        .map(|i| Component::new(kelvins(tc[i]), pascals(pc[i]), omega[i]).expect("valid"))
        .collect();
    Mixture::new(components, kij).expect("a valid mixture")
}

fn methane_butane() -> Mixture {
    mix(
        &[190.56, 425.12],
        &[4_599_200.0, 3_796_000.0],
        &[0.01142, 0.2002],
        vec![0.0, 0.05, 0.05, 0.0],
    )
}

/// At one component the mixture's departure functions *are* `eos.pr_departure`.
///
/// The strongest check available, and an exact one rather than a tolerance: the
/// mixture form is the pure form with `psi` replaced by a weighted average, and at
/// `N = 1` that average is the single component's own `psi`. `h_dep_rt` is asserted
/// **bit-identical**; `s_dep_r` gets one ulp, because it is formed by subtracting the
/// Gibbs term and the pure calc forms it by adding two others - a different summation
/// order and nothing more.
#[test]
fn the_departures_reduce_to_pr_departure_at_one_component() {
    for (tc, pc, omega, t, p) in [
        (369.83, 4_248_000.0, 0.1523, 300.0, 1_000_000.0),
        (369.83, 4_248_000.0, 0.1523, 350.0, 2_000_000.0),
        (425.12, 3_796_000.0, 0.2002, 350.0, 1_000_000.0),
        (190.56, 4_599_200.0, 0.01142, 200.0, 3_000_000.0),
    ] {
        let mixture = mix(&[tc], &[pc], &[omega], vec![0.0]);
        let reduced = mixture
            .reduced_parameters(kelvins(t), pascals(p))
            .expect("a state");
        let state = mixture
            .phase_state(&reduced, &[1.0], RootSide::Vapour)
            .expect("a phase");

        let kappa = azoth_eos::pr_kappa(omega).expect("a coefficient").kappa;
        let pure = azoth_eos::pr_departure(reduced.a[0], reduced.b[0], state.z, kappa, t / tc)
            .expect("a departure");

        assert_eq!(
            state.psi_bar, reduced.psi[0],
            "Tc={tc}, T={t}: psi_bar should be the component's own psi"
        );
        assert_eq!(
            state.h_dep_rt, pure.h_dep_rt,
            "Tc={tc}, T={t}: the mixture departure enthalpy should be bit-identical"
        );
        assert!(
            (state.s_dep_r - pure.s_dep_r).abs() < 1e-15,
            "Tc={tc}, T={t}: the entropy departures differ by {:e}, more than the \
             one ulp a different summation order explains",
            state.s_dep_r - pure.s_dep_r
        );
    }
}

/// The Gibbs identity holds for a mixture too, and it is what defines `s_dep_r`.
///
/// `h_dep_rt - s_dep_r = sum_i z_i ln phi_i` is exact for a mixture as it is for a
/// pure component - `G = H - TS` makes the departure Gibbs energy the composition-
/// weighted log of the fugacity coefficients - so a disagreement is a defect rather
/// than a residual to be tolerated. It is asserted at `N = 2` and `N = 3` because at
/// `N = 1` the reduction test above already covers it.
#[test]
fn the_gibbs_identity_holds_for_a_mixture() {
    let ternary = mix(
        &[190.56, 369.83, 425.12],
        &[4_599_200.0, 4_248_000.0, 3_796_000.0],
        &[0.01142, 0.1523, 0.2002],
        vec![0.0; 9],
    );
    for (mixture, t, p, z) in [
        (methane_butane(), 330.0, 2_500_000.0, vec![0.6, 0.4]),
        (methane_butane(), 300.0, 3_000_000.0, vec![0.1, 0.9]),
        (ternary, 320.0, 2_000_000.0, vec![0.5, 0.3, 0.2]),
    ] {
        for side in [RootSide::Liquid, RootSide::Vapour] {
            let reduced = mixture
                .reduced_parameters(kelvins(t), pascals(p))
                .expect("a state");
            let state = mixture.phase_state(&reduced, &z, side).expect("a phase");
            let g_dep_rt: f64 = z
                .iter()
                .zip(&state.ln_phi)
                .map(|(&z_i, &lp)| z_i * lp)
                .sum();
            assert!(
                (state.h_dep_rt - state.s_dep_r - g_dep_rt).abs() < 1e-14,
                "T={t}, P={p}, z={z:?}: h_dep_rt - s_dep_r is {}, but sum z ln phi is \
                 {g_dep_rt}",
                state.h_dep_rt - state.s_dep_r
            );
        }
    }
}

/// `psi_bar` is a weighted average, so it lies between the extremes of the `psi_i`.
///
/// A property rather than a value, and it is what would fail if the weights were
/// dropped, transposed, or summed as `i` alone. The weights are the `A_ij`, and a
/// one-component mixture is the degenerate case where the bound is attained - which is
/// the same statement as the reduction test above, read as an inequality.
#[test]
fn psi_bar_lies_between_the_components_psi() {
    let ternary = mix(
        &[190.56, 369.83, 425.12],
        &[4_599_200.0, 4_248_000.0, 3_796_000.0],
        &[0.01142, 0.1523, 0.2002],
        vec![0.0; 9],
    );
    for (t, p, z) in [
        (320.0, 2_000_000.0, vec![0.5, 0.3, 0.2]),
        (350.0, 5_000_000.0, vec![0.2, 0.3, 0.5]),
    ] {
        let reduced = ternary
            .reduced_parameters(kelvins(t), pascals(p))
            .expect("a state");
        let state = ternary
            .phase_state(&reduced, &z, RootSide::Vapour)
            .expect("a phase");
        let low = reduced.psi.iter().copied().fold(f64::INFINITY, f64::min);
        let high = reduced
            .psi
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            state.psi_bar >= low - 1e-15 && state.psi_bar <= high + 1e-15,
            "T={t}: psi_bar is {} but the components' psi lie in [{low}, {high}]",
            state.psi_bar
        );
    }
}

/// Both departures are proportional to pressure in the low-pressure limit.
///
/// At low pressure a fluid approaches the ideal gas and its departure functions go to
/// zero - and the *shape* of that approach is the checkable part, so this asserts the
/// shape rather than an arbitrary smallness. Measured on methane/n-butane at 330 K,
/// `h_dep_rt` is `-1.926095e-04` at 1 kPa and `-1.926771e-03` at 10 kPa, a ratio of
/// **10.0035**; at 100 kPa and 1 MPa the ratios are 10.035 and 10.382, drifting up as
/// the ideal-gas limit is left behind. The first ratio is the one asserted, at one per
/// cent, because it is the one that is still in the linear regime.
///
/// A sign error or a stray constant survives every point-value case one might choose
/// at ordinary pressures and shows up here, where the answer is nearly zero and its
/// *scaling* is all there is to see.
#[test]
fn the_departures_are_proportional_to_pressure_at_low_pressure() {
    let mixture = methane_butane();
    let at = |p: f64| {
        let reduced = mixture
            .reduced_parameters(kelvins(330.0), pascals(p))
            .expect("a state");
        mixture
            .phase_state(&reduced, &[0.6, 0.4], RootSide::Vapour)
            .expect("a phase")
    };
    let low = at(1_000.0);
    let mid = at(10_000.0);

    for (name, a, b) in [
        ("h_dep_rt", low.h_dep_rt, mid.h_dep_rt),
        ("s_dep_r", low.s_dep_r, mid.s_dep_r),
    ] {
        let ratio = b / a;
        assert!(
            (ratio - 10.0).abs() < 0.1,
            "{name}: a tenfold pressure rise should multiply the departure by ten, \
             but {a:e} became {b:e} - a ratio of {ratio}"
        );
    }
}
