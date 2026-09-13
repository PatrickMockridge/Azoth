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

/// The Helmholtz energy's gradient *is* the fugacity coefficient, by a different route.
///
/// `d(A^R/RT)/dn_i` must equal `ln phi_i + ln Z`. The right-hand side comes from
/// [`Mixture::phase_state_at`], which reaches it by differentiating a departure
/// function; the left is a central difference of the energy directly. Two routes to
/// one quantity, so the agreement is evidence rather than a restatement - and it is
/// the identity that ties the critical point's machinery to the flash's.
#[test]
fn the_energy_differentiates_to_the_fugacity_coefficient() {
    let mixture = methane_butane();
    let (t, p, n) = (330.0, 2_500_000.0, [0.6, 0.4]);
    let reduced = mixture
        .reduced_parameters(kelvins(t), pascals(p))
        .expect("a state");
    let z = 0.827_448_258_840_078_9;
    let state = mixture.phase_state_at(&reduced, &n, z).expect("a phase");

    let h = 1e-6;
    for i in 0..2 {
        let mut up = n;
        let mut down = n;
        up[i] += h;
        down[i] -= h;
        let gradient = (mixture
            .helmholtz_energy(&reduced, &up, z)
            .expect("an energy")
            - mixture
                .helmholtz_energy(&reduced, &down, z)
                .expect("an energy"))
            / (2.0 * h);
        let wanted = state.ln_phi[i] + z.ln();
        assert!(
            (gradient - wanted).abs() < 1e-9,
            "component {i}: the gradient is {gradient}, but `ln phi + ln Z` is {wanted}"
        );
    }
}

/// The Hessian is the second derivative of the energy, checked by finite difference.
///
/// The step is `1e-3` and the tolerance is loose because a central second difference
/// of a function evaluated in `f64` cannot do better: the truncation error falls as
/// `h**2` and the round-off rises as `eps/h**2`, and the measured minimum of the sum
/// is around `3e-8` at `h = 1e-3`. A tighter tolerance would be a claim about the
/// difference quotient rather than about the Hessian - and it would still catch any
/// error in a term, which moves an entry by order `0.1`.
#[test]
fn the_hessian_is_the_second_derivative_of_the_energy() {
    let mixture = methane_butane();
    let (t, p, n) = (330.0, 2_500_000.0, [0.6, 0.4]);
    let reduced = mixture
        .reduced_parameters(kelvins(t), pascals(p))
        .expect("a state");
    let z = 0.827_448_258_840_078_9;
    let hessian = mixture
        .helmholtz_hessian(&reduced, &n, z)
        .expect("a hessian");

    let h = 1e-3;
    for i in 0..2 {
        for j in 0..2 {
            let corner = |di: f64, dj: f64| {
                let mut shifted = n;
                shifted[i] += di * h;
                shifted[j] += dj * h;
                mixture
                    .helmholtz_energy(&reduced, &shifted, z)
                    .expect("an energy")
            };
            let fd = (corner(1.0, 1.0) - corner(1.0, -1.0) - corner(-1.0, 1.0)
                + corner(-1.0, -1.0))
                / (4.0 * h * h);
            assert!(
                (hessian[i][j] - fd).abs() < 1e-6,
                "H[{i}][{j}] is {} but the finite difference gives {fd}",
                hessian[i][j]
            );
        }
    }
}

/// `H_ij == H_ji` exactly, not merely closely.
///
/// The expression is symmetric term by term, so the two must agree bit for bit. The
/// check matters because the eigenvalues are only real if the matrix is symmetric,
/// and a `f64` asymmetry that appeared here would be invisible until it produced a
/// complex eigenvector somewhere downstream.
#[test]
fn the_hessian_and_the_criticality_matrix_are_exactly_symmetric() {
    let mixture = methane_butane();
    let reduced = mixture
        .reduced_parameters(kelvins(330.0), pascals(2_500_000.0))
        .expect("a state");
    let n = [0.6, 0.4];
    let z = 0.827_448_258_840_078_9;

    let hessian = mixture
        .helmholtz_hessian(&reduced, &n, z)
        .expect("a hessian");
    assert_eq!(hessian[0][1], hessian[1][0], "the Hessian is not symmetric");

    let q = mixture
        .criticality_matrix(&reduced, &n, z)
        .expect("a matrix");
    assert_eq!(q[0][1], q[1][0], "Q is not symmetric");
}

/// `Q` vanishes at a pure component's critical point - the check a port cannot inherit.
///
/// Heidemann & Khalil's first condition is that the smallest eigenvalue of `Q` reach
/// zero. For one component `Q` is `1 x 1` and equals `H_11 + 1`, and PR's critical
/// point is known in closed form: `Tr = Pr = 1, Z_c = (1 - omega_b)/3`. So the whole
/// construction - the constant-volume Hessian, the ideal part, the scaling - is
/// checked against an analytic answer rather than against another implementation.
///
/// **Every component gives the same number**, because at `Tr = Pr = 1` they all have
/// the same reduced state `(A, B) = (omega_a, omega_b)`. That is worth asserting
/// separately: it is the statement that this construction depends on `(A, B, Z)` and
/// on nothing else, which is what "the eos core carries no dimensioned quantity"
/// means in practice.
#[test]
fn the_criticality_matrix_vanishes_at_a_pure_components_critical_point() {
    let z_critical = (1.0 - azoth_eos::OMEGA_B) / 3.0;
    let mut seen: Option<f64> = None;
    for (tc, pc, omega) in [
        (369.83, 4_248_000.0, 0.1523),
        (190.56, 4_599_200.0, 0.01142),
        (425.12, 3_796_000.0, 0.2002),
        (304.13, 7_377_000.0, 0.2239),
    ] {
        let mixture = mix(&[tc], &[pc], &[omega], vec![0.0]);
        // `Tr = Pr = 1` is the critical state, so the reduced parameters are the
        // universal `(omega_a, omega_b)` and are the same for every component here.
        let reduced = mixture
            .reduced_parameters(kelvins(tc), pascals(pc))
            .expect("a state");
        assert!((reduced.a[0] - azoth_eos::OMEGA_A).abs() < 1e-15);
        assert!((reduced.b[0] - azoth_eos::OMEGA_B).abs() < 1e-15);

        let q = mixture
            .criticality_matrix(&reduced, &[1.0], z_critical)
            .expect("a matrix");
        assert!(
            q[0][0].abs() < 1e-14,
            "Tc={tc}: Q is {} at the critical point, and it must be zero",
            q[0][0]
        );
        match seen {
            None => seen = Some(q[0][0]),
            Some(first) => assert_eq!(first, q[0][0], "Tc={tc}: Q differs between components"),
        }
    }
}

/// And the zero is a minimum, not a small number that happens to be near one.
///
/// Along `P = Pc` the critical point is the temperature at which `Q` reaches zero,
/// and it is positive on both sides. Without this, `Q = 1e-16` at one state would be
/// equally consistent with a construction that is uniformly near zero everywhere -
/// which is what a missing diagonal term would give.
#[test]
fn the_criticality_matrix_is_minimised_at_the_critical_temperature() {
    let z_critical = (1.0 - azoth_eos::OMEGA_B) / 3.0;
    let (tc, pc, omega) = (369.83, 4_248_000.0, 0.1523);
    let mixture = mix(&[tc], &[pc], &[omega], vec![0.0]);

    let at_critical = {
        let reduced = mixture
            .reduced_parameters(kelvins(tc), pascals(pc))
            .expect("a state");
        mixture
            .criticality_matrix(&reduced, &[1.0], z_critical)
            .expect("a matrix")[0][0]
    };
    assert!(at_critical.abs() < 1e-14);

    for offset in [0.98, 0.99, 0.999, 1.001, 1.01, 1.02] {
        let reduced = mixture
            .reduced_parameters(kelvins(tc * offset), pascals(pc))
            .expect("a state");
        let root = azoth_eos::pr_z_factor(reduced.a[0], reduced.b[0]).expect("a cubic");
        let z = if offset < 1.0 { root.z_min } else { root.z_max };
        let q = mixture
            .criticality_matrix(&reduced, &[1.0], z)
            .expect("a matrix")[0][0];
        assert!(
            q > 1e-4,
            "T/Tc = {offset}: Q is {q}, but away from the critical point it must be \
             clearly positive"
        );
    }
}

/// The two implementations agree on the Helmholtz layer, to the last few digits.
///
/// The values are pinned here and in `python/tests/eos/test_mixture_layer.py`, and
/// that is the whole check: two languages, the same algebra written out separately,
/// one set of numbers. The tolerance is `1e-12` rather than bit-equality because
/// `ln` is not correctly rounded in either language - the same reason the departure
/// tests are not bit-exact - and `+ - * / sqrt` alone would have permitted `to_bits`.
#[test]
fn both_implementations_agree_on_the_helmholtz_layer() {
    let mixture = methane_butane();
    let reduced = mixture
        .reduced_parameters(kelvins(330.0), pascals(2_500_000.0))
        .expect("a state");
    let n = [0.6, 0.4];
    let z = 0.827_448_258_840_078_9;

    let energy = mixture
        .helmholtz_energy(&reduced, &n, z)
        .expect("an energy");
    assert!(
        (energy + 0.184_316_374_846_176_38).abs() < 1e-12,
        "the energy is {energy}"
    );

    let hessian = mixture
        .helmholtz_hessian(&reduced, &n, z)
        .expect("a hessian");
    let wanted = [
        [-0.070_616_465_542_426_04, -0.266_557_518_655_936_9],
        [-0.266_557_518_655_936_9, -1.060_519_445_307_295_5],
    ];
    for i in 0..2 {
        for j in 0..2 {
            assert!(
                (hessian[i][j] - wanted[i][j]).abs() < 1e-12,
                "H[{i}][{j}] is {} but the other implementation gives {}",
                hessian[i][j],
                wanted[i][j]
            );
        }
    }

    let q = mixture
        .criticality_matrix(&reduced, &n, z)
        .expect("a matrix");
    let wanted = [
        [0.957_630_120_674_544_4, -0.130_585_981_561_890_6],
        [-0.130_585_981_561_890_6, 0.575_792_221_877_081_8],
    ];
    for i in 0..2 {
        for j in 0..2 {
            assert!(
                (q[i][j] - wanted[i][j]).abs() < 1e-12,
                "Q[{i}][{j}] is {} but the other implementation gives {}",
                q[i][j],
                wanted[i][j]
            );
        }
    }
}
