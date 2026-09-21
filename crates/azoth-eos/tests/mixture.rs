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
use azoth_eos::databank;
use azoth_eos::mixture::{Mixture, ReducedParameters, RootSide};
use azoth_eos::{Alpha, Cubic, MixingRule, SoreideWhitsonRole, pr_z_factor};

/// The pair, with the interaction parameter NeqSim fits for it.
///
/// Resolved through the databank rather than typed here, like every other fixture in
/// this tree. The literals this function used to carry had drifted from NeqSim's
/// `COMP.csv`: methane is 0.0115 and 4 599 000 Pa in the table, not 0.01142 and
/// 4 599 200, and the `kij` was an illustrative 0.05 where `INTER.csv` fits 0.01289789.
fn methane_butane() -> Mixture {
    databank::mixture_of(&["methane", "n-butane"], Cubic::Pr, None)
        .expect("the pair resolves")
        .0
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
    for (name, t, p) in [
        ("propane", 300.0, 1_000_000.0),
        ("propane", 350.0, 2_000_000.0),
        ("n-butane", 350.0, 1_000_000.0),
        ("methane", 200.0, 3_000_000.0),
    ] {
        let entry = databank::entry(name, None).expect("the databank has it");
        let (tc, omega) = (entry.tc, entry.omega);
        let mixture = databank::mixture_of(&[name], Cubic::Pr, None)
            .expect("the substance resolves")
            .0;
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
            "{name}, T={t}: psi_bar should be the component's own psi"
        );
        assert_eq!(
            state.h_dep_rt, pure.h_dep_rt,
            "{name}, T={t}: the mixture departure enthalpy should be bit-identical"
        );
        assert!(
            (state.s_dep_r - pure.s_dep_r).abs() < 1e-15,
            "{name}, T={t}: the entropy departures differ by {:e}, more than the \
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
    let ternary = databank::mixture_of(&["methane", "propane", "n-butane"], Cubic::Pr, None)
        .expect("the three resolve")
        .0;
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
    let ternary = databank::mixture_of(&["methane", "propane", "n-butane"], Cubic::Pr, None)
        .expect("the three resolve")
        .0;
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
/// `h_dep_rt` is `-1.961381e-04` at 1 kPa and `-1.962095e-03` at 10 kPa, a ratio of
/// **10.0036**; at 100 kPa and 1 MPa the ratios are 10.037 and 10.397, drifting up as
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

/// The mixture's own root at a state and composition.
///
/// Computed rather than written down. Two tests below took a literal `Z`, and it was a
/// literal for a fluid that no longer exists: `0.8274482588400789` belonged to the
/// hand-typed methane/n-butane with an illustrative `kij` of 0.05, and the fixture now
/// resolves that pair through the databank. The energy identity is a statement about
/// the *state*, so a `Z` that is not that state's makes it fail by percent for a reason
/// that looks like a defect.
fn compressibility_of(mixture: &Mixture, t: f64, p: f64, n: &[f64]) -> f64 {
    let reduced = mixture
        .reduced_parameters(kelvins(t), pascals(p))
        .expect("a state");
    let (a_mix, b_mix) = mixture.mixture_parameters(&reduced, n);
    azoth_eos::pr_z_factor(a_mix, b_mix).expect("a cubic").z_max
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
    let z = compressibility_of(&mixture, t, p, &n);
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

/// The same identity for an **associating** mixture, where the energy has two parts.
///
/// `helmholtz_energy` carried the cubic alone until it carried the association too, and
/// the two roots of an associating mixture at 356 K and 1 bar differ by `0.157 RT`. The
/// identity is what notices: with the cubic-only energy the two sides of this disagree by
/// a factor of three, while the *root choice* they feed still comes out the same way - so
/// a test of the verdict does not see it and this one does.
#[test]
fn the_energy_differentiates_to_the_fugacity_for_an_associating_mixture() {
    let (mixture, _) = databank::associating_mixture_of(&["water", "methanol"], Cubic::Srk, None)
        .expect("water and methanol bond");
    let (t, p, n) = (356.0, 1.0e5, [0.6, 0.4]);
    let reduced = mixture
        .reduced_parameters(kelvins(t), pascals(p))
        .expect("a state");

    let h = 1.0e-6;
    for liquid in [true, false] {
        let state = mixture
            .phase_state(
                &reduced,
                &n,
                if liquid {
                    RootSide::Liquid
                } else {
                    RootSide::Vapour
                },
            )
            .expect("a root");
        for i in 0..2 {
            let mut up = n;
            let mut down = n;
            up[i] += h;
            down[i] -= h;
            // `n` at a fixed compressibility is a perturbation at fixed volume: this
            // function's volume is `Z R T/P` whatever the moles.
            let gradient = (mixture
                .helmholtz_energy(&reduced, &up, state.z)
                .expect("an energy")
                - mixture
                    .helmholtz_energy(&reduced, &down, state.z)
                    .expect("an energy"))
                / (2.0 * h);
            let wanted = state.ln_phi[i] + state.z.ln();
            assert!(
                (gradient / wanted - 1.0).abs() < 1.0e-6,
                "{} root, component {i}: the gradient is {gradient} but `ln phi + ln Z` is {wanted}",
                if liquid { "liquid" } else { "vapour" }
            );
        }
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
    let z = compressibility_of(&mixture, 330.0, 2_500_000.0, &n);

    let hessian = mixture
        .helmholtz_hessian(&reduced, &n, z)
        .expect("a hessian");
    assert_eq!(hessian[0][1], hessian[1][0], "the Hessian is not symmetric");

    let q = mixture
        .criticality_matrix(&reduced, &n, z)
        .expect("a matrix");
    assert_eq!(q[0][1], q[1][0], "Q is not symmetric");
}

/// The triple-root solution of the Peng-Robinson cubic, which this library does *not*
/// ship - see `eos.pr_alpha_ab`'s assumptions and the constants' own comments.
///
/// Named here because the tests below are statements about the cubic rather than about
/// the constants: `Q` vanishes, and is minimised, at the critical point of a cubic whose
/// `Omega` pair satisfies the triple-root condition. NeqSim's pair does not, so the
/// shipped cubic's critical point is 4.8e-5 away from `Tr = Pr = 1` and its `Q` there is
/// 1.37e-4 rather than zero.
const TRIPLE_ROOT_A: f64 = 0.457_235_528_921_382_2;
const TRIPLE_ROOT_B: f64 = 0.077_796_073_903_888_46;

/// A one-component reduced state at a chosen `(A, B)`, built by hand.
///
/// Built directly rather than reached through a substance: the claim is about
/// `(A, B, Z)`, which is the point the rest of this file makes about the construction
/// carrying no dimensioned quantity. The `psi` fields are zero because nothing here
/// differentiates the alpha function - `Q` is a statement about `(A, B, Z)` alone.
///
/// The mixture a caller passes alongside this to `criticality_matrix` is a carrier and
/// nothing more: at one component there is no interaction parameter for it to hold.
fn one_component_state(a: f64, b: f64) -> ReducedParameters {
    ReducedParameters {
        a: vec![a],
        b: vec![b],
        psi: vec![0.0],
        psi_t: vec![0.0],
        reduced_temperatures: vec![1.0],
        t_kelvin: 300.0,
        pressure: 1.0e5,
        kij: vec![0.0],
        warnings: Vec::new(),
    }
}

/// The check a port cannot inherit, against an answer known in closed form.
///
/// Heidemann & Khalil's first condition is that the smallest eigenvalue of `Q` reach
/// zero. For one component `Q` is `1 x 1` and equals `H_11 + 1`, and the critical point
/// of a Peng-Robinson cubic is analytic: `Tr = Pr = 1`, `A = Omega_a`, `B = Omega_b`,
/// `Z_c = (1 - Omega_b)/3`. So the whole construction - the constant-volume Hessian,
/// the ideal part, the scaling - is checked against a closed form rather than against
/// another implementation.
///
/// **At the triple-root pair, which is what makes the closed form hold.** The state is
/// built directly rather than reached through a substance at its own tabulated `Tc` and
/// `Pc`, because the substance route lands on the cubic's critical point only for as
/// long as the library ships the triple-root `Omega` pair. It does not - it ships
/// NeqSim's, deliberately, and `eos.pr_alpha_ab`'s assumptions record why - so a
/// substance at `Tr = Pr = 1` is 4.8e-5 off it and `Q` there is 1.37e-4. Asserting the
/// closed form requires the pair the closed form is about, and building the reduced
/// state directly is the honest way to ask for it.
#[test]
fn the_criticality_matrix_vanishes_at_the_cubics_critical_point() {
    let mixture = databank::mixture_of(&["propane"], Cubic::Pr, None)
        .expect("propane resolves")
        .0;
    let q = mixture
        .criticality_matrix(
            &one_component_state(TRIPLE_ROOT_A, TRIPLE_ROOT_B),
            &[1.0],
            (1.0 - TRIPLE_ROOT_B) / 3.0,
        )
        .expect("a matrix");
    assert!(
        q[0][0].abs() < 1e-14,
        "Q is {} at the critical point, and it must be zero",
        q[0][0]
    );
}

/// How far the shipped pair puts the cubic's critical point from `Tr = Pr = 1`.
///
/// Recorded rather than asserted loosely, because it is the one number that says what
/// carrying NeqSim's literals costs: at the shipped pair the same construction gives
/// `Q` on the order of 1e-4 at `(A, B) = (Omega_a, Omega_b)`, and the cubic is at a
/// triple root at `Tr = 1 + 4.8e-5` instead. Bounded on both sides so a change to the
/// constants is visible here as well as in `eos.pr_alpha_ab`.
#[test]
fn the_shipped_omegas_leave_a_critical_point_residue() {
    let mixture = databank::mixture_of(&["propane"], Cubic::Pr, None)
        .expect("propane resolves")
        .0;
    let q = mixture
        .criticality_matrix(
            &one_component_state(azoth_eos::OMEGA_A, azoth_eos::OMEGA_B),
            &[1.0],
            (1.0 - azoth_eos::OMEGA_B) / 3.0,
        )
        .expect("a matrix");
    assert!(
        1e-5 < q[0][0] && q[0][0] < 1e-3,
        "Q at the shipped pair's Tr = Pr = 1 is {}",
        q[0][0]
    );
}

/// The residue is a minimum, and clearly positive on both sides of it.
///
/// Along `P = Pc` the critical point is the temperature at which `Q` reaches its least
/// value, and it is orders of magnitude larger either side. Without this, a `Q` that
/// was uniformly near zero everywhere would pass the test above - which is what a
/// missing diagonal term would give.
///
/// **The minimum is not zero, because the shipped `Omega` pair is not the triple-root
/// one** - measured at 2.5e-3 for propane at `T/Tc = 1`. What survives the port is
/// where the minimum sits (`T/Tc = 1`, to five figures) and how sharply it rises away
/// from it.
#[test]
fn the_criticality_matrix_is_minimised_at_the_critical_temperature() {
    let entry = databank::entry("propane", None).expect("the databank has propane");
    let (tc, pc) = (entry.tc, entry.pc);
    let mixture = databank::mixture_of(&["propane"], Cubic::Pr, None)
        .expect("propane resolves")
        .0;

    let q_at = |offset: f64| {
        let reduced = mixture
            .reduced_parameters(kelvins(tc * offset), pascals(pc))
            .expect("a state");
        let root = azoth_eos::pr_z_factor(reduced.a[0], reduced.b[0]).expect("a cubic");
        let z = if offset < 1.0 { root.z_min } else { root.z_max };
        mixture
            .criticality_matrix(&reduced, &[1.0], z)
            .expect("a matrix")[0][0]
    };

    let at_critical = q_at(1.0);
    assert!(
        1e-3 < at_critical && at_critical < 1e-2,
        "Q at T/Tc = 1 is {at_critical}"
    );

    for offset in [0.98, 0.99, 0.999, 1.001, 1.01, 1.02] {
        let value = q_at(offset);
        assert!(
            value > 1e-2,
            "T/Tc = {offset}: Q is {value}, but away from the critical point it must \
             be clearly positive"
        );
        assert!(
            value > at_critical,
            "T/Tc = {offset}: Q is {value}, below the value at the critical point - \
             the critical temperature is not the minimum"
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
    let z = compressibility_of(&mixture, 330.0, 2_500_000.0, &n);

    let energy = mixture
        .helmholtz_energy(&reduced, &n, z)
        .expect("an energy");
    assert!(
        (energy + 0.189_257_994_678_617_87).abs() < 1e-12,
        "the energy is {energy}, z is {z}, a_mix/b_mix check follows"
    );

    let hessian = mixture
        .helmholtz_hessian(&reduced, &n, z)
        .expect("a hessian");
    let wanted = [
        [-0.070_535_600_455_860_16, -0.283_559_373_215_204_75],
        [-0.283_559_373_215_204_75, -1.064_193_906_442_352_5],
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
        [0.957_678_639_726_483_9, -0.138_915_155_232_134_2],
        [-0.138_915_155_232_134_2, 0.574_322_437_423_059],
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

/// The Twu-Sim-Tassone cubic: Soave's constants, Peng-Robinson's geometry, Twu's alpha.
///
/// `ComponentTST` is the one cubic in NeqSim's set that is not a setting of the three
/// shapes `Cubic` already carried - it pairs Soave's `omega` pair (rounded to six
/// decimals) with Peng-Robinson's `delta`. The reduction checks all three axes at once:
/// the geometry literals, the alpha the builder selects, and the reduced `a`/`b` the
/// mixture actually computes.
#[test]
fn the_tst_cubic_reduces_to_its_constants() {
    assert_eq!(Cubic::Tst.omega_a(), 0.427481);
    assert_eq!(Cubic::Tst.omega_b(), 0.086641);
    assert_eq!(Cubic::Tst.delta1(), 1.0 + std::f64::consts::SQRT_2);
    assert_eq!(Cubic::Tst.delta2(), 1.0 - std::f64::consts::SQRT_2);

    let entry = databank::entry("methane", None).expect("the databank has it");
    let mixture = databank::mixture_of(&["methane"], Cubic::Pr, None)
        .expect("methane resolves")
        .0
        .with_cubic(Cubic::Tst);
    assert_eq!(mixture.alpha(), Alpha::Twu);

    let (t, p) = (300.0, 1_000_000.0);
    let reduced = mixture
        .reduced_parameters(kelvins(t), pascals(p))
        .expect("a state");
    let tr = t / entry.tc;
    let pr = p / entry.pc;
    let kappa = azoth_eos::twu_kappa(entry.omega)
        .expect("a coefficient")
        .kappa;
    let alpha = (1.0 + kappa * (1.0 - tr.sqrt())).powi(2);

    assert!(
        (reduced.a[0] - Cubic::Tst.omega_a() * alpha * pr / (tr * tr)).abs() < 1e-12,
        "a_reduced {} against the TST constant",
        reduced.a[0]
    );
    assert!(
        (reduced.b[0] - Cubic::Tst.omega_b() * pr / tr).abs() < 1e-12,
        "b_reduced {} against the TST constant",
        reduced.b[0]
    );

    let psi_expected = -kappa * tr.sqrt() / (1.0 + kappa * (1.0 - tr.sqrt()));
    assert!(
        (reduced.psi[0] - psi_expected).abs() < 1e-12,
        "psi {} against Twu's Soave form",
        reduced.psi[0]
    );

    // Peng-Robinson geometry means the vapour root is `pr_z_factor`'s, not `srk_z_factor`'s.
    let state = mixture
        .phase_state(&reduced, &[1.0], RootSide::Vapour)
        .expect("a phase");
    let roots = pr_z_factor(reduced.a[0], reduced.b[0]).expect("roots");
    assert_eq!(state.z, roots.z_max);
}

/// The volume translation mixes linearly and subtracts from the untranslated volume.
///
/// `ComponentPRvolcor`/`ComponentSrkvolcor` are PR/SRK plus a per-component Peneloux
/// shift; the shift mixes like the co-volume and is subtracted from `pr_molar_volume`'s
/// `v = z R T/P`. The check is the two-component linear mix, then the corrected volume
/// formed from the untranslated one.
#[test]
fn the_volume_translation_mixes_linearly_and_subtracts() {
    let (base, _) =
        databank::mixture_of(&["methane", "n-butane"], Cubic::Pr, None).expect("the pair resolves");
    let k = base.kij(0, 1);
    let methane = base.components()[0].clone().with_volume_shift(1.0e-6);
    let butane = base.components()[1].clone().with_volume_shift(3.0e-6);
    let mixture = Mixture::new(vec![methane, butane], vec![0.0, k, k, 0.0]).expect("a mixture");

    let x = [0.25, 0.75];
    let c_mix = 0.25 * 1.0e-6 + 0.75 * 3.0e-6;
    assert!(
        (mixture.volume_shift(&x) - c_mix).abs() < 1e-18,
        "the shift is the composition-weighted sum, not {}",
        mixture.volume_shift(&x)
    );

    let reduced = mixture
        .reduced_parameters(kelvins(330.0), pascals(2_500_000.0))
        .expect("a state");
    let state = mixture
        .phase_state(&reduced, &x, RootSide::Vapour)
        .expect("a phase");
    let v = azoth_eos::pr_molar_volume(state.z, kelvins(330.0), pascals(2_500_000.0))
        .expect("a volume")
        .v
        .value;
    assert!(
        (v - mixture.volume_shift(&x)) < v,
        "the translation subtracts from the untranslated molar volume"
    );
}

/// The Soreide-Whitson rule changes the aqueous `a_mix` through the mixture layer.
///
/// `MixingRule::phase_kij` is oracle-checked in `mixing_rule.rs`; this checks the other
/// half of the port - that `reduced_parameters` carries the reduced temperatures the rule
/// reads, and that `mixture_parameters` actually resolves against the phase composition
/// rather than the base matrix.
#[test]
fn the_soreide_whitson_rule_changes_the_aqueous_a_mix() {
    let (base, _) =
        databank::mixture_of(&["water", "methane"], Cubic::Pr, None).expect("the pair resolves");
    let k = base.kij(0, 1);
    let mixture = base.with_mixing_rule(MixingRule::SoreideWhitson {
        kij: vec![0.0, k, k, 0.0],
        roles: vec![SoreideWhitsonRole::Water, SoreideWhitsonRole::Hydrocarbon],
        salinity: 2.0,
    });
    let reduced = mixture
        .reduced_parameters(kelvins(300.0), pascals(10_000_000.0))
        .expect("a state");

    for (i, component) in mixture.components().iter().enumerate() {
        assert!(
            (reduced.reduced_temperatures[i] - 300.0 / component.tc.value).abs() < 1e-15,
            "reduced temperature {i} should be T / Tc"
        );
    }

    // The aqueous a_mix differs from the base-kij a_mix, because the (methane, water)
    // interaction is replaced by the salinity correlation.
    let aqueous = mixture.mixture_parameters(&reduced, &[0.9, 0.1]);
    let (a0, a1) = (reduced.a[0], reduced.a[1]);
    let base_a_mix =
        0.9 * 0.9 * a0 + 2.0 * 0.9 * 0.1 * (a0 * a1).sqrt() * (1.0 - k) + 0.1 * 0.1 * a1;
    assert!(
        (aqueous.0 - base_a_mix).abs() > 1e-9,
        "the aqueous a_mix must reflect the salinity correlation, not the base kij"
    );
}

/// The Huron-Vidal rule reproduces NeqSim's CLASSIC_HV water/ethanol liquid phase.
///
/// `a_mix` and the fugacity coefficients are checked against a NeqSim 3.20.0 TP flash
/// at T = 350 K, P = 1 bar, x = 0.5/0.5, with the fitted NRTL parameters its database
/// carries. The tolerance is the databank's, not the port's.
#[test]
fn the_huron_vidal_rule_matches_neqsims_water_ethanol_phase() {
    let (base, _) =
        databank::mixture_of(&["water", "ethanol"], Cubic::Srk, None).expect("the pair resolves");
    let k = base.kij(0, 1);
    let mixture = base
        .with_cubic(Cubic::Srk)
        .with_mixing_rule(MixingRule::HuronVidal {
            kij: vec![0.0, k, k, 0.0],
            hv_gij: vec![0.0, -2612.51, 2207.03, 0.0],
            hv_gij_t: vec![0.0, 7.3, -4.6, 0.0],
            hv_alpha: vec![0.0, 0.2245, 0.2245, 0.0],
            hv_pairs: vec![false, true, true, false],
        });
    let reduced = mixture
        .reduced_parameters(kelvins(350.0), pascals(100_000.0))
        .expect("a state");
    let state = mixture
        .phase_state(&reduced, &[0.5, 0.5], RootSide::Liquid)
        .expect("a phase");

    assert!(
        (state.a_mix - 0.017_011).abs() < 1e-3,
        "a_mix {}",
        state.a_mix
    );
    assert!(
        (state.ln_phi[0] - (-0.553_545)).abs() < 1e-3,
        "ln phi water {}",
        state.ln_phi[0]
    );
    assert!(
        (state.ln_phi[1] - 0.167_826).abs() < 1e-3,
        "ln phi ethanol {}",
        state.ln_phi[1]
    );
}

/// The Wong-Sandler rule reproduces NeqSim's water/ethanol liquid phase.
///
/// Same GE model as Huron-Vidal, but a GE-dependent `b_mix`, the rule's own `kij` and
/// its own temperature coefficient. Checked against a NeqSim 3.20.0 flash at
/// T = 350 K, P = 1 bar, x = 0.5/0.5.
///
/// Every matrix comes from the databank now rather than from this file: the parameters
/// were typed here until `parse_kij` learned to read the columns they come from, which
/// is what made them a caller's problem in the first place.
#[test]
fn the_wong_sandler_rule_matches_neqsims_water_ethanol_phase() {
    let (base, _) =
        databank::mixture_of(&["water", "ethanol"], Cubic::Pr, None).expect("the pair resolves");
    let ws =
        databank::wong_sandler_parameters(&["water", "ethanol"], None).expect("the pair resolves");
    assert!(
        ws.hv_pairs == vec![false, true, true, false],
        "water/ethanol is an `WS` pair: {:?}",
        ws.hv_pairs
    );
    let k = base.kij(0, 1);
    let mixture = base
        .with_cubic(Cubic::Srk)
        .with_mixing_rule(MixingRule::WongSandler {
            kij: vec![0.0, k, k, 0.0],
            hv_gij: ws.hv_gij,
            hv_gij_t: ws.hv_gij_t,
            hv_alpha: ws.hv_alpha,
            hv_pairs: ws.hv_pairs,
        });
    let reduced = mixture
        .reduced_parameters(kelvins(350.0), pascals(100_000.0))
        .expect("a state");
    let state = mixture
        .phase_state(&reduced, &[0.5, 0.5], RootSide::Liquid)
        .expect("a phase");

    assert!((state.z - 0.001_348_2).abs() < 1e-3, "z {}", state.z);
    assert!(
        (state.ln_phi[0] - (-1.605_35)).abs() < 1e-3,
        "ln phi water {}",
        state.ln_phi[0]
    );
    assert!(
        (state.ln_phi[1] - (-1.332_47)).abs() < 1e-3,
        "ln phi ethanol {}",
        state.ln_phi[1]
    );
}

/// `ln phi_i` as a function of mole numbers, which is what the surface differentiates.
///
/// The analytic derivative is the one of the *intensive* function `ln phi_i(T, P, x)`
/// with `x = n / sum(n)`, so the difference quotient has to be taken of the same
/// function. Calling [`Mixture::phase_state`] with a raw perturbed vector would
/// differentiate a different function - one whose `A` is `sum sum n_i n_j A_ij` rather
/// than the same over `N**2` - and the two agree only on the sum-to-one surface.
fn ln_phi_of(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    n: &[f64],
    side: RootSide,
) -> Vec<f64> {
    let total: f64 = n.iter().sum();
    let x: Vec<f64> = n.iter().map(|value| value / total).collect();
    mixture
        .phase_state(reduced, &x, side)
        .expect("a phase")
        .ln_phi
}

/// The composition derivative is the one the fugacity coefficient actually has.
///
/// A central difference of `ln phi_i(n)` at the `N = 1` the surface is stated at, with
/// the second-order truncation error falling as `h**2`; `h = 1e-6` puts the round-off
/// term below `1e-9` and the tolerance is an order looser than the identity the
/// Helmholtz layer is checked to.
#[test]
fn the_composition_derivative_is_the_fugacity_coefficients_own() {
    let mixture = methane_butane();
    let (t, p, n) = (330.0, 2_500_000.0, [0.6, 0.4]);
    let reduced = mixture
        .reduced_parameters(kelvins(t), pascals(p))
        .expect("a state");
    let side = RootSide::Liquid;
    // The root the difference is taken on, which is the one the surface is asked for.
    let (a_mix, b_mix) = mixture.mixture_parameters(&reduced, &n);
    let roots = azoth_eos::pr_z_factor(a_mix, b_mix).expect("a cubic");
    let z = roots.z_min;
    let derivatives = mixture
        .phase_derivatives(&reduced, &n, z)
        .expect("the derivative surface");

    let h = 1e-6;
    for j in 0..2 {
        let mut up = n;
        let mut down = n;
        up[j] += h;
        down[j] -= h;
        let above = ln_phi_of(&mixture, &reduced, &up, side);
        let below = ln_phi_of(&mixture, &reduced, &down, side);
        for i in 0..2 {
            let difference = (above[i] - below[i]) / (2.0 * h);
            assert!(
                (difference - derivatives.d_ln_phi_dn[i][j]).abs() < 1e-8,
                "d ln phi_{i} / d n_{j}: analytic {}, difference {difference}",
                derivatives.d_ln_phi_dn[i][j]
            );
        }
    }
}

/// The temperature and pressure derivatives are the fugacity coefficient's own too.
///
/// Both at constant composition, which is what the surface states them at and what a
/// flash's Newton step is written in.
#[test]
fn the_state_derivatives_are_the_fugacity_coefficients_own() {
    let mixture = methane_butane();
    let (t, p, n) = (330.0, 2_500_000.0, [0.6, 0.4]);
    let side = RootSide::Vapour;
    let at = |temperature: f64, pressure: f64| {
        let reduced = mixture
            .reduced_parameters(kelvins(temperature), pascals(pressure))
            .expect("a state");
        ln_phi_of(&mixture, &reduced, &n, side)
    };

    let reduced = mixture
        .reduced_parameters(kelvins(t), pascals(p))
        .expect("a state");
    let (a_mix, b_mix) = mixture.mixture_parameters(&reduced, &n);
    let z = azoth_eos::pr_z_factor(a_mix, b_mix).expect("a cubic").z_max;
    let derivatives = mixture
        .phase_derivatives(&reduced, &n, z)
        .expect("the derivative surface");

    // `h = 1e-3 K` and `1e-3 Pa` keep the round-off term of `1/h` small against the
    // values themselves, which are of order `1e-2` and `1e-6`.
    let (h_t, h_p) = (1e-3, 1e-2);
    let above = at(t + h_t, p);
    let below = at(t - h_t, p);
    for i in 0..2 {
        let difference = (above[i] - below[i]) / (2.0 * h_t);
        assert!(
            (difference - derivatives.d_ln_phi_dt[i]).abs() < 1e-7,
            "d ln phi_{i} / dT: analytic {}, difference {difference}",
            derivatives.d_ln_phi_dt[i]
        );
    }
    let above = at(t, p + h_p);
    let below = at(t, p - h_p);
    for i in 0..2 {
        let difference = (above[i] - below[i]) / (2.0 * h_p);
        assert!(
            (difference - derivatives.d_ln_phi_dp[i]).abs() < 1e-9,
            "d ln phi_{i} / dP: analytic {}, difference {difference}",
            derivatives.d_ln_phi_dp[i]
        );
    }
}

/// The composition derivative obeys Gibbs-Duhem, which no difference quotient can show.
///
/// At constant temperature and pressure `sum_i n_i d ln phi_i = 0`, so every column of
/// the matrix sums to zero against the composition. It is a *structural* property of
/// the whole matrix - it couples all `N**2` entries - so a single wrong term cannot
/// satisfy it, and it holds at any state rather than at the one a difference is taken
/// at. It is also the test a matrix built by differencing the raw vector fails.
/// **The fitted Soave coefficient is the one a pseudo-component carries.**
///
/// NeqSim's `AttractiveTermSrk.setm` puts a cut's own `m` on its attractive term, and its
/// TBP machinery does that for every pseudo-component. **It also overwrites the component's
/// acentric factor** with the root of the polynomial `m` belongs to, so the two are not
/// alternatives: this variant is what `setm` means, and a caller that supplies `m` should
/// supply the acentric factor it implies rather than the one `calcAcentricFactor` returns.
///
/// The reduction this is checked by: at one component the mixture's attraction must be the
/// registered `eos.srk_alpha_ab` at the coefficient the component carries, whatever that
/// coefficient is.
#[test]
fn the_fitted_soave_coefficient_is_the_components_own() {
    use azoth_eos::mixture::Component;
    use azoth_eos::{srk_alpha_ab, tbp_fraction_properties};

    // A `C19`-scale cut, whose fitted coefficient the correlation gives.
    let cut = tbp_fraction_properties(0.200, 800.0).expect("the cut computes");
    let fitted = cut.attraction_exponent;
    assert!((fitted - 1.564_540_560).abs() < 1e-9, "m = {fitted}");

    let component = Component::new(cut.tc, cut.pc, cut.acentric_factor)
        .expect("the cut's constants are positive")
        .with_molar_mass(Some(0.200))
        .with_alpha_params(vec![fitted]);
    let mixture = Mixture::new(vec![component.clone()], vec![0.0])
        .expect("one component and a 1x1 matrix")
        .with_cubic(Cubic::Srk)
        .with_alpha(Alpha::SrkFitted);

    let (t, p) = (450.0, 1.0e7);
    let reduced = mixture
        .reduced_parameters(kelvins(t), pascals(p))
        .expect("reduces");
    let expected = srk_alpha_ab(fitted, t / cut.tc.value, p / cut.pc.value).expect("alpha");
    assert!((reduced.a[0] - expected.a_reduced).abs() < 1e-14);
    assert!((reduced.a[0] / expected.a_reduced - 1.0).abs() < 1e-14);

    // **This fixture's acentric factor is `calcAcentricFactor`'s**, which is what NeqSim
    // computes and then discards - so here the two variants *do* disagree, and the check
    // below is that the variant reads the component's own coefficient rather than deriving
    // one. A caller following `setm` supplies the root instead, and then they agree.
    let correlated = Mixture::new(
        vec![component.clone().with_alpha_params(Vec::new())],
        vec![0.0],
    )
    .expect("one component")
    .with_cubic(Cubic::Srk)
    .with_alpha(Alpha::Srk);
    let other = correlated
        .reduced_parameters(kelvins(t), pascals(p))
        .expect("reduces");
    let relative = (reduced.a[0] - other.a[0]).abs() / other.a[0];
    assert!(
        relative > 0.01,
        "the two alphas differ by only {relative}, so this is not the seam it claims"
    );

    // A component that takes the fitted alpha and carries no coefficient is refused rather
    // than given Soave's default for an acentric factor nobody stated.
    let bare = Mixture::new(
        vec![Component::new(cut.tc, cut.pc, cut.acentric_factor).expect("positive")],
        vec![0.0],
    )
    .expect("one component")
    .with_cubic(Cubic::Srk)
    .with_alpha(Alpha::SrkFitted);
    let error = bare
        .reduced_parameters(kelvins(t), pascals(p))
        .expect_err("no coefficient, no alpha");
    assert!(
        matches!(error, azoth_core::AzothError::InvalidInput { .. }),
        "{error:?}"
    );
}

#[test]
fn the_composition_derivative_obeys_gibbs_duhem() {
    let mixture = methane_butane();
    for (t, p, n, side) in [
        (330.0, 2_500_000.0, [0.6, 0.4], RootSide::Liquid),
        (330.0, 2_500_000.0, [0.6, 0.4], RootSide::Vapour),
        (300.0, 3_000_000.0, [0.1, 0.9], RootSide::Liquid),
        (400.0, 1_000_000.0, [0.5, 0.5], RootSide::Vapour),
    ] {
        let reduced = mixture
            .reduced_parameters(kelvins(t), pascals(p))
            .expect("a state");
        let z = {
            let (a_mix, b_mix) = mixture.mixture_parameters(&reduced, &n);
            let roots = azoth_eos::pr_z_factor(a_mix, b_mix).expect("a cubic");
            match side {
                RootSide::Liquid => roots.z_min,
                RootSide::Vapour => roots.z_max,
            }
        };
        let derivatives = mixture
            .phase_derivatives(&reduced, &n, z)
            .expect("the derivative surface");
        for j in 0..n.len() {
            let sum: f64 = (0..n.len())
                .map(|i| n[i] * derivatives.d_ln_phi_dn[i][j])
                .sum();
            assert!(
                sum.abs() < 1e-9,
                "T={t}, P={p}, z={n:?}: column {j} sums to {sum}, not to zero"
            );
        }
    }
}

/// The surface refuses the rules it does not differentiate, by name.
///
/// A silent return of the classical matrix for an activity-coefficient rule would be a
/// wrong answer that looks like a right one, which is the failure mode this crate's
/// error types exist to make impossible.
#[test]
fn the_activity_rules_are_refused_rather_than_approximated() {
    let mixture = databank::mixture_of(&["water", "ethanol"], Cubic::Pr, None)
        .expect("the pair resolves")
        .0
        .with_mixing_rule(MixingRule::HuronVidal {
            kij: vec![0.0; 4],
            hv_gij: vec![0.0; 4],
            hv_gij_t: vec![0.0; 4],
            hv_alpha: vec![0.0; 4],
            hv_pairs: vec![false; 4],
        });
    let reduced = mixture
        .reduced_parameters(kelvins(350.0), pascals(1.0e5))
        .expect("a state");
    let error = mixture
        .phase_derivatives(&reduced, &[0.5, 0.5], 0.9)
        .expect_err("the activity rules are not differentiated");
    assert!(
        matches!(error, azoth_core::AzothError::InvalidInput { .. }),
        "got {error:?}"
    );
}

/// The composition-derivative matrix is symmetric, which couples every pair of entries.
///
/// `d ln phi_i / d n_j` is `(1/RT) d2(G^R)/dn_i dn_j` at constant `T` and `P`, and a
/// Gibbs energy's second derivative is symmetric. It is a different statement from
/// Gibbs-Duhem - that one couples a column to itself, this one couples a column to
/// another - so a matrix can satisfy either and fail the other, and the two together
/// pin the whole `N x N` block rather than a sum of it.
#[test]
fn the_composition_derivative_is_symmetric() {
    let mixture = methane_butane();
    for (t, p, n, side) in [
        (330.0, 2_500_000.0, [0.6, 0.4], RootSide::Liquid),
        (330.0, 2_500_000.0, [0.6, 0.4], RootSide::Vapour),
        (300.0, 3_000_000.0, [0.1, 0.9], RootSide::Vapour),
    ] {
        let reduced = mixture
            .reduced_parameters(kelvins(t), pascals(p))
            .expect("a state");
        let (a_mix, b_mix) = mixture.mixture_parameters(&reduced, &n);
        let roots = azoth_eos::pr_z_factor(a_mix, b_mix).expect("a cubic");
        let z = match side {
            RootSide::Liquid => roots.z_min,
            RootSide::Vapour => roots.z_max,
        };
        let d = mixture
            .phase_derivatives(&reduced, &n, z)
            .expect("the derivative surface");
        for i in 0..n.len() {
            for j in 0..n.len() {
                let gap = (d.d_ln_phi_dn[i][j] - d.d_ln_phi_dn[j][i]).abs();
                assert!(
                    gap < 1e-12,
                    "T={t}, P={p}, n={n:?}: d ln phi_{i}/dn_{j} and d ln phi_{j}/dn_{i} \
                     differ by {gap}"
                );
            }
        }
    }
}

/// The derivative surface against NeqSim's own, which is a third derivation of it.
///
/// `ComponentEos.dFdNdN`, `dFdNdT` and `dFdNdV` are NeqSim's analytic surface, reached
/// through `SystemThermo.init(3)` - `init(1)` leaves every one of its derivative arrays
/// at a placeholder, which is what the first run of the probe printed. The values below
/// are what `validation/neqsim/FugacityDerivativeProbe.java` printed for
/// methane/n-butane at 330 K and 25 bar, under Peng-Robinson with the classic mixing
/// rule, at the two compositions and roots NeqSim's own `TPflash` converges to:
///
/// ```text
///   gas, x = 0.69459826444/0.30540173556, Z = 0.86861125971
///     dfugdt [9.4100184113568750e-05, 0.0039357013111080480]
///     dfugdp [2.4218250434082655e-05, -0.017263724964853800]
///     dfugdx [[-0.068448243246064850, 0.15567701629505004],
///             [ 0.15567701629504996, -0.35406801187585346]]
///   liquid, x = 0.09509515775/0.90490484225, Z = 0.09315049464
///     dfugdt [0.00089198433650824160, 0.022295574688458470]
///     dfugdp [-0.036524606661905790, -0.036247642238949045]
///     dfugdx [[-0.78944362582463430, 0.082961503382128090],
///             [ 0.082961503382126390, -0.0087183059287235100]]
/// ```
///
/// Two conventions were measured rather than assumed. NeqSim's `dfugdx` is
/// `dfugdn * numberOfMolesInPhase`, and it is the constant-*pressure* composition
/// derivative: its columns satisfy `sum_i x_i dfugdx[i][j] = 0`, which is Gibbs-Duhem at
/// constant `T` and `P` and is not what a constant-volume frame gives. That is the same
/// matrix this returns, at unit total mole number. `dfugdp` is per bar where this is per
/// pascal, a factor of `1e5` and nothing else - the comparison below converts it.
///
/// The probe's fluid matters as much: `setMixingRule("classic")` and `setAttractiveTerm(1)`
/// are what `FlashTp.java` builds, and the integer overload of `setMixingRule` is a
/// different rule that splits the feed to a different composition. The first run did
/// that and disagreed with this surface by 1%.
#[test]
fn the_derivative_surface_matches_neqsims() {
    let mixture = methane_butane();
    let (t, p) = (330.0, 2_500_000.0);
    let reduced = mixture
        .reduced_parameters(kelvins(t), pascals(p))
        .expect("a state");

    for (name, x, side, z, wanted_dt, wanted_dp_bar, wanted_dn) in [
        (
            "gas",
            [0.694_598_264_442_796_2, 0.305_401_735_557_203_96],
            RootSide::Vapour,
            0.868_611_259_706_769_8,
            [9.410_018_411_356_875e-05, 3.935_701_311_108_048e-03],
            [2.421_825_043_408_265_5e-05, -1.726_372_496_485_38e-02],
            [
                [-0.068_448_243_246_064_85, 0.155_677_016_295_050_04],
                [0.155_677_016_295_049_96, -0.354_068_011_875_853_46],
            ],
        ),
        (
            "liquid",
            [0.095_095_157_748_050_7, 0.904_904_842_251_949_3],
            RootSide::Liquid,
            0.093_150_494_638_982_19,
            [8.919_843_365_082_416e-04, 2.229_557_468_845_847e-02],
            [-3.652_460_666_190_579e-02, -3.624_764_223_894_904e-02],
            [
                [-0.789_443_625_824_634_3, 0.082_961_503_382_128_09],
                [0.082_961_503_382_126_39, -0.008_718_305_928_723_51],
            ],
        ),
    ] {
        // The root is checked too: a derivative at a different root is a derivative at
        // a different phase, and the probe printed NeqSim's.
        let state = mixture.phase_state(&reduced, &x, side).expect("a phase");
        assert!(
            (state.z - z).abs() < 1e-12,
            "{name}: the root is {} but NeqSim's is {z}",
            state.z
        );
        let d = mixture
            .phase_derivatives(&reduced, &x, state.z)
            .expect("the derivative surface");

        for i in 0..2 {
            assert!(
                (d.d_ln_phi_dt[i] - wanted_dt[i]).abs() < 1e-13,
                "{name}: dfugdt[{i}] is {} but NeqSim's is {}",
                d.d_ln_phi_dt[i],
                wanted_dt[i]
            );
            // NeqSim reports per bar and this per pascal.
            assert!(
                (d.d_ln_phi_dp[i] * 1.0e5 - wanted_dp_bar[i]).abs() < 1e-14,
                "{name}: dfugdp[{i}] is {} per bar but NeqSim's is {}",
                d.d_ln_phi_dp[i] * 1.0e5,
                wanted_dp_bar[i]
            );
            for (j, &wanted) in wanted_dn[i].iter().enumerate() {
                assert!(
                    (d.d_ln_phi_dn[i][j] - wanted).abs() < 1e-13,
                    "{name}: dfugdx[{i}][{j}] is {} but NeqSim's is {wanted}",
                    d.d_ln_phi_dn[i][j]
                );
            }
        }
    }
}

/// **The two pressure-carrying terms are refused together rather than one winning.**
///
/// `PhaseElectrolyteCPA` is the model that has both, and it is carried rather than ported -
/// so a mixture naming the Wertheim association and the Fürst electrolyte at once is a
/// request for a model that does not exist here. Silently summing the two pressures would
/// give a root for neither.
#[test]
fn a_mixture_carrying_two_pressure_terms_is_refused() {
    use azoth_eos::furst_dielectric::MixingRule;
    use azoth_eos::furst_electrolyte::{FurstElectrolyte, FurstSpecies};

    let species = vec![
        FurstSpecies {
            name: "methane".into(),
            charge: 0.0,
            diameter_m: 2.52e-10,
            table_diameter: 2.52,
            dielectric_coefficients: [2.0, 0.0, 0.0, 0.0, 0.0],
            critical_volume: 9.9e-5,
            dielectric_at_reference: 2.0,
        },
        FurstSpecies {
            name: "water".into(),
            charge: 0.0,
            diameter_m: 2.52e-10,
            table_diameter: 2.52,
            dielectric_coefficients: [-19.2905, 29814.5, -0.019678, 0.000132, -3.11e-07],
            critical_volume: 5.6e-5,
            dielectric_at_reference: 78.4,
        },
    ];
    let term = FurstElectrolyte::new(species, MixingRule::default_for_the_model())
        .expect("the two build a table");

    // A plain mixture and a term is a statement about the model, and it is accepted.
    let (base, _) = azoth_eos::databank::mixture_of(&["methane", "water"], Cubic::Pr, None)
        .expect("the databank has both");
    let furst = base.clone().with_furst(term.clone());
    assert!(furst.furst().is_some());
    assert!(
        base.furst().is_none(),
        "the term is opt-in, like the association: the same substances are a classical \
         mixture under another model"
    );

    // Both at once is refused, at the phase rather than at construction, because which
    // terms a mixture runs is not a shape error.
    let both = furst.with_association();
    if let Ok(both) = both {
        let reduced = both
            .reduced_parameters(kelvins(298.15), pascals(1.0e5))
            .expect("reduces");
        let error = both
            .phase_state(&reduced, &[0.5, 0.5], RootSide::Vapour)
            .expect_err("the two terms cannot both carry the pressure");
        assert_eq!(error.field(), Some("mixture"), "{error:?}");
    }
}
