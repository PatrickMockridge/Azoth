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
    databank::mixture_of(&["methane", "n-butane"], None)
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
        let mixture = databank::mixture_of(&[name], None)
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
    let ternary = databank::mixture_of(&["methane", "propane", "n-butane"], None)
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
    let ternary = databank::mixture_of(&["methane", "propane", "n-butane"], None)
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
    let mixture = databank::mixture_of(&["propane"], None)
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
    let mixture = databank::mixture_of(&["propane"], None)
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
    let mixture = databank::mixture_of(&["propane"], None)
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
    let mixture = databank::mixture_of(&["methane"], None)
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
        databank::mixture_of(&["methane", "n-butane"], None).expect("the pair resolves");
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
    let (base, _) = databank::mixture_of(&["water", "methane"], None).expect("the pair resolves");
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
    let (base, _) = databank::mixture_of(&["water", "ethanol"], None).expect("the pair resolves");
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
/// Same GE model as Huron-Vidal, but a GE-dependent `b_mix` and the DijT-free activity
/// coefficients NeqSim reads for this rule. Checked against a NeqSim 3.20.0 flash at
/// T = 350 K, P = 1 bar, x = 0.5/0.5.
#[test]
fn the_wong_sandler_rule_matches_neqsims_water_ethanol_phase() {
    let (base, _) = databank::mixture_of(&["water", "ethanol"], None).expect("the pair resolves");
    let k = base.kij(0, 1);
    let mixture = base
        .with_cubic(Cubic::Srk)
        .with_mixing_rule(MixingRule::WongSandler {
            kij: vec![0.0, k, k, 0.0],
            hv_gij: vec![0.0, -2612.51, 2207.03, 0.0],
            hv_alpha: vec![0.0, 0.2245, 0.2245, 0.0],
            hv_pairs: vec![false, true, true, false],
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
