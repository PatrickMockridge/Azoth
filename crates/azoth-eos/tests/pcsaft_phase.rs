//! PC-SAFT's volume solve, against the volumes NeqSim converges to.
//!
//! The oracle is `validation/neqsim/PcsaftProbe.java`, run at each state below:
//! `volumeSAFT` is the molar volume its Newton loop stopped at, `Z` the compressibility
//! factor recorded on the same phase, and `nSAFT` the packing fraction, which the probe
//! prints and this module also reports. Four states, two of them on the vapour branch
//! and two on the liquid one, because the branches are different searches and a solve
//! that works on one says nothing about the other.
//!
//! **Those four are not enough on their own.** At each of them the isotherm has one root
//! in the region the phase is in, so a liquid branch that quietly returned the vapour
//! root would pass - which is what one draft did. The fifth state is inside the
//! isotherm's loop, where the two roots are different numbers, and it is the one that
//! holds the branches apart.

use azoth_eos::mixture::RootSide;
use azoth_eos::pcsaft::PcsaftComponent;
use azoth_eos::pcsaft_phase::molar_volume;

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

fn propane() -> PcsaftComponent {
    PcsaftComponent {
        m: 2.002,
        sigma: 3.6184e-10,
        epsik: 208.11,
    }
}

/// **`1e-8`, and the limit is NeqSim's rather than this crate's.** Its Newton loop stops
/// on a relative *step* of `1e-10`, so the volume it prints is its own answer to about
/// that; and the packing fraction it prints is looser still, because `volInit` builds it
/// from a volume the phase does not report elsewhere - the `1.1e-9` its own bookkeeping
/// costs, which `tests/pcsaft.rs` records. Measured, the volumes here agree to `1.3e-9`:
/// inside what the oracle can offer, and not a licence to ask for less.
fn matches(actual: f64, expected: f64, context: &str) {
    let relative = (actual / expected - 1.0).abs();
    assert!(
        relative < 1.0e-8,
        "{context}: {actual} against the probe's {expected}, a relative {relative:e}"
    );
}

#[test]
fn methane_and_the_binary_converge_to_neqsims_vapour_volume() {
    let pure = molar_volume(&[methane()], &[0.0], &[1.0], 300.0, 5.0e6, RootSide::Vapour)
        .expect("a volume");
    matches(pure.v, 4.55032070500990e-4, "methane v");
    matches(pure.z, 0.912129702491155, "methane Z");
    matches(pure.eta, 0.0324636218320417, "methane eta");

    let binary = molar_volume(
        &[methane(), n_butane()],
        &[0.0, 0.022, 0.022, 0.0],
        &[0.6, 0.4],
        350.0,
        3.0e6,
        RootSide::Vapour,
    )
    .expect("a volume");
    matches(binary.v, 8.14109734626316e-4, "binary v");
    matches(binary.z, 0.839270581279057, "binary Z");
    matches(binary.eta, 0.0281366301313126, "binary eta");
}

/// The liquid branch, where the root is found by walking up from the packing floor rather
/// than by seeding. Both states are dense enough that the phase NeqSim converges to is
/// the liquid one.
#[test]
fn methane_and_propane_converge_to_neqsims_liquid_volume() {
    let cold = molar_volume(&[methane()], &[0.0], &[1.0], 150.0, 5.0e6, RootSide::Liquid)
        .expect("a volume");
    matches(cold.v, 4.34500702054451e-5, "methane v");
    matches(cold.z, 0.174194753210883, "methane Z");
    matches(cold.eta, 0.362239782531743, "methane eta");

    let dense = molar_volume(&[propane()], &[0.0], &[1.0], 300.0, 1.0e7, RootSide::Liquid)
        .expect("a volume");
    matches(dense.v, 8.57488611218831e-5, "propane v");
    matches(dense.z, 0.343773937094727, "propane Z");
    matches(dense.eta, 0.333379192828988, "propane eta");
}

/// Above the isotherm's loop there is one root, so **both sides are the same answer** -
/// which is what makes the branches a statement about the phase rather than about the
/// fluid. Methane at 400 K is past the loop and the probe records one volume.
#[test]
fn an_isotherm_with_one_root_gives_it_to_both_sides() {
    let vapour = molar_volume(&[methane()], &[0.0], &[1.0], 400.0, 5.0e6, RootSide::Vapour)
        .expect("a volume");
    let liquid = molar_volume(&[methane()], &[0.0], &[1.0], 400.0, 5.0e6, RootSide::Liquid)
        .expect("a volume");
    matches(vapour.v, 6.49279638776721e-4, "methane v");
    matches(vapour.z, 0.976129951289097, "methane Z");
    // **To a tolerance, not exactly.** The two branches reach this root by different
    // routes - Newton from the ideal gas and a bisection of the walk - so they meet in
    // the twelfth digit rather than in the last one.
    matches(vapour.v, liquid.v, "the two branches at one root");
}

/// **The vapour branch's Newton converges, rather than the walk answering for it.**
///
/// Every state above reproduces NeqSim whether or not the Newton works, because a vapour
/// solve that does not converge falls back to the walk and the walk finds the same root.
/// That is exactly what hid a sign error in the step: the first version of this solve
/// reported `iterations = 0` at every state and matched the oracle to ten digits through
/// the fallback alone. So the count is asserted here rather than the answer.
#[test]
fn the_vapour_branch_converges_by_newton() {
    for (t, p) in [(300.0, 5.0e6), (400.0, 5.0e6)] {
        let solved =
            molar_volume(&[methane()], &[0.0], &[1.0], t, p, RootSide::Vapour).expect("a volume");
        assert!(
            solved.iterations > 0,
            "at {t} K and {p} Pa the vapour branch took the walk's answer, so the Newton \
             step is not converging"
        );
    }
    // The liquid branch is a bisection by construction and reports no steps, which is how
    // a caller can tell the two paths apart in the result.
    let liquid = molar_volume(&[methane()], &[0.0], &[1.0], 150.0, 5.0e6, RootSide::Liquid)
        .expect("a volume");
    assert_eq!(
        liquid.iterations, 0,
        "the liquid branch is bracketed, not Newtoned"
    );
}

/// A pressure at or below zero is refused rather than solved for. The solve divides by
/// it, and `R T/P` is the seed.
#[test]
fn a_pressure_that_is_not_positive_is_refused() {
    let error = molar_volume(&[methane()], &[0.0], &[1.0], 300.0, 0.0, RootSide::Vapour)
        .expect_err("zero pressure");
    assert!(
        error.to_string().contains("P"),
        "the failure should name the pressure: {error}"
    );
}

/// The residual is what the solve is *for*, so it is checked at the answer rather than
/// only the answer being checked against the oracle: a volume that reproduces NeqSim but
/// does not satisfy `P_calc = P` would be a coincidence of the oracle's own tolerance.
#[test]
fn the_answer_satisfies_its_own_equation() {
    let solved = molar_volume(
        &[methane(), n_butane()],
        &[0.0, 0.022, 0.022, 0.0],
        &[0.6, 0.4],
        350.0,
        3.0e6,
        RootSide::Vapour,
    )
    .expect("a volume");
    let residual = azoth_eos::pcsaft::pressure_over_rt(
        &[methane(), n_butane()],
        &[0.0, 0.022, 0.022, 0.0],
        &[0.6, 0.4],
        350.0,
        solved.v,
    )
    .expect("a pressure")
        * 8.3144621
        * 350.0
        - 3.0e6;
    assert!(
        residual.abs() / 3.0e6 < 1.0e-10,
        "the residual at the converged volume is {residual} Pa"
    );
}

/// **Two roots, and the two branches must part.** At 150 K and 50 bara the isotherm has
/// run past its loop and only the liquid root is left, so that state cannot tell the
/// branches apart - which is why it is not enough. Here the pressure is inside the loop,
/// so there is a dilute root and a dense one, and the sides have to disagree.
///
/// The check is self-contained because NeqSim's flash reports one phase per call and the
/// two-phase states are not what the probe reads; what it can assert is that each branch
/// satisfies the equation it was solved for, and that they are on the sides they claim.
#[test]
fn a_pressure_inside_the_loop_gives_the_two_branches_different_volumes() {
    let components = [methane()];
    let kij = [0.0];
    let x = [1.0];
    let (t, p) = (150.0, 5.0e5);

    let vapour = molar_volume(&components, &kij, &x, t, p, RootSide::Vapour).expect("a volume");
    let liquid = molar_volume(&components, &kij, &x, t, p, RootSide::Liquid).expect("a volume");

    assert!(
        vapour.v > liquid.v,
        "the vapour root {} is not above the liquid one {}",
        vapour.v,
        liquid.v
    );
    for (side, v) in [("vapour", vapour.v), ("liquid", liquid.v)] {
        let residual = azoth_eos::pcsaft::pressure_over_rt(&components, &kij, &x, t, v)
            .expect("a pressure")
            * 8.3144621
            * t
            - p;
        assert!(
            residual.abs() / p < 1.0e-8,
            "the {side} volume {v} has a residual of {residual} Pa"
        );
    }
}
