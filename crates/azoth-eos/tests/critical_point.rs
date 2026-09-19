//! Spec-driven tests for the `eos.critical_point` model.
//!
//! The spec's cases pin the two implementations to each other. These tests are for
//! what a case cannot say, and one of them matters more than the rest: **a pure
//! component's critical point is known in closed form**, so the whole construction -
//! the constant-volume Hessian, the ideal part, the scaling, the cubic form and the
//! nesting - is checked against an analytic answer rather than against a second
//! reading of the same method.

use azoth_core::{AzothError, CalcResult};
use azoth_eos::Cubic;
use azoth_eos::critical_point::symmetric_eigen;
use azoth_eos::mixture::Mixture;
use azoth_eos::{critical_point, databank, model_gen};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.critical_point";

/// The methane/n-butane pair the spec's cases use, resolved through the databank.
fn methane_butane() -> Mixture {
    databank::mixture_of(&["methane", "n-butane"], Cubic::Pr, None)
        .expect("the pair resolves")
        .0
}

fn mixture_from_case(case: &azoth_core::spec::TestCase) -> Mixture {
    let names = case
        .list("components")
        .expect("the case declares components");
    databank::mixture_of(names, Cubic::Pr, None)
        .expect("the case's components resolve")
        .0
}

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::CriticalPointResult {
    let mixture = mixture_from_case(case);
    critical_point(&mixture, case.vector("z").expect("the case declares z"))
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        common::assert_close(
            result.tc.value,
            common::expected(case, "tc"),
            case.tolerance,
            &format!("{}::{} (tc)", spec.id, case.id),
        );
        common::assert_close(
            result.pc.value,
            common::expected(case, "pc"),
            case.tolerance,
            &format!("{}::{} (pc)", spec.id, case.id),
        );
        common::assert_close(
            result.vc.value,
            common::expected(case, "vc"),
            case.tolerance,
            &format!("{}::{} (vc)", spec.id, case.id),
        );
        common::assert_close(
            result.z_c,
            common::expected(case, "z_c"),
            case.tolerance,
            &format!("{}::{} (z_c)", spec.id, case.id),
        );
        assert!(
            result.iterations > 0,
            "{}::{}: the solve took no steps",
            spec.id,
            case.id
        );
        assert!(
            result.residual <= spec.algorithm.expect("a procedure").tolerance,
            "{}::{}: residual {:e} exceeds the declared tolerance",
            spec.id,
            case.id,
            result.residual
        );
        common::assert_consistent(&result, case.id);
    }
}

/// A pure component's critical point is analytic, and this is where the `Z_c` that
/// every test in this repository asserts comes from: `(1 - omega_b)/3`.
///
/// The check is on the whole construction rather than on the arithmetic, because
/// every piece of it - the Hessian at constant volume, the ideal part
/// `delta_ij/n_i`, the `sqrt(n_i n_j)` scaling, the ideal third derivative in the
/// cubic form, the nesting - has to be right for this to come out.
#[test]
fn a_pure_component_reproduces_the_analytic_critical_point() {
    let expected_z_c = (1.0 - azoth_eos::OMEGA_B) / 3.0;
    for name in ["propane", "methane", "n-butane", "co2"] {
        let entry = databank::entry(name, None).expect("a databank entry");
        let (tc, pc) = (entry.tc, entry.pc);
        let mixture = databank::mixture_of(&[name], Cubic::Pr, None)
            .expect("resolves")
            .0;
        let r = critical_point(&mixture, &[1.0]).expect("a critical point");

        // Measured, on both implementations: Tc to 5.5e-05 relative, Pc to 1.5e-04,
        // Z_c to 2.4e-06 absolute - and the same four figures for all four components,
        // which is the signature of a systematic departure rather than of noise.
        //
        // The departure is the port's, and it is understood: `(1 - omega_b)/3` is the
        // critical compressibility of a cubic whose Omega pair satisfies the
        // triple-root condition, and NeqSim's does not - it is 5.6e-06 off it. So
        // NeqSim's cubic reaches its triple root at `Tr = 1 + 4.8e-05` rather than at
        // `Tr = 1`, and the construction lands there exactly as it should. The
        // tolerances below were 1e-10 and 1e-08 before the port adopted NeqSim's
        // constants. See `eos.pr_alpha_ab`'s assumptions and
        // `validation/eos/methane_butane_flash_against_neqsim.json`.
        common::assert_close(r.tc.value, tc, 2e-4, &format!("{name}: Tc"));
        common::assert_close(r.pc.value, pc, 2e-4, &format!("{name}: Pc"));
        assert!(
            (r.z_c - expected_z_c).abs() < 5e-6,
            "{name}: Z_c is {} against the analytic {expected_z_c}",
            r.z_c
        );
    }
}

/// A mixture's `Z_c` varies with composition, and that is the whole discriminating test.
///
/// The mechanical conditions - solving `dP/dV = d2P/dV2 = 0` at fixed composition -
/// reproduce a pure component's critical point exactly and return `(1 - omega_b)/3`
/// for **every** mixture, because in reduced variables they have a single universal
/// root. A pure-component check cannot tell the two routes apart. This can: the
/// spread below is 0.14, against a constant.
#[test]
fn a_mixtures_critical_compressibility_varies_with_composition() {
    let mut seen = Vec::new();
    for methane_fraction in [0.2, 0.4, 0.6, 0.8] {
        let mixture = methane_butane();
        let r = critical_point(&mixture, &[methane_fraction, 1.0 - methane_fraction])
            .expect("a critical point");
        seen.push(r.z_c);
    }
    let spread = seen.iter().cloned().fold(f64::MIN, f64::max)
        - seen.iter().cloned().fold(f64::MAX, f64::min);
    assert!(
        spread > 0.1,
        "Z_c varied by only {spread} across the composition range: {seen:?}. A spread \
         near zero is what the mechanical conditions give, and is not a mixture \
         critical point."
    );
}

/// The `Tc` locus of a binary falls between its pure endpoints, and monotonically.
///
/// A shape check rather than a value check, and it is the one that would catch an
/// iteration that converges on the wrong root: the critical locus of a binary is a
/// continuous curve from one pure component to the other, so a point outside the
/// bracket or out of order is wrong whatever it is near.
#[test]
fn a_binarys_critical_locus_falls_between_its_pure_endpoints() {
    let mixture = methane_butane();

    let light = databank::entry("methane", None).expect("methane").tc;
    let heavy = databank::entry("n-butane", None).expect("n-butane").tc;
    let mut previous = heavy;
    for methane_fraction in [0.2, 0.4, 0.6, 0.8] {
        let r = critical_point(&mixture, &[methane_fraction, 1.0 - methane_fraction])
            .expect("a critical point");
        assert!(
            r.tc.value > light && r.tc.value < heavy,
            "Tc = {} at z = {methane_fraction} is outside the pure endpoints",
            r.tc.value
        );
        assert!(
            r.tc.value < previous,
            "Tc = {} at z = {methane_fraction} did not fall below {previous}",
            r.tc.value
        );
        previous = r.tc.value;
    }
}

/// The Jacobi solver, against eigenvalues that are known without it.
#[test]
fn the_eigensolver_agrees_with_the_closed_form() {
    // [[2, 1], [1, 2]] has eigenvalues 1 and 3, and eigenvectors along (1, -1) and
    // (1, 1).
    let e = symmetric_eigen(&[vec![2.0, 1.0], vec![1.0, 2.0]], 1e-15, 100);
    common::assert_close(e.values[0], 1.0, 1e-14, "the smaller eigenvalue");
    common::assert_close(e.values[1], 3.0, 1e-14, "the larger eigenvalue");
    let root_half = std::f64::consts::FRAC_1_SQRT_2;
    common::assert_close(e.vectors[0][0].abs(), root_half, 1e-14, "the eigenvector");

    // A diagonal matrix must come back exactly, and ascending.
    let e = symmetric_eigen(&[vec![5.0, 0.0], vec![0.0, 3.0]], 1e-15, 100);
    assert_eq!(e.values, vec![3.0, 5.0]);

    // An already-diagonal 1x1 is its own answer.
    let e = symmetric_eigen(&[vec![-2.0]], 1e-15, 100);
    assert_eq!(e.values, vec![-2.0]);
    assert_eq!(e.vectors, vec![vec![1.0]]);
}

/// A composition that is not a composition is refused rather than renormalised.
#[test]
fn a_malformed_composition_is_refused() {
    let mixture = methane_butane();
    for (label, z) in [
        ("too short", vec![0.5]),
        ("too long", vec![0.5, 0.25, 0.25]),
        ("negative", vec![1.5, -0.5]),
        ("does not sum to one", vec![0.5, 0.6]),
    ] {
        let error = critical_point(&mixture, &z).expect_err(label);
        assert!(
            matches!(error, AzothError::InvalidInput { .. }),
            "{label}: expected an invalid-input error"
        );
    }
}

/// The result is clean at an ordinary composition, and reports what the spec says.
#[test]
fn the_result_carries_the_declared_fields() {
    let mixture = methane_butane();
    let r = critical_point(&mixture, &[0.4, 0.6]).expect("a critical point");
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
    assert_eq!(
        <azoth_eos::CriticalPointResult as CalcResult>::CALC_ID,
        MODEL_ID
    );
    // `Z_c` is `Pc Vc/(R Tc)` by definition, so the two ways of forming it agree.
    let recomputed = r.pc.value * r.vc.value / (8.31446261815324 * r.tc.value);
    common::assert_close(recomputed, r.z_c, 1e-12, "Z_c from the returned state");
}
