//! Spec-driven tests for `eos.pr_z_factor`.

use azoth_core::AzothError;
use azoth_eos::spec_gen;
use azoth_eos::{RootStructure, pr_z_factor};
use azoth_test_support as common;

const CALC_ID: &str = "eos.pr_z_factor";

/// The cubic the calc solves, as `f(z)`, from the two parameters.
///
/// Written out here rather than reached for from the crate so the test can check
/// the implementation against the *published* equation rather than against
/// itself. A shared helper would make a wrong coefficient in the crate wrong in
/// the test too, which is the failure this file exists to catch.
fn cubic(a: f64, b: f64, z: f64) -> f64 {
    let (c2, c1, c0) = (
        -(1.0 - b),
        a - 3.0 * b * b - 2.0 * b,
        -(a * b - b * b - b * b * b),
    );
    ((z + c2) * z + c1) * z + c0
}

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::PrZFactorResult {
    pr_z_factor(
        common::input(case, "a_reduced"),
        common::input(case, "b_reduced"),
    )
    .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    common::assert_skips_are_explained(spec);

    let mut executed = 0;
    for case in spec.all_tests() {
        if !case.is_active() {
            continue;
        }
        match case.kind {
            "worked_example" | "reference" => {
                let result = call(case);
                let a = common::input(case, "a_reduced");
                let b = common::input(case, "b_reduced");
                for field in ["z_min", "z_max"] {
                    let actual = if field == "z_min" {
                        result.z_min
                    } else {
                        result.z_max
                    };
                    common::assert_close(
                        actual,
                        common::expected(case, field),
                        case.tolerance,
                        &format!("{}::{} ({field})", spec.id, case.id),
                    );
                    // The residual check comes free with every case and is the one
                    // check that does not depend on the expected value being right -
                    // which matters most at the critical point, where the spec
                    // asserts loosely because the answer genuinely is not determined.
                    let residual = cubic(a, b, actual);
                    assert!(
                        residual.abs() < 1e-12,
                        "{}::{}: {field} = {actual} does not satisfy the cubic: f = {residual:e}",
                        spec.id,
                        case.id
                    );
                }
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "a_reduced" => Some(a),
                        "b_reduced" => Some(b),
                        "z_min" => Some(result.z_min),
                        "z_max" => Some(result.z_max),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
            }
            "property" => match case.property {
                // Both consistency tests and the monotonic one are hand-written
                // below; the spec's `note` on each says what it checks.
                Some("consistency_with") | Some("monotonic") => {}
                other => panic!("{}::{}: unhandled property {other:?}", spec.id, case.id),
            },
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 4,
        "expected several active cases, ran {executed}"
    );
}

/// Every returned root satisfies the cubic, across a sweep of parameters.
///
/// The check that does not need to know the answer, which is what makes it usable
/// where the answer is not determined. A polish that has been removed, an analytic
/// branch taken on the wrong side of the discriminant, or a root returned
/// unpolished all produce plausible values that fail here.
#[test]
fn every_returned_root_satisfies_the_cubic() {
    for a in [0.05, 0.1, 0.20206500174625697, 0.3, 0.4572355289213822] {
        for b in [0.005, 0.012, 0.02431127309496514, 0.05, 0.07779607390388846] {
            let r = pr_z_factor(a, b).unwrap();
            for z in [r.z_min, r.z_max] {
                let residual = cubic(a, b, z);
                assert!(
                    residual.abs() < 1e-12,
                    "A = {a}, B = {b}: z = {z} leaves f = {residual:e}"
                );
            }
        }
    }
}

/// When three roots are admissible, the middle one is not among the returned pair.
///
/// `z_min` and `z_max` must be the *outermost* two, which means the polynomial has
/// exactly one further root strictly between them. Counted by sign change on a
/// grid, because the endpoints are themselves roots and a sign test at them would
/// be comparing rounding noise to zero.
///
/// An implementation returning the two *smallest* roots, or the two it happened to
/// find, passes the residual check above and fails here.
#[test]
fn the_middle_root_is_not_returned() {
    let (a, b) = (0.20206500174625697, 0.02431127309496514);
    let r = pr_z_factor(a, b).unwrap();
    assert_eq!(r.root_structure, RootStructure::Three);

    // The cubic is -,+,-,+ across r1, r2, r3, so crossing the middle root is the
    // only sign change strictly inside (z_min, z_max). Two would mean the outer
    // pair is not the outer pair; zero would mean there is no third root at all.
    let steps = 4000;
    let mut changes = 0;
    let mut previous = cubic(a, b, r.z_min + (r.z_max - r.z_min) * 0.5 / steps as f64);
    for i in 2..steps {
        let fraction = i as f64 / steps as f64;
        let z = r.z_min + (r.z_max - r.z_min) * fraction;
        let value = cubic(a, b, z);
        if previous * value < 0.0 {
            changes += 1;
        }
        previous = value;
    }
    assert_eq!(
        changes, 1,
        "exactly one root should lie strictly between z_min and z_max, found {changes}"
    );
}

/// The admissible root count is one or three, and never two.
///
/// The test for the theorem in `RootStructure`: `f(b_reduced)` is exactly
/// `-2*b_reduced**2`, so `B` lies below all three roots or between the middle and
/// largest ones, and there is no in-between case.
///
/// Swept rather than argued, because the argument rests on the cubic's coefficients
/// being exactly what the spec says. A sign error in the constant term would break
/// `f(B) = -2*B**2` and could make a count of two reachable, at which point the enum
/// would be missing a value it needs - and that should fail here rather than surface
/// as a caller's `AttributeError` or a silently mislabelled structure.
#[test]
fn the_admissible_count_is_never_two() {
    let mut seen_one = 0_u32;
    let mut seen_three = 0_u32;
    for a_step in 0..=120 {
        for b_step in 1..=98 {
            let a = 0.6 * f64::from(a_step) / 120.0;
            let b = 0.002 + 0.098 * f64::from(b_step) / 98.0;
            let r = pr_z_factor(a, b).unwrap();
            match r.root_structure {
                RootStructure::One => seen_one += 1,
                RootStructure::Three => seen_three += 1,
            }
            // The theorem, stated directly on the polynomial as well: f(B) is
            // negative, which is what forces the count.
            assert!(
                cubic(a, b, b) < 0.0,
                "f(B) should be -2*B**2 < 0 for A = {a}, B = {b}"
            );
        }
    }
    assert!(
        seen_one > 0 && seen_three > 0,
        "the sweep should reach both cases"
    );
}

/// One admissible root means the two outputs are equal, and three means they are not.
#[test]
fn root_structure_matches_what_the_outputs_say() {
    let cases = [
        (
            0.20206500174625697,
            0.02431127309496514,
            RootStructure::Three,
        ),
        (0.25, 0.03125, RootStructure::Three),
        (
            0.08448417260831159,
            0.025932024634629486,
            RootStructure::One,
        ),
        (0.4572355289213822, 0.07779607390388846, RootStructure::One),
    ];
    for (a, b, expected) in cases {
        let r = pr_z_factor(a, b).unwrap();
        assert_eq!(r.root_structure, expected, "for A = {a}, B = {b}");
        if expected == RootStructure::One {
            assert_eq!(
                r.z_min.to_bits(),
                r.z_max.to_bits(),
                "a single admissible root must be reported as one, not as two equal ones"
            );
        } else {
            assert!(r.z_max > r.z_min, "for A = {a}, B = {b}");
        }
    }
}

/// Both extreme roots fall as `A` rises and rise as `B` rises.
///
/// More attraction pulls the fluid together, more repulsion pushes it apart. A
/// sign error in either coefficient of the cubic, or a transposed one, cannot
/// survive both directions.
#[test]
fn monotonic() {
    let z = |a: f64, b: f64| pr_z_factor(a, b).unwrap();

    for pair in [(0.10, 0.15), (0.15, 0.25), (0.25, 0.30)] {
        let (low, high) = pair;
        assert!(
            z(high, 0.02431127309496514).z_max < z(low, 0.02431127309496514).z_max,
            "z_max must fall as A rises: {low} -> {high}"
        );
    }

    for pair in [(0.005, 0.012), (0.012, 0.024), (0.024, 0.06)] {
        let (low, high) = pair;
        assert!(
            z(0.20206500174625697, high).z_max > z(0.20206500174625697, low).z_max,
            "z_max must rise as B rises: {low} -> {high}"
        );
    }
}

#[test]
fn a_non_positive_b_reduced_is_an_error() {
    // B is proportional to pressure and vanishes only at zero pressure, where the
    // cubic degenerates to z**2 (z - 1) - a trivial double root at zero, which
    // would come back as a phase split with nothing to say it was not one.
    for b in [0.0, -0.01] {
        let err = pr_z_factor(0.2, b).unwrap_err();
        assert!(matches!(err, AzothError::OutOfRange { .. }));
        assert_eq!(err.field(), Some("b_reduced"), "for B = {b}");
    }
}

#[test]
fn a_negative_a_reduced_is_an_error() {
    // A is Omega_a * alpha * Pr / Tr**2, a squared term times a positive
    // coefficient, so no Peng-Robinson state produces a negative one.
    let err = pr_z_factor(-0.01, 0.024).unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("a_reduced"));
}

#[test]
fn zero_a_reduced_is_the_hard_sphere_limit_and_is_allowed() {
    // The bound is inclusive at zero on purpose: it is a real thing to ask about,
    // and refusing it would be refusing a limit rather than a mistake.
    let b = 0.02431127309496514;
    let r = pr_z_factor(0.0, b).unwrap();

    // With no attraction the cubic has roots -0.0587, 0.0101 and 1.0243, and the
    // first two are at or below B - so they are discarded as unphysical and only
    // the third survives. That is worth asserting precisely because it is
    // counter-intuitive: "three real roots" does not mean "three admissible ones",
    // and this is the case that shows the difference.
    assert_eq!(r.root_structure, RootStructure::One);
    assert!(r.z_max > b, "the surviving root must be above B");
    assert!(cubic(0.0, b, r.z_max).abs() < 1e-12);

    // And the two discarded roots really are below B, so the filter is doing
    // something rather than being vacuous here.
    assert!(cubic(0.0, b, 0.0) > 0.0, "a root lies between 0 and B");
    assert!(cubic(0.0, b, b) < 0.0, "and B is not one of them");
}

/// At the critical point the answer is not determined, and the test says how little.
///
/// This is the measurement behind the spec's `assumptions` entry. The polynomial
/// stays within `1e-12` of zero across a window about `2e-4` wide in `z`, so no
/// implementation can pin the last four digits - and rather than assert a value it
/// cannot know, this asserts the *window*, which it can.
#[test]
fn the_critical_point_is_only_determined_to_about_one_part_in_ten_thousand() {
    // The pair that places a triple root at the cubic's own critical point, which is
    // what this test is about. **Not the pair the crate ships**: NeqSim's literals are
    // 5.6e-6 off it, deliberately, so `(1 - B)/3` is no longer the triple root of the
    // cubic `OMEGA_A`/`OMEGA_B` build - it is 4.8e-5 away, which is thirty times this
    // test's tolerance. See `eos.pr_alpha_ab`'s assumptions.
    let (a, b) = (0.4572355289213822, 0.07779607390388846);
    // c2 = -(1 - B), and the triple root of `w**3 + c2*w**2 + c1*w + c0` sits at
    // -c2/3 = (1 - B)/3 when the constants are the ones that place it there.
    let analytic_triple_root = (1.0 - b) / 3.0;

    let r = pr_z_factor(a, b).unwrap();
    assert!(
        (r.z_max - analytic_triple_root).abs() < 1e-4,
        "got {}, expected to be near {analytic_triple_root}",
        r.z_max
    );

    // And the flatness itself, so the tolerance above is justified by a
    // measurement rather than by caution: a quarter of the way across the window
    // the polynomial is still indistinguishable from zero.
    let offset = 5e-5;
    assert!(
        cubic(a, b, analytic_triple_root + offset).abs() < 1e-12,
        "the cubic should be flat near the critical point"
    );
}
