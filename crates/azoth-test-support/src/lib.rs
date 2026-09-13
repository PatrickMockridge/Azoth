//! Shared helpers for the spec-driven integration tests, across every namespace.
//!
//! The tests are not hand-written per calc. They walk the spec's `tests` list and
//! execute whatever it declares, so a test added to a spec YAML runs in both
//! languages with no new test code. That is the whole point of generating
//! `spec_gen.rs`: the specs are the test plan.
//!
//! # Why this is a crate and not a `tests/common/mod.rs`
//!
//! It was `crates/azoth-hydraulics/tests/common/mod.rs`, which works for one
//! namespace and fails for two: Rust has no way to share a test module between
//! crates, so the second namespace would have copied the file. Two copies of
//! `assert_warnings_agree_with_spec` is the worst possible thing to duplicate -
//! it is the check that a warning fires exactly when a spec bound says it should,
//! and a copy that drifted would silently weaken the guarantee in one namespace
//! while the other stayed green.
//!
//! The only thing in here that ever knew about a namespace was the spec lookup.
//! That is now a parameter: each test binary passes its own `spec_gen::specs()`,
//! which is the one line of a spec-driven test that is namespace-specific.

use azoth_core::spec::{CalcSpec, TestCase};
use azoth_core::{CalcResult, Warning, WarningCode};

/// Fetch a spec from a namespace's generated tables, failing loudly if the id is
/// wrong. A missing spec means the test file and the registry disagree, which is
/// a test bug, not a runtime condition.
///
/// The tables are passed in rather than looked up here because this crate cannot
/// depend on the namespace crates - they depend on it.
#[track_caller]
pub fn spec(all: &'static [&'static CalcSpec], id: &str) -> &'static CalcSpec {
    all.iter().copied().find(|s| s.id == id).unwrap_or_else(|| {
        let known: Vec<&str> = all.iter().map(|s| s.id).collect();
        panic!("no spec for `{id}`; registry has {known:?}")
    })
}

/// Relative comparison with a readable failure message.
///
/// Relative rather than absolute, because the quantities here span nine orders
/// of magnitude (Re ~ 1e5, f ~ 1e-2) and a single absolute tolerance would be
/// meaningless at both ends.
#[track_caller]
pub fn assert_close(actual: f64, expected: f64, tolerance: f64, context: &str) {
    assert!(
        actual.is_finite(),
        "{context}: got {actual}, which is not finite (expected {expected})"
    );
    let scale = expected.abs().max(f64::MIN_POSITIVE);
    let rel = (actual - expected).abs() / scale;
    assert!(
        rel <= tolerance,
        "{context}: got {actual}, expected {expected}\n  \
         relative error {rel:.3e} exceeds tolerance {tolerance:.3e}"
    );
}

/// A scalar input from a test case, with a clear message if the spec omitted it.
#[track_caller]
pub fn input(case: &TestCase, name: &str) -> f64 {
    case.input(name).unwrap_or_else(|| {
        panic!(
            "test `{}` does not supply input `{name}`; it has {:?}",
            case.id,
            case.numbers.iter().map(|(k, _)| *k).collect::<Vec<_>>()
        )
    })
}

/// An expected output from a test case.
#[track_caller]
pub fn expected(case: &TestCase, name: &str) -> f64 {
    case.expected_value(name).unwrap_or_else(|| {
        panic!(
            "test `{}` does not assert `{name}`; it asserts {:?}",
            case.id,
            case.expected.iter().map(|(k, _)| *k).collect::<Vec<_>>()
        )
    })
}

/// A list-valued input from a test case.
#[track_caller]
pub fn list_input(case: &TestCase, name: &str) -> &'static [&'static str] {
    case.list(name)
        .unwrap_or_else(|| panic!("test `{}` does not supply list `{name}`", case.id))
}

/// Assert that a skipped test says why it is skipped.
///
/// A skip with no reason is indistinguishable from an oversight, and the whole
/// convention of marking a source `TODO: source needed` depends on the reason
/// travelling with the test.
#[track_caller]
pub fn assert_skips_are_explained(spec: &CalcSpec) {
    for case in spec.all_tests() {
        if !case.is_active() {
            assert!(
                case.skip_reason.is_some_and(|r| !r.trim().is_empty()),
                "{}::{} is skipped with no reason",
                spec.id,
                case.id
            );
        }
    }
}

/// Assert the warnings a result carries agree exactly with the spec's own range
/// checks.
///
/// This is the strongest available check on the warning mechanism, and it is
/// fully generic: it reads the bounds out of the spec and confirms that the
/// implementation warned precisely when a warning-severity bound was violated -
/// no more, and no fewer.
///
/// It catches the two failure modes that matter. Warning when nothing is wrong
/// is warning fatigue, which teaches callers to ignore the field. Failing to
/// warn when something is wrong is the dangerous one: an out-of-range result
/// indistinguishable from a validated one.
#[track_caller]
pub fn assert_warnings_agree_with_spec(
    spec: &CalcSpec,
    warnings: &[Warning],
    resolve: impl Fn(&str) -> Option<f64>,
    context: &str,
) {
    use azoth_core::Severity;
    use std::collections::{BTreeMap, BTreeSet};

    let has = |code: WarningCode, field: &str| {
        warnings
            .iter()
            .any(|w| w.code == code && w.field.as_deref() == Some(field))
    };

    // Checks are grouped by (quantity, code), not examined one at a time. A
    // spec may legitimately carry several bounds on one quantity sharing a code
    // - swamee_jain has `relative_roughness < 0` as an error and
    // `relative_roughness <= 1e-6` as a warning, both OUT_OF_VALID_RANGE. Since a
    // warning is identified by its code and field, "bound A was satisfied" only
    // implies "no warning" if no other bound in the group fired. Getting that
    // wrong makes this helper report false failures, which is exactly what
    // happened before it was grouped.
    let mut expected: BTreeMap<(&str, WarningCode), bool> = BTreeMap::new();
    let mut unresolvable: BTreeSet<&str> = BTreeSet::new();

    for check in spec.all_checks() {
        match resolve(check.quantity) {
            None => {
                unresolvable.insert(check.quantity);
            }
            Some(value) => {
                let fired = check.severity == Severity::Warning && check.violated(value);
                let entry = expected
                    .entry((check.quantity, check.code))
                    .or_insert(false);
                *entry = *entry || fired;
            }
        }
    }

    for quantity in unresolvable {
        assert!(
            has(WarningCode::RangeCheckSkipped, quantity),
            "{context}: `{quantity}` could not be resolved, so its check did not run, \
             but no RANGE_CHECK_SKIPPED warning was emitted - an unchecked value must \
             not look like a checked one"
        );
    }

    for ((quantity, code), should_warn) in expected {
        let warned = has(code, quantity);
        if should_warn {
            assert!(
                warned,
                "{context}: `{quantity}` violates a warning-severity spec bound yet no \
                 {code} warning was emitted"
            );
        } else {
            assert!(
                !warned,
                "{context}: `{quantity}` satisfies every warning-severity spec bound for \
                 {code}, yet a {code} warning was emitted"
            );
        }
    }
}

/// Assert a result is internally consistent for reporting purposes.
///
/// Cheap invariant, but it catches a class of bug where a warning list is
/// rebuilt by hand at the end of a calculation and loses its earlier entries.
#[track_caller]
pub fn assert_consistent<T: CalcResult>(result: &T, context: &str) {
    assert_eq!(
        result.is_clean(),
        result.warnings().is_empty(),
        "{context}: is_clean disagrees with the warning list"
    );
    for w in result.warnings() {
        assert!(
            !w.message.trim().is_empty(),
            "{context}: warning {:?} has an empty message",
            w.code
        );
    }
}

/// Run the `unit_round_trip` property: the same physical state described in two
/// unit systems must give bit-identical results.
///
/// This is the sharpest available check that unit handling is correct, because
/// any accidental conversion factor shows up immediately as a discrepancy. The
/// caller supplies both computations; the helper only compares them.
#[track_caller]
pub fn assert_same_state(si_value: f64, other_unit_value: f64, tolerance: f64, context: &str) {
    assert_close(other_unit_value, si_value, tolerance, context);
}
