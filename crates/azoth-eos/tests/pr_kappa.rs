//! Spec-driven tests for `eos.pr_kappa`.

use azoth_core::CalcResult;
use azoth_eos::pr_kappa;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.pr_kappa";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::PrKappaResult {
    pr_kappa(common::input(case, "omega"))
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
                common::assert_close(
                    result.kappa,
                    common::expected(case, "kappa"),
                    case.tolerance,
                    &format!("{}::{} (kappa)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "omega" => Some(common::input(case, "omega")),
                        "kappa" => Some(result.kappa),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
            }
            "property" => match case.property {
                Some("monotonic") => monotonic(),
                other => panic!("{}::{}: unhandled property {other:?}", spec.id, case.id),
            },
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 3,
        "expected several active cases, ran {executed}"
    );
}

/// `kappa` must increase with `omega`.
///
/// Its derivative, `1.54226 - 0.53984*omega`, is positive for every `omega` below
/// 2.856883521043272, so over the whole range a real fluid occupies - roughly
/// -0.4 to 1.1 - the coefficient is strictly increasing in the acentric factor.
///
/// This is what a sign error in either coefficient breaks, and neither point-value
/// case would necessarily catch it: a wrong formula can still be right at the two
/// sets of inputs the spec happens to name.
#[test]
fn monotonic() {
    let at = |omega: f64| pr_kappa(omega).unwrap().kappa;

    // Walk the acentric factor across the range real fluids occupy, rather than
    // asserting one adjacent pair - a coefficient wrong in its curvature can be
    // locally increasing and globally wrong.
    let samples: Vec<f64> = (0..=40).map(|i| -0.4 + 1.5 * f64::from(i) / 40.0).collect();
    for pair in samples.windows(2) {
        let (low, high) = (pair[0], pair[1]);
        assert!(
            at(high) > at(low),
            "kappa must increase with omega: kappa({high}) = {} is not greater than \
             kappa({low}) = {}",
            at(high),
            at(low)
        );
    }
}

#[test]
fn the_worked_example_is_exact_in_binary() {
    // The spec claims every intermediate is the correctly rounded double for its
    // exact value, which is what justifies a 1e-12 tolerance on an arithmetic
    // whose inputs are decimal. Assert the end of that claim directly: the f64
    // result must equal the f64 nearest the exact decimal 0.60282728832.
    //
    // If this ever fails, the spec's derivation claim has become false and the
    // tolerance is no longer justified by it - so the failure is about
    // documentation, not about a wrong number.
    let r = pr_kappa(0.152).unwrap();
    assert_eq!(
        r.kappa.to_bits(),
        0.60282728832_f64.to_bits(),
        "the worked example is no longer bit-exact; the spec's derivation needs revising"
    );
}

#[test]
fn a_negative_kappa_is_returned_with_a_warning_rather_than_refused() {
    // Helium's acentric factor is below the polynomial's root at -0.233383499424,
    // so the coefficient comes back negative. That is a real answer to a real
    // question - "what does this correlation give for helium" - and the honest
    // response is to say the value is out of range, not to refuse to compute it.
    let r = pr_kappa(-0.385).unwrap();
    assert!((r.kappa + 0.259138992).abs() < 1e-12, "got {}", r.kappa);
    assert!(
        !r.is_clean(),
        "a negative kappa means the alpha function is not the one PR models, and \
         the caller has to be able to see that without inferring it from the sign"
    );
}

#[test]
fn a_positive_kappa_at_the_same_magnitude_is_silent() {
    // The bound has to discriminate rather than merely fire near zero: hydrogen
    // sits just the other side of the root and must come back clean.
    let r = pr_kappa(-0.216).unwrap();
    assert!((r.kappa - 0.02891845248).abs() < 1e-12, "got {}", r.kappa);
    assert!(
        r.is_clean(),
        "hydrogen's kappa is positive and in range; got {:?}",
        r.warnings
    );
}
