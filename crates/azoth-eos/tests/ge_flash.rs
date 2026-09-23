//! Spec-driven tests for the `eos.ge_flash` model.
//!
//! The oracle is `validation/neqsim/VanLaarGammaPhiProbe.java` and its capture
//! `captures/van_laar_gamma_phi_probe.tsv`: the acid fluid's split, both phases' `ln phi`,
//! the K-values and the vapour root, plus the bubble and dew pressures the same route gives.

use azoth_core::spec::TestCase;
use azoth_core::units::{kelvins, pascals};
use azoth_eos::ge_flash::ge_flash;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.ge_flash";

fn call(case: &TestCase) -> azoth_eos::GeFlashResult {
    ge_flash(
        case.list("components").expect("components"),
        common::input_str(case, "cubic"),
        common::input_str(case, "liquid_model"),
        kelvins(common::input(case, "T")),
        pascals(common::input(case, "P")),
        case.vector("z").expect("z"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = azoth_eos::model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);

        for (name, actual) in [
            ("x", &result.x),
            ("y", &result.y),
            ("k", &result.k),
            ("ln_phi_liquid", &result.ln_phi_liquid),
            ("ln_phi_vapour", &result.ln_phi_vapour),
        ] {
            // A case pins the fields it is about and leaves the others to the case that
            // is about them: this fluid's liquid and its vapour are pinned separately,
            // because the system's tuning reaches one and not the other.
            let Some(expected) = case.expected_vector(name) else {
                continue;
            };
            for (i, (&got, &want)) in actual.iter().zip(expected).enumerate() {
                common::assert_close(
                    got,
                    want,
                    case.tolerance,
                    &format!("{context} ({name}[{i}])"),
                );
            }
        }
        for (name, got) in [
            ("z_vapour", result.z_vapour),
            ("iterations", f64::from(result.iterations)),
        ] {
            if let Some(want) = case.expected_value(name) {
                common::assert_close(got, want, case.tolerance, &format!("{context} ({name})"));
            }
        }
        if let (Some(beta), Some(want)) = (result.beta, case.expected_value("beta")) {
            common::assert_close(beta, want, case.tolerance, &format!("{context} (beta)"));
        }
        common::assert_consistent(&result, context);
    }
}

/// **The acid liquid's uncovered-species penalty decides the carrier gas's fate.** CO2 is
/// not one of the three substances the Van Laar acid model covers, so its liquid fugacity
/// coefficient is the model's `1e12` - the capture's `ln_phi_liquid[0]` is `ln(1e12)` to
/// the last digit - and its K-value is `1.0065493354718438e12`. A carrier gas in this
/// liquid is a vapour by construction, not by the flash.
#[test]
fn the_uncovered_carrier_gas_is_a_vapour_by_construction() {
    let spec = azoth_eos::model_gen::model(MODEL_ID).expect("the model is in its table");
    let case = &spec.cases[0];
    let result = call(case);
    assert!(
        (result.ln_phi_liquid[0] - 1.0e12_f64.ln()).abs() < 1.0e-12,
        "the penalty is exact: {}",
        result.ln_phi_liquid[0]
    );
    assert!(
        result.k[0] > 1.0e11,
        "and the K-value follows it: {}",
        result.k[0]
    );
    assert!(
        result.y[0] > 0.99,
        "so CO2 is in the vapour: {}",
        result.y[0]
    );
}

/// **The gap to NeqSim is real and bounded, and this is what says so.** The untuned loop
/// reproduces the acid fluid's liquid to `2.0e-5` worst and its vapour to `3.8e-3` worst on
/// the two species the system tunes - and *not* to the `1e-10` the loop converges its
/// K-values to,
/// because the system it is compared against overrides three hooks this model does not
/// carry: a fitted vapour fugacity coefficient for an activity component in a CO2-rich
/// vapour, a `0.45` cap on nitric acid's K, and a `0.25` (or `0.05`) damping of the update.
/// If this ever closes to the loop's own tolerance, one of those has been ported and the
/// case's reason has been superseded - and this test is what would say so.
#[test]
fn the_untuned_loop_is_close_but_not_equal_to_the_tuned_system() {
    let spec = azoth_eos::model_gen::model(MODEL_ID).expect("the model is in its table");
    let case = spec
        .cases
        .iter()
        .find(|case| case.id == "the_same_fluid_vapour_at_its_own_tolerance")
        .expect("the vapour case is in the spec");
    let result = call(case);
    let expected = case.expected_vector("y").expect("the vapour case states y");
    let worst = result
        .y
        .iter()
        .zip(expected)
        .map(|(got, want)| (got - want).abs() / want.abs().max(f64::MIN_POSITIVE))
        .fold(0.0_f64, f64::max);
    assert!(
        worst > 1.0e-4,
        "the tuned system's three hooks still account for {worst:.3e} of the difference"
    );
    assert!(
        worst < 1.0e-2,
        "and the untuned loop is still close: {worst:.3e}"
    );
}

/// A liquid this library has no phase model for is refused rather than defaulted.
#[test]
fn an_unknown_liquid_is_refused() {
    assert!(
        ge_flash(
            &["CO2", "water"],
            "srk",
            "uniquac",
            kelvins(273.15),
            pascals(1.0e5),
            &[0.5, 0.5]
        )
        .is_err(),
        "NeqSim has no UNIQUAC system, so this is not a liquid the route reaches"
    );
}
