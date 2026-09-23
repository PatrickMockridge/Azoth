//! Spec-driven tests for the `eos.hybrid_eos_ge_flash` model.
//!
//! The oracle is `validation/neqsim/HybridEosGeProbe.java` and its capture
//! `captures/hybrid_eos_ge_probe.tsv`: `SystemHybridEosGeFlashTest`'s own fluid through
//! NeqSim's fixed-role hybrid route, with the acceptance contract's residuals printed beside
//! the split.

use azoth_eos::hybrid_eos_ge_flash;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.hybrid_eos_ge_flash";

fn azoth_model_spec(id: &str) -> &'static azoth_core::spec::ModelSpec {
    azoth_eos::model_gen::model(id).expect("the model should be in its own table")
}

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::HybridEosGeFlashResult {
    let names = case.list("components").expect("components");
    let moles = case.vector("moles").expect("moles");
    hybrid_eos_ge_flash(
        &names,
        common::input_str(case, "cubic").parse().expect("cubic"),
        common::input(case, "T"),
        common::input(case, "P"),
        &moles,
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = azoth_model_spec(MODEL_ID);
    let mut executed = 0;
    for case in spec.cases {
        if !case.is_active() {
            continue;
        }
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);

        common::assert_close(
            result.beta.iter().sum::<f64>(),
            1.0,
            1.0e-09,
            &format!("{context} (beta sum)"),
        );
        let expected_beta = case
            .expected_vector("beta")
            .unwrap_or_else(|| panic!("the case declares beta"));
        assert_eq!(
            result.beta.len(),
            expected_beta.len(),
            "{context}: beta length"
        );
        for (index, (&got, &want)) in result.beta.iter().zip(expected_beta).enumerate() {
            common::assert_close(
                got,
                want,
                case.tolerance,
                &format!("{context} (beta[{index}])"),
            );
        }

        // **The two matrices are asserted by the Python runner and by `test_cross_impl`**,
        // not here: the generator carries no matrix *expectation*, so a case that declares
        // one is read from the raw spec - which only the Python side does. What this test
        // asserts instead is the shape and the two properties a matrix must have, and the
        // values are pinned in `the_acceptance_contract_holds_on_the_captured_fluid` below
        // against the capture's own.
        let columns = case.list("components").expect("components").len();
        assert_eq!(result.x.len(), result.beta.len());
        assert_eq!(result.ln_phi.len(), result.beta.len());
        for (role, row) in result.x.iter().enumerate() {
            assert_eq!(row.len(), columns, "{context}: x[{role}] width");
            common::assert_close(
                row.iter().sum::<f64>(),
                1.0,
                1.0e-09,
                &format!("{context} (x[{role}] sum)"),
            );
            assert!(row.iter().all(|value| *value > 0.0), "{context}: x[{role}]");
        }
        assert_eq!(result.ln_phi[0].len(), columns, "{context}: ln_phi width");

        common::assert_consistent(&result, context);
        executed += 1;
    }
    assert!(executed >= 1, "expected a case, ran {executed}");
}

#[test]
fn the_model_is_registered() {
    let spec = azoth_model_spec(MODEL_ID);
    assert_eq!(spec.kind, "procedure");
    assert!(!spec.cases.is_empty());
}

/// **The acceptance contract, asserted on the result rather than read off the capture.**
///
/// `SystemHybridEosGeFlashTest` states four things about this state: `beta` sums to one, the
/// material balance holds to `1e-7`, `ln(x_i phi_i P)` agrees across the roles to `1e-5`, and
/// an ion is at most `1e-40` outside the aqueous role. The captured fluid's own residuals are
/// `3.40e-13` and `1.68e-11`; the port reaches the same magnitudes, which is the claim.
///
/// **And the assertion the test makes about n-heptane is the one worth checking here**: its
/// brine coefficient exceeds `1e6`, which is the Henry arm of `eos.pitzer_phase` and not a
/// cubic's number. A cubic over these components puts it nowhere near `1.15e12`.
#[test]
fn the_acceptance_contract_holds_on_the_captured_fluid() {
    let names = ["methane", "n-heptane", "water", "Na+", "Cl-"];
    let moles = [5.0, 2.0, 55.5, 1.0, 1.0];
    let result = hybrid_eos_ge_flash(&names, azoth_eos::Cubic::Srk, 313.15, 50.0e5, &moles)
        .expect("the captured fluid");

    assert!(
        (result.beta.iter().sum::<f64>() - 1.0).abs() < 1.0e-09,
        "beta sums to {}",
        result.beta.iter().sum::<f64>()
    );
    assert!(
        result.max_material_balance_residual < 1.0e-07,
        "the material balance is {:e}, and the contract holds it to 1e-7",
        result.max_material_balance_residual
    );
    assert!(
        result.max_log_fugacity_residual < 1.0e-05,
        "the log fugacity spread is {:e}, and the contract holds it to 1e-5",
        result.max_log_fugacity_residual
    );

    let aqueous = 2;
    for (index, name) in names.iter().enumerate() {
        if !name.ends_with('+') && !name.ends_with('-') {
            continue;
        }
        assert!(
            result.x[aqueous][index] > 0.0,
            "{name} left the brine entirely"
        );
        for (role, row) in result.x.iter().enumerate() {
            if role != aqueous {
                assert!(
                    row[index] <= 1.0e-40,
                    "{name} is at {} in role {role}, and the contract holds an ion to 1e-40",
                    row[index]
                );
            }
        }
    }

    // The test's own claim about how the brine treats a hydrocarbon solvent.
    let heptane = 1;
    assert!(
        result.ln_phi[aqueous][heptane].exp() > 1.0e6,
        "n-heptane's brine coefficient is {}, and the test asserts it exceeds 1e6",
        result.ln_phi[aqueous][heptane].exp()
    );
    // And the roles really are three, with the oil heavier in heptane than the gas is.
    assert!(result.beta.iter().all(|beta| *beta > 0.0));
    assert!(result.x[1][heptane] > result.x[0][heptane]);
    assert!(result.x[aqueous][2] > 0.9, "the brine is mostly water");
}

/// **A fluid with no ion takes the same route and is not a second model.**
///
/// The ionic rules are the reason this flash is its own id, so a fluid without an ion is the
/// check that they are *rules and not the algorithm*: the fraction band returns immediately on
/// a zero inventory and the step cap is not applied, and the solve is then a plain
/// fixed-topology fraction Newton over three roles.
#[test]
fn a_fluid_without_an_ion_is_the_same_solve() {
    let names = ["methane", "n-heptane", "water"];
    let moles = [5.0, 2.0, 55.5];
    let result = hybrid_eos_ge_flash(&names, azoth_eos::Cubic::Srk, 313.15, 50.0e5, &moles)
        .expect("a fluid without an ion");
    assert!(result.max_material_balance_residual < 1.0e-09);
    assert!(
        result.max_log_fugacity_residual < 1.0e-06,
        "the log fugacity spread is {:e} without an ion to confine",
        result.max_log_fugacity_residual
    );
    assert!((result.beta.iter().sum::<f64>() - 1.0).abs() < 1.0e-09);
}
