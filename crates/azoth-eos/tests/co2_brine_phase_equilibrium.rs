//! Spec-driven tests for the `eos.co2_brine_phase_equilibrium` model.

use azoth_core::AzothError;
use azoth_core::units::{kelvins, pascals};
use azoth_eos::co2_brine_phase_equilibrium::co2_brine_phase_equilibrium;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.co2_brine_phase_equilibrium";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::Co2BrinePhaseEquilibriumResult {
    let components: Vec<String> = case
        .list("components")
        .expect("components")
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    co2_brine_phase_equilibrium(
        &components,
        kelvins(common::input(case, "T")),
        pascals(common::input(case, "P")),
        case.vector("moles").expect("moles"),
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
        common::assert_close(
            f64::from(result.phase_count),
            case.expected_value("phase_count").expect("phase_count"),
            case.tolerance,
            &format!("{context} (phase_count)"),
        );
        // **The root, which a split on the wrong one can still hide.** Both trials are run
        // and the lower-Gibbs one wins, so a port that kept whichever split converged first
        // closes its own balances and lands on a defensible pair of fractions - the numbers
        // beside this cannot tell the two apart, and this can.
        //
        // A case may leave it unstated, and the one that does is the 200 bara state, where
        // the two roots are the same root. Its `source` gives the measurement; the check is
        // not skipped, it is not applicable there.
        if let Some(wanted) = case.expected_string("accepted_root") {
            assert_eq!(
                result.accepted_root, wanted,
                "{context}: the root the first phase was solved on"
            );
        }
        for (i, (got, want)) in result
            .phase_fractions
            .iter()
            .zip(
                case.expected_vector("phase_fractions")
                    .expect("phase_fractions"),
            )
            .enumerate()
        {
            common::assert_close(
                *got,
                *want,
                case.tolerance,
                &format!("{context} (phase_fractions[{i}])"),
            );
        }
        // The compositions are a matrix, which `TestCase` carries for an input but not for an
        // expectation - the Python runner reads the case file directly and compares the rows,
        // and `test_cross_impl` holds the two kernels to each other there. What is asserted
        // *here* is the shape - one row per phase, one entry per component - so that a result
        // that dropped a phase cannot pass on the scalars alone.
        let components = case.list("components").expect("components").len();
        assert_eq!(
            result.x.len(),
            result.phase_count as usize,
            "{context}: one composition row per phase"
        );
        for row in &result.x {
            assert_eq!(
                row.len(),
                components,
                "{context}: one composition entry per component"
            );
        }
        common::assert_consistent(&result, context);
    }
}

/// **A feed outside the class's scope is refused rather than answered.**
///
/// `hasSupportedComposition` asks for CO2, water and ions alone with more water than CO2 - a
/// five-component brine with a hydrocarbon in it is a different question, and NeqSim throws
/// rather than running. The three refusals below are the three ways that is reached.
#[test]
fn a_feed_outside_the_scope_is_refused() {
    let brine = |names: &[&str], moles: &[f64]| {
        let components: Vec<String> = names.iter().map(|name| (*name).to_string()).collect();
        co2_brine_phase_equilibrium(&components, kelvins(283.15), pascals(4.0e6), moles)
    };

    // A substance that is neither CO2, water nor an ion.
    let error = brine(&["CO2", "water", "methane"], &[1.0, 10.0, 0.01])
        .expect_err("methane is not in the class's scope");
    assert!(
        matches!(error, AzothError::InvalidInput { .. }),
        "{error:?}"
    );
    assert_eq!(error.field(), Some("components"));

    // CO2 in excess of the water: a gas with a little water in it, not a brine.
    let error =
        brine(&["CO2", "water"], &[10.0, 1.0]).expect_err("water must be in excess of the CO2");
    assert!(
        matches!(error, AzothError::InvalidInput { .. }),
        "{error:?}"
    );

    // No water at all.
    let error = brine(&["CO2", "Na+"], &[1.0, 1.0]).expect_err("there is no water");
    assert!(
        matches!(error, AzothError::InvalidInput { .. }),
        "{error:?}"
    );
}
