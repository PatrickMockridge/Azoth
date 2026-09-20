//! Spec-driven tests for the `eos.kent_eisenberg_phase` model.

use azoth_core::AzothError;
use azoth_eos::kent_eisenberg_phase::kent_eisenberg_phase;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.kent_eisenberg_phase";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::KentEisenbergPhaseResult {
    kent_eisenberg_phase(
        case.list("components").expect("components"),
        common::input(case, "T"),
        common::input(case, "P"),
        case.vector("x").expect("x"),
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
        for name in ["gamma", "ln_gamma", "ln_phi"] {
            let actual = match name {
                "gamma" => &result.gamma,
                "ln_gamma" => &result.ln_gamma,
                _ => &result.ln_phi,
            };
            let expected = case
                .expected_vector(name)
                .unwrap_or_else(|| panic!("the case declares {name}"));
            for (i, (&got, &want)) in actual.iter().zip(expected).enumerate() {
                common::assert_close(
                    got,
                    want,
                    case.tolerance,
                    &format!("{context} ({name}[{i}])"),
                );
            }
        }
        common::assert_consistent(&result, context);
    }
}

/// **The activity coefficient is identically one**, which is the model's defining property
/// and not a computed column: `PhaseKentEisenberg.getActivityCoefficient` returns it.
///
/// Asserted on a phase whose *reference states* differ - a solvent, a neutral solute and an
/// ion at once - because if `gamma` were being computed from anything, that is where it
/// would show.
#[test]
fn the_activity_coefficient_is_one_for_every_branch() {
    let r = kent_eisenberg_phase(
        &["water", "Na+", "Cl-", "CO2"],
        313.15,
        500_000.0,
        &[0.89, 0.04, 0.04, 0.03],
    )
    .unwrap();
    assert_eq!(r.gamma, vec![1.0; 4]);
    assert_eq!(r.ln_gamma, vec![0.0; 4]);
}

/// The three branches, told apart by what they read.
///
/// An ion's coefficient is the constant `1e8` at every temperature, because that branch
/// takes no argument. The other two move, and they move differently - so a port that had
/// swapped `P0` for `H` would agree at no temperature, and one that had scaled `P` wrongly
/// would disagree at all three.
#[test]
fn the_three_branches_read_different_things() {
    let ion = 1.0e8f64.ln();
    let mut water = Vec::new();
    let mut carbon_dioxide = Vec::new();
    for temperature in [298.15, 313.15, 373.15] {
        let r = kent_eisenberg_phase(
            &["water", "Na+", "Cl-", "CO2"],
            temperature,
            500_000.0,
            &[0.89, 0.04, 0.04, 0.03],
        )
        .unwrap();
        assert!(
            (r.ln_phi[1] - ion).abs() < 1.0e-12,
            "an ion's ln_phi is {} at {temperature} K, not ln(1e8) = {ion}",
            r.ln_phi[1]
        );
        assert_eq!(r.ln_phi[1], r.ln_phi[2], "both ions take the same branch");
        water.push(r.ln_phi[0]);
        carbon_dioxide.push(r.ln_phi[3]);
    }
    assert!(water[0] < water[1] && water[1] < water[2], "{water:?}");
    assert!(
        carbon_dioxide[0] < carbon_dioxide[1] && carbon_dioxide[1] < carbon_dioxide[2],
        "{carbon_dioxide:?}"
    );
    // **And they are not the same curve scaled.** `P0` rises faster than `H` over this
    // range, so the gap between them narrows - which is what says two correlations were
    // read rather than one.
    let first_gap = carbon_dioxide[0] - water[0];
    let last_gap = carbon_dioxide[2] - water[2];
    assert!(first_gap > last_gap, "{first_gap} against {last_gap}");
}

/// A `solvent` component with no vapour pressure is refused rather than given the filler
/// number NeqSim returns.
///
/// `MDEA` is tagged `solvent` and carries no Antoine row, so `ComponentKentEisenberg` gives
/// it `NaN / P` - a fifteen-digit NaN that reaches `ln_phi` without comment. This model
/// refuses it, which is the one place it diverges from NeqSim.
#[test]
fn a_solvent_with_no_vapour_pressure_is_refused() {
    let err = kent_eisenberg_phase(
        &["water", "MDEA", "Na+", "Cl-"],
        313.15,
        500_000.0,
        &[0.80, 0.08, 0.06, 0.06],
    )
    .unwrap_err();
    assert!(
        matches!(err, AzothError::PropertyUnavailable { .. }),
        "MDEA has no vapour-pressure correlation: {err:?}"
    );
    assert!(
        err.to_string().contains("Antoine"),
        "the refusal should name what is missing: {err}"
    );
}

/// **The `1e8` is this model's own.** Pitzer gives an ion `1e12` and Desmukh-Mather gives
/// it `1e-15`, so a port that shared one constant between the three would be wrong at every
/// state - and on a two-component mixture it would be wrong by a *plausible* amount.
#[test]
fn the_insoluble_ion_constant_is_not_shared() {
    let r = kent_eisenberg_phase(
        &["water", "Na+", "Cl-"],
        313.15,
        500_000.0,
        &[0.9, 0.05, 0.05],
    )
    .unwrap();
    assert!((r.ln_phi[1] - 1.0e8f64.ln()).abs() < 1.0e-12);
    assert!(
        (r.ln_phi[1] - 1.0e12f64.ln()).abs() > 1.0,
        "Pitzer's constant is not this model's"
    );
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let err = kent_eisenberg_phase(
        &["water", "Na+", "Cl-"],
        313.15,
        500_000.0,
        &[0.9, 0.05, 0.2],
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}

#[test]
fn a_pressure_that_is_not_positive_is_refused() {
    let err = kent_eisenberg_phase(&["water", "Na+", "Cl-"], 313.15, 0.0, &[0.9, 0.05, 0.05])
        .unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }), "{err:?}");
}
