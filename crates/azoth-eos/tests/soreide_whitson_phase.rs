//! Spec-driven tests for the `eos.soreide_whitson_phase` model.

use azoth_core::AzothError;
use azoth_core::units::{kelvins, pascals};
use azoth_eos::soreide_whitson_phase::soreide_whitson_phase;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.soreide_whitson_phase";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::SoreideWhitsonPhaseResult {
    let components: Vec<String> = case
        .list("components")
        .expect("components")
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    soreide_whitson_phase(
        &components,
        kelvins(common::input(case, "T")),
        pascals(common::input(case, "P")),
        case.vector("x").expect("x"),
        common::input(case, "salinity"),
        common::input_str(case, "compressed_phase"),
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
            result.z_factor,
            case.expected_value("z_factor").expect("z_factor"),
            case.tolerance,
            &format!("{context} (z_factor)"),
        );
        let expected = case.expected_vector("ln_phi").expect("ln_phi");
        for (i, (&got, &want)) in result.ln_phi.iter().zip(expected).enumerate() {
            common::assert_close(
                got,
                want,
                case.tolerance,
                &format!("{context} (ln_phi[{i}])"),
            );
        }
        common::assert_consistent(&result, context);
    }
}

/// **A brine reaches both of the model's salinity paths, and the rule's gate decides which
/// of them a phase takes.**
///
/// Water's alpha reads the salinity directly; the aqueous interaction matrix reads it
/// through the rule. A port that wired only one of the two would agree with NeqSim at zero
/// salinity - which is every state the shipped test uses - and be wrong on every real brine.
///
/// The gate is NeqSim's `water.x > 0.8`, and the two phases here sit either side of it:
/// water is `0.9976` in the liquid and `0.0028` in the gas, so the correlation fires for
/// one and not the other. That is what the magnitudes below measure. The gas moves through
/// water's alpha alone - its interaction matrix keeps the base values - and the liquid takes
/// the correlation on top, an order of magnitude more.
#[test]
fn a_brine_reaches_both_paths_and_the_gate_decides_which() {
    let names = ["nitrogen", "CO2", "methane", "ethane", "water"];
    let vapour = [
        0.110818402794,
        0.221422141728,
        0.332455899714,
        0.332455590773,
        0.00284796499037,
    ];
    let aqueous = [
        4.43466786008e-05,
        0.00207206064762,
        0.000126652518641,
        0.00012950695283,
        0.997627433202,
    ];

    let names: Vec<String> = names.iter().map(|name| (*name).to_string()).collect();
    let at = |x: &[f64], phase: &str, salinity: f64| {
        soreide_whitson_phase(&names, kelvins(318.15), pascals(4.0e6), x, salinity, phase)
            .expect("the state should compute")
    };

    let fresh = at(&vapour, "vapour", 0.0);
    let salty = at(&vapour, "vapour", 4.0);
    assert!(
        (salty.ln_phi[4] - fresh.ln_phi[4]).abs() > 1.0e-3,
        "water's own alpha is salinity-dependent in every phase, and here it moves by {}",
        (salty.ln_phi[4] - fresh.ln_phi[4]).abs()
    );
    for i in 0..4 {
        assert!(
            (salty.ln_phi[i] - fresh.ln_phi[i]).abs() < 1.0e-4,
            "the gate is closed at x_water = 0.0028, so the gas matrix keeps its base \
             values and component {i} moves by {}",
            (salty.ln_phi[i] - fresh.ln_phi[i]).abs()
        );
    }

    let fresh = at(&aqueous, "liquid", 0.0);
    let salty = at(&aqueous, "liquid", 4.0);
    for i in 0..5 {
        assert!(
            (salty.ln_phi[i] - fresh.ln_phi[i]).abs() > 1.0e-2,
            "the liquid takes the correlation as well, so component {i} moves by {}",
            (salty.ln_phi[i] - fresh.ln_phi[i]).abs()
        );
    }
}

/// **A negative salinity is refused rather than evaluated.** The correlation's `s**0.75`
/// and `s**1.1` are real for a negative `s`, so the expression would return a number - a
/// fitted value for a brine that cannot exist.
#[test]
fn a_negative_salinity_is_refused() {
    let err = soreide_whitson_phase(
        &["methane".to_string(), "water".to_string()],
        kelvins(318.15),
        pascals(4.0e6),
        &[0.5, 0.5],
        -1.0,
        "vapour",
    )
    .expect_err("a molality cannot be negative");
    assert_eq!(err.field(), Some("salinity"), "{err:?}");
}

/// **An ion is refused, with the reason the mixture layer gives.** The Soreide-Whitson
/// model is a cubic with a brine in it, so the salt is the scalar the rule reads and not a
/// component - and a caller who added one has asked for a different model.
#[test]
fn an_ion_is_refused() {
    let err = soreide_whitson_phase(
        &["water".to_string(), "Na+".to_string(), "Cl-".to_string()],
        kelvins(318.15),
        pascals(4.0e6),
        &[0.9, 0.05, 0.05],
        0.0,
        "liquid",
    )
    .expect_err("a cubic has no notion of an ion");
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
}

/// The root a caller names has two spellings, and neither of them is a default.
#[test]
fn an_unknown_root_is_refused() {
    for bad in ["gas", "vapor ", "Liquid"] {
        let err = soreide_whitson_phase(
            &["methane".to_string(), "water".to_string()],
            kelvins(318.15),
            pascals(4.0e6),
            &[0.5, 0.5],
            0.0,
            bad,
        )
        .expect_err("every spelling but the two is refused");
        assert!(
            matches!(err, AzothError::InvalidInput { .. }),
            "{bad}: {err:?}"
        );
    }
}
