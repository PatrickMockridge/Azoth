//! Spec-driven tests for the `eos.desmukh_mather_phase` model.

use azoth_core::AzothError;
use azoth_eos::desmukh_mather_phase::desmukh_mather_phase;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.desmukh_mather_phase";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::DesmukhMatherPhaseResult {
    desmukh_mather_phase(
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
        for name in ["gamma", "ln_gamma", "molality", "ln_phi"] {
            let actual = match name {
                "gamma" => &result.gamma,
                "ln_gamma" => &result.ln_gamma,
                "molality" => &result.molality,
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
        for (name, want) in [
            ("ionic_strength", result.ionic_strength),
            ("solvent_molar_mass", result.solvent_molar_mass),
        ] {
            let expected = case
                .expected_value(name)
                .unwrap_or_else(|| panic!("the case declares {name}"));
            common::assert_close(
                want,
                expected,
                case.tolerance,
                &format!("{context} ({name})"),
            );
        }
        common::assert_consistent(&result, context);
    }
}

/// **The pair sum is live, and it is the only thing that makes two like ions differ.**
///
/// `MDEA+` and `Cl-` are both monovalent, at the same molality, in a brine with no
/// same-sign pair - so the Debye-Huckel term is identical for them and any difference is
/// the pair sum. Measured, the difference is exactly `2 aij m_CO2`.
#[test]
fn the_pair_sum_makes_two_like_ions_differ() {
    let r = desmukh_mather_phase(
        &["water", "MDEA+", "Cl-", "CO2"],
        313.15,
        500_000.0,
        &[0.89, 0.04, 0.04, 0.03],
    )
    .unwrap();
    let expected = 2.0 * 0.0678732 * r.molality[3];
    let difference = r.ln_gamma[1] - r.ln_gamma[2];
    assert!(
        (difference - expected).abs() < 1.0e-12,
        "the two ions differ by {difference}, and `2 aij m(CO2)` is {expected}"
    );

    // And on a brine whose pairs are all zero, the same two ions agree exactly.
    let plain = desmukh_mather_phase(
        &["water", "Na+", "Cl-"],
        313.15,
        500_000.0,
        &[0.90, 0.05, 0.05],
    )
    .unwrap();
    assert_eq!(plain.ln_gamma[1], plain.ln_gamma[2]);
}

/// **The solvent is selected by `REFERENCESTATETYPE`, not by the name `water`.**
///
/// Methanol is tagged `solvent`, so it joins water in `getSolventWeight` and the mean molar
/// mass is `0.0187943` - neither component's. Water's activity coefficient is then
/// `1.11111` rather than the `1.17647` a water-only solvent gives.
#[test]
fn the_solvent_is_whatever_the_reference_state_says() {
    let mixed = desmukh_mather_phase(
        &["water", "methanol", "Na+", "Cl-"],
        313.15,
        500_000.0,
        &[0.85, 0.05, 0.05, 0.05],
    )
    .unwrap();
    assert!(
        (mixed.solvent_molar_mass - 0.018_794_277_777_777_774).abs() < 1.0e-15,
        "the mean molar mass is {}",
        mixed.solvent_molar_mass
    );
    assert_ne!(mixed.solvent_molar_mass, 0.018015, "not water's");

    // Water alone gives water's own molar mass and `1 / x_water`, because its pair sum is
    // zero and the conversion is `m M / x` with `m M = 1`.
    let alone = desmukh_mather_phase(
        &["water", "Na+", "Cl-"],
        313.15,
        500_000.0,
        &[0.90, 0.05, 0.05],
    )
    .unwrap();
    assert_eq!(alone.solvent_molar_mass, 0.018015);
    assert!((alone.gamma[0] - 1.0 / 0.90).abs() < 1.0e-12);
    assert!(
        (mixed.gamma[0] - 1.0 / 0.85).abs() > 0.05,
        "the mixed solvent should not give the water-only answer"
    );
}

/// **An ion gets `1e-15`**, the smallest of the tranche's three insoluble-ion constants.
///
/// Pitzer's is `1e12` and Kent-Eisenberg's `1e8`, and the three say different things: this
/// one says the ion is *absent* from the vapour rather than sparingly present.
#[test]
fn the_insoluble_ion_constant_is_not_shared() {
    let r = desmukh_mather_phase(
        &["water", "Na+", "Cl-"],
        313.15,
        500_000.0,
        &[0.90, 0.05, 0.05],
    )
    .unwrap();
    assert!((r.ln_phi[1] - 1.0e-15f64.ln()).abs() < 1.0e-12);
    assert!(
        (r.ln_phi[1] - 1.0e8f64.ln()).abs() > 1.0,
        "Kent-Eisenberg's constant is not this model's"
    );
}

/// A mixture with no `solvent`-reference component is refused rather than divided by zero.
#[test]
fn a_mixture_with_no_solvent_is_refused() {
    let err =
        desmukh_mather_phase(&["methane", "CO2"], 313.15, 500_000.0, &[0.5, 0.5]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let err = desmukh_mather_phase(
        &["water", "Na+", "Cl-"],
        313.15,
        500_000.0,
        &[0.9, 0.05, 0.2],
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}
