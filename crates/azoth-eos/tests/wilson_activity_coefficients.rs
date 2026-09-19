//! Spec-driven tests for the `eos.wilson_activity_coefficients` model.

use azoth_core::{AzothError, CalcResult};
use azoth_eos::Cubic;
use azoth_eos::databank;
use azoth_eos::mixture::{Component, Mixture};
use azoth_eos::{model_gen, wilson_activity_coefficients};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.wilson_activity_coefficients";

/// The n-butane/nc12 mixture the databank resolves for the case, which is also the
/// doc-test's and the refusal tests' stand-in.
fn butane_nc12() -> Mixture {
    databank::mixture_of(&["n-butane", "nc12"], Cubic::Pr, None)
        .expect("the case's components resolve")
        .0
}

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::WilsonActivityCoefficientsResult {
    let (mixture, _) = databank::mixture_of(
        case.list("components").expect("components"),
        Cubic::Pr,
        None,
    )
    .unwrap_or_else(|e| panic!("case `{}` should resolve but failed: {e}", case.id));
    wilson_activity_coefficients(
        &mixture,
        common::input(case, "T"),
        case.vector("x").expect("x"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        for (name, actual) in [("ln_gamma", &result.ln_gamma), ("gamma", &result.gamma)] {
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

#[test]
fn the_result_is_clean_at_an_ordinary_state() {
    let r = wilson_activity_coefficients(&butane_nc12(), 298.15, &[0.5, 0.5]).unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

/// The correlation reads the molar mass and the critical temperature, and the case's
/// recorded values are those two databank columns. This pins the resolution rather
/// than the numbers, which every case above would follow anywhere.
#[test]
fn the_databank_resolves_the_molar_mass_and_critical_temperature() {
    let (mixture, _) = databank::mixture_of(&["n-butane", "nc12"], Cubic::Pr, None).unwrap();
    let components = mixture.components();
    let mass: Vec<f64> = components.iter().filter_map(|c| c.molar_mass).collect();
    let tc: Vec<f64> = components.iter().map(|c| c.tc.value).collect();
    assert_eq!(mass.len(), 2, "both components carry a molar mass");
    assert!((mass[0] - 0.058_123).abs() < 1e-12, "{}", mass[0]);
    assert!((mass[1] - 0.1703).abs() < 1e-12, "{}", mass[1]);
    assert!((tc[0] - 425.12).abs() < 1e-12, "{}", tc[0]);
    assert!((tc[1] - 658.0).abs() < 1e-12, "{}", tc[1]);
}

/// The molar mass is what the correlation turns into a carbon number and a fusion
/// temperature, so a component without one is refused rather than given a zero - which
/// would also flip the `M_i > M_j` branch that decides `Lambda_ij`.
#[test]
fn a_component_without_a_molar_mass_is_refused() {
    let components = vec![
        Component::new(
            azoth_core::units::kelvins(425.12),
            azoth_core::units::pascals(3_796_000.0),
            0.2,
        )
        .unwrap(),
    ];
    let mixture = Mixture::new(components, vec![0.0]).unwrap();
    let err = wilson_activity_coefficients(&mixture, 298.15, &[1.0]).unwrap_err();
    assert!(
        matches!(err, AzothError::PropertyUnavailable { .. }),
        "{err:?}"
    );
}

#[test]
fn a_length_mismatch_is_refused() {
    let err = wilson_activity_coefficients(&butane_nc12(), 298.15, &[1.0]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("components"));
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let err = wilson_activity_coefficients(&butane_nc12(), 298.15, &[0.6, 0.6]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}
