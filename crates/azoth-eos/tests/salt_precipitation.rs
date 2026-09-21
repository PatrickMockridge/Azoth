//! Spec-driven tests for the `eos.salt_precipitation` model.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::{model_gen, salt_precipitation};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.salt_precipitation";

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty());

    for case in spec.cases {
        let names = case.list("components").expect("components");
        let result = salt_precipitation(
            names,
            common::input_str(case, "salt"),
            kelvins(common::input(case, "T")),
            pascals(common::input(case, "P")),
            case.vector("z").expect("z"),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        let context = &format!("{}::{}", spec.id, case.id);
        for (field, actual, expected) in [
            (
                "precipitated_moles",
                result.precipitated_moles,
                common::expected(case, "precipitated_moles"),
            ),
            (
                "initial_saturation_ratio",
                result.initial_saturation_ratio,
                common::expected(case, "initial_saturation_ratio"),
            ),
        ] {
            common::assert_close(
                actual,
                expected,
                case.tolerance,
                &format!("{context} ({field})"),
            );
        }
        assert!(
            (result.final_saturation_ratio - 1.0).abs() < 1.0e-6
                || (result.final_saturation_ratio - result.initial_saturation_ratio).abs() < 1.0e-9,
            "{context}: the final ratio is {} and the initial {}",
            result.final_saturation_ratio,
            result.initial_saturation_ratio
        );
    }
}

/// **The extent is monotone in the ratio, and the answer is where it reaches one.**
///
/// `SR(e)` falls as the mineral takes ions out - removing them can only lower the ion activity
/// product - so the bracket is `[0, min_i n_i/stoc_i]` and the answer is the crossing. This
/// checks the shape rather than the number: the ratio at zero is the initial one, and the
/// extent that came back puts it at one. `NaCl`'s product is about `38` and this brine is a
/// quarter over it, so the crossing is well inside what the sodium and the chloride can give -
/// which is what makes it a crossing rather than an exhaustion at the bracket's end.
#[test]
fn the_ratio_falls_with_the_extent() {
    let names = ["water", "Na+", "Cl-", "CO3--", "HCO3-"];
    let z = [
        0.802_568_218_298_555,
        0.096_308_186_195_826_6,
        0.096_308_186_195_826_6,
        0.003_210_272_873_194_22,
        0.001_605_136_436_597_11,
    ];
    let result =
        salt_precipitation(&names, "NaCl", kelvins(298.15), pascals(1.0e6), &z).expect("computes");

    // Halite is over its product by a quarter, so the extent is well inside what the sodium
    // and the chloride can give - which is what makes this a crossing rather than an
    // exhaustion.
    assert!(result.precipitated_moles > 0.0);
    assert!(result.precipitated_moles < z[1]);
    common::assert_close(result.final_saturation_ratio, 1.0, 1.0e-6, "at the answer");
    // The extent is bounded by P8's Pitzer agreement and not by this solve - see the case.
    common::assert_close(
        result.initial_saturation_ratio,
        1.272_199_198_956_26,
        1.0e-4,
        "the initial ratio",
    );
}

/// A mineral the brine cannot form is refused rather than reported as a zero.
#[test]
fn a_mineral_the_brine_cannot_form_is_refused() {
    let names = ["water", "Na+", "Cl-"];
    let error = salt_precipitation(
        &names,
        "BaSO4",
        kelvins(298.15),
        pascals(1.0e6),
        &[0.98, 0.01, 0.01],
    )
    .expect_err("no barium, no barite");
    assert!(
        matches!(error, azoth_core::AzothError::InvalidInput { .. }),
        "{error:?}"
    );

    let error = salt_precipitation(
        &names,
        "unobtainium",
        kelvins(298.15),
        pascals(1.0e6),
        &[0.98, 0.01, 0.01],
    )
    .expect_err("not a row");
    assert!(
        matches!(error, azoth_core::AzothError::InvalidInput { .. }),
        "{error:?}"
    );
}
