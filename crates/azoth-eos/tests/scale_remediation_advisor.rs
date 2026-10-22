//! Spec-driven tests for `eos.scale_remediation_advisor`.

use azoth_core::spec::RecordField;
use azoth_eos::model_gen;
use azoth_eos::results::{RemediationOptionRow, ScaleRemediationAdvisorResult};
use azoth_eos::scale_remediation_advisor::{canonical_key, scale_remediation_advisor};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.scale_remediation_advisor";

fn call(case: &azoth_core::spec::TestCase) -> ScaleRemediationAdvisorResult {
    scale_remediation_advisor(common::input_str(case, "deposit"))
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

fn assert_row(expected: &[RecordField], row: &RemediationOptionRow, at: &str, tol: f64) {
    for field in [
        "scale_type",
        "category",
        "dissolver",
        "concentration",
        "method",
        "effectiveness",
        "cautions",
        "standard",
    ] {
        let declared = common::record_name(expected, field)
            .unwrap_or_else(|| panic!("{at}: the case does not state `{field}`"));
        let actual = match field {
            "scale_type" => &row.scale_type,
            "category" => &row.category,
            "dissolver" => &row.dissolver,
            "concentration" => &row.concentration,
            "method" => &row.method,
            "effectiveness" => &row.effectiveness,
            "cautions" => &row.cautions,
            _ => &row.standard,
        };
        assert_eq!(actual, declared, "{at}::{field}");
    }
    common::assert_close(
        row.min_temperature.value,
        common::record_number(expected, "min_temperature").expect("the case states the window"),
        tol,
        &format!("{at}::min_temperature"),
    );
    common::assert_close(
        row.max_temperature.value,
        common::record_number(expected, "max_temperature").expect("the case states the window"),
        tol,
        &format!("{at}::max_temperature"),
    );
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        if let Some(expected) = case.expected_string("canonical_key") {
            assert_eq!(
                result.canonical_key, expected,
                "{context}: the canonical key"
            );
        }
        let expected = common::expected_record(case, "options");
        assert_eq!(
            result.options.len(),
            expected.len(),
            "{context}: the option count"
        );
        for (index, (row, declared)) in result.options.iter().zip(expected).enumerate() {
            assert_row(
                declared,
                row,
                &format!("{context}::options[{index}]"),
                case.tolerance,
            );
        }
        common::assert_consistent(&result, context);
    }
}

/// **The alias table is case-insensitive and the unmatched name is not folded** - the
/// class's own asymmetry, which is why `canonicalKey` splits trimming from lower-casing.
#[test]
fn canonical_key_maps_aliases_and_preserves_the_rest() {
    assert_eq!(canonical_key("calcite"), "CaCO3");
    assert_eq!(canonical_key("Barite"), "BaSO4");
    assert_eq!(canonical_key("  mackinawite "), "FeS");
    assert_eq!(canonical_key("scavenger solids"), "Dithiazine");
    assert_eq!(canonical_key("Gypsumx"), "Gypsumx");
    assert_eq!(canonical_key("  "), "UNKNOWN");
}

/// **The unknown fallback is the knowledge base's UNKNOWN row**, not an empty list, and
/// it is the only row returned for a name no alias reaches.
#[test]
fn an_unknown_deposit_gets_the_sample_and_identify_row() {
    let result = scale_remediation_advisor("some-mystery-deposit").expect("the lookup runs");
    assert_eq!(result.canonical_key, "some-mystery-deposit");
    assert_eq!(result.options.len(), 1);
    assert_eq!(result.options[0].scale_type, "UNKNOWN");
    assert_eq!(result.options[0].category, "UNSPECIFIED");
}

/// **The order is descending effectiveness and stable within a rank**, checked over a
/// deposit with two rows of different ranks.
#[test]
fn options_are_ordered_by_effectiveness() {
    let result = scale_remediation_advisor("dithiazine").expect("the lookup runs");
    assert_eq!(result.options[0].effectiveness, "MEDIUM");
    assert_eq!(result.options[1].effectiveness, "LOW");
}
