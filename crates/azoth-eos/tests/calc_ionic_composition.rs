//! Spec-driven tests for the `eos.calc_ionic_composition` model.

use azoth_core::AzothError;
use azoth_core::units::{kelvins, pascals};
use azoth_eos::calc_ionic_composition::calc_ionic_composition;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.calc_ionic_composition";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::CalcIonicCompositionResult {
    let components: Vec<String> = case
        .list("components")
        .expect("components")
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    calc_ionic_composition(
        &components,
        kelvins(common::input(case, "T")),
        pascals(common::input(case, "P")),
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

        // **The composite test is the solvent's own row.** Water's molality is `1/M_water`
        // and its mass ratio is one, whatever the state - so a port that divided by the
        // wrong component's mass, or that dropped the solvent from its own table, fails on
        // the entry that needs no property of the mixture at all.
        let water = case
            .list("components")
            .expect("components")
            .iter()
            .position(|name| name.eq_ignore_ascii_case("water"))
            .expect("the class reports per kilogram of water");
        common::assert_close(
            result.mass_ratio[water],
            1.0,
            1.0e-12,
            &format!("{context} (water's mass ratio is one)"),
        );

        for (field, got) in [
            ("molality", &result.molality),
            ("mass_ratio", &result.mass_ratio),
            ("activity_coefficient", &result.activity_coefficient),
        ] {
            let wanted = case
                .expected_vector(field)
                .unwrap_or_else(|| panic!("case `{}` does not state {field}", case.id));
            assert_eq!(got.len(), wanted.len(), "{context}: {field} length");
            for (i, (value, want)) in got.iter().zip(wanted).enumerate() {
                common::assert_close(
                    *value,
                    *want,
                    case.tolerance,
                    &format!("{context} ({field}[{i}])"),
                );
            }
        }
        common::assert_consistent(&result, context);
    }
}

/// **A phase with no water cannot be reported on.** The two concentration columns are per
/// kilogram of it, and the class asks the phase for its `water` component rather than for one
/// whose reference state is `solvent` - so a phase of `CO2` and an ion is refused rather than
/// silently divided by something that is not there.
#[test]
fn a_phase_with_no_water_is_refused() {
    let components: Vec<String> = ["CO2", "Na+"]
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    let error = calc_ionic_composition(&components, kelvins(283.15), pascals(4.0e6), &[0.9, 0.1])
        .expect_err("there is no water to report per kilogram of");
    assert!(
        matches!(error, AzothError::InvalidInput { .. }),
        "{error:?}"
    );
    assert_eq!(error.field(), Some("components"));
}
