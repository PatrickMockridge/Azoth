//! Spec-driven tests for `eos.production_chemical_scale_scenario`.

use azoth_core::spec::{RecordField, TestCase};
use azoth_core::units::{kelvins, kilograms_per_cubic_meter, pascals};
use azoth_eos::model_gen;
use azoth_eos::production_chemical_scale_scenario::{
    ProductionChemicalDose, production_chemical_scale_scenario,
};
use azoth_eos::results::ProductionChemicalScaleScenarioResult;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.production_chemical_scale_scenario";

fn chemical(row: &[RecordField]) -> ProductionChemicalDose {
    let number = |name: &str| common::record_number(row, name);
    ProductionChemicalDose {
        chemical_type: common::record_name(row, "chemical_type")
            .expect("every row states its type")
            .to_owned(),
        name: common::record_name(row, "name")
            .expect("every row states its name")
            .to_owned(),
        active_ingredient: common::record_name(row, "active_ingredient")
            .expect("every row states its active")
            .to_owned(),
        active_wt_pct: number("active_wt_pct").expect("every row states its active wt%"),
        dosage: kilograms_per_cubic_meter(number("dosage").expect("every row states its dose")),
        alkalinity_capacity: number("alkalinity_capacity"),
        h2s_capacity: number("h2s_capacity"),
    }
}

fn call(case: &TestCase) -> ProductionChemicalScaleScenarioResult {
    let mass = |name: &str| kilograms_per_cubic_meter(common::input(case, name));
    let chemicals: Vec<ProductionChemicalDose> = common::record_input(case, "chemicals")
        .iter()
        .map(|row| chemical(row))
        .collect();
    production_chemical_scale_scenario(
        mass("calcium"),
        mass("barium"),
        mass("strontium"),
        mass("iron"),
        mass("magnesium"),
        mass("sodium"),
        mass("sulphate"),
        mass("bicarbonate"),
        mass("total_dissolved_solids"),
        common::input(case, "ph"),
        mass("dissolved_h2s"),
        kelvins(common::input(case, "temperature")),
        pascals(common::input(case, "pressure")),
        common::input_str(case, "activity_model"),
        &chemicals,
        case.input("pitzer_ionic_strength_molal"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        for (field, got) in [
            ("treated_ph", result.treated_ph),
            ("treated_bicarbonate", result.treated_bicarbonate.value),
            ("base_dose_equivalents", result.base_dose_equivalents.value),
            ("scavenger_capacity", result.scavenger_capacity.value),
            (
                "residual_dissolved_h2s",
                result.residual_dissolved_h2s.value,
            ),
        ] {
            if let Some(want) = case.expected_value(field) {
                common::assert_close(got, want, case.tolerance, &format!("{context} ({field})"));
            }
        }

        for (name, rows) in [
            ("baseline_si", &result.baseline_si),
            ("treated_si", &result.treated_si),
        ] {
            let expected = case.record(name).unwrap_or(&[]);
            if case.record(name).is_some() {
                assert_eq!(rows.len(), expected.len(), "{context}: {name} rows");
            }
            for (index, (row, declared)) in rows.iter().zip(expected).enumerate() {
                if let Some(mineral) = common::record_name(declared, "mineral") {
                    assert_eq!(row.mineral, mineral, "{context}: {name}[{index}] mineral");
                }
                if let Some(want) = common::record_number(declared, "saturation_index") {
                    common::assert_close(
                        row.saturation_index.expect("the case pins an index"),
                        want,
                        case.tolerance,
                        &format!("{context}: {name}[{index}]"),
                    );
                }
            }
        }

        let changes = case.record("si_change").unwrap_or(&[]);
        if case.record("si_change").is_some() {
            assert_eq!(
                result.si_change.len(),
                changes.len(),
                "{context}: change rows"
            );
        }
        for (index, (row, declared)) in result.si_change.iter().zip(changes).enumerate() {
            if let Some(mineral) = common::record_name(declared, "mineral") {
                assert_eq!(row.mineral, mineral, "{context}: change[{index}] mineral");
            }
            if let Some(want) = common::record_number(declared, "change") {
                common::assert_close(
                    row.change.expect("the case pins a change"),
                    want,
                    case.tolerance,
                    &format!("{context}: change[{index}]"),
                );
            }
        }

        let warnings = case.record("model_warnings").unwrap_or(&[]);
        if case.record("model_warnings").is_some() {
            assert_eq!(
                result.model_warnings.len(),
                warnings.len(),
                "{context}: warning rows"
            );
        }
        for (index, (row, declared)) in result.model_warnings.iter().zip(warnings).enumerate() {
            if let Some(code) = common::record_name(declared, "code") {
                assert_eq!(row.code, code, "{context}: warning[{index}] code");
            }
            if let Some(message) = common::record_name(declared, "message") {
                assert_eq!(row.message, message, "{context}: warning[{index}] message");
            }
        }
        common::assert_consistent(&result, context);
    }
}

fn caustic(capacity: Option<f64>) -> ProductionChemicalDose {
    ProductionChemicalDose {
        chemical_type: "ph_adjuster".to_owned(),
        name: "caustic".to_owned(),
        active_ingredient: "NaOH".to_owned(),
        active_wt_pct: 50.0,
        dosage: kilograms_per_cubic_meter(0.1),
        alkalinity_capacity: capacity,
        h2s_capacity: None,
    }
}

fn run(chemicals: &[ProductionChemicalDose], bicarbonate: f64) -> ProductionChemicalScaleScenarioResult {
    production_chemical_scale_scenario(
        kilograms_per_cubic_meter(2.0),
        kilograms_per_cubic_meter(0.05),
        kilograms_per_cubic_meter(0.1),
        kilograms_per_cubic_meter(0.005),
        kilograms_per_cubic_meter(0.5),
        kilograms_per_cubic_meter(20.0),
        kilograms_per_cubic_meter(0.5),
        kilograms_per_cubic_meter(bicarbonate),
        kilograms_per_cubic_meter(80.0),
        6.5,
        kilograms_per_cubic_meter(0.02),
        kelvins(333.15),
        pascals(1.0e7),
        "davies",
        chemicals,
        None,
    )
    .expect("the state runs")
}

/// **A stated capacity overrides the name-based fit**, as the class's non-NaN branch does.
#[test]
fn a_stated_capacity_overrides_the_fit() {
    let fitted = run(&[caustic(None)], 0.8);
    let stated = run(&[caustic(Some(50.0))], 0.8);
    // 50 mol/kg * 50 mg active/L = 2.5e-3 mol/L, twice the name-based fit.
    assert!((stated.base_dose_equivalents.value - 2.5).abs() < 1e-9);
    assert!(stated.treated_ph > fitted.treated_ph);
}

/// **With no bicarbonate the carbonate minerals are absent**, both from the equilibrium and
/// from the change list, rather than a zero row.
#[test]
fn no_bicarbonate_leaves_carbonate_absent() {
    let result = run(&[caustic(None)], 0.0);
    assert!(
        result
            .baseline_si
            .iter()
            .all(|row| !matches!(row.mineral.as_str(), "CaCO3" | "FeCO3"))
    );
    assert!(
        result
            .si_change
            .iter()
            .all(|row| !matches!(row.mineral.as_str(), "CaCO3" | "FeCO3"))
    );
    assert_eq!(result.treated_bicarbonate.value, 0.0);
}

/// **An unrecognised active is a warning, not an effect**: the dose contributes nothing.
#[test]
fn an_unknown_active_changes_nothing() {
    let mut unknown = caustic(None);
    unknown.active_ingredient = "KOH".to_owned();
    unknown.name = "bespoke".to_owned();
    let result = run(&[unknown], 0.8);
    assert_eq!(result.base_dose_equivalents.value, 0.0);
    assert_eq!(result.treated_ph, 6.5);
    assert_eq!(result.model_warnings.len(), 1);
    assert_eq!(
        result.model_warnings[0].code,
        "unsupported_ph_adjuster_bespoke"
    );
}
