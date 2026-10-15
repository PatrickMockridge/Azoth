//! Spec-driven tests for the `characterization.oil_assay_derivation` model.
//!
//! Every expected value is a row of `validation/neqsim/captures/oil_assay_derivation_probe.tsv`.

use azoth_characterization::{
    AssayBasis, AssayCutDerivation, OilAssayDerivationResult, model_gen, oil_assay_derivation,
};
use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_core::{AzothError, spec::TestCase};
use azoth_test_support as common;

const MODEL_ID: &str = "characterization.oil_assay_derivation";

/// One cut, with only the fields the test is about stated.
fn cut(name: &str, declared_fraction: f64) -> AssayCutDerivation {
    AssayCutDerivation {
        name: name.to_string(),
        declared_fraction,
        density: None,
        api_gravity: None,
        molar_mass: None,
        average_boiling_point: None,
        watson_characterization_factor: None,
        standard_component_name: None,
    }
}

fn with_density(cut: AssayCutDerivation, density: f64) -> AssayCutDerivation {
    AssayCutDerivation {
        density: Some(kilograms_per_cubic_meter(density)),
        ..cut
    }
}

fn with_boiling_point(cut: AssayCutDerivation, boiling_point: f64) -> AssayCutDerivation {
    AssayCutDerivation {
        average_boiling_point: Some(boiling_point),
        ..cut
    }
}

fn with_molar_mass(cut: AssayCutDerivation, molar_mass: f64) -> AssayCutDerivation {
    AssayCutDerivation {
        molar_mass: Some(kilograms_per_mole(molar_mass)),
        ..cut
    }
}

fn with_gravity(cut: AssayCutDerivation, gravity: f64) -> AssayCutDerivation {
    AssayCutDerivation {
        api_gravity: Some(gravity),
        ..cut
    }
}

fn with_watson(cut: AssayCutDerivation, watson: f64) -> AssayCutDerivation {
    AssayCutDerivation {
        watson_characterization_factor: Some(watson),
        ..cut
    }
}

fn with_standard(cut: AssayCutDerivation, name: &str) -> AssayCutDerivation {
    AssayCutDerivation {
        standard_component_name: Some(name.to_string()),
        ..cut
    }
}

/// The cuts a case states, as the kernel takes them.
fn cuts_of(case: &TestCase) -> Vec<AssayCutDerivation> {
    common::record_input(case, "cuts")
        .iter()
        .map(|row| AssayCutDerivation {
            name: common::record_name(row, "name")
                .expect("every cut is named")
                .to_string(),
            declared_fraction: common::record_number(row, "declared_fraction")
                .expect("every cut states a yield"),
            density: common::record_number(row, "density").map(kilograms_per_cubic_meter),
            api_gravity: common::record_number(row, "api_gravity"),
            molar_mass: common::record_number(row, "molar_mass").map(kilograms_per_mole),
            average_boiling_point: common::record_number(row, "average_boiling_point"),
            watson_characterization_factor: common::record_number(
                row,
                "watson_characterization_factor",
            ),
            standard_component_name: common::record_name(row, "standard_component_name")
                .map(str::to_string),
        })
        .collect()
}

fn call(case: &TestCase) -> OilAssayDerivationResult {
    let basis: AssayBasis = case
        .string("basis")
        .expect("the case names a basis")
        .parse()
        .expect("valid");
    oil_assay_derivation(
        basis,
        &cuts_of(case),
        common::input(case, "total_assay_mass"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = format!("{}::{}", spec.id, case.id);
        for (name, actual, expected) in [
            (
                "mass_fraction",
                result.mass_fraction.clone(),
                case.expected_vector("mass_fraction"),
            ),
            (
                "molar_mass",
                result.molar_mass.iter().map(|q| q.value).collect(),
                case.expected_vector("molar_mass"),
            ),
            ("moles", result.moles.clone(), case.expected_vector("moles")),
        ] {
            let Some(expected) = expected else { continue };
            assert_eq!(
                actual.len(),
                expected.len(),
                "{context}::{name}: the component count"
            );
            for (index, (got, want)) in actual.iter().zip(expected).enumerate() {
                common::assert_close(
                    *got,
                    *want,
                    case.tolerance,
                    &format!("{context}::{name}[{index}]"),
                );
            }
        }
        if let Some(want) = case.expected_value("total_moles") {
            common::assert_close(result.total_moles, want, case.tolerance, &context);
        }
        // **The balance is a contract rather than a value**: the class gates it at 1e-10 and the
        // numbers close by construction, so what a test can assert is that the gate held.
        assert!(
            result.mass_balance_error <= 1e-10,
            "{context}: the mass balance is {}",
            result.mass_balance_error
        );
        common::assert_consistent(&result, &context);
    }
}

#[test]
fn the_three_routes_are_taken_in_the_class_s_order() {
    // The capture's first fluid: a correlated cut, a standard component and a stated molar mass.
    // The standard component's 0.100205 kg/mol is the database's `n-heptane` and not a
    // correlation of anything in the slate, which is what says that route ran.
    let cuts = [
        with_boiling_point(with_density(cut("cut_a", 0.4), 780.0), 400.0),
        with_standard(cut("cut_b", 0.3), "n-heptane"),
        with_molar_mass(with_density(cut("cut_c", 0.3), 850.0), 0.25),
    ];
    let result = oil_assay_derivation(AssayBasis::Mass, &cuts, 1000.0).expect("computes");
    assert_eq!(
        result.component_name,
        vec!["cut_a_PC", "n-heptane", "cut_c_PC"]
    );
    common::assert_close(
        result.molar_mass[0].value,
        0.112_610_323_807_509,
        1e-12,
        "cut_a",
    );
    common::assert_close(result.molar_mass[1].value, 0.100_205, 1e-12, "n-heptane");
    common::assert_close(result.molar_mass[2].value, 0.25, 1e-12, "cut_c");
    assert!((result.mass_balance_error - 0.0).abs() < 1e-15);
}

#[test]
fn a_gravity_and_a_boiling_point_each_fall_back_to_the_cut_s_other_field() {
    // The capture's second and third fluids, in one place: gravity from an API gravity, boiling
    // point from a Watson factor.
    let by_api = oil_assay_derivation(
        AssayBasis::Mass,
        &[with_boiling_point(
            with_gravity(cut("light", 1.0), 55.0),
            340.0,
        )],
        100.0,
    )
    .expect("computes");
    common::assert_close(
        by_api.molar_mass[0].value,
        0.078_528_223_866_995_3,
        1e-12,
        "api",
    );

    let by_watson = oil_assay_derivation(
        AssayBasis::Mass,
        &[with_watson(with_density(cut("w1", 1.0), 800.0), 11.5)],
        100.0,
    )
    .expect("computes");
    common::assert_close(
        by_watson.molar_mass[0].value,
        0.132_491_038_438_475,
        1e-12,
        "watson",
    );
}

#[test]
fn a_cut_with_no_positive_yield_produces_no_component() {
    // The class skips rather than deriving a zero - so the answer is shorter than the slate, and
    // the names are what re-align it.
    let cuts = [
        with_boiling_point(with_density(cut("kept", 1.0), 780.0), 400.0),
        with_boiling_point(with_density(cut("empty", 0.0), 900.0), 500.0),
    ];
    let result = oil_assay_derivation(AssayBasis::Mass, &cuts, 100.0).expect("computes");
    assert_eq!(result.component_name, vec!["kept_PC"]);
}

#[test]
fn the_class_s_own_validations_are_refusals_here() {
    // Two cuts with one name. `validateUniqueCutNames` is case-sensitive, so `A` and `a` are two
    // cuts - which is what makes this a comparison and not a normalisation.
    let named_twice = [
        with_boiling_point(with_density(cut("cut", 0.5), 780.0), 400.0),
        with_boiling_point(with_density(cut("cut", 0.5), 800.0), 420.0),
    ];
    assert!(matches!(
        oil_assay_derivation(AssayBasis::Mass, &named_twice, 100.0),
        Err(AzothError::InvalidInput { .. })
    ));

    // Two cuts naming one standard component, which `validateUniqueStandardComponentNames`
    // compares lower-cased.
    let same_substance = [
        with_standard(cut("a", 0.5), "n-heptane"),
        with_standard(cut("b", 0.5), "N-Heptane"),
    ];
    assert!(matches!(
        oil_assay_derivation(AssayBasis::Mass, &same_substance, 100.0),
        Err(AzothError::InvalidInput { .. })
    ));

    // A standard component on a volume basis, and one that also states pseudo-component
    // properties: the two refusals inside the class's standard branch.
    assert!(matches!(
        oil_assay_derivation(
            AssayBasis::Volume,
            &[with_standard(cut("a", 1.0), "n-heptane")],
            100.0
        ),
        Err(AzothError::InvalidInput { .. })
    ));
    assert!(matches!(
        oil_assay_derivation(
            AssayBasis::Mass,
            &[with_density(
                with_standard(cut("a", 1.0), "n-heptane"),
                800.0
            )],
            100.0
        ),
        Err(AzothError::InvalidInput { .. })
    ));

    // A substance the database does not carry, which the class reports as an unknown component.
    assert!(matches!(
        oil_assay_derivation(
            AssayBasis::Mass,
            &[with_standard(cut("a", 1.0), "unobtainium")],
            100.0
        ),
        Err(AzothError::InvalidInput { .. })
    ));
}

#[test]
fn a_cut_that_can_reach_no_molar_mass_is_refused_by_the_field_it_lacks() {
    // The correlation needs a boiling point and a gravity, and each has a fallback - so a cut
    // missing both halves is refused, and the refusal names the field that would have supplied it.
    let no_gravity = oil_assay_derivation(
        AssayBasis::Mass,
        &[with_boiling_point(cut("a", 1.0), 400.0)],
        100.0,
    );
    assert!(matches!(no_gravity, Err(AzothError::InvalidInput { .. })));

    let no_boiling_point = oil_assay_derivation(
        AssayBasis::Mass,
        &[with_density(cut("a", 1.0), 780.0)],
        100.0,
    );
    assert!(matches!(
        no_boiling_point,
        Err(AzothError::InvalidInput { .. })
    ));

    // **A stated molar mass still needs a gravity** - and not for the correlation, which this
    // route never reaches: the class reads `resolveDensity` for every cut it derives, before it
    // asks whether the cut has a molar mass at all. So the third route is taken with no boiling
    // point, and not with no density.
    assert!(matches!(
        oil_assay_derivation(
            AssayBasis::Mass,
            &[with_molar_mass(cut("a", 1.0), 0.25)],
            100.0
        ),
        Err(AzothError::InvalidInput { .. })
    ));
    let stated = oil_assay_derivation(
        AssayBasis::Mass,
        &[with_molar_mass(with_density(cut("a", 1.0), 850.0), 0.25)],
        100.0,
    )
    .expect("computes");
    common::assert_close(stated.molar_mass[0].value, 0.25, 1e-15, "stated");
}
