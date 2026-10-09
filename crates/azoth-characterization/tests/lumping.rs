//! Spec-driven tests for the `characterization.lumping` model.
//!
//! Every expected value is the lumping model's own `fractionOfHeavyEnd` and the added
//! components' molar masses and gravities, in
//! `validation/neqsim/captures/plus_fraction_probe.tsv`.

use azoth_characterization::{LumpingResult, lumping, model_gen};
use azoth_core::units::{MassDensity, MolarMass, kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_core::{AzothError, spec::TestCase};
use azoth_test_support as common;

const MODEL_ID: &str = "characterization.lumping";

fn call(case: &TestCase) -> LumpingResult {
    let masses: Vec<MolarMass> = case
        .vector("cut_molar_mass")
        .expect("the case states a cut table")
        .iter()
        .map(|value| kilograms_per_mole(*value))
        .collect();
    let densities: Vec<MassDensity> = case
        .vector("cut_density")
        .expect("the case states a cut table")
        .iter()
        .map(|value| kilograms_per_cubic_meter(*value))
        .collect();
    lumping(
        kilograms_per_mole(common::input(case, "molar_mass")),
        common::input(case, "mole_fraction"),
        case.vector("cut_z").expect("the case states a cut table"),
        &masses,
        &densities,
        case.input("number_of_lumps").map(|n| n as usize),
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
                "fraction_of_heavy_end",
                result.fraction_of_heavy_end.clone(),
                case.expected_vector("fraction_of_heavy_end"),
            ),
            (
                "lump_molar_mass",
                result.lump_molar_mass.iter().map(|q| q.value).collect(),
                case.expected_vector("lump_molar_mass"),
            ),
            (
                "lump_density",
                result.lump_density.iter().map(|q| q.value).collect(),
                case.expected_vector("lump_density"),
            ),
        ] {
            let Some(expected) = expected else { continue };
            assert_eq!(actual.len(), expected.len(), "{context}::{name}: the lump count");
            for (index, (got, want)) in actual.iter().zip(expected).enumerate() {
                common::assert_close(
                    *got,
                    *want,
                    case.tolerance,
                    &format!("{context}::{name}[{index}]"),
                );
            }
        }
        common::assert_consistent(&result, &context);
    }
}

#[test]
fn the_default_is_the_classs_own_seven() {
    // `LumpingModel.numberOfPseudocomponents` is 7, and every captured row is seven lumps
    // because nothing sets it. An absent count and a stated seven are the same run.
    let case = &model_gen::model(MODEL_ID).expect("in its own table").cases[0];
    let defaulted = call(case);
    let stated = lumping(
        kilograms_per_mole(common::input(case, "molar_mass")),
        common::input(case, "mole_fraction"),
        case.vector("cut_z").expect("a cut table"),
        &case
            .vector("cut_molar_mass")
            .expect("a cut table")
            .iter()
            .map(|v| kilograms_per_mole(*v))
            .collect::<Vec<_>>(),
        &case
            .vector("cut_density")
            .expect("a cut table")
            .iter()
            .map(|v| kilograms_per_cubic_meter(*v))
            .collect::<Vec<_>>(),
        Some(7),
    )
    .expect("computes");
    // Compared field by field rather than as a whole result: an absent optional input makes the
    // spec's range check on it *skipped*, and the skipped check is a warning, so the two runs
    // differ by that and by nothing else.
    assert_eq!(defaulted.fraction_of_heavy_end, stated.fraction_of_heavy_end);
    assert_eq!(defaulted.lump_mole_fraction, stated.lump_mole_fraction);
    assert_eq!(defaulted.lump_molar_mass, stated.lump_molar_mass);
    assert_eq!(defaulted.lump_density, stated.lump_density);
    assert_eq!(defaulted.fraction_of_heavy_end.len(), 7);
}

#[test]
fn a_lump_count_above_the_cut_count_is_refused() {
    // The class can close one lump per cut and sizes its array to the count asked for, so more
    // lumps than cuts leaves entries it never fills and never reports as unfilled.
    let refused = lumping(
        kilograms_per_mole(0.4),
        0.1,
        &[0.06, 0.04],
        &[kilograms_per_mole(0.3), kilograms_per_mole(0.5)],
        &[
            kilograms_per_cubic_meter(800.0),
            kilograms_per_cubic_meter(900.0),
        ],
        Some(3),
    );
    let error = refused.expect_err("three lumps from two cuts is refused");
    assert!(error.to_string().contains('3'), "the refusal names the count: {error}");
    assert!(matches!(error, AzothError::InvalidInput { .. }));
}

#[test]
fn the_three_vectors_have_to_be_one_table() {
    let refused = lumping(
        kilograms_per_mole(0.4),
        0.1,
        &[0.05, 0.05],
        &[kilograms_per_mole(0.3)],
        &[
            kilograms_per_cubic_meter(800.0),
            kilograms_per_cubic_meter(900.0),
        ],
        None,
    );
    assert!(matches!(refused, Err(AzothError::InvalidInput { .. })));
    let empty = lumping(
        kilograms_per_mole(0.4),
        0.1,
        &[],
        &[],
        &[],
        None,
    );
    assert!(matches!(empty, Err(AzothError::InvalidInput { .. })));
}

#[test]
fn the_lump_mole_fractions_are_the_fractions_scaled_and_the_cuts_respent() {
    // `lump_mole_fraction` is an output with no row of its own - the capture prints
    // `fractionOfHeavyEnd` and the added components' moles - so what is assertable is the two
    // relations it stands in: each entry is its fraction times the plus fraction's mole fraction,
    // and the row sums to the cut table's own total.
    let case = &model_gen::model(MODEL_ID).expect("in its own table").cases[0];
    let result = call(case);
    for (index, (fraction, mole)) in result
        .fraction_of_heavy_end
        .iter()
        .zip(&result.lump_mole_fraction)
        .enumerate()
    {
        common::assert_close(
            *mole,
            fraction * common::input(case, "mole_fraction"),
            1e-15,
            &format!("lump_mole_fraction[{index}] against its fraction"),
        );
    }
    let total: f64 = case.vector("cut_z").expect("a cut table").iter().sum();
    common::assert_close(
        result.lump_mole_fraction.iter().sum::<f64>(),
        total,
        1e-15,
        "the lump mole fractions against the cut table's total",
    );
}
