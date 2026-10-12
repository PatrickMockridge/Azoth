//! Spec-driven tests for the `characterization.delumping` model.
//!
//! Every expected value is a row of `validation/neqsim/captures/pseudo_component_delumping_probe.tsv`.

use azoth_characterization::{DelumpingResult, Lump, SubFraction, delumping, model_gen};
use azoth_core::spec::TestCase;
use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_test_support as common;

const MODEL_ID: &str = "characterization.delumping";

/// One lump, with only the fields the test is about stated.
fn lump(name: &str, moles: f64, molar_mass: f64) -> Lump {
    Lump {
        name: name.to_string(),
        moles,
        molar_mass: kilograms_per_mole(molar_mass),
        density: None,
        boiling_point: None,
    }
}

/// The lumps a case declares, as the kernel takes them.
fn lumps(case: &TestCase) -> Vec<Lump> {
    common::record_input(case, "lumps")
        .iter()
        .map(|row| Lump {
            name: common::record_name(row, "name")
                .expect("every lump is named")
                .to_string(),
            moles: common::record_number(row, "moles").expect("every lump states its moles"),
            molar_mass: kilograms_per_mole(
                common::record_number(row, "molar_mass").expect("every lump states a molar mass"),
            ),
            density: common::record_number(row, "density").map(kilograms_per_cubic_meter),
            boiling_point: common::record_number(row, "boiling_point"),
        })
        .collect()
}

fn call(case: &TestCase) -> DelumpingResult {
    let resolution = case
        .input("resolution")
        .expect("every case states a resolution");
    delumping(&lumps(case), resolution as i64)
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

/// One row's field, compared against the case's own, absence included.
fn assert_row(expected: &[azoth_core::spec::RecordField], fraction: &SubFraction, at: &str, tol: f64) {
    assert_eq!(
        fraction.name,
        common::record_name(expected, "name").expect("every row is named"),
        "{at}: the name"
    );
    common::assert_close(
        fraction.moles,
        common::record_number(expected, "moles").expect("every row states moles"),
        tol,
        &format!("{at}::moles"),
    );
    common::assert_close(
        fraction.molar_mass.value,
        common::record_number(expected, "molar_mass").expect("every row states a molar mass"),
        tol,
        &format!("{at}::molar_mass"),
    );
    // **Both of these are absent per row and not per vector**, which is the whole reason the
    // answer is a record list rather than four parallel vectors: a lump that states no boiling
    // point yields sub-fractions that state none, and one that states none on one side and a
    // number on the other is a port that filled an absence in.
    for (label, got, name) in [
        ("density", fraction.density.map(|q| q.value), "density"),
        ("boiling_point", fraction.boiling_point, "boiling_point"),
    ] {
        match (got, common::record_number(expected, name)) {
            (Some(got), Some(want)) => common::assert_close(got, want, tol, &format!("{at}::{label}")),
            (None, None) => {}
            (got, want) => panic!("{at}::{label}: {got:?} against {want:?}"),
        }
    }
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = format!("{}::{}", spec.id, case.id);
        let expected = common::expected_record(case, "fractions");
        assert_eq!(
            result.fractions.len(),
            expected.len(),
            "{context}: the sub-fraction count"
        );
        for (index, (fraction, row)) in result.fractions.iter().zip(expected).enumerate() {
            assert_row(row, fraction, &format!("{context}::fractions[{index}]"), case.tolerance);
        }
        common::assert_consistent(&result, &context);
    }
}

#[test]
fn the_moles_and_the_mass_are_each_conserved_by_a_different_mechanism() {
    // **Two conservations, and a port that got one by the other's means fails here.** The mole
    // weights sum to one, so the parent's moles come back with no rescale at all; the molar-mass
    // grid is then scaled by one factor so the mass closes exactly. Scaling the *moles* instead
    // would hold the mass and break the mole count, so asserting them on one result is a claim
    // neither mechanism can satisfy alone.
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    let case = &spec.cases[0];
    let result = call(case);

    for parent in lumps(case) {
        let prefix = format!("{}_d", parent.name);
        let children: Vec<&SubFraction> = result
            .fractions
            .iter()
            .filter(|f| f.name.starts_with(&prefix))
            .collect();
        assert!(
            !children.is_empty(),
            "{}: a split lump's sub-fractions are named after it",
            parent.name
        );
        let moles: f64 = children.iter().map(|f| f.moles).sum();
        let mass: f64 = children.iter().map(|f| f.moles * f.molar_mass.value).sum();
        common::assert_close(moles, parent.moles, 1e-12, &format!("{}::moles", parent.name));
        common::assert_close(
            mass,
            parent.moles * parent.molar_mass.value,
            1e-12,
            &format!("{}::mass", parent.name),
        );
    }
}

#[test]
fn a_resolution_of_one_is_the_identity_and_not_an_error() {
    // The class's own early return: asking for one sub-fraction per lump is asking for the lumps.
    // It keeps the order it was given, which is the one thing that differs from a split answer.
    let lumps = vec![lump("heavy", 2.0, 0.30), lump("light", 1.0, 0.10)];
    let result = delumping(&lumps, 1).expect("computes");
    assert_eq!(
        result.fractions.iter().map(|f| f.name.as_str()).collect::<Vec<_>>(),
        vec!["heavy", "light"],
        "the identity keeps the input's order"
    );
    assert_eq!(result.fractions[0].moles, 2.0);
}

#[test]
fn a_lump_the_method_cannot_split_is_passed_through_rather_than_refused() {
    // `parentMoles > 0` and `parentMolarMass > 0`, the class's own guards. Kept rather than
    // refused: the method defines what it does with one, and it is unsplit and unsorted-by-split.
    let lumps = vec![lump("empty", 0.0, 0.20), lump("real", 1.0, 0.20)];
    let result = delumping(&lumps, 2).expect("computes");
    assert_eq!(result.fractions.len(), 3, "one passed through, one split in two");
    assert!(
        result.fractions.iter().any(|f| f.name == "empty"),
        "the lump it could not split is still in the answer"
    );
}
