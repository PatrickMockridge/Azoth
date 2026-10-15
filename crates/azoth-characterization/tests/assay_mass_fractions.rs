//! Spec-driven tests for the `characterization.assay_mass_fractions` model.
//!
//! Every expected value is a row of `validation/neqsim/captures/oil_assay_probe.tsv`.

use azoth_characterization::{
    AssayBasis, AssayCut, AssayMassFractionsResult, assay_mass_fractions, model_gen,
};
use azoth_core::units::kilograms_per_cubic_meter;
use azoth_core::{AzothError, spec::TestCase};
use azoth_test_support as common;

const MODEL_ID: &str = "characterization.assay_mass_fractions";

/// The case's cuts, as the kernel takes them.
///
/// A field the case omits is a cut that states none, which is the shape the record exists for:
/// `record_number` returns `None` rather than failing, and the `Option` carries it through.
fn cuts_of(case: &TestCase) -> Vec<AssayCut> {
    common::record_input(case, "cuts")
        .iter()
        .map(|row| AssayCut {
            declared_fraction: common::record_number(row, "declared_fraction")
                .expect("every cut states its declared fraction"),
            density: common::record_number(row, "density").map(kilograms_per_cubic_meter),
        })
        .collect()
}

fn cut(declared_fraction: f64, density: Option<f64>) -> AssayCut {
    AssayCut {
        declared_fraction,
        density: density.map(kilograms_per_cubic_meter),
    }
}

fn call(case: &TestCase) -> AssayMassFractionsResult {
    let basis: AssayBasis = case
        .string("basis")
        .expect("the case names a basis")
        .parse()
        .expect("valid");
    assay_mass_fractions(basis, &cuts_of(case))
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = format!("{}::{}", spec.id, case.id);
        let expected = case
            .expected_vector("mass_fraction")
            .expect("the case expects a row of mass fractions");
        assert_eq!(
            result.mass_fraction.len(),
            expected.len(),
            "{context}: the cut count"
        );
        for (index, (got, want)) in result.mass_fraction.iter().zip(expected).enumerate() {
            common::assert_close(
                *got,
                *want,
                case.tolerance,
                &format!("{context}::w[{index}]"),
            );
        }
        if let Some(want) = case.expected_value("total_declared_fraction") {
            common::assert_close(
                result.total_declared_fraction,
                want,
                case.tolerance,
                &context,
            );
        }
        // **A scalar output, so `expected_value`.** This read `expected_vector` while the case
        // states a number, which meant `None` and a comparison that never ran - every case here
        // carried a bulk density nothing checked.
        if let Some(want) = case.expected_value("bulk_density") {
            let got = result
                .bulk_density
                .expect("the case gave every cut a density, so the bulk density is answerable");
            common::assert_close(
                got.value,
                want,
                case.tolerance,
                &format!("{context}::bulk_density"),
            );
        }
        common::assert_consistent(&result, &context);
    }
}

#[test]
fn the_mass_basis_resolves_to_itself() {
    // The capture's second row. `[0.4, 0.35, 0.25]` is already a mass basis, so the answer is the
    // input - which is what says the gravity conversion is skipped and not merely inert here.
    let result = assay_mass_fractions(
        AssayBasis::Mass,
        &[
            cut(0.4, Some(700.0)),
            cut(0.35, Some(820.0)),
            cut(0.25, Some(900.0)),
        ],
    )
    .expect("computes");
    assert_eq!(result.mass_fraction, vec![0.4, 0.35, 0.25]);
    // And with no densities at all, which a mass basis is allowed: only the bulk density is lost.
    let without = assay_mass_fractions(
        AssayBasis::Mass,
        &[cut(0.4, None), cut(0.35, None), cut(0.25, None)],
    )
    .expect("computes");
    assert_eq!(without.mass_fraction, result.mass_fraction);
    assert_eq!(without.bulk_density, None);
}

#[test]
fn a_mass_basis_may_state_a_density_for_some_cuts_and_not_others() {
    // **The shape a parallel vector could not have.** `density` used to be one vector aligned with
    // the fractions, so a mass-basis assay stated one for every cut or for none; here the middle
    // cut does not, and the fractions are unaffected while the bulk density goes.
    let partial = assay_mass_fractions(
        AssayBasis::Mass,
        &[
            cut(0.4, Some(700.0)),
            cut(0.35, None),
            cut(0.25, Some(900.0)),
        ],
    )
    .expect("computes");
    assert_eq!(partial.mass_fraction, vec![0.4, 0.35, 0.25]);
    assert_eq!(
        partial.bulk_density, None,
        "one cut without a density is not a bulk density averaged over the others"
    );
}

#[test]
fn a_volume_basis_is_the_gravity_weighted_renormalisation() {
    // One volume fraction times its gravity over the sum of both, which is not the mean of the two
    // gravities and not a uniform shift of the row.
    let result = assay_mass_fractions(
        AssayBasis::Volume,
        &[cut(0.5, Some(700.0)), cut(0.5, Some(900.0))],
    )
    .expect("computes");
    assert!((result.mass_fraction[0] - 0.4375).abs() < 1e-15);
    assert!((result.mass_fraction[1] - 0.5625).abs() < 1e-15);
}

#[test]
fn the_closure_is_a_hard_error_and_its_tolerance_is_a_thousandth() {
    // `1.0005` is inside and is normalised; `1.003` is outside and raises with the sum it saw.
    let inside = assay_mass_fractions(AssayBasis::Mass, &[cut(0.5, None), cut(0.5005, None)])
        .expect("computes");
    assert!((inside.total_declared_fraction - 1.0005).abs() < 1e-15);
    assert!((inside.mass_fraction.iter().sum::<f64>() - 1.0).abs() < 1e-15);
    assert!((inside.mass_fraction[0] - 0.499_750_124_937_531).abs() < 1e-15);

    let outside = assay_mass_fractions(AssayBasis::Mass, &[cut(0.5, None), cut(0.51, None)]);
    let error = outside.expect_err("1.01 is outside the closure tolerance");
    assert!(
        error.to_string().contains("0.001"),
        "the refusal names the tolerance: {error}"
    );
    assert!(matches!(error, AzothError::InvalidInput { .. }));
}

#[test]
fn a_volume_basis_needs_a_density_on_every_cut_and_names_the_one_that_lacks_it() {
    // The refusal used to be about the vector's length - a count against a count. It is now about
    // a particular cut, and the field names which one.
    assert!(matches!(
        assay_mass_fractions(AssayBasis::Volume, &[cut(0.5, None), cut(0.5, None)]),
        Err(AzothError::InvalidInput { .. })
    ));
    let hole = assay_mass_fractions(AssayBasis::Volume, &[cut(0.5, Some(700.0)), cut(0.5, None)]);
    let error = hole.expect_err("the second cut states no density");
    assert!(
        error.to_string().contains("cuts[1].density"),
        "the refusal names the cut: {error}"
    );
    assert!(matches!(
        assay_mass_fractions(AssayBasis::Mass, &[]),
        Err(AzothError::InvalidInput { .. })
    ));
}
