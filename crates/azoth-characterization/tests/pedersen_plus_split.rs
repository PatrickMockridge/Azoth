//! Spec-driven tests for the `characterization.pedersen_plus_split` model.
//!
//! Every expected value is a row of `validation/neqsim/captures/plus_fraction_probe.tsv`.

use azoth_characterization::{PedersenPlusSplitResult, model_gen, pedersen_plus_split};
use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_core::{AzothError, spec::TestCase};
use azoth_test_support as common;

const MODEL_ID: &str = "characterization.pedersen_plus_split";

/// The plus fraction's specific gravity in g/cm3, which is the scale the class solves in.
const ANCHOR_ROW_C6: f64 = 0.663_999_915_122_986;

fn call(case: &TestCase) -> PedersenPlusSplitResult {
    pedersen_plus_split(
        kilograms_per_mole(common::input(case, "molar_mass")),
        kilograms_per_cubic_meter(common::input(case, "density")),
        common::input(case, "mole_fraction"),
        common::input(case, "first_carbon_number") as usize,
        common::input(case, "last_carbon_number") as usize,
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
        for (name, actual) in [
            ("z_intercept", result.z_intercept),
            ("z_slope", result.z_slope),
            ("density_intercept", result.density_intercept),
            ("density_slope", result.density_slope),
        ] {
            if let Some(expected) = case.expected_value(name) {
                common::assert_close(actual, expected, case.tolerance, &format!("{context}::{name}"));
            }
        }
        if let Some(expected) = case.expected_vector("cut_z") {
            assert_eq!(result.cut_z.len(), expected.len(), "{context}: the cut count");
            for (index, (actual, want)) in result.cut_z.iter().zip(expected).enumerate() {
                common::assert_close(
                    *actual,
                    *want,
                    case.tolerance,
                    &format!("{context}::cut_z[{index}]"),
                );
            }
        }
        for (name, actual, expected) in [
            (
                "cut_molar_mass",
                result
                    .cut_molar_mass
                    .iter()
                    .map(|q| q.value)
                    .collect::<Vec<_>>(),
                case.expected_vector("cut_molar_mass"),
            ),
            (
                "cut_density",
                result.cut_density.iter().map(|q| q.value).collect::<Vec<_>>(),
                case.expected_vector("cut_density"),
            ),
        ] {
            let Some(expected) = expected else { continue };
            assert_eq!(
                actual.len(),
                expected.len(),
                "{context}::{name}: the cut count"
            );
            for (index, (got, want)) in actual.iter().zip(expected).enumerate() {
                common::assert_close(*got, *want, case.tolerance, &format!("{context}::{name}[{index}]"));
            }
        }
        common::assert_consistent(&result, &context);
    }
}

#[test]
fn the_gravity_line_is_anchored_one_carbon_below_the_first_cut() {
    // The class writes `c + d*ln(first - 1) == rho_table[first]`, so the anchor expression is
    // evaluated at the carbon number *below* the first cut while the gravity it is set to is the
    // first cut's own. A port that anchored at `ln(first)` would satisfy its residual at a
    // different `(c, d)` and still look self-consistent, which is why this is asserted here
    // rather than left to the case rows.
    let result = pedersen_plus_split(
        kilograms_per_mole(0.4),
        kilograms_per_cubic_meter(850.0),
        0.1,
        20,
        80,
    )
    .expect("the capture's first row computes");
    // `PVTsimDensities[20 - 6]`, the table's fifteenth row.
    let anchor = 0.862_999_975_681_305;
    let residual =
        result.density_intercept + result.density_slope * 19.0_f64.ln() - anchor;
    assert!(residual.abs() < 1.0e-12, "the anchor residual is {residual}");
}

#[test]
fn the_cuts_molar_masses_are_the_table_and_nothing_else() {
    // The abundance and gravity lines are solved; the molar masses are read. C20 is the table's
    // fifteenth row, 275 g/mol, and C79 (the last cut of the default range) is 2796 g/mol.
    let result = pedersen_plus_split(
        kilograms_per_mole(0.4),
        kilograms_per_cubic_meter(850.0),
        0.1,
        20,
        80,
    )
    .expect("the capture's first row computes");
    assert!((result.cut_molar_mass[0].value - 0.275).abs() < 1.0e-15);
    assert!((result.cut_molar_mass[59].value - 1.102).abs() < 1.0e-15);
    assert_eq!(result.cut_molar_mass.len(), 60);
    assert!(
        result.cut_molar_mass[0].value > ANCHOR_ROW_C6 / 1000.0,
        "the range starts well above the table's first row"
    );
}

#[test]
fn the_two_models_are_one_solve_over_two_ranges() {
    // `PedersenHeavyOilPlusModel` overrides only its constructor. The capture shows it by giving
    // the same coefficients as the explicitly-requested heavy row, and this reads the same
    // property off the port: the first cut of the longer range is the same number.
    let short = pedersen_plus_split(
        kilograms_per_mole(0.4),
        kilograms_per_cubic_meter(850.0),
        0.1,
        20,
        80,
    )
    .expect("computes");
    let long = pedersen_plus_split(
        kilograms_per_mole(0.4),
        kilograms_per_cubic_meter(850.0),
        0.1,
        20,
        200,
    )
    .expect("computes");
    assert_eq!(long.cut_molar_mass.len(), 180);
    assert_eq!(short.cut_molar_mass[0], long.cut_molar_mass[0]);
    assert!(
        (short.z_intercept - long.z_intercept).abs() > 1.0e-3,
        "the two ranges solve to different coefficients, which is what makes them one solve"
    );
}

#[test]
fn a_plus_fraction_lighter_than_the_tables_first_cut_is_refused() {
    // The capture's `m_plus = 0.2` rows: `PVTsimMolarMass[14] = 275` g/mol against `0.2 * 1000`,
    // where NeqSim logs and returns without a split. A port that carried on would return cuts.
    let refused = pedersen_plus_split(
        kilograms_per_mole(0.2),
        kilograms_per_cubic_meter(780.0),
        0.1,
        20,
        80,
    );
    let error = refused.expect_err("a lighter plus fraction than C20 is refused");
    assert!(
        error.to_string().contains("275"),
        "the refusal should name the table's cut, got {error}"
    );
    assert!(matches!(error, AzothError::InvalidInput { .. }));
}

#[test]
fn an_empty_range_is_refused() {
    let refused = pedersen_plus_split(
        kilograms_per_mole(0.4),
        kilograms_per_cubic_meter(850.0),
        0.1,
        20,
        20,
    );
    assert!(matches!(refused, Err(AzothError::InvalidInput { .. })));
}

#[test]
fn a_carbon_number_outside_the_tables_is_refused() {
    // Below C6 the molar-mass table has no row; above C80 the density table has none. NeqSim
    // indexes both unchecked, so a port that kept its `usize` subtraction would panic.
    for (first, last) in [(5_usize, 80_usize), (81, 200), (20, 201)] {
        let refused = pedersen_plus_split(
            kilograms_per_mole(5.0),
            kilograms_per_cubic_meter(850.0),
            0.1,
            first,
            last,
        );
        assert!(
            matches!(
                refused,
                Err(AzothError::OutOfRange { .. } | AzothError::InvalidInput { .. })
            ),
            "C{first}..C{last} should be refused, got {refused:?}"
        );
    }
}
