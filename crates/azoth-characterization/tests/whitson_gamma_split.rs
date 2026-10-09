//! Spec-driven tests for the `characterization.whitson_gamma_split` model.
//!
//! Every expected value is a `# gamma variant = ...` row of
//! `validation/neqsim/captures/plus_fraction_probe.tsv`.

use azoth_characterization::{
    WhitsonDensityModel, WhitsonGammaSplitResult, gamma, model_gen, p0_p1, whitson_gamma_split,
};
use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_core::{AzothError, spec::TestCase};
use azoth_test_support as common;

const MODEL_ID: &str = "characterization.whitson_gamma_split";

fn call(case: &TestCase) -> WhitsonGammaSplitResult {
    whitson_gamma_split(
        kilograms_per_mole(common::input(case, "molar_mass")),
        kilograms_per_cubic_meter(common::input(case, "density")),
        common::input(case, "mole_fraction"),
        common::input(case, "first_carbon_number") as usize,
        common::input(case, "last_carbon_number") as usize,
        case.input("alpha"),
        case.input("eta").map(kilograms_per_mole),
        case.string("density_model")
            .unwrap_or("uop")
            .parse()
            .expect("the case names a density model"),
        case.flag("auto_estimate_shape").unwrap_or(false),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

/// The capture's own row: 0.4 kg/mol, 850 kg/m3, every default.
fn defaults() -> WhitsonGammaSplitResult {
    whitson_gamma_split(
        kilograms_per_mole(0.4),
        kilograms_per_cubic_meter(850.0),
        0.1,
        1,
        80,
        None,
        None,
        WhitsonDensityModel::Uop,
        false,
    )
    .expect("the capture's default row computes")
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = format!("{}::{}", spec.id, case.id);
        for (name, actual) in [("shape", result.shape)] {
            if let Some(expected) = case.expected_value(name) {
                common::assert_close(actual, expected, case.tolerance, &format!("{context}::{name}"));
            }
        }
        for (name, actual, expected) in [
            (
                "cut_z",
                result.cut_z.clone(),
                case.expected_vector("cut_z"),
            ),
            (
                "cut_molar_mass",
                result.cut_molar_mass.iter().map(|q| q.value).collect(),
                case.expected_vector("cut_molar_mass"),
            ),
            (
                "cut_density",
                result.cut_density.iter().map(|q| q.value).collect(),
                case.expected_vector("cut_density"),
            ),
        ] {
            let Some(expected) = expected else { continue };
            assert_eq!(actual.len(), expected.len(), "{context}::{name}: the cut count");
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
fn the_carbon_numbers_label_the_row_and_do_not_move_the_windows() {
    // The windows step from `eta` in molar mass, so `first_carbon_number` cannot reach them:
    // the same model over C1..C79 and over C20..C79 must agree cut for cut on the overlap. A
    // port that derived the window from the carbon number - `eta + 14*CN` - would not.
    let from_one = defaults();
    let from_twenty = whitson_gamma_split(
        kilograms_per_mole(0.4),
        kilograms_per_cubic_meter(850.0),
        0.1,
        20,
        80,
        None,
        None,
        WhitsonDensityModel::Uop,
        false,
    )
    .expect("computes");
    assert_eq!(from_twenty.cut_z.len(), 60);
    // Every window but the tail is fixed by `eta` and the increment alone, so the molar masses
    // agree exactly - the last window of each run is widened to 10 000 and so is the one that
    // differs.
    for index in 0..59 {
        assert_eq!(
            from_one.cut_molar_mass[index],
            from_twenty.cut_molar_mass[index]
        );
    }
    // The abundances are the same windows renormalised to the plus fraction's own mole fraction,
    // and the two runs have different totals because their tail windows are at different cuts.
    // A ratio is what the windowing fixes and what a port deriving the window from the carbon
    // number would not reproduce.
    let one = from_one.cut_z[5] / from_one.cut_z[0];
    let twenty = from_twenty.cut_z[5] / from_twenty.cut_z[0];
    common::assert_close(one, twenty, 1e-12, "the abundance ratio over the first six cuts");
}

#[test]
fn the_two_gravity_correlations_are_separable_from_the_split() {
    // The capture's `density=soreide` row. Only the densities move, which is what says the
    // correlation is fed the split rather than the split being fitted to it.
    let uop = defaults();
    let soreide = whitson_gamma_split(
        kilograms_per_mole(0.4),
        kilograms_per_cubic_meter(850.0),
        0.1,
        1,
        80,
        None,
        None,
        WhitsonDensityModel::Soreide,
        false,
    )
    .expect("computes");
    assert_eq!(uop.cut_z, soreide.cut_z);
    assert_eq!(uop.cut_molar_mass, soreide.cut_molar_mass);
    assert_ne!(uop.cut_density, soreide.cut_density);
    // The Søreide row's own first cut, from the capture.
    assert!((soreide.cut_density[0].value - 699.840_980_914_405).abs() < 1e-9);
}

#[test]
fn a_plus_fraction_at_or_below_eta_is_refused() {
    // `beta = (M_plus - eta)/alpha`, refused when it is not positive. At exactly `eta` the
    // class's own log line fires and it returns without a split.
    for molar_mass in [0.09, 0.05] {
        let refused = whitson_gamma_split(
            kilograms_per_mole(molar_mass),
            kilograms_per_cubic_meter(850.0),
            0.1,
            1,
            80,
            None,
            None,
            WhitsonDensityModel::Uop,
            false,
        );
        assert!(
            matches!(refused, Err(AzothError::InvalidInput { .. })),
            "a plus fraction at {molar_mass} kg/mol should be refused, got {refused:?}"
        );
    }
}

#[test]
fn the_scale_is_the_shape_it_used_that_derives_it() {
    // `scale` is an output with no oracle of its own - the capture prints `getCoefs()`'s two
    // values and not `betta` - so what is assertable is the relation, on a shape the capture
    // *did* measure: the auto-estimated `0.621390782386720`.
    let estimated = whitson_gamma_split(
        kilograms_per_mole(0.4),
        kilograms_per_cubic_meter(850.0),
        0.1,
        1,
        80,
        None,
        None,
        WhitsonDensityModel::Uop,
        true,
    )
    .expect("computes");
    common::assert_close(estimated.shape, 0.621_390_782_386_720, 1e-12, "the estimated shape");
    let expected = (0.4 - estimated.minimum_molar_mass.value) / estimated.shape;
    common::assert_close(estimated.scale.value, expected, 1e-15, "the derived scale");
}

#[test]
fn the_gamma_fit_is_the_one_the_class_carries() {
    // An eight-term polynomial for `Gamma(1+z)` on `[0, 1)` with a reflection below one. At one
    // it is exact; at a half it lands on `sqrt(pi)` to `1.4e-7`, which is the fit's own accuracy
    // and the check that says the coefficients are the right eight rather than a library call
    // wearing their names. **That error is part of every `P0` and `P1` the split is made of**, so
    // the tolerance is the fit's and not a tighter one this could be tuned to pass.
    assert_eq!(gamma(1.0), 1.0);
    assert!((gamma(0.5) - std::f64::consts::PI.sqrt()).abs() < 3e-7);
    // A window exactly at `eta` returns zero without evaluating, which is what makes the first
    // window's lower edge free rather than a limit.
    assert_eq!(p0_p1(90.0, 90.0, 1.0, 310.0), [0.0, 0.0]);
}
