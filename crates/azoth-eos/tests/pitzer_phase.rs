//! Spec-driven tests for the `eos.pitzer_phase` model.

use azoth_core::AzothError;
use azoth_eos::pitzer_phase::pitzer_phase;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.pitzer_phase";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::PitzerPhaseResult {
    pitzer_phase(
        case.list("components").expect("components"),
        common::input(case, "T"),
        common::input(case, "P"),
        case.vector("x").expect("x"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = azoth_model_spec(MODEL_ID);
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);

        for name in ["gamma", "ln_gamma", "molality"] {
            let actual = match name {
                "gamma" => &result.gamma,
                "ln_gamma" => &result.ln_gamma,
                _ => &result.molality,
            };
            let expected = case
                .expected_vector(name)
                .unwrap_or_else(|| panic!("the case declares {name}"));
            assert_eq!(actual.len(), expected.len(), "{context}: {name} length");
            for (i, (&got, &want)) in actual.iter().zip(expected).enumerate() {
                common::assert_close(
                    got,
                    want,
                    case.tolerance,
                    &format!("{context} ({name}[{i}])"),
                );
            }
        }

        for (name, want) in [
            ("ionic_strength", result.ionic_strength),
            ("osmotic_coefficient", result.osmotic_coefficient),
            ("water_activity", result.water_activity),
        ] {
            let expected = case
                .expected_value(name)
                .unwrap_or_else(|| panic!("the case declares {name}"));
            common::assert_close(
                want,
                expected,
                case.tolerance,
                &format!("{context} ({name})"),
            );
        }

        // **The dataset is an output and is asserted as one.** A coefficient without it is
        // not reproducible: the two datasets disagree on pairs they share.
        assert_eq!(
            result.dataset.name(),
            case.expected_string("dataset")
                .unwrap_or_else(|| panic!("the case declares dataset")),
            "{context}: which dataset answered"
        );

        common::assert_consistent(&result, context);
    }
}

/// Which dataset applies is a property of the **topology**, and one ion moves it.
///
/// This is the sabotage check the tranche's plan names: the catalogue covers
/// `water + Na+ + Cl-` and does not cover `water + Na+ + HCO3-`, and the two datasets give
/// different numbers for the pairs they share. Removing a component from the second brine
/// would take it back to the first dataset, so the answer has to move with the ions rather
/// than with the build.
#[test]
fn the_dataset_moves_with_the_topology() {
    let catalogue = pitzer_phase(
        &["water", "na+", "cl-"],
        298.15,
        100000.0,
        &[0.88, 0.06, 0.06],
    )
    .unwrap();
    let legacy = pitzer_phase(
        &["water", "na+", "hco3-"],
        298.15,
        100000.0,
        &[0.88, 0.06, 0.06],
    )
    .unwrap();
    assert_eq!(catalogue.dataset.name(), "phreeqc");
    assert_eq!(legacy.dataset.name(), "legacy");

    // **The same ionic strength on the two datasets**, and the sodium's activity
    // coefficient differs by more than a third - so a brine that took the wrong dataset
    // would still return a plausible-looking number.
    assert!(
        (catalogue.ionic_strength - legacy.ionic_strength).abs() < 1.0e-12,
        "the two brines have the same ionic strength: {} against {}",
        catalogue.ionic_strength,
        legacy.ionic_strength
    );
    assert!(
        (catalogue.ln_gamma[1] - legacy.ln_gamma[1]).abs() > 0.4,
        "ln gamma(Na+) is {} under the catalogue and {} under the CSV",
        catalogue.ln_gamma[1],
        legacy.ln_gamma[1]
    );
    assert!(
        catalogue.osmotic_coefficient > 1.0 && legacy.osmotic_coefficient < 1.0,
        "phi is {} on the catalogue and {} on the CSV",
        catalogue.osmotic_coefficient,
        legacy.osmotic_coefficient
    );
}

/// **A topology the loaded dataset cannot cover is refused, not evaluated at zero.**
///
/// `NH4+`/`Cl-` is a pair neither dataset carries, and the audit is what says so. NeqSim
/// *would* initialize this brine - `validateParameterCoverageOncePerState` only enforces
/// the audit for a mixture with more than one cation or anion - and evaluate the pair at
/// zero, which is the answer this library refuses to give.
#[test]
fn an_uncovered_pair_is_refused() {
    let err = pitzer_phase(
        &["water", "nh4+", "cl-"],
        298.15,
        100000.0,
        &[0.98, 0.01, 0.01],
    )
    .unwrap_err();
    assert!(
        matches!(err, AzothError::InvalidInput { .. }),
        "a pair in neither dataset should be refused: {err:?}"
    );
    let message = err.to_string();
    assert!(
        message.contains("missingBinary=[cl-|nh4+]"),
        "the refusal should name the pair that is absent: {message}"
    );

    // A mixed brine forced onto the legacy dataset is refused too, for the same reason but
    // a different family: that dataset carries no same-sign rows at all.
    let err = pitzer_phase(
        &["water", "na+", "k+", "cl-", "hco3-"],
        298.15,
        100000.0,
        &[0.96, 0.01, 0.01, 0.01, 0.01],
    )
    .unwrap_err();
    let message = err.to_string();
    assert!(
        message.contains("missingTheta=[k+|na+]"),
        "the refusal should name the absent same-sign pair: {message}"
    );
}

/// The neutral layer is a *branch*, not a contribution: a component with no charge takes
/// it, and on a brine the catalogue covers for it, it is the whole of that component's
/// answer.
#[test]
fn a_neutral_solute_takes_the_neutral_branch() {
    // The CO2/sulphate brine is the catalogue's own neutral topology; CO2/chloride is not,
    // which the tranche measured and which a separate test covers.
    let brine = pitzer_phase(
        &["water", "na+", "so4--", "co2"],
        298.15,
        100000.0,
        &[0.86, 0.06, 0.03, 0.05],
    )
    .unwrap();
    assert_eq!(brine.dataset.name(), "phreeqc");
    assert!(
        brine.ln_gamma[3] != 0.0,
        "the catalogue carries CO2's rows, so its activity coefficient is not ideal"
    );
    // Water is on the osmotic route and an ion on the extended expression, so no two of
    // the three branches agree - which is what says they are three.
    assert_ne!(brine.ln_gamma[0], brine.ln_gamma[1]);
    assert_ne!(brine.ln_gamma[1], brine.ln_gamma[3]);
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let err = pitzer_phase(
        &["water", "na+", "cl-"],
        298.15,
        100000.0,
        &[0.88, 0.06, 0.12],
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}

#[test]
fn a_brine_with_no_water_is_refused() {
    let err = pitzer_phase(&["na+", "cl-"], 298.15, 100000.0, &[0.5, 0.5]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
}

/// The model table is the one the spec generates, and it is reachable from the crate root.
#[test]
fn the_model_is_registered() {
    let spec = azoth_model_spec(MODEL_ID);
    assert_eq!(spec.kind, "direct");
    assert!(!spec.cases.is_empty());
}

fn azoth_model_spec(id: &str) -> &'static azoth_core::spec::ModelSpec {
    azoth_eos::model_gen::model(id).expect("the model should be in its own table")
}

/// **The one captured fluid this model refuses**, and the reason is a number rather than a
/// shortcut.
///
/// `water-methanol-NaCl` is the probe's fourth fluid, and NeqSim answers it: methanol has no
/// Henry correlation, the compiled table's four zeros evaluate to `1.802 bar`, and that
/// value goes into `gamma H (m/x) / P` to give `23.606550097141273`. `1.802` is water's own
/// molar-mass product and not a property of methanol, so it is a unit-conversion artefact
/// being read as a Henry constant - and this library refuses it in `eos.kent_eisenberg_phase`
/// and `eos.desmukh_mather_phase` for the same reason.
///
/// **The ion branch does not refuse the same rows**, which is the reason this is a refusal
/// of one branch and not of the model: an ion is capped by its class before the row is read,
/// so `1e12` is what an ion's entry gives whatever its correlation says. That is why the
/// captured ion entries pin and this one cannot.
#[test]
fn a_neutral_with_no_henry_correlation_is_refused_where_neqsim_answers() {
    let err = pitzer_phase(
        &["water", "methanol", "na+", "cl-"],
        313.15,
        5.0e5,
        &[
            0.847457627118644,
            0.0847457627118644,
            0.03389830508474576,
            0.03389830508474576,
        ],
    )
    .unwrap_err();
    let message = err.to_string();
    assert!(message.contains("methanol"), "{message}");
    assert!(
        message.contains("Henry coefficient"),
        "the refusal should name what it cannot build: {message}"
    );

    // The same fluid's ions *do* answer, so the refusal is the neutral arm's alone - and
    // the fluid is otherwise the captured one, at the captured state.
    let brine = pitzer_phase(
        &["water", "na+", "cl-"],
        313.15,
        5.0e5,
        &[
            0.9090909090909091,
            0.045454545454545456,
            0.045454545454545456,
        ],
    )
    .expect("a brine of ions and water has no neutral arm to take");
    assert!(brine.ln_phi[1] > 0.0, "the ionic arm answered");
}
