//! Spec-driven tests for the `eos.ge_nrtl_phase` model.

use azoth_core::{AzothError, CalcResult};
use azoth_eos::databank::{self, ge_nrtl_phase_parameters};
use azoth_eos::ge_nrtl_phase;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.ge_nrtl_phase";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::GeNrtlPhaseResult {
    let params = ge_nrtl_phase_parameters(case.list("components").expect("components"), None)
        .unwrap_or_else(|e| panic!("case `{}` should resolve but failed: {e}", case.id));
    ge_nrtl_phase(
        &params,
        common::input(case, "T"),
        common::input(case, "P"),
        case.vector("x").expect("x"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = azoth_eos::model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        for name in ["gamma", "ln_gamma", "ln_phi"] {
            let actual = match name {
                "gamma" => &result.gamma,
                "ln_gamma" => &result.ln_gamma,
                _ => &result.ln_phi,
            };
            let expected = case
                .expected_vector(name)
                .unwrap_or_else(|| panic!("the case declares {name}"));
            for (i, (&got, &want)) in actual.iter().zip(expected).enumerate() {
                common::assert_close(
                    got,
                    want,
                    case.tolerance,
                    &format!("{context} ({name}[{i}])"),
                );
            }
        }
        let expected = case
            .expected_vector("p_sat")
            .expect("the case declares p_sat");
        for (i, (&got, &want)) in result.p_sat.iter().zip(expected).enumerate() {
            common::assert_close(
                got.value,
                want,
                case.tolerance,
                &format!("{context} (p_sat[{i}])"),
            );
        }
        if let Some(expected) = case.expected_vector("henry") {
            for (i, (&got, &want)) in result.henry.iter().zip(expected).enumerate() {
                common::assert_close(
                    got.value,
                    want,
                    case.tolerance,
                    &format!("{context} (henry[{i}])"),
                );
            }
        }
        if let Some(expected) = case.expected_vector("gamma_inf") {
            for (i, (&got, &want)) in result.gamma_inf.iter().zip(expected).enumerate() {
                common::assert_close(
                    got,
                    want,
                    case.tolerance,
                    &format!("{context} (gamma_inf[{i}])"),
                );
            }
        }
        common::assert_consistent(&result, context);
    }
}

#[test]
fn the_result_is_clean_at_an_ordinary_state() {
    let params = ge_nrtl_phase_parameters(&["methanol", "water"], None).unwrap();
    let r = ge_nrtl_phase(&params, 298.15, 100_000.0, &[0.5, 0.5]).unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

/// The phase is `gamma_i P0_i / P` and nothing else.
///
/// The identity is the whole of `ComponentGE.fugcoef`, so asserting it against the two
/// factors the result reports is what says the composition happened rather than some
/// other arrangement of the same three numbers.
#[test]
fn the_fugacity_coefficient_is_gamma_times_p0_over_p() {
    let params = ge_nrtl_phase_parameters(&["methanol", "water"], None).unwrap();
    let (t, p) = (298.15, 100_000.0);
    let r = ge_nrtl_phase(&params, t, p, &[0.5, 0.5]).unwrap();
    for i in 0..2 {
        let composed = r.gamma[i].ln() + (r.p_sat[i].value / p).ln();
        assert!(
            (r.ln_phi[i] - composed).abs() < 1e-15,
            "component {i}: ln_phi is {} but ln(gamma P0 / P) is {composed}",
            r.ln_phi[i]
        );
    }
}

/// The activity coefficients are `eos.nrtl_activity_coefficients`', not a second copy.
///
/// The phase resolves its NRTL matrices through the same function that model does, and
/// this is what says so: a phase that built its own matrices could drift from the model
/// whose spec the matrices belong to.
#[test]
fn the_activity_coefficients_are_the_nrtl_models() {
    let names = ["methanol", "water"];
    let phase = ge_nrtl_phase_parameters(&names, None).unwrap();
    let model = databank::nrtl_parameters(&names, None).unwrap();
    assert_eq!(phase.alpha, model.alpha);
    assert_eq!(phase.dij, model.dij);

    let r = ge_nrtl_phase(&phase, 298.15, 100_000.0, &[0.5, 0.5]).unwrap();
    let pure = azoth_eos::nrtl_activity_coefficients(&model, 298.15, &[0.5, 0.5]).unwrap();
    assert_eq!(r.gamma, pure.gamma);
}

/// A component the databank carries no Antoine correlation for is refused.
///
/// The phase's fugacity coefficient is `gamma_i P0_i / P`, and a component without a
/// correlation has no `P0` - a card-added substance is exactly that case, since a card
/// states the parameters a cubic reads.
#[test]
fn a_component_without_a_vapour_pressure_correlation_is_refused() {
    let overlay = {
        let mut card = azoth_eos::databank::Overlay::new();
        card.set_component(
            "unobtainium",
            azoth_eos::databank::ComponentOverride {
                tc: Some(500.0),
                pc: Some(4_000_000.0),
                omega: Some(0.2),
                cp: None,
                association: None,
                ..Default::default()
            },
        );
        card
    };
    let err = ge_nrtl_phase_parameters(&["methanol", "unobtainium"], Some(&overlay)).unwrap_err();
    assert!(
        matches!(err, AzothError::PropertyUnavailable { .. }),
        "{err:?}"
    );
}

/// A solute in water takes the IAPWS aqueous arm, and the table shadows its databank row.
#[test]
fn the_aqueous_arm_reads_the_iapws_table() {
    use azoth_eos::ge_phase::{HenryState, ReferenceState};

    let params = ge_nrtl_phase_parameters(&["water", "CO2"], None).unwrap();
    let ReferenceState::Henry(HenryState::Iapws(_)) = params.references[1] else {
        panic!(
            "CO2 in water should take the IAPWS arm: {:?}",
            params.references[1]
        );
    };
    let r = ge_nrtl_phase(&params, 298.15, 100_000.0, &[0.98, 0.02]).unwrap();
    // IAPWS G7-04 at 298.15 K; the databank row would give 51.65 bar.
    assert!((r.henry[1].value - 165_644_589.179_205_27).abs() < 1e-3);
    assert!((r.gamma_inf[1] - 0.012_921_263_405_098_09).abs() < 1e-15);
    assert_eq!(
        r.henry[0].value, 0.0,
        "the solvent arm reads no Henry coefficient"
    );
    assert_eq!(
        r.gamma_inf[0], 1.0,
        "the solvent arm does not divide by gamma_inf"
    );
    assert_eq!(r.p_sat[1].value, 0.0, "the Henry arm reads no P0");
}

/// Without water the database row answers, and `gamma^inf` is the diluted full composition.
#[test]
fn the_dry_arm_reads_the_database_row() {
    use azoth_eos::ge_phase::{HenryState, ReferenceState};

    let params = ge_nrtl_phase_parameters(&["methanol", "CO2"], None).unwrap();
    let ReferenceState::Henry(HenryState::Database { ion: false, .. }) = params.references[1]
    else {
        panic!(
            "CO2 without water should take the database arm: {:?}",
            params.references[1]
        );
    };
    let r = ge_nrtl_phase(&params, 298.15, 1_000_000.0, &[0.98, 0.02]).unwrap();
    assert!((r.henry[1].value - 2_895_000.000_000_065).abs() < 1e-6);
    assert!((r.gamma_inf[1] - 124_733.561_757_904_7).abs() < 1e-9);
}

/// **Every solute row is fitted now**, where two of these pairs used to be refused. The
/// refusal survives for a row the table fits nothing for, and `henry.rs` tests it on `mdea`.
#[test]
fn every_solute_row_is_fitted_now() {
    use azoth_eos::ge_phase::ReferenceState;

    for (a, b) in [("water", "n-butane"), ("methanol", "n-hexane")] {
        let params = ge_nrtl_phase_parameters(&[a, b], None)
            .unwrap_or_else(|error| panic!("{a}/{b} used to be refused: {error}"));
        assert!(matches!(params.references[1], ReferenceState::Henry(_)));
    }
}

/// **The old `900` sentinel is now a fitted row.** Upstream PR #4045 re-sourced the table,
/// so methanol/methane's Henry is methane's own `723.75 bar` rather than the cap.
#[test]
fn the_old_sentinel_is_now_a_fitted_row() {
    let params = ge_nrtl_phase_parameters(&["methanol", "methane"], None).unwrap();
    let r = ge_nrtl_phase(&params, 298.15, 1_000_000.0, &[0.98, 0.02]).unwrap();
    assert!((r.henry[1].value - 72_374_999.999_999_51).abs() < 1e-3);
}

/// An ion is capped whatever its row says.
#[test]
fn an_ion_is_capped() {
    let params = ge_nrtl_phase_parameters(&["water", "Na+", "Cl-"], None).unwrap();
    let r = ge_nrtl_phase(&params, 298.15, 100_000.0, &[0.96, 0.02, 0.02]).unwrap();
    assert_eq!(r.henry[1].value, 1.0e17);
    assert_eq!(r.henry[2].value, 1.0e17);
    assert!(r.ln_phi.iter().all(|value| value.is_finite()));
}

/// A solvent-tagged IAPWS species in water is refused, because upstream throws.
#[test]
fn a_solvent_tagged_iapws_species_is_refused_in_water() {
    for name in ["argon", "hydrogen", "helium"] {
        let err = ge_nrtl_phase_parameters(&[name, "water"], None).unwrap_err();
        assert!(
            matches!(err, AzothError::PropertyUnavailable { .. }),
            "{name}: {err:?}"
        );
    }
    // The same species without water takes the Raoult arm, as any solvent does.
    ge_nrtl_phase_parameters(&["argon", "methanol"], None).unwrap();
}

/// The substances the database calls `solvent` take the Raoult arm.
#[test]
fn the_solvent_tagged_substances_take_the_raoult_arm() {
    use azoth_eos::ge_phase::ReferenceState;

    for name in ["water", "methanol", "ethanol", "MEG", "acetone"] {
        let params = ge_nrtl_phase_parameters(&[name, "water"], None).unwrap_or_else(|e| {
            panic!("{name} is tagged solvent, so the phase should resolve it: {e}")
        });
        assert!(
            matches!(params.references[0], ReferenceState::Solvent(_)),
            "{name}"
        );
        assert_eq!(params.references.len(), 2, "{name}");
    }
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let params = ge_nrtl_phase_parameters(&["methanol", "water"], None).unwrap();
    let err = ge_nrtl_phase(&params, 298.15, 100_000.0, &[0.6, 0.6]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}

#[test]
fn a_length_mismatch_is_refused() {
    let params = ge_nrtl_phase_parameters(&["methanol", "water"], None).unwrap();
    let err = ge_nrtl_phase(&params, 298.15, 100_000.0, &[1.0]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("components"));
}

/// **A one-component GE phase resolves** (upstream #3837): `PhaseGE.getActivityCoefficientInfDil`
/// passed a literal `2` as the component count, so a one-component phase read a component that
/// does not exist. This port evaluates the reference over the components it has.
#[test]
fn a_one_component_ge_phase_resolves() {
    let params = ge_nrtl_phase_parameters(&["water"], None).expect("water resolves");
    let result = ge_nrtl_phase(&params, 298.15, 1.0e5, &[1.0]).expect("one component");
    assert!(result.ln_phi[0].is_finite(), "{:?}", result.ln_phi);
}
