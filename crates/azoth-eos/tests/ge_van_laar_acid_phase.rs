//! Spec-driven tests for the `eos.ge_van_laar_acid_phase` model.

use azoth_core::AzothError;
use azoth_eos::databank::ge_van_laar_acid_phase_parameters;
use azoth_eos::ge_van_laar_acid_phase::ge_van_laar_acid_phase;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.ge_van_laar_acid_phase";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::GeVanLaarAcidPhaseResult {
    let params =
        ge_van_laar_acid_phase_parameters(case.list("components").expect("components"), None)
            .unwrap_or_else(|e| panic!("case `{}` should resolve but failed: {e}", case.id));
    ge_van_laar_acid_phase(
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
        common::assert_consistent(&result, context);
    }
}

/// The three modelled species satisfy `phi = gamma P0 / P`; the phase is that and nothing
/// else, which is the Raoult identity `ComponentGEVanLaarAcid.fugcoef` exists to enforce.
#[test]
fn the_fugacity_coefficient_is_gamma_times_p0_over_p() {
    let names = ["water", "nitric acid", "sulfuric acid"];
    let params = ge_van_laar_acid_phase_parameters(&names, None).unwrap();
    let (t, p) = (250.0, 100_000.0);
    let r = ge_van_laar_acid_phase(&params, t, p, &[0.5, 0.3, 0.2]).unwrap();
    for i in 0..3 {
        let composed = r.gamma[i].ln() + (r.p_sat[i].value / p).ln();
        assert!(
            (r.ln_phi[i] - composed).abs() < 1e-15,
            "component {i}: ln_phi is {} but ln(gamma P0 / P) is {composed}",
            r.ln_phi[i]
        );
    }
}

/// **This phase does not refuse a Henry's-law solute, and its four siblings do.**
///
/// `ComponentGEVanLaarAcid.fugcoef` overrides the branch that reads `referenceStateType`,
/// and it has to: both acids are tagged `solute` in NeqSim's database, so the inherited
/// method would give them a Henry's-law coefficient instead of the Raoult one the model is
/// written for. A resolver that applied the shared guard here would refuse the model's own
/// subject - which is what this test says, by resolving the pair that guard rejects.
#[test]
fn the_henrys_law_refusal_does_not_apply_here() {
    let params = ge_van_laar_acid_phase_parameters(&["water", "nitric acid"], None).unwrap();
    assert_eq!(
        params.acid_index,
        vec![1, 2],
        "water and nitric acid resolve"
    );

    // The same two names are refused by every other GE phase, because there the tag *is*
    // read. Stated here so the inconsistency is a decision rather than an oversight.
    let refused = azoth_eos::databank::ge_nrtl_phase_parameters(&["water", "nitric acid"], None);
    assert!(
        matches!(refused, Err(AzothError::InvalidInput { .. })),
        "the shared guard should refuse this pair: {refused:?}"
    );
}

/// A component the model does not cover takes the 1.0e12 penalty, whatever its
/// composition.
///
/// NeqSim returns the penalty *before* computing anything else, so its `ln_phi` is
/// `ln(1e12)` rather than `ln(gamma P0 / P)`. It is a penalty rather than physics -
/// driving the species out of the liquid is the intent - and finite rather than infinite
/// so a flash over the mixture can still iterate while it does so.
#[test]
fn an_uncovered_component_takes_the_penalty() {
    let names = ["water", "nitric acid", "nitrogen"];
    let params = ge_van_laar_acid_phase_parameters(&names, None).unwrap();
    assert_eq!(params.acid_index, vec![1, 2, 0]);

    let r = ge_van_laar_acid_phase(&params, 250.0, 100_000.0, &[0.5, 0.3, 0.2]).unwrap();
    let penalty = 1.0e12_f64.ln();
    assert!(
        (r.ln_phi[2] - penalty).abs() < 1e-12,
        "an uncovered component's ln_phi is {}, not ln(1e12) = {penalty}",
        r.ln_phi[2]
    );
    // And it is not the composed expression, which is the thing the penalty replaces.
    let composed = r.gamma[2].ln() + (r.p_sat[2].value / 100_000.0).ln();
    assert_ne!(r.ln_phi[2], composed);
}

/// The three acids take their `P0` from the correlation, not from the databank Antoine.
///
/// This is what makes the phase a composition of two ported pieces rather than one: the
/// acids' stored Antoine rows are the all-zero ones, so a phase reading them would give
/// zero for all three.
#[test]
fn the_acids_vapour_pressure_is_the_correlation_not_antoine() {
    let names = ["water", "nitric acid", "sulfuric acid"];
    let params = ge_van_laar_acid_phase_parameters(&names, None).unwrap();
    // The stored rows are zero, which is *why* the phase has to override `P0`.
    for index in [1usize, 2] {
        assert!(
            params.antoine[index].coefficients.iter().all(|v| *v == 0.0),
            "the databank carries an all-zero Antoine row for {names:?}[{index}], and the \
             phase must not read it"
        );
    }

    let acids = azoth_eos::nitric_sulfuric_acid_vapor_pressure(kelvins(250.0)).unwrap();
    let r = ge_van_laar_acid_phase(&params, 250.0, 100_000.0, &[0.5, 0.3, 0.2]).unwrap();
    assert_eq!(r.p_sat[0].value, acids.p_water.value);
    assert_eq!(r.p_sat[1].value, acids.p_nitric_acid.value);
    assert_eq!(r.p_sat[2].value, acids.p_sulfuric_acid.value);
}

use azoth_core::units::kelvins;

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let params = ge_van_laar_acid_phase_parameters(&["water", "nitric acid"], None).unwrap();
    let err = ge_van_laar_acid_phase(&params, 250.0, 100_000.0, &[0.6, 0.6]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}
