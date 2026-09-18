//! Spec-driven tests for the `eos.ge_wilson_phase` model.

use azoth_core::AzothError;
use azoth_eos::databank::{ge_wilson_phase_parameters, mixture_of};
use azoth_eos::ge_wilson_phase::ge_wilson_phase;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.ge_wilson_phase";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::GeWilsonPhaseResult {
    let names = case.list("components").expect("components");
    let params = ge_wilson_phase_parameters(names, None)
        .unwrap_or_else(|e| panic!("case `{}` should resolve but failed: {e}", case.id));
    let (mixture, _) = mixture_of(names, None).expect("the components resolve");
    ge_wilson_phase(
        &params,
        &mixture,
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

/// The phase is `gamma_i P0_i / P` and nothing else - the same identity every GE phase
/// satisfies, because `ComponentGE.fugcoef` is one method and none of them overrides it.
#[test]
fn the_fugacity_coefficient_is_gamma_times_p0_over_p() {
    let names = ["nc10", "n-octane"];
    let params = ge_wilson_phase_parameters(&names, None).unwrap();
    let (mixture, _) = mixture_of(&names, None).unwrap();
    let (t, p) = (298.15, 100_000.0);
    let r = ge_wilson_phase(&params, &mixture, t, p, &[0.5, 0.5]).unwrap();
    for i in 0..2 {
        let composed = r.gamma[i].ln() + (r.p_sat[i].value / p).ln();
        assert!(
            (r.ln_phi[i] - composed).abs() < 1e-15,
            "component {i}: ln_phi is {} but ln(gamma P0 / P) is {composed}",
            r.ln_phi[i]
        );
    }
}

/// The activity coefficients are `eos.wilson_activity_coefficients`', not a second copy.
#[test]
fn the_activity_coefficients_are_the_wilson_models() {
    let names = ["nc10", "n-octane"];
    let params = ge_wilson_phase_parameters(&names, None).unwrap();
    let (mixture, _) = mixture_of(&names, None).unwrap();
    let r = ge_wilson_phase(&params, &mixture, 298.15, 100_000.0, &[0.5, 0.5]).unwrap();

    let activity = azoth_eos::wilson_activity_coefficients::wilson_activity_coefficients(
        &mixture,
        298.15,
        &[0.5, 0.5],
    )
    .unwrap();
    assert_eq!(r.gamma, activity.gamma);
}

/// The three phases that share `ComponentGE.fugcoef` agree on `P0` and disagree on
/// `gamma`.
///
/// NRTL and UNIFAC were the first two; Wilson is the third. That all three satisfy one
/// identity at one state is what says the shared arithmetic is shared, and that their
/// activity coefficients differ is what says they are three models.
#[test]
fn three_ge_phases_share_the_vapour_pressure_and_not_the_activity() {
    let names = ["methanol", "water"];
    let (t, p) = (298.15, 100_000.0);
    let x = [0.5, 0.5];
    let unifac = azoth_eos::ge_unifac_phase::ge_unifac_phase(
        &azoth_eos::databank::ge_unifac_phase_parameters(&names, None).unwrap(),
        t,
        p,
        &x,
    )
    .unwrap();
    let nrtl = azoth_eos::ge_nrtl_phase::ge_nrtl_phase(
        &azoth_eos::databank::ge_nrtl_phase_parameters(&names, None).unwrap(),
        t,
        p,
        &x,
    )
    .unwrap();
    let (mixture, _) = mixture_of(&names, None).unwrap();
    let wilson = ge_wilson_phase(
        &ge_wilson_phase_parameters(&names, None).unwrap(),
        &mixture,
        t,
        p,
        &x,
    )
    .unwrap();

    assert_eq!(unifac.p_sat, nrtl.p_sat);
    assert_eq!(unifac.p_sat, wilson.p_sat);
    assert_ne!(unifac.gamma, nrtl.gamma);
    assert_ne!(unifac.gamma, wilson.gamma);
    assert_ne!(nrtl.gamma, wilson.gamma);
}

/// A component NeqSim's database tags a Henry's-law solute is refused.
///
/// `n-butane` is one, which is why the pair the Wilson *activity* model is validated on -
/// `n-butane`/`nc12` - is not a pair this *phase* can be stated against.
#[test]
fn a_henrys_law_component_is_refused_rather_than_computed() {
    for name in ["n-butane", "CO2", "methane", "n-hexane"] {
        let err = ge_wilson_phase_parameters(&[name, "n-octane"], None).unwrap_err();
        assert!(
            matches!(err, AzothError::InvalidInput { .. }),
            "{name}: {err:?}"
        );
    }
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let names = ["nc10", "n-octane"];
    let params = ge_wilson_phase_parameters(&names, None).unwrap();
    let (mixture, _) = mixture_of(&names, None).unwrap();
    let err = ge_wilson_phase(&params, &mixture, 298.15, 100_000.0, &[0.6, 0.6]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}
