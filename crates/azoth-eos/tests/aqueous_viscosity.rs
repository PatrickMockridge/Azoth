//! Spec-driven tests for the `eos.aqueous_viscosity` model.
//!
//! **The port's own numbers, against `validation/neqsim/captures/aqueous_viscosity_probe.tsv`.**
//! The correlation and the pressure correction are closed forms, so the expected values are
//! NeqSim's to every digit its probe printed, and the cases carry `1e-12`.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::Cubic;
use azoth_eos::aqueous_viscosity;
use azoth_eos::databank;
use azoth_eos::model_gen;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.aqueous_viscosity";

fn from_case(case: &azoth_core::spec::TestCase) -> azoth_eos::mixture::Mixture {
    let names = case
        .list("components")
        .expect("the case declares components");
    databank::mixture_of(names, Cubic::Pr, None)
        .expect("the case's components resolve")
        .0
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let mixture = from_case(case);
        let result = aqueous_viscosity(
            &mixture,
            kelvins(common::input(case, "T")),
            pascals(common::input(case, "P")),
            case.vector("z").expect("z"),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        common::assert_close(
            result.viscosity.value,
            common::expected(case, "viscosity"),
            case.tolerance,
            &format!("{}::{} (viscosity)", spec.id, case.id),
        );
        common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
    }
}

/// **The refusal the spec's own assumption is about.** NeqSim's `LIQVISC` rows for `na+`
/// and `cl-` are the same four numbers as methanol's, so a brine computed from them comes
/// out less viscous than pure water - measured on the probe, `8.4026e-4` against `8.5510e-4`
/// at 300 K and 5 bar. This library refuses the state rather than reproducing an answer
/// that rests on another substance's correlation.
///
/// **The mixture is built here rather than resolved by name**, because `databank::mixture_of`
/// refuses an ion one level up - a cubic has no notion of one - so this guard is reachable
/// only for a caller holding the components directly. It is kept for the same reason
/// `reactive_tp_flash` keeps its own: it is the model's contract, and the spec states it.
#[test]
fn a_brine_is_refused_with_the_reason() {
    use azoth_eos::databank::ION;
    use azoth_eos::{Component, Mixture};

    let water = Component::new(kelvins(647.3), pascals(22.089e6), 0.344)
        .expect("water's constants are a state")
        .with_molar_mass(Some(0.018015))
        .with_liquid_viscosity([-27.952757828, 4665.22592993, 0.052323342, -3.8356e-05], 3);
    let salt = Component::new(kelvins(3000.0), pascals(1.0e8), 0.0)
        .expect("the ion's placeholder constants are positive")
        .with_molar_mass(Some(0.02299))
        .with_class(ION.to_string());
    let mixture = Mixture::new(vec![water, salt], vec![0.0; 4])
        .expect("a two-component mixture")
        .with_names(vec!["water".to_string(), "na+".to_string()]);

    let error = aqueous_viscosity(&mixture, kelvins(300.0), pascals(5.0e5), &[0.98, 0.02])
        .expect_err("an ion's viscosity row is not its own");

    let message = format!("{error}");
    assert!(
        message.contains("na+"),
        "the refusal names the ion: {message}"
    );
    assert!(
        message.contains("methanol"),
        "and says whose numbers the row carries: {message}"
    );
}
