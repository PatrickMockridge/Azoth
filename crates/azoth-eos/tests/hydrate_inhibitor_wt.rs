//! Spec-driven tests for the `eos.hydrate_inhibitor_wt` model.
//!
//! Two oracles, both from `validation/neqsim/HydrateInhibitorProbe.java`: the dose, from
//! `captures/hydrate_inhibitor_probe.tsv`, and the **phase label** its inner step needs, from
//! `captures/phase_type_probe.tsv`.
//!
//! **The dose agrees to machine precision** - `3.3e-16`, `8.9e-16` and `2.2e-16` relative at
//! targets of 0.30, 0.50 and 0.70 - because what is left of this model once the labelling is
//! right is an ordinary cubic flash, and the two libraries already agree on those. Its
//! `...ConcentrationFlash` sibling manages `2e-8`, its inner step being a hydrate equilibrium
//! with tolerances of its own.

// The oracle values are NeqSim's ground truth, so their full digits are the point.
#![allow(clippy::excessive_precision)]

use azoth_core::AzothError;
use azoth_core::units::{kelvins, pascals};
use azoth_eos::hydrate_inhibitor_wt::{PhaseLabel, hydrate_inhibitor_wt, label};
use azoth_eos::{Cubic, model_gen, pt_flash};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.hydrate_inhibitor_wt";

/// NeqSim's own `main` for the sibling class: methane, ethane, propane, i-butane, MEG, water.
const NAMES: [&str; 6] = ["methane", "ethane", "propane", "i-butane", "MEG", "water"];

const MOLES: [f64; 6] = [1.0, 0.10, 0.050, 0.0050, 0.1, 1.0];

/// `(target mass fraction, the capture's plain-SRK inhibitor moles)`.
const DOSES: [(f64, f64); 3] = [
    (0.30, 0.124_375_874_492_093),
    (0.50, 0.290_210_522_059_004),
    (0.70, 0.677_156_064_897_633),
];

/// `(label, components, moles, T in K, P in bara, the probe's types per phase)`, from
/// `PhaseTypeProbe`. The fluids are chosen so each branch of NeqSim's rule is exercised.
type State = (
    &'static str,
    &'static [&'static str],
    &'static [f64],
    f64,
    f64,
    &'static [PhaseLabel],
);

const STATES: &[State] = &[
    (
        "methane gas",
        &["methane"],
        &[1.0],
        300.0,
        50.0,
        &[PhaseLabel::Gas],
    ),
    (
        "lng liquid",
        &["methane", "ethane", "propane"],
        &[0.5, 0.3, 0.2],
        150.0,
        50.0,
        &[PhaseLabel::Oil],
    ),
    (
        "meg and water",
        &["MEG", "water"],
        &[0.1, 1.0],
        300.0,
        1.0,
        &[PhaseLabel::Aqueous],
    ),
    (
        "inhibitor feed",
        &NAMES,
        &MOLES,
        273.15,
        100.0,
        &[PhaseLabel::Gas, PhaseLabel::Aqueous],
    ),
    (
        "co2 liquid",
        &["CO2"],
        &[1.0],
        260.0,
        60.0,
        &[PhaseLabel::Oil],
    ),
];

fn mixture_of(names: &[&str]) -> azoth_eos::mixture::Mixture {
    let (mixture, _) = azoth_eos::databank::mixture_of(names, Cubic::Srk, None).expect("resolves");
    mixture
}

#[test]
fn the_dose_reproduces_the_capture() {
    let fluid = mixture_of(&NAMES);
    for (target, want) in DOSES {
        let result = hydrate_inhibitor_wt(
            &fluid,
            "MEG",
            &MOLES,
            target,
            kelvins(273.15),
            pascals(1.0e7),
        )
        .unwrap_or_else(|e| panic!("{target}: {e:?}"));
        assert!(
            (result.inhibitor_moles / want - 1.0).abs() < 1.0e-12,
            "{target}: {} mol against NeqSim's {want}",
            result.inhibitor_moles
        );
        // **The target is met in the aqueous phase**, which is what the residual is on.
        assert!(
            (result.weight_fraction - target).abs() <= 1.0e-5,
            "{target}: the aqueous phase reached {}",
            result.weight_fraction
        );
        assert_eq!(
            result.phases, 2,
            "{target}: the capture's state has two phases"
        );
    }
}

/// **The label, which is the new machinery.** Each state in `PhaseTypeProbe` lands in one of
/// NeqSim's three branches and the capture records which; this reproduces every one.
#[test]
fn the_phase_labels_reproduce_the_probe() {
    for (name, names, moles, t, p_bar, want) in STATES {
        let fluid = mixture_of(names);
        let total: f64 = moles.iter().sum();
        let z: Vec<f64> = moles.iter().map(|value| value / total).collect();
        let t_q = kelvins(*t);
        let p_q = pascals(p_bar * 1.0e5);
        let flash = pt_flash(&fluid, t_q, p_q, &z).expect("the state flashes");
        let reduced = fluid
            .reduced_parameters(t_q, p_q)
            .expect("reduced parameters");

        let phases: Vec<(Vec<f64>, f64)> = match flash.phase {
            // **Vapour first**, which is the order the probe prints and the order the model's
            // own helper uses.
            azoth_eos::results::Phase::TwoPhase => vec![
                (flash.y.clone(), flash.z_vapour),
                (flash.x.clone(), flash.z_liquid),
            ],
            azoth_eos::results::Phase::AllLiquid => vec![(flash.x.clone(), flash.z_liquid)],
            _ => vec![(flash.y.clone(), flash.z_vapour)],
        };

        let mut labels = Vec::new();
        for (composition, z_factor) in &phases {
            labels.push(
                label(&fluid, &reduced, composition, *z_factor)
                    .unwrap_or_else(|e| panic!("{name}: {e:?}")),
            );
        }
        assert_eq!(
            labels.as_slice(),
            *want,
            "{name}: the labels the probe prints"
        );
    }
}

/// A component with no `COMPTYPE` is refused rather than assumed aqueous, because the branch
/// it decides is between an oil and an aqueous phase.
#[test]
fn a_component_with_no_class_is_refused() {
    let fluid = azoth_eos::mixture::Mixture::new(
        vec![
            azoth_eos::mixture::Component::new(kelvins(190.56), pascals(4.599e6), 0.011)
                .expect("methane's constants")
                .with_molar_mass(Some(0.016043)),
            azoth_eos::mixture::Component::new(kelvins(647.1), pascals(2.2064e7), 0.344)
                .expect("water's constants")
                .with_molar_mass(Some(0.018015)),
        ],
        vec![0.0; 4],
    )
    .expect("the pair resolves");
    let reduced = fluid
        .reduced_parameters(kelvins(300.0), pascals(1.0e5))
        .expect("reduced parameters");
    let error = label(&fluid, &reduced, &[0.5, 0.5], 1.0e-3)
        .expect_err("neither component carries a COMPTYPE");
    assert!(
        matches!(error, AzothError::InvalidInput { .. }),
        "{error:?}"
    );
}

/// An inhibitor the feed does not carry, and a feed with no water, are both refused.
#[test]
fn a_feed_this_cannot_dose_is_refused() {
    let fluid = mixture_of(&NAMES);
    let error = hydrate_inhibitor_wt(
        &fluid,
        "methanol",
        &MOLES,
        0.30,
        kelvins(273.15),
        pascals(1.0e7),
    )
    .expect_err("the feed has no methanol");
    assert!(
        matches!(error, AzothError::InvalidInput { .. }),
        "{error:?}"
    );

    let dry = mixture_of(&["methane", "ethane"]);
    let error = hydrate_inhibitor_wt(
        &dry,
        "methane",
        &[0.5, 0.5],
        0.30,
        kelvins(273.15),
        pascals(1.0e7),
    )
    .expect_err("no water, no aqueous phase");
    assert!(
        matches!(error, AzothError::InvalidInput { .. }),
        "{error:?}"
    );
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty());
    for case in spec.cases {
        let names = case
            .list("components")
            .expect("the case declares components");
        let fluid = mixture_of(names);
        let result = hydrate_inhibitor_wt(
            &fluid,
            common::input_str(case, "inhibitor"),
            case.vector("moles").expect("moles"),
            common::input(case, "wt_target"),
            kelvins(common::input(case, "T")),
            pascals(common::input(case, "P")),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));
        common::assert_close(
            result.inhibitor_moles,
            common::expected(case, "inhibitor_moles"),
            case.tolerance,
            &format!("{}::{} (inhibitor moles)", spec.id, case.id),
        );
        common::assert_close(
            result.weight_fraction,
            common::expected(case, "weight_fraction"),
            case.tolerance,
            &format!("{}::{} (aqueous mass fraction)", spec.id, case.id),
        );
    }
}
