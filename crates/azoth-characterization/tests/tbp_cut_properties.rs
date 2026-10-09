//! `characterization.tbp_cut_properties` against the NeqSim capture.
//!
//! Every expected value is a row of `validation/neqsim/captures/characterization_probe.tsv`,
//! which `CharacterizationProbe` printed from the pinned jar. `pc` is the capture's bar column
//! times `1e5`, because the spec declares pascals and the class computes bar.
//!
//! The cuts are the probe's own: cut[4] is `M = 500` g/mol, `d = 0.88`, which every model runs;
//! cut[3] is `M = 300`, `d = 0.85`, which is where `RiaziDaubert`'s `> 300` switch and
//! `Cavett`'s API correction both fall on the other side.

use azoth_characterization::{TbpModel, tbp_cut_properties};
use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};

/// One row of the capture, as the test states it.
struct Expected {
    model: TbpModel,
    tc: f64,
    pc: f64,
    tb: f64,
    acentric: f64,
    exponent: Option<f64>,
    watson_k: f64,
}

const CUT4: [Expected; 10] = [
    Expected {
        model: TbpModel::PedersenSrk,
        tc: 891.945256085426,
        pc: 1196293.59599742,
        tb: 746.89,
        acentric: 1.36586663262266,
        exponent: Some(2.228110216),
        watson_k: 12.5418146375894,
    },
    Expected {
        model: TbpModel::PedersenSrkHeavyOil,
        tc: 891.945256085426,
        pc: 1196293.59599742,
        tb: 746.89,
        acentric: 1.36586663262266,
        exponent: Some(2.228110216),
        watson_k: 12.5418146375894,
    },
    Expected {
        model: TbpModel::PedersenPr,
        tc: 974.879772951610,
        pc: 1239174.46749789,
        tb: 746.89,
        acentric: 0.526722654332295,
        exponent: Some(1.897865692),
        watson_k: 12.5418146375894,
    },
    Expected {
        model: TbpModel::PedersenPr2,
        tc: 974.879772951610,
        pc: 1239174.46749789,
        tb: 777.612898051266,
        acentric: 0.837080393811892,
        exponent: Some(1.897865692),
        watson_k: 12.7114759871488,
    },
    Expected {
        model: TbpModel::PedersenPrHeavyOil,
        tc: 862.021877542583,
        pc: 1686247.92785472,
        tb: 746.89,
        acentric: 2.39525799802058,
        exponent: Some(-1598.7338588),
        watson_k: 12.5418146375894,
    },
    Expected {
        model: TbpModel::RiaziDaubert,
        tc: 891.945256085426,
        pc: 1196293.59599742,
        tb: 765.014651094374,
        acentric: 1.37706289398442,
        exponent: None,
        watson_k: 12.6424549207867,
    },
    Expected {
        model: TbpModel::LeeKesler,
        tc: 906.211218915845,
        pc: 663946.297092384,
        tb: 785.228404434496,
        acentric: 1.43229568844577,
        exponent: None,
        watson_k: 12.752837620579,
    },
    Expected {
        model: TbpModel::Twu,
        tc: 917.899204856315,
        pc: 724270.609562486,
        tb: 785.228404434496,
        acentric: 1.16668324925980,
        exponent: None,
        watson_k: 12.752837620579,
    },
    Expected {
        model: TbpModel::Cavett,
        tc: 894.305850190858,
        pc: 735341.797507176,
        tb: 765.014651094374,
        acentric: 1.18279252714060,
        exponent: None,
        watson_k: 12.6424549207867,
    },
    Expected {
        model: TbpModel::Standing,
        tc: 942.447643628006,
        pc: 522599.858109319,
        tb: 785.228404434496,
        acentric: 1.21363914242666,
        exponent: None,
        watson_k: 12.752837620579,
    },
];

fn check(molar_mass: f64, density: f64, expected: &Expected) {
    let r = tbp_cut_properties(
        Some(expected.model),
        kilograms_per_mole(molar_mass),
        kilograms_per_cubic_meter(density),
        None,
    )
    .expect("a cut with a positive molar mass and density is in range");

    let tol = 1e-9;
    let name = format!("{:?}", expected.model);
    assert!(
        (r.tc.value - expected.tc).abs() < tol,
        "{name}: tc {} against {}",
        r.tc.value,
        expected.tc
    );
    assert!(
        (r.pc.value - expected.pc).abs() < tol * expected.pc.abs(),
        "{name}: pc {} against {}",
        r.pc.value,
        expected.pc
    );
    assert!(
        (r.boiling_temperature.value - expected.tb).abs() < tol,
        "{name}: tb {} against {}",
        r.boiling_temperature.value,
        expected.tb
    );
    assert!(
        (r.acentric_factor - expected.acentric).abs() < tol * expected.acentric.abs(),
        "{name}: acentric {} against {}",
        r.acentric_factor,
        expected.acentric
    );
    assert!(
        (r.watson_k - expected.watson_k).abs() < tol * expected.watson_k.abs(),
        "{name}: watson_k {} against {}",
        r.watson_k,
        expected.watson_k
    );
    match (r.attraction_exponent, expected.exponent) {
        (None, None) => {}
        (Some(got), Some(want)) => assert!(
            (got - want).abs() < tol * want.abs().max(1.0),
            "{name}: attraction_exponent {got} against {want}"
        ),
        (got, want) => panic!("{name}: attraction_exponent {got:?} against {want:?}"),
    }
}

#[test]
fn every_model_matches_the_capture_on_cut_4() {
    for expected in &CUT4 {
        check(0.5, 880.0, expected);
    }
}

#[test]
fn pedersen_srk_heavy_oil_is_a_no_op() {
    // The class shadows its parent's coefficient fields, so this model must give
    // `PedersenSRK`'s numbers and not the heavy set's. A port that "fixed" the shadowing
    // would fail here and nowhere else.
    let srk = tbp_cut_properties(
        Some(TbpModel::PedersenSrk),
        kilograms_per_mole(0.5),
        kilograms_per_cubic_meter(880.0),
        None,
    )
    .expect("in range");
    let heavy = tbp_cut_properties(
        Some(TbpModel::PedersenSrkHeavyOil),
        kilograms_per_mole(0.5),
        kilograms_per_cubic_meter(880.0),
        None,
    )
    .expect("in range");
    assert_eq!(srk, heavy);
}

#[test]
fn riazi_daubert_takes_its_own_form_at_exactly_300() {
    // The class's test is `molarMass > 300`, so 300 g/mol is *not* over it and the
    // Riazi-Daubert pair applies. Reading the threshold as `>=` agrees at 500 and differs here.
    let r = tbp_cut_properties(
        Some(TbpModel::RiaziDaubert),
        kilograms_per_mole(0.3),
        kilograms_per_cubic_meter(850.0),
        None,
    )
    .expect("in range");
    assert!(
        (r.tc.value - 815.802701858322).abs() < 1e-9,
        "tc {}",
        r.tc.value
    );
    assert!((r.boiling_temperature.value - 644.547792011568).abs() < 1e-9);
    assert!((r.acentric_factor - 0.697599640892113).abs() < 1e-12);
}

#[test]
fn cavett_applies_no_correction_above_api_30() {
    // `d = 0.85` gives API 34.9, over 30, so both corrections are skipped - where cut[4]'s
    // `d = 0.88` gives 29.295 and both apply.
    let r = tbp_cut_properties(
        Some(TbpModel::Cavett),
        kilograms_per_mole(0.3),
        kilograms_per_cubic_meter(850.0),
        None,
    )
    .expect("in range");
    assert!(
        (r.tc.value - 800.246923082386).abs() < 1e-9,
        "tc {}",
        r.tc.value
    );
    assert!(
        (r.pc.value - 1165584.40264508).abs() < 1e-3,
        "pc {}",
        r.pc.value
    );
    assert!((r.acentric_factor - 0.882074292558076).abs() < 1e-12);
}

#[test]
fn a_supplied_boiling_point_overrides_the_correlation() {
    // Every model short-circuits on a positive stored boiling point, which `addTBPfraction`
    // never sets because it leaves the field at zero.
    let r = tbp_cut_properties(
        Some(TbpModel::PedersenSrk),
        kilograms_per_mole(0.5),
        kilograms_per_cubic_meter(880.0),
        Some(azoth_core::units::kelvins(700.0)),
    )
    .expect("in range");
    assert!((r.boiling_temperature.value - 700.0).abs() < 1e-12);
    // `watson_k` follows the boiling point, so it moves with it.
    assert!((r.watson_k - (1.8_f64 * 700.0).cbrt() / 0.88).abs() < 1e-12);
}

#[test]
fn a_skipped_check_is_reported_and_not_a_pass() {
    // The spec's `boiling_point` bound is optional, so with no boiling point the check cannot
    // run and `apply_checks` emits RANGE_CHECK_SKIPPED.
    let r = tbp_cut_properties(
        Some(TbpModel::LeeKesler),
        kilograms_per_mole(0.5),
        kilograms_per_cubic_meter(880.0),
        None,
    )
    .expect("in range");
    assert!(
        r.warnings
            .iter()
            .any(|w| w.code.as_str() == "RANGE_CHECK_SKIPPED"
                && w.field.as_deref() == Some("boiling_point")),
        "expected a RANGE_CHECK_SKIPPED on boiling_point, got {:?}",
        r.warnings
    );
}

#[test]
fn a_non_positive_molar_mass_is_refused() {
    let refused = tbp_cut_properties(
        Some(TbpModel::PedersenSrk),
        kilograms_per_mole(0.0),
        kilograms_per_cubic_meter(880.0),
        None,
    );
    assert!(refused.is_err());
}
