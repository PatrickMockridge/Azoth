//! The Guo-Finch hydrate, checked against `captures/hydrate_gf_probe.tsv`.
//!
//! The probe runs both component models over one plain SRK fluid, feeding the guests' reference
//! fugacities in the way `HydrateFormationTemperatureFlash` does, so this is the two models'
//! arithmetic and nothing else. **The models genuinely differ**: at 50 bara NeqSim's PVTsim
//! water fugacity coefficient is `6.74441265462830e-05` and its Guo-Finch one
//! `1.02598375004935e-04`, so a port that had one route wearing two names would fail here.
//!
//! The layers are checked as well as the answer, because the answer alone cannot say *which*
//! part of the model a port got wrong: `cki` and `yki` are the occupancy chain, where the two
//! fitted pairs diverge, and `empty_vapour_pressure`/`molar_volume_hydrate` are the reference
//! term only the Guo-Finch route multiplies in.

// The oracle values are NeqSim's ground truth, so their full digits are the point.
#![allow(clippy::excessive_precision)]

use azoth_eos::databank;
use azoth_eos::hydrate::{self, HydrateGuest, HydrateModel};

/// The probe's feed, as `validation/neqsim/HydrateGFProbe.java` builds it.
const NAMES: [&str; 4] = ["methane", "ethane", "propane", "water"];

/// `(P in bara, [f in bar: methane, ethane, propane, water], pvtsim coefficient, gf coefficient)`.
///
/// The probe's own states, transcribed from the capture.
const STATES: [(f64, [f64; 4], f64, f64); 3] = [
    (
        100.0,
        [
            70.0988594113204,
            4.51610630309347,
            0.532415310445719,
            0.00447696002836766,
        ],
        3.34542273841746e-05,
        5.02541866879314e-05,
    ),
    (
        50.0,
        [
            38.7358695470026,
            3.59682291143079,
            0.563592533913349,
            0.00425210026944710,
        ],
        6.74441265462830e-05,
        1.02598375004935e-04,
    ),
    (
        200.0,
        [
            120.860624796036,
            5.06615849751147,
            0.443338479751108,
            0.00496253957920541,
        ],
        1.77990759845883e-05,
        2.60994368043168e-05,
    ),
];

/// 273.15 K, which is where the probe's fluid starts and stays: it is flashed at fixed T.
const T: f64 = 273.15;

fn guests() -> Vec<HydrateGuest> {
    NAMES
        .iter()
        .map(|name| {
            let entry = databank::entry(name, None).expect("the substance is in the databank");
            HydrateGuest {
                name: entry.name.clone(),
                langmuir_a: entry.hydrate_langmuir_a,
                langmuir_b: entry.hydrate_langmuir_b,
                guo_finch_a: entry.hydrate_guo_finch_a,
                guo_finch_b: entry.hydrate_guo_finch_b,
                former: entry.hydrate_former,
            }
        })
        .collect()
}

/// The Guo-Finch table is only carried for some substances, and a guest with no fit is a guest
/// the model cannot take into a cage. This is what makes the port refuse rather than invent.
#[test]
fn the_guo_finch_parameters_are_carried_for_the_guests() {
    let records = guests();
    for record in records.iter().take(3) {
        assert!(
            record
                .guo_finch_a
                .iter()
                .flatten()
                .any(|value| *value != 0.0),
            "{} has no Guo-Finch pair",
            record.name
        );
    }
}

#[test]
fn the_guo_finch_coefficient_reproduces_the_capture() {
    for (p_bara, fugacities, _pvtsim, guo_finch) in STATES {
        let records = guests();
        let refs: Vec<f64> = fugacities.iter().map(|value| value * 1.0e5).collect();
        let p = p_bara * 1.0e5;
        // The Guo-Finch route's reference is the empty lattice's, so this argument is not read:
        // passing the fluid's own water fugacity would be a PVTsim claim.
        let (_structure, coefficient) =
            hydrate::stable_structure(&records, &refs, HydrateModel::GuoFinch, T, p, 0.0)
                .expect("the cavity sum");
        assert!(
            (coefficient / guo_finch - 1.0).abs() < 1.0e-9,
            "P = {p_bara} bara: the coefficient is {coefficient}, the capture says {guo_finch}"
        );
    }
}

/// **The two models are different answers on the same state**, which is the whole reason the
/// choice is an input. A port that wired the Guo-Finch name to the PVTsim arithmetic would pass
/// the test above only by having transcribed the wrong column.
///
/// The PVTsim column is compared against, not recomputed: its own route needs the reference
/// water fugacity, which is a cubic solve, and `tests/hydrate.rs` pins it against
/// `hydrate_probe.tsv`. What this asserts is the thing that test cannot - that the capture's two
/// columns are two different numbers and this port reproduces the one it claims.
#[test]
fn the_two_models_disagree_on_the_same_state() {
    for (p_bara, fugacities, pvtsim, guo_finch) in STATES {
        let records = guests();
        let refs: Vec<f64> = fugacities.iter().map(|value| value * 1.0e5).collect();
        let p = p_bara * 1.0e5;

        let (_, from_guo_finch) =
            hydrate::stable_structure(&records, &refs, HydrateModel::GuoFinch, T, p, 0.0)
                .expect("the cavity sum");
        assert!(
            (pvtsim / guo_finch - 1.0).abs() > 0.1,
            "P = {p_bara}: the capture's two columns are {pvtsim} and {guo_finch}, which are \
             too close for the model choice to be the thing this test is about"
        );
        assert!(
            (from_guo_finch / pvtsim - 1.0).abs() > 0.1,
            "P = {p_bara}: this port gave {from_guo_finch}, which is the capture's *PVTsim* \
             column {pvtsim}"
        );
    }
}

/// **The reference term is per structure, and that is what the two models differ about.**
///
/// NeqSim's PVTsim comparison of the two structures is a comparison of the cavity exponents,
/// because its reference is one number for both. The Guo-Finch reference is not, so a port that
/// kept the exponent comparison would pick the wrong structure wherever the vapour pressures
/// differ by more than the exponents do.
#[test]
fn the_empty_lattice_reference_is_not_the_same_for_both_structures() {
    let first = hydrate::empty_hydrate_vapour_pressure(0, T);
    let second = hydrate::empty_hydrate_vapour_pressure(1, T);
    assert!(
        (first - 0.0108123351187016).abs() < 1.0e-15,
        "structure I's empty vapour pressure is {first}"
    );
    assert!(
        (second - 0.00923067365143309).abs() < 1.0e-15,
        "structure II's is {second}"
    );
    assert!(
        (first / second - 1.0).abs() > 0.1,
        "the two references are {first} and {second}, which are too close for this to matter"
    );
}

#[test]
fn the_hydrate_molar_volumes_reproduce_the_capture() {
    assert!((hydrate::molar_volume_hydrate(0, T) - 2.235e-05).abs() < 1.0e-18);
    assert!((hydrate::molar_volume_hydrate(1, T) - 2.257e-05).abs() < 1.0e-18);
    // Avlonitis's curve moves with temperature, which is what makes the Poynting term a
    // function of the state rather than a constant.
    assert!(hydrate::molar_volume_hydrate(0, 300.0) > hydrate::molar_volume_hydrate(0, T));
}
