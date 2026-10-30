//! The phase capacities `eos.sound_speed` and the leak blowdown share.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::cubic::Cubic;
use azoth_eos::databank;
use azoth_eos::heat_capacity::phase_capacity;
use azoth_eos::mixture::RootSide;
use azoth_eos::pt_flash::pt_flash;
use azoth_eos::results::Phase;

fn capacity(components: &[&str], t: f64, p: f64, composition: &[f64], side: RootSide) {
    let (mixture, ideal_gas) =
        databank::mixture_of(components, Cubic::Pr, None).expect("the fluid");
    let reduced = mixture
        .reduced_parameters(kelvins(t), pascals(p))
        .expect("the state");
    let root = mixture
        .phase_state(&reduced, composition, side)
        .expect("the root")
        .z;
    let capacity = phase_capacity(
        &mixture,
        &ideal_gas,
        &reduced,
        kelvins(t),
        pascals(p),
        composition,
        root,
    )
    .expect("the capacities");
    assert!(
        capacity.cp.is_finite() && capacity.cp > 0.0,
        "cp={}",
        capacity.cp
    );
    assert!(
        capacity.cv.is_finite() && capacity.cv > 0.0,
        "cv={}",
        capacity.cv
    );
    let gamma = capacity.cp / capacity.cv;
    assert!((1.0..=2.0).contains(&gamma), "gamma={gamma}");
    assert!(capacity.d_p_d_v < 0.0, "dP/dV={}", capacity.d_p_d_v);
}

#[test]
fn methane_vapour_has_a_physical_capacity_ratio() {
    capacity(&["methane"], 300.0, 500_000.0, &[1.0], RootSide::Vapour);
}

#[test]
fn propane_liquid_has_a_physical_capacity_ratio() {
    capacity(&["propane"], 300.0, 2_000_000.0, &[1.0], RootSide::Liquid);
}

#[test]
fn a_two_phase_split_has_a_capacity_on_each_side() {
    let (mixture, ideal_gas) =
        databank::mixture_of(&["methane", "ethane"], Cubic::Pr, None).expect("the fluid");
    let t = kelvins(200.0);
    let p = pascals(1_000_000.0);
    let z = [0.5, 0.5];
    let flash = pt_flash(&mixture, t, p, &z).expect("the flash");
    assert_eq!(
        flash.phase,
        Phase::TwoPhase,
        "the state is on the saturation line"
    );
    let reduced = mixture.reduced_parameters(t, p).expect("the state");
    for (composition, root) in [(&flash.y, flash.z_vapour), (&flash.x, flash.z_liquid)] {
        let capacity = phase_capacity(&mixture, &ideal_gas, &reduced, t, p, composition, root)
            .expect("the capacities");
        let gamma = capacity.cp / capacity.cv;
        assert!(
            capacity.cp > 0.0 && capacity.cv > 0.0,
            "capacities are positive"
        );
        assert!((1.0..=2.0).contains(&gamma), "gamma={gamma}");
    }
}
