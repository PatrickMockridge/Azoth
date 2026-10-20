//! `Mixture::helmholtz_gradient`: the first composition derivative the gradient-theory
//! solvers need, held to two independent things.
//!
//! **The identity.** `d(A^R/RT)/dn_i = ln phi_i + ln Z` is stated in `helmholtz_energy`'s
//! own documentation; here it is checked from the other side, against `phase_state_at`.
//!
//! **The finite difference.** The analytic gradient is also compared with a central
//! difference of `helmholtz_energy` at fixed `Z`, so a slip in the algebra that the
//! identity happens to share - both are built from the same `ln phi` - still fails.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::Cubic;
use azoth_eos::{databank, pt_flash};

fn methane_butane() -> azoth_eos::mixture::Mixture {
    databank::mixture_of(&["methane", "n-butane"], Cubic::Pr, None)
        .expect("the pair resolves")
        .0
}

#[test]
fn the_gradient_is_ln_phi_plus_ln_z() {
    let mixture = methane_butane();
    let t = kelvins(300.0);
    let p = pascals(2.0e6);
    let reduced = mixture.reduced_parameters(t, p).expect("the state reduces");
    let flash = pt_flash::pt_flash(&mixture, t, p, &[0.5, 0.5]).expect("the flash runs");
    for (x, z) in [(&flash.x, flash.z_liquid), (&flash.y, flash.z_vapour)] {
        let gradient = mixture
            .helmholtz_gradient(&reduced, x, z)
            .expect("the gradient runs");
        let state = mixture
            .phase_state_at(&reduced, x, z)
            .expect("the state evaluates");
        for (i, (value, ln_phi)) in gradient.iter().zip(&state.ln_phi).enumerate() {
            let expected = ln_phi + z.ln();
            assert!(
                (value - expected).abs() < 1e-9,
                "component {i}: gradient {value} against ln phi + ln Z {expected}"
            );
        }
    }
}

#[test]
fn the_gradient_matches_a_finite_difference() {
    let mixture = methane_butane();
    let t = kelvins(300.0);
    let p = pascals(2.0e6);
    let reduced = mixture.reduced_parameters(t, p).expect("the state reduces");
    let flash = pt_flash::pt_flash(&mixture, t, p, &[0.5, 0.5]).expect("the flash runs");
    let n = flash.x.clone();
    let z = flash.z_liquid;
    let gradient = mixture
        .helmholtz_gradient(&reduced, &n, z)
        .expect("the gradient runs");
    let step = 1.0e-6;
    for (i, value) in gradient.iter().enumerate() {
        let mut plus = n.clone();
        plus[i] += step;
        let mut minus = n.clone();
        minus[i] -= step;
        let plus_energy = mixture
            .helmholtz_energy(&reduced, &plus, z)
            .expect("the perturbed state evaluates");
        let minus_energy = mixture
            .helmholtz_energy(&reduced, &minus, z)
            .expect("the perturbed state evaluates");
        let numerical = (plus_energy - minus_energy) / (2.0 * step);
        assert!(
            (value - numerical).abs() < 1e-5,
            "component {i}: analytic {value} against finite difference {numerical}"
        );
    }
}
