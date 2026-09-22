//! The kernels' balance invariants: moles and enthalpy are conserved.

use azoth_core::units::{kelvins, pascals, watts_per_kelvin};
use azoth_process::Stream;
use azoth_process::kernels::{heat_exchanger, mixer, pump, separator, splitter, throttling_valve};

fn close(a: f64, b: f64) {
    let scale = 1.0 + a.abs() + b.abs();
    assert!((a - b).abs() < 1e-9 * scale, "{a} vs {b}");
}

fn binary(z_methane: f64, n: f64, p: f64, t: f64) -> Stream {
    Stream::from_pt(
        vec!["methane".into(), "n-butane".into()],
        vec![z_methane, 1.0 - z_methane],
        n,
        pascals(p),
        kelvins(t),
    )
    .expect("methane / n-butane resolves and flashes")
}

#[test]
fn a_splitter_conserves_moles_and_state() {
    let feed = binary(0.5, 100.0, 1e5, 300.0);
    let outs = splitter(&feed, &[0.3, 0.7]).expect("splitter");

    assert_eq!(outs.len(), 2);
    close(outs[0].n + outs[1].n, feed.n);
    for out in &outs {
        assert_eq!(out.z, feed.z);
        assert_eq!(out.t, feed.t);
        assert_eq!(out.p, feed.p);
        assert_eq!(out.h, feed.h);
    }
}

#[test]
fn a_mixer_conserves_moles_and_enthalpy() {
    let a = binary(1.0, 50.0, 1e5, 300.0);
    let b = binary(0.0, 50.0, 1e5, 300.0);
    let m = mixer(&[a.clone(), b.clone()], None).expect("mixer");

    close(m.n, 100.0);
    close(m.z[0], 0.5);
    close(m.z[1], 0.5);
    close(m.h.value * m.n, a.h.value * a.n + b.h.value * b.n);
}

#[test]
fn a_separator_conserves_moles_and_energy() {
    let feed = binary(0.5, 100.0, 5e5, 270.0);
    let (vapour, liquid) = separator(&feed, pascals(0.0), 0.0, None).expect("separator");

    close(vapour.n + liquid.n, feed.n);
    // **Energy too, and only because nothing was asked of the vessel.** The flash is at
    // the feed's own temperature and pressure, so the split is the equilibrium one and the
    // sums match whatever the fluid is. A pressure drop or a heat input would move it, and
    // that is a case's business rather than an invariant's.
    close(
        vapour.h.value * vapour.n + liquid.h.value * liquid.n,
        feed.h.value * feed.n,
    );
    if vapour.n > 0.0 && liquid.n > 0.0 {
        for i in 0..2 {
            close(
                vapour.z[i] * vapour.n + liquid.z[i] * liquid.n,
                feed.z[i] * feed.n,
            );
        }
    }
}

#[test]
fn a_separator_carries_vapour_into_the_liquid() {
    let feed = binary(0.5, 100.0, 5e5, 270.0);
    let (plain, _) = separator(&feed, pascals(0.0), 0.0, None).expect("separator");
    let (carried, liquid) = separator(&feed, pascals(0.0), 0.1, None).expect("separator");

    // A tenth of the vapour's moles move, and they move as material: the vapour keeps its
    // composition and the liquid's becomes the mixture of what it had and what arrived.
    close(carried.n, plain.n * 0.9);
    close(carried.z[0], plain.z[0]);
    close(carried.n + liquid.n, feed.n);
    close(
        carried.z[0] * carried.n + liquid.z[0] * liquid.n,
        feed.z[0] * feed.n,
    );
}

#[test]
fn a_valve_is_isenthalpic() {
    let feed = binary(0.5, 100.0, 5e5, 300.0);
    let out = throttling_valve(&feed, pascals(1e5)).expect("valve");

    close(out.h.value, feed.h.value);
    close(out.p.value, 1e5);
    close(out.n, feed.n);
}

#[test]
fn a_heat_exchanger_conserves_energy() {
    let hot = binary(1.0, 50.0, 1e5, 350.0);
    let cold = binary(0.0, 50.0, 1e5, 300.0);
    let (hot_out, cold_out) = heat_exchanger(
        &hot,
        &cold,
        Some(watts_per_kelvin(50.0)),
        "counterflow",
        None,
        None,
    )
    .expect("heat_exchanger");

    // Whatever the rating decides, the two sides' duties are one duty.
    close(
        (hot.h.value - hot_out.h.value) * hot.n,
        (cold_out.h.value - cold.h.value) * cold.n,
    );
    close(hot_out.n, hot.n);
    close(cold_out.n, cold.n);
    close(hot_out.p.value, hot.p.value);
    close(cold_out.p.value, cold.p.value);
}

#[test]
fn a_pinned_outlet_temperature_is_honoured() {
    let hot = binary(1.0, 50.0, 1e5, 350.0);
    let cold = binary(0.0, 50.0, 1e5, 300.0);
    let (hot_out, cold_out) =
        heat_exchanger(&hot, &cold, None, "counterflow", Some(kelvins(330.0)), None)
            .expect("heat_exchanger");

    close(hot_out.t.value, 330.0);
    close(
        (hot.h.value - hot_out.h.value) * hot.n,
        (cold_out.h.value - cold.h.value) * cold.n,
    );
}

#[test]
fn an_exchanger_refuses_the_shapes_it_cannot_answer() {
    let hot = binary(1.0, 50.0, 1e5, 350.0);
    let cold = binary(0.0, 50.0, 1e5, 300.0);
    // Both outlets pinned leaves the other side nothing to solve for.
    assert!(
        heat_exchanger(
            &hot,
            &cold,
            None,
            "counterflow",
            Some(kelvins(330.0)),
            Some(kelvins(340.0)),
        )
        .is_err()
    );
    // The rating without a conductance is NeqSim's 500 W/K, which this will not guess.
    assert!(heat_exchanger(&hot, &cold, None, "counterflow", None, None).is_err());
    // An arrangement NeqSim does not know falls through to counterflow; this refuses.
    assert!(
        heat_exchanger(
            &hot,
            &cold,
            Some(watts_per_kelvin(50.0)),
            "concentric tube counterflow",
            None,
            None,
        )
        .is_err()
    );
}

#[test]
fn a_pump_raises_pressure_and_adds_work() {
    let feed = binary(0.0, 100.0, 5e5, 250.0);
    let out = pump(&feed, pascals(2e6), 0.75).expect("pump");

    close(out.p.value, 2e6);
    close(out.n, feed.n);
    assert!(out.h.value > feed.h.value, "the pump adds work as enthalpy");
}
