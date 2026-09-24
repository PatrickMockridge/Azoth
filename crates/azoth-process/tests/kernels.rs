//! The kernels' balance invariants: moles and enthalpy are conserved.

use azoth_core::units::{kelvins, pascals, watts_per_kelvin};
use azoth_process::Stream;
use azoth_process::kernels::{
    filter, heat_exchanger, heater, mixer, pump, separator, splitter, throttling_valve,
};

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

#[test]
fn a_heater_holds_the_temperature_it_is_given() {
    let feed = binary(0.9, 100.0, 3e6, 320.0);
    let out = heater(&feed, Some(kelvins(380.0)), None, None).expect("heater");

    close(out.outlet.t.value, 380.0);
    close(out.outlet.n, feed.n);
    assert_eq!(out.outlet.z, feed.z);
    assert_eq!(out.outlet.p, feed.p);
    // **The duty is the state's, not a number the caller gave.** It is `n * dh` over the
    // same two enthalpies the outlet carries, so it is exactly what a reader of the outlet
    // record would compute - which is `Heater.run`'s own recomputation after its flash.
    close(
        out.duty.value,
        feed.n * (out.outlet.h.value - feed.h.value),
    );
    assert!(out.duty.value > 0.0, "380 K is above the inlet's 320");
}

#[test]
fn an_unstated_heater_is_an_isothermal_drop_and_not_a_throttling() {
    // Neither a temperature nor a duty: `run`'s else branch is `T_in + dT` with `dT` zero,
    // so what is left is the pressure drop - the branch a kernel written from the class's
    // documentation alone would miss.
    //
    // **It holds the temperature and not the enthalpy**, which is the whole difference
    // between this and `throttling_valve`: a real fluid's `h` depends on `P` at fixed `T`,
    // so the outlet enthalpy moves by about `50` J/mol over these two bar and the reported
    // duty is that movement - nonzero, and correctly so, because `run` reports `newH - oldH`
    // for whatever its flash did. An isenthalpic reading of this branch is the mistake.
    let feed = binary(0.9, 100.0, 3e6, 320.0);
    let out = heater(&feed, None, None, Some(pascals(2e5))).expect("heater");

    close(out.outlet.p.value, 2.8e6);
    close(out.outlet.t.value, feed.t.value);
    close(
        out.duty.value,
        feed.n * (out.outlet.h.value - feed.h.value),
    );
}

#[test]
fn a_heater_refuses_a_temperature_and_a_duty_together() {
    // `Heater`'s two setters clear each other's flags, so the class's answer to both is the
    // order they were called in. A model that picked one would be inventing a rule the
    // class does not have, and a model that silently preferred the temperature would report
    // a duty the caller never asked for.
    let feed = binary(0.9, 100.0, 3e6, 320.0);
    assert!(
        heater(&feed, Some(kelvins(380.0)), Some(azoth_core::units::watts(5000.0)), None)
            .is_err()
    );
}

#[test]
fn a_heater_refuses_a_drop_past_zero_pressure() {
    let feed = binary(0.9, 100.0, 3e5, 320.0);
    assert!(heater(&feed, None, None, Some(pascals(3e5))).is_err());
    assert!(heater(&feed, None, None, Some(pascals(4e5))).is_err());
}

#[test]
fn a_filter_drops_the_pressure_and_holds_the_temperature() {
    let feed = binary(0.9, 100.0, 3e6, 320.0);
    let outcome = filter(&feed, pascals(1e5)).expect("filter");

    close(outcome.outlet.p.value, 2.9e6);
    close(outcome.outlet.t.value, feed.t.value);
    close(outcome.applied_drop.value, 1e5);
    close(outcome.outlet.n, feed.n);
    // **The enthalpy moves, and that is the whole difference from `throttling_valve`.** A
    // filter's drop is isothermal, so a real fluid's pressure dependence shows up here; the
    // valve holds `h` and lets `T` fall. Measured on this fluid the two differ by `24.8`
    // J/mol over one bar, which is far above any tolerance this suite compares at.
    assert!(
        (outcome.outlet.h.value - feed.h.value).abs() > 1.0,
        "an isothermal drop moves the enthalpy: {} against {}",
        outcome.outlet.h.value,
        feed.h.value
    );
}

#[test]
fn a_filter_clamps_a_drop_past_the_inlet_pressure_into_a_state_the_cubic_refuses() {
    // The class clamps: `min(max(0, dP), max(0, P_in - 1e-6 bar))`, and logs a warning rather
    // than failing. So the outlet lands a millionth of a bar above vacuum - and **that is a
    // state no cubic can evaluate**: NeqSim extrapolates to `1968.3` J/mol there, and this
    // library's solver does not converge, so the port reproduces the clamp and then refuses
    // the flash it implies.
    //
    // The refusal is the assertion. A kernel that returned a number here would be returning
    // an extrapolation, and the case records the row as uncased for exactly this reason.
    let feed = binary(0.9, 100.0, 3e6, 320.0);
    assert!(
        filter(&feed, pascals(35e5)).is_err(),
        "the clamped state is below what the cubic converges at, and refusing it is the port"
    );
    // A drop that leaves a state the equation of state can be evaluated at is answered, so
    // the refusal above is about the state and not about large numbers.
    assert!(filter(&feed, pascals(29.9e5)).is_ok());
}
