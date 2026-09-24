//! The kernels' balance invariants: moles and enthalpy are conserved.

use azoth_core::units::{kelvins, meters, pascals, watts_per_kelvin};
use azoth_eos::RootSide;
use azoth_process::Stream;
use azoth_process::kernels::ejector::EjectorSetup;
use azoth_process::kernels::three_phase_separator::Entrainment;
use azoth_process::kernels::{
    FlowRegime, compressor, expander, filter, heat_exchanger, heater, mixer, pipe, pump, separator,
    shortcut_distillation_column, splitter, throttling_valve,
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
    let (plain, settled) = separator(&feed, pascals(0.0), 0.0, None).expect("separator");
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

    // **An entrainment also re-runs that outlet, and without one it is left as a phase.**
    // `Separator.run` calls `liquidOutStream.run(id)` under a non-zero fraction and only
    // there, and a stream `run` is a `TPflash` of the carried composition - so the two
    // branches answer with different states, and the difference is the rule rather than
    // noise. The capture's entrainment row is the measurement: `-13999.03` J/mol reported
    // against the liquid root's `-14916.40`.
    let root = Stream::from_side(
        liquid.components.clone(),
        liquid.z.clone(),
        liquid.n,
        liquid.p,
        liquid.t,
        RootSide::Liquid,
    )
    .expect("the carried composition has a liquid root");
    assert!(
        (root.h.value - liquid.h.value).abs() > 100.0,
        "the re-flashed liquid is {} J/mol against the liquid root's {}, which is the \
         distinction this branch is",
        liquid.h.value,
        root.h.value
    );
    // With no entrainment the outlet *is* the phase, so the root is what it reports.
    close(settled.h.value, {
        let side = Stream::from_side(
            settled.components.clone(),
            settled.z.clone(),
            settled.n,
            settled.p,
            settled.t,
            RootSide::Liquid,
        )
        .expect("the flash's liquid has a root");
        side.h.value
    });
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
    close(out.duty.value, feed.n * (out.outlet.h.value - feed.h.value));
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
    close(out.duty.value, feed.n * (out.outlet.h.value - feed.h.value));
}

#[test]
fn a_heater_refuses_a_temperature_and_a_duty_together() {
    // `Heater`'s two setters clear each other's flags, so the class's answer to both is the
    // order they were called in. A model that picked one would be inventing a rule the
    // class does not have, and a model that silently preferred the temperature would report
    // a duty the caller never asked for.
    let feed = binary(0.9, 100.0, 3e6, 320.0);
    assert!(
        heater(
            &feed,
            Some(kelvins(380.0)),
            Some(azoth_core::units::watts(5000.0)),
            None
        )
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

#[test]
fn a_compressor_raises_the_pressure_along_the_inlet_entropy() {
    let feed = binary(0.9, 100.0, 3e6, 320.0);
    let reversible = compressor(&feed, pascals(6e6), 1.0).expect("compressor");
    let real = compressor(&feed, pascals(6e6), 0.75).expect("compressor");

    close(reversible.p.value, 6e6);
    close(reversible.n, feed.n);
    assert!(reversible.t.value > feed.t.value, "compression heats");
    assert!(
        reversible.h.value > feed.h.value,
        "compression adds enthalpy"
    );

    // **The efficiency divides the step, and the two runs pin it without the oracle.**
    // `h_out = h_in + (h_is - h_in) / eta` with the reversible run's answer *being* `h_is`,
    // so the real run's step is the reversible one over 0.75 to the last digit.
    let isentropic_step = reversible.h.value - feed.h.value;
    close(real.h.value - feed.h.value, isentropic_step / 0.75);
}

#[test]
fn an_efficiency_of_one_leaves_the_entropy_alone() {
    // **The isentropic claim, as an invariant rather than a number** - and the reason the
    // record needs no `s` field: the entropy is derived from `(T, P, z)` at both ends, so a
    // reversible step is one whose computed entropies agree.
    let feed = binary(0.9, 100.0, 3e6, 320.0);
    for machine in [
        compressor(&feed, pascals(6e6), 1.0),
        expander(&feed, pascals(1.5e6), 1.0),
    ] {
        let out = machine.expect("the machine runs");
        close(
            out.entropy().expect("the outlet's entropy"),
            feed.entropy().expect("the inlet's"),
        );
    }
    // And an irreversible one does not, which is what makes the invariant worth asserting:
    // at 0.75 the compressor's outlet carries more entropy than it came in with.
    let real = compressor(&feed, pascals(6e6), 0.75).expect("compressor");
    assert!(
        real.entropy().expect("entropy") > feed.entropy().expect("entropy"),
        "an efficiency below one produces entropy"
    );
}

#[test]
fn an_expander_multiplies_the_step_where_a_compressor_divides_it() {
    // `Expander.run` multiplies by the efficiency where `Compressor.run` divides, and the
    // reason is the sign: an expansion's isentropic difference is negative, so dividing by an
    // efficiency below one would recover *more* work than the reversible machine.
    let feed = binary(0.9, 100.0, 6e6, 320.0);
    let reversible = expander(&feed, pascals(3e6), 1.0).expect("expander");
    let real = expander(&feed, pascals(3e6), 0.75).expect("expander");

    assert!(
        reversible.h.value < feed.h.value,
        "an expansion removes enthalpy"
    );
    assert!(reversible.t.value < feed.t.value, "an expansion cools");

    let isentropic_step = reversible.h.value - feed.h.value;
    assert!(
        isentropic_step < 0.0,
        "the step is negative: {isentropic_step}"
    );
    close(real.h.value - feed.h.value, isentropic_step * 0.75);
    // The direction the division would take it, stated as the assertion it fails: a divided
    // step is *below* the reversible one, i.e. more work out than reversible.
    assert!(
        isentropic_step / 0.75 < isentropic_step,
        "dividing would beat the reversible machine"
    );
}

/// The line's own three rows, against `validation/neqsim/captures/process_pipe.tsv`.
///
/// **The gas row is a different fluid from the other captures in this tier, and the reason
/// is a measurement**: 0.9/0.1 methane/n-butane at 320 K and 30 bara is the one state where
/// azoth's cubic lands on a different root from NeqSim's (`Z = 0.87296` against `0.91976`,
/// recorded in `tests/stream.rs`), and a pipe reads `Z` for its velocity and its
/// `P1^2 - P2^2` term - so a gas row on that fluid would measure `eos.pt_flash` rather than
/// the hydraulics.
#[test]
fn a_pipe_solves_a_gas_line() {
    let feed = Stream::from_pt(
        vec!["methane".into(), "CO2".into()],
        vec![0.7, 0.3],
        1.0,
        pascals(50.0e5),
        kelvins(300.0),
    )
    .expect("methane/CO2 resolves");
    let out = pipe(&feed, meters(1000.0), meters(0.1), meters(1.0e-5)).expect("pipe");

    assert!((out.outlet.p.value - 49.99978227836236e5).abs() / 50.0e5 < 1e-6);
    assert!((out.velocity - 0.05491620934533985).abs() / 0.05491620934533985 < 1e-4);
    assert!((out.reynolds - 21377.405123798144).abs() / 21377.405123798144 < 1e-4);
    assert!((out.friction_factor - 0.025488330998681204).abs() < 1e-6);
    assert_eq!(out.regime, FlowRegime::Turbulent);
}

/// **The liquid row, and the two densities it takes.** This reproduces the capture's state
/// to fourteen digits, and it is the row that found the quirk: NeqSim's *velocity* takes the
/// physical-properties density (`567.33`, translated) and its *Reynolds number* takes
/// `getKinematicViscosity()`, which divides by the untranslated cubic (`603.86`). A port
/// that used one density for both is 6% out on the Reynolds number.
#[test]
fn a_pipe_solves_a_liquid_line_with_two_densities() {
    let feed = Stream::from_pt(
        vec!["n-butane".into()],
        vec![1.0],
        1.0,
        pascals(20.0e5),
        kelvins(300.0),
    )
    .expect("n-butane resolves");
    let out = pipe(&feed, meters(1000.0), meters(0.1), meters(1.0e-5)).expect("pipe");

    assert!((out.reynolds - 4912.595137899346).abs() / 4912.595137899346 < 1e-9);
    assert!((out.friction_factor - 0.038002601910602105).abs() / 0.038002601910602105 < 1e-9);
    assert!((out.outlet.p.value - 1999981.6572976355).abs() / 20.0e5 < 1e-9);
    assert_eq!(out.regime, FlowRegime::Turbulent);
}

/// **The water row, which was the gap and is now closed.** NeqSim's aqueous phase takes
/// `WaterPhysicalProperties` and the liquid `Viscosity` correlation (`eos.aqueous_viscosity`);
/// before that id existed this row diverged by `1.61` on the Reynolds number, because
/// `eos.viscosity` is the PFCT form the *gas* and *oil* branches use and is 38% low on water.
/// Now it reproduces the capture: `Re = 231.24476346337397` and the laminar `64 / Re`.
#[test]
fn a_water_line_takes_the_aqueous_viscosity() {
    let feed = Stream::from_pt(
        vec!["water".into()],
        vec![1.0],
        1.0,
        pascals(5.0e5),
        kelvins(300.0),
    )
    .expect("water resolves");
    let out = pipe(&feed, meters(1000.0), meters(0.1), meters(1.0e-5)).expect("pipe");

    assert!((out.velocity - 0.002331069853012549).abs() / 0.002331069853012549 < 1e-6);
    assert!(
        (out.reynolds - 231.24476346337397).abs() / 231.24476346337397 < 1e-6,
        "the aqueous branch, against the capture's own Reynolds number: {}",
        out.reynolds
    );
    assert!((out.friction_factor - 0.2767630239122657).abs() < 1e-6);
    assert_eq!(out.regime, FlowRegime::Laminar);
    assert!((out.outlet.p.value - 499992.6009196372).abs() / 5.0e5 < 1e-9);
}

/// The manifold's rows, printed so the case can be written from azoth's own numbers.
#[test]
fn a_manifold_joins_and_divides() {
    let first = Stream::from_pt(
        vec!["methane".into(), "n-butane".into()],
        vec![0.9, 0.1],
        1.0,
        pascals(30.0e5),
        kelvins(320.0),
    )
    .expect("first feed");
    let second = Stream::from_pt(
        vec!["methane".into(), "n-butane".into()],
        vec![0.3, 0.7],
        2.0,
        pascals(10.0e5),
        kelvins(300.0),
    )
    .expect("second feed");

    let outs =
        azoth_process::kernels::manifold::manifold(&[first.clone(), second.clone()], &[0.25, 0.75])
            .expect("manifold");
    println!("outs {}", outs.len());
    for (i, out) in outs.iter().enumerate() {
        println!(
            "out{i}: n={} p={} t={} h={} z={:?}",
            out.n, out.p.value, out.t.value, out.h.value, out.z
        );
    }
    println!("neqsim products0_n=0.75 products1_n=2.25 T=292.53258807480535 P=1000000 z=0.5/0.5");
    println!("neqsim product_n=3.0 h=-6804.587611589724 (its branch h: -42407.18, -8650.65)");

    // The zero-flow row: the first feed is stated at nothing, so the mixture is the second.
    let mut empty = first.clone();
    empty.n = 0.0;
    let outs = azoth_process::kernels::manifold::manifold(&[empty, second.clone()], &[0.5, 0.5])
        .expect("manifold");
    println!(
        "zero-flow row: out0 n={} t={} h={} z={:?}",
        outs[0].n, outs[0].t.value, outs[0].h.value, outs[0].z
    );
    println!("neqsim: n=1.0, T=299.9999999999398, h=-10824.760384299883, z=0.3/0.7");
}

// ==================== the shortcut distillation column ====================

/// A four-component PR feed at 20 bara, the shape every row below shares.
fn column_feed(t: f64) -> Stream {
    Stream::from_pt(
        vec![
            "methane".into(),
            "ethane".into(),
            "propane".into(),
            "n-butane".into(),
        ],
        vec![0.1, 0.3, 0.4, 0.2],
        1.0,
        pascals(20.0e5),
        kelvins(t),
    )
    .expect("the fluid resolves")
}

fn relative(a: f64, b: f64, tolerance: f64, what: &str) {
    let scale = 1.0 + a.abs().max(b.abs());
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{what}: {a} vs {b}, {} relative",
        (a - b).abs() / scale
    );
}

/// An absolute comparison, for a quantity that is a difference and can sit near zero.
///
/// **The two libraries' enthalpies differ by a systematic offset that is absolute, not
/// relative** - `tests/stream.rs` measures it at -0.0634 J/mol at 300 K, linear in
/// temperature - so a relative bound is the wrong measure wherever the value is small. A
/// product whose enthalpy is -13 J/mol would fail at one part in ten thousand while being
/// 0.06 J/mol out.
fn absolute(a: f64, b: f64, tolerance: f64, what: &str) {
    assert!(
        (a - b).abs() <= tolerance,
        "{what}: {a} vs {b}, {} absolute",
        (a - b).abs()
    );
}

/// **The capture's first row, and the whole FUG chain on it.** A two-phase feed, so the
/// K-values are the flash's `y/x` and every number below is a rearrangement of them:
/// measured against NeqSim, the relative volatility through to the duties agree to `3e-12`.
/// That is the point of porting this class first - it is closed form, so an agreement this
/// tight is plumbing and a disagreement would be a plumbing bug rather than a physics one.
///
/// **The two enthalpies are the one place the two libraries part**, and by the systematic
/// ideal-gas offset `tests/stream.rs` records rather than by anything here: azoth is 0.0636
/// J/mol below NeqSim on the distillate at 300 K, which is that offset to three figures.
/// `reboiler_duty` inherits a thousandth of it, because the class's own estimate adds one per
/// cent of the feed's enthalpy.
#[test]
fn a_shortcut_column_splits_a_four_component_feed() {
    let out = shortcut_distillation_column(
        &column_feed(300.0),
        "propane",
        "n-butane",
        0.98,
        0.98,
        1.2,
        None,
        None,
    )
    .expect("the column solves");

    relative(
        out.relative_volatility,
        2.7769364052320538,
        1e-6,
        "relative_volatility",
    );
    relative(
        out.minimum_stages,
        7.620946291100026,
        1e-6,
        "minimum_stages",
    );
    relative(
        out.minimum_reflux_ratio,
        0.44008894984866465,
        1e-6,
        "minimum_reflux_ratio",
    );
    relative(
        out.actual_reflux_ratio,
        0.5281067398183975,
        1e-6,
        "actual_reflux_ratio",
    );
    relative(out.actual_stages, 20.510953860048968, 1e-6, "actual_stages");
    assert_eq!(out.feed_tray_number, 14, "the feed tray, which is a count");
    relative(
        out.condenser_duty.value,
        -36472.85166598551,
        1e-6,
        "condenser_duty",
    );
    relative(
        out.reboiler_duty.value,
        36378.53354647986,
        1e-6,
        "reboiler_duty",
    );

    relative(out.distillate.n, 0.7956000000005571, 1e-6, "distillate_n");
    let expected_z = [
        0.12556561086035364,
        0.3766968325795219,
        0.49270990447366647,
        0.005027652086458033,
    ];
    for (value, expected) in out.distillate.z.iter().zip(expected_z) {
        relative(*value, expected, 1e-6, "distillate_z");
    }
    assert_eq!(out.distillate.p.value, 20.0e5, "the feed's pressure");
    assert_eq!(out.distillate.t.value, 300.0, "the feed's temperature");
    relative(
        out.distillate.h.value,
        -2243.812676863934,
        1e-4,
        "distillate_h",
    );

    relative(out.bottoms.n, 0.20439999999944294, 1e-6, "bottoms_n");
    let expected_bottoms_z = [
        0.0004892367906107846,
        0.0014677103718262435,
        0.03913894324859146,
        0.9589041095889715,
    ];
    for (value, expected) in out.bottoms.z.iter().zip(expected_bottoms_z) {
        relative(*value, expected, 1e-6, "bottoms_z");
    }
    relative(out.bottoms.h.value, -18697.707831107058, 1e-4, "bottoms_h");

    // The balance the class's own split fractions close, component by component.
    for i in 0..4 {
        relative(
            out.distillate.n * out.distillate.z[i] + out.bottoms.n * out.bottoms.z[i],
            column_feed(300.0).n * column_feed(300.0).z[i],
            1e-9,
            "the component balance",
        );
    }
}

/// The same feed and keys with both product pressures stated.
///
/// **Every FUG number is identical**, because the class reads the two pressures only where it
/// creates the products. What moves is the products' state, and the distillate's enthalpy
/// goes from -2243.81 to -13.32 J/mol - not two bar of compression, but a phase change: at 18
/// bara the distillate flashes to a single vapour, at the feed's own temperature.
#[test]
fn a_stated_product_pressure_moves_the_product_phase() {
    let out = shortcut_distillation_column(
        &column_feed(300.0),
        "propane",
        "n-butane",
        0.98,
        0.98,
        1.2,
        Some(pascals(18.0e5)),
        Some(pascals(21.0e5)),
    )
    .expect("the column solves");

    relative(
        out.minimum_stages,
        7.620946291100026,
        1e-6,
        "minimum_stages",
    );
    assert_eq!(out.distillate.p.value, 18.0e5);
    assert_eq!(out.bottoms.p.value, 21.0e5);
    assert_eq!(
        out.distillate.t.value, 300.0,
        "the flash holds the temperature"
    );
    // An absolute bound, because the distillate is nearly a gas here and its enthalpy is
    // close to zero: azoth is -13.385 against NeqSim's -13.322, which is the libraries'
    // ideal-gas offset and a *relative* error of four parts in a thousand.
    absolute(
        out.distillate.h.value,
        -13.322032609178178,
        0.1,
        "distillate_h",
    );
    relative(out.bottoms.h.value, -18694.58756409396, 1e-4, "bottoms_h");
}

/// **The Wilson branch.** At 450 K the flash finds one phase, so the class has no `y` to
/// divide by an `x` and estimates the K-values instead - `(Pc/P) exp(5.373 (1 + omega) (1 -
/// Tc/T))`, with the class's own `5.373`. The feed quality is zero, because the one phase is
/// a gas.
#[test]
fn a_single_phase_feed_takes_the_wilson_estimates() {
    let out = shortcut_distillation_column(
        &column_feed(450.0),
        "propane",
        "n-butane",
        0.98,
        0.98,
        1.2,
        None,
        None,
    )
    .expect("the column solves");

    relative(
        out.relative_volatility,
        2.360742641862776,
        1e-6,
        "relative_volatility",
    );
    relative(
        out.minimum_stages,
        9.061531808008356,
        1e-6,
        "minimum_stages",
    );
    relative(
        out.minimum_reflux_ratio,
        0.8355991286740769,
        1e-6,
        "minimum_reflux_ratio",
    );
    relative(out.actual_stages, 22.4423070346342, 1e-6, "actual_stages");
    assert_eq!(out.feed_tray_number, 15);
    assert_eq!(out.distillate.t.value, 450.0);
}

/// A binary, which is the row a reader can check by hand.
#[test]
fn a_shortcut_column_solves_a_binary() {
    let feed = Stream::from_pt(
        vec!["methane".into(), "n-butane".into()],
        vec![0.5, 0.5],
        1.0,
        pascals(20.0e5),
        kelvins(300.0),
    )
    .expect("the fluid resolves");
    let out =
        shortcut_distillation_column(&feed, "methane", "n-butane", 0.99, 0.99, 1.2, None, None)
            .expect("the column solves");

    relative(
        out.relative_volatility,
        46.648767420893286,
        1e-6,
        "relative_volatility",
    );
    relative(
        out.minimum_stages,
        2.391643282325193,
        1e-6,
        "minimum_stages",
    );
    relative(out.actual_stages, 8.193966146245264, 1e-6, "actual_stages");
    assert_eq!(out.feed_tray_number, 5);
    relative(out.distillate.z[0], 0.99, 1e-9, "the light key's recovery");
    relative(out.bottoms.z[1], 0.99, 1e-9, "the heavy key's recovery");
}

/// **The class's refusal, and the two degeneracies it does not refuse.** A light key less
/// volatile than the heavy key is `solved = false` with every answer at its field
/// initialiser, which is a refusal and not an answer. A reflux multiplier of exactly one
/// leaves Gilliland's `X` at zero and the class returns `actual_stages = Infinity` with a
/// feed tray of zero - `(int) Math.round(Infinity) + 1` wrapping through `Integer.MIN_VALUE` -
/// and below one it silently takes its `Y = 0.5` fallback. Neither is a column, so both are
/// refused here, and both rows stay in the capture as evidence.
#[test]
fn a_shortcut_column_refuses_what_the_class_cannot_answer() {
    let swapped = shortcut_distillation_column(
        &column_feed(300.0),
        "n-butane",
        "propane",
        0.98,
        0.98,
        1.2,
        None,
        None,
    );
    assert!(swapped.is_err(), "a light key below the heavy key");

    let multiplier_one = shortcut_distillation_column(
        &column_feed(300.0),
        "propane",
        "n-butane",
        0.98,
        0.98,
        1.0,
        None,
        None,
    );
    assert!(
        multiplier_one.is_err(),
        "a reflux multiplier at the minimum"
    );

    let unknown_key = shortcut_distillation_column(
        &column_feed(300.0),
        "butane",
        "n-butane",
        0.98,
        0.98,
        1.2,
        None,
        None,
    );
    assert!(unknown_key.is_err(), "a key that is not a component");
}

/// The component splitter's rows, printed so the case can be written from azoth's numbers.
#[test]
fn a_component_splitter_routes_each_component() {
    let feed = Stream::from_pt(
        vec!["methane".into(), "n-butane".into(), "n-pentane".into()],
        vec![0.5, 0.3, 0.2],
        1.0,
        pascals(20.0e5),
        kelvins(300.0),
    )
    .expect("the three components resolve");
    let (overhead, bottoms) =
        azoth_process::kernels::component_splitter::component_splitter(&feed, &[0.98, 0.05, 0.02])
            .expect("component splitter");

    println!(
        "overhead n={} t={} h={} z={:?}",
        overhead.n, overhead.t.value, overhead.h.value, overhead.z
    );
    println!(
        "bottoms n={} t={} h={} z={:?}",
        bottoms.n, bottoms.t.value, bottoms.h.value, bottoms.z
    );
    println!(
        "neqsim overhead n=0.509 z=0.962671905697446/0.029469548133595282/0.007858546168958742 h=551.7022676649342"
    );
    println!(
        "neqsim bottoms n=0.491 z=0.02036659877800409/0.5804480651731161/0.39918533604887985 h=-14431.31737097011"
    );
    assert_eq!(overhead.z.len(), 3);
    assert_eq!(bottoms.z.len(), 3);
    // ---- and the other two rows, which the case pins too ----
    let (even_h, even_b) =
        azoth_process::kernels::component_splitter::component_splitter(&feed, &[0.5, 0.5, 0.5])
            .expect("component splitter");
    println!(
        "even: n={} h={} z={:?} | n={} h={} z={:?}",
        even_h.n, even_h.h.value, even_h.z, even_b.n, even_b.h.value, even_b.z
    );
    let (all_h, all_b) =
        azoth_process::kernels::component_splitter::component_splitter(&feed, &[1.0, 0.5, 0.0])
            .expect("component splitter");
    println!(
        "all-of-one: n={} h={} z={:?} | n={} h={} z={:?}",
        all_h.n, all_h.h.value, all_h.z, all_b.n, all_b.h.value, all_b.z
    );
}

/// **The tank's two rows, and the two things the kernel is.**
///
/// The first is the claim it rests on: `VUflash` at a fluid's own volume and internal energy
/// returns that fluid's state, so a tank with no drop and no heat input *is* the flash it
/// holds - and the captured two-phase row is `process.separator`'s first row to the last
/// digit. The second is the row that found `Stream::from_side`.
#[test]
fn a_tank_is_the_flash_it_holds() {
    let feed = binary(0.7, 1.0, 20.0e5, 300.0);
    let joined = mixer(std::slice::from_ref(&feed), None).expect("mixer");
    let (gas, liquid) =
        azoth_process::kernels::tank::tank(std::slice::from_ref(&feed)).expect("tank");
    let (vapour, condensate) = separator(&joined, pascals(0.0), 0.0, None).expect("separator");

    println!(
        "gas n={} t={} h={} z={:?}",
        gas.n, gas.t.value, gas.h.value, gas.z
    );
    println!(
        "liquid n={} t={} h={} z={:?}",
        liquid.n, liquid.t.value, liquid.h.value, liquid.z
    );
    println!("neqsim gas n=0.8182211906421442 h=530.1529163915924");
    println!("neqsim liquid n=0.18177880935785584 h=-17268.958496924704");

    // The tank is `process.separator` at zero drop and no entrainment, field for field.
    close(gas.n, vapour.n);
    assert_eq!(gas.z, vapour.z);
    close(gas.h.value, vapour.h.value);
    close(liquid.n, condensate.n);
    assert_eq!(liquid.z, condensate.z);

    // Moles are exact. The energy is closed to the phases' own roots and not to the feed's
    // flash, which is what `Stream::from_side` is: the residual is measured, not assumed.
    close(gas.n + liquid.n, feed.n);
    let balance = gas.n * gas.h.value + liquid.n * liquid.h.value;
    println!(
        "energy residual = {} J/mol of feed",
        balance - feed.n * feed.h.value
    );
    absolute(
        balance,
        feed.n * feed.h.value,
        1.0,
        "the tank's energy balance",
    );
}

/// **A tank's outlet is a *phase*, and re-flashing its composition is a different question.**
///
/// On the water-bearing feed the two part: the gas's composition, flashed on its own, splits
/// again at a vapour fraction of `0.767` and reports `-10043.91` J/mol, where the phase's own
/// root is `441.77` - against NeqSim's `441.83` from `setThermoSystemFromPhase`. The last
/// assertion is the one that keeps this: an outlet built by a re-flash fails it and says why.
#[test]
fn a_tanks_outlet_is_the_phase_it_was_split_into() {
    let feed = Stream::from_pt(
        vec!["methane".into(), "n-butane".into(), "water".into()],
        vec![0.5, 0.3, 0.2],
        1.0,
        pascals(20.0e5),
        kelvins(300.0),
    )
    .expect("the three components resolve");
    let (gas, liquid) =
        azoth_process::kernels::tank::tank(std::slice::from_ref(&feed)).expect("tank");

    println!("gas n={} h={} z={:?}", gas.n, gas.h.value, gas.z);
    println!(
        "liquid n={} h={} z={:?}",
        liquid.n, liquid.h.value, liquid.z
    );
    println!("neqsim gas n=0.8036298485392942 h=441.8331830403947 w=0.2343073049287422");
    println!("neqsim liquid n=0.1963701514607058 h=-16971.673011528655");

    // The aqueous phase is not an outlet: a tank never asks for a third one, so the water
    // leaves in the gas at `0.234` where a `ThreePhaseSeparator` leaves it at `0.0016`.
    close(gas.n + liquid.n, feed.n);
    relative(gas.z[2], 0.2343073049287422, 1e-9, "the gas's water");
    relative(gas.n, 0.8036298485392942, 1e-9, "the gas flow");
    relative(liquid.n, 0.1963701514607058, 1e-9, "the liquid flow");

    // The same composition, re-flashed on its own, is a different state - the divergence
    // this port does not reproduce.
    let reflashed = Stream::from_pt(gas.components.clone(), gas.z.clone(), gas.n, gas.p, gas.t)
        .expect("the gas's composition resolves on its own");
    println!("re-flashed h={}", reflashed.h.value);
    assert!(
        (reflashed.h.value - gas.h.value).abs() > 1000.0,
        "a re-flash of the outlet's composition is {} J/mol against the phase's own {}",
        reflashed.h.value,
        gas.h.value
    );
}

/// **A three-phase separator splits a feed three ways, and it is the only kernel here that
/// does.** The capture's first row is the state to reproduce; the invariants are that the
/// three outlets carry the feed's moles and that the two entrainment directions move what
/// they say they move.
#[test]
fn a_three_phase_separator_splits_a_feed_three_ways() {
    let feed = Stream::from_pt(
        vec!["methane".into(), "n-butane".into(), "water".into()],
        vec![0.5, 0.3, 0.2],
        1.0,
        pascals(20.0e5),
        kelvins(300.0),
    )
    .expect("the three components resolve");
    let none = Entrainment::default();
    let (vapour, oil, aqueous) =
        azoth_process::kernels::three_phase_separator::three_phase_separator(
            &feed,
            pascals(0.0),
            None,
            none,
        )
        .expect("the separator splits three ways");

    println!(
        "vapour n={} h={} z={:?}",
        vapour.n, vapour.h.value, vapour.z
    );
    println!("oil n={} h={} z={:?}", oil.n, oil.h.value, oil.z);
    println!(
        "aqueous n={} h={} z={:?}",
        aqueous.n, aqueous.h.value, aqueous.z
    );
    println!("neqsim vapour n=0.5742666398533198 h=530.2999342332546");
    println!("neqsim oil n=0.22670512468681528 h=-17267.597049152402");
    println!("neqsim aqueous n=0.199028235459865 h=-44702.48623278945");

    close(vapour.n + oil.n + aqueous.n, feed.n);
    for i in 0..3 {
        close(
            vapour.z[i] * vapour.n + oil.z[i] * oil.n + aqueous.z[i] * aqueous.n,
            feed.z[i] * feed.n,
        );
    }
    relative(vapour.n, 0.5742666398533198, 1e-9, "the vapour flow");
    relative(oil.n, 0.22670512468681528, 1e-9, "the oil flow");
    relative(aqueous.n, 0.199028235459865, 1e-9, "the aqueous flow");
    // The three phases are what they are: a gas, a hydrocarbon liquid and a water-rich one.
    relative(
        vapour.z[0],
        0.8323862772745304,
        1e-9,
        "the vapour's methane",
    );
    relative(oil.z[1], 0.9026669734727175, 1e-9, "the oil's n-butane");
    relative(aqueous.z[2], 0.9999999860500847, 1e-12, "the aqueous water");
}

/// **An entrainment moves a share of one phase's moles into another, and the two directions
/// are not symmetric.**
#[test]
fn an_entrainment_moves_moles_between_the_phases() {
    let feed = Stream::from_pt(
        vec!["methane".into(), "n-butane".into(), "water".into()],
        vec![0.5, 0.3, 0.2],
        1.0,
        pascals(20.0e5),
        kelvins(300.0),
    )
    .expect("the three components resolve");
    let plain = azoth_process::kernels::three_phase_separator::three_phase_separator(
        &feed,
        pascals(0.0),
        None,
        Entrainment::default(),
    )
    .expect("the equilibrium split");
    let into_oil = azoth_process::kernels::three_phase_separator::three_phase_separator(
        &feed,
        pascals(0.0),
        None,
        Entrainment {
            gas_in_oil: 0.05,
            ..Entrainment::default()
        },
    )
    .expect("gas into oil");
    let out_of_oil = azoth_process::kernels::three_phase_separator::three_phase_separator(
        &feed,
        pascals(0.0),
        None,
        Entrainment {
            oil_in_gas: 0.05,
            ..Entrainment::default()
        },
    )
    .expect("oil into gas");

    // A twentieth of the vapour's moles move: the vapour keeps its composition, its flow
    // falls by that share, and the oil takes both.
    close(into_oil.0.n, plain.0.n * 0.95);
    // The composition is held to the last bit of the normalisation rather than exactly: the
    // transfer is proportional component by component, so the share that comes back out of
    // `amounts / beta` is the same number up to the division's rounding.
    for i in 0..3 {
        close(into_oil.0.z[i], plain.0.z[i]);
    }
    close(into_oil.1.n, plain.1.n + plain.0.n * 0.05);
    close(into_oil.0.n + into_oil.1.n + into_oil.2.n, feed.n);
    // And the other way: the vapour gains, the oil loses, and the water-rich phase is
    // untouched by either - the fractions are per-pair, not per-vessel.
    close(out_of_oil.0.n, plain.0.n + plain.1.n * 0.05);
    close(out_of_oil.1.n, plain.1.n * 0.95);
    close(out_of_oil.2.n, plain.2.n);
}

/// **An ejector entrains a low-pressure stream with a high-pressure one, and discharges
/// between them.** The invariants are the machine's own: the moles are the two inlets' sum,
/// the outlet is above the suction and below the motive, and a worse nozzle discharges
/// hotter.
#[test]
fn an_ejector_discharges_between_its_inlets() {
    let motive = binary(0.9, 1.0, 30.0e5, 400.0);
    let suction = binary(0.9, 0.5, 5.0e5, 300.0);
    let defaults = EjectorSetup {
        discharge_pressure: pascals(10.0e5),
        motive_nozzle_efficiency: 0.75,
        suction_nozzle_efficiency: 0.90,
        mixing_efficiency: 0.85,
        diffuser_efficiency: 0.80,
    };
    let out =
        azoth_process::kernels::ejector::ejector(&motive, &suction, defaults).expect("ejector");

    println!(
        "outlet n={} p={} T={} h={} z={:?}",
        out.n, out.p.value, out.t.value, out.h.value, out.z
    );
    println!("neqsim n=1.5 p=1000000.0 T=357.55282081845144 h=3446.5843564496936");

    close(out.n, motive.n + suction.n);
    close(out.p.value, 10.0e5);
    close(out.z[0], 0.9);
    assert!(
        out.p.value > suction.p.value && out.p.value < motive.p.value,
        "an ejector discharges between its inlets, and {} Pa is not between {} and {}",
        out.p.value,
        suction.p.value,
        motive.p.value
    );
    relative(out.n, 1.5, 1e-12, "the outlet flow");
    relative(
        out.t.value,
        357.55282081845144,
        1e-4,
        "the outlet temperature",
    );

    // **The parameter the palette entry did not declare moves the answer.** A poorer motive
    // nozzle keeps less of the expansion's enthalpy drop, so the jet is slower, the mixing
    // velocity and the diffuser's recovery fall with it, and the outlet comes back hotter.
    let poorer = azoth_process::kernels::ejector::ejector(
        &motive,
        &suction,
        EjectorSetup {
            motive_nozzle_efficiency: 0.5,
            ..defaults
        },
    )
    .expect("ejector");
    println!("poorer nozzle: T={} h={}", poorer.t.value, poorer.h.value);
    println!("neqsim poorer: T=359.5821295331475 h=3541.935633243748");
    assert!(
        poorer.t.value > out.t.value,
        "a worse nozzle leaves the outlet hotter: {} against {}",
        poorer.t.value,
        out.t.value
    );
    relative(
        poorer.t.value,
        359.5821295331475,
        1e-4,
        "the poorer nozzle's temperature",
    );
}

/// **A flare is a pass-through with a report.** The record does not move; the duty and the
/// emission are the class's own two numbers, and the duty multiplies an energy density at
/// one reference by a volumetric flow at another.
#[test]
fn a_flare_passes_its_record_through_and_reports_two_numbers() {
    let inlet = Stream::from_pt(
        vec!["methane".into(), "n-butane".into()],
        vec![0.9, 0.1],
        1.0,
        pascals(1.01325e5),
        kelvins(288.15),
    )
    .expect("the pair resolves");
    let (product, numbers) =
        azoth_process::kernels::flare::flare(&inlet).expect("the standard has both rows");

    println!(
        "product n={} p={} t={} h={} duty={} co2={}",
        product.n,
        product.p.value,
        product.t.value,
        product.h.value,
        numbers.heat_duty.value,
        numbers.co2_emission.value
    );
    println!("neqsim duty=1046833.0648978995 co2=0.057213");

    // The record is the inlet's, field for field.
    close(product.n, inlet.n);
    assert_eq!(product.z, inlet.z);
    close(product.p.value, inlet.p.value);
    close(product.t.value, inlet.t.value);
    close(product.h.value, inlet.h.value);

    // The carbon is the gas's: 0.9 methane (one atom) and 0.1 n-butane (four) is 1.3 mol/s of
    // carbon, times 44.01e-3 kg/mol.
    close(numbers.co2_emission.value, 1.3 * 44.01e-3);
    relative(
        numbers.heat_duty.value,
        1_046_833.064_897_899_5,
        1e-12,
        "the duty",
    );

    // **The two reference states are not the same one.** The duty is the calorific value per
    // normal cubic metre at 0 C times a volume at 15 C, so it is 288.15/273.15 above the
    // consistently-referenced one - which is the class's arithmetic and not a defect this
    // port may quietly fix.
    let consistently = numbers.heat_duty.value * 273.15 / 288.15;
    close(numbers.heat_duty.value / consistently, 288.15 / 273.15);
    assert!(
        (numbers.heat_duty.value - consistently) / consistently > 0.05,
        "the two references differ by more than five per cent: {} against {}",
        numbers.heat_duty.value,
        consistently
    );
}
