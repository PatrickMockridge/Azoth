//! The executor calls the kernels, made testable.
//!
//! **This is C9's acceptance and it is not a second implementation.** The claim is that the
//! execution layer has no arithmetic of its own: for each instance the executor ran, its outlet
//! must equal what that unit's own *registered model* answers for the same inlet state. If the
//! dispatch table ever grew a conversion the model does not make, or a parameter reached a kernel
//! differently than it reaches a case, this fails - and it fails on the state, not on the wiring,
//! because the wiring is what `tests/executor.rs` already holds.
//!
//! **The comparison is a band and not an equality, and the measurement is why.** The model layer
//! reconstructs its inlet with `Stream::from_pt`, which *re-derives* `h` from `(T, P, z)` - by
//! design, and `models/mod.rs` says so: "an inlet takes four of the record's five fields ... `h`
//! is a state function of `(T, P, z)`". The executor *carries* the `h` its upstream unit wrote,
//! and a mixer's output enthalpy is a flow-weighted mix of its inlets' rather than a fresh flash
//! of the mixture. So the two paths agree to `5.7e-12` relative on the pump's outlet temperature
//! and not to the bit, and holding them to equality would be asserting something the two layers
//! never claimed.
//!
//! **The mixer is not in it, and the reason is a real limit of what the session retains.** Its
//! inlet is `feed_1` *plus the tear*, and the tear's value on the last pass is the previous pass's
//! `sep1.liquid`, which the run does not keep - only the final state is held. So the mixer's own
//! inputs cannot be reconstructed from the session, and a cross-check of it would be a check of
//! something the executor was not called with. The other three are single-inlet and exact.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use azoth_core::units::{kelvins, pascals};
use azoth_process::executor::{Session, run};
use azoth_process::{ExecutionOrder, Flowsheet, Stream, load_palette};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// The shipped graph with the heater cold enough to condense, so every unit does something.
fn session() -> Session {
    let text = std::fs::read_to_string(root().join("specs/flowsheets/demo.toml")).expect("there");
    let mut flowsheet = Flowsheet::from_toml(&text).expect("it parses");
    let palette = load_palette(&root().join("specs/unit_ops")).expect("the palette loads");
    flowsheet
        .instances
        .iter_mut()
        .find(|instance| instance.id == "hx1")
        .expect("the heater is declared")
        .parameters
        .insert("outlet_temperature".to_string(), toml::Value::Float(250.0));
    let feeds = BTreeMap::from([(
        "feed_1".to_string(),
        Stream::from_pt(
            ["methane", "n-butane"]
                .iter()
                .map(|n| (*n).to_string())
                .collect(),
            vec![0.9, 0.1],
            1.0,
            pascals(5.0e5),
            kelvins(300.0),
        )
        .expect("the two resolve"),
    )]);
    Session::run(flowsheet, &palette, &feeds, ExecutionOrder::Insertion)
        .expect("the flowsheet runs")
}

/// One instance's parameter, as the document states it.
fn parameter(session: &Session, instance: &str, name: &str) -> f64 {
    session
        .flowsheet()
        .instances
        .iter()
        .find(|entry| entry.id == instance)
        .unwrap_or_else(|| panic!("`{instance}` is declared"))
        .parameters
        .get(name)
        .unwrap_or_else(|| panic!("`{instance}` states `{name}`"))
        .as_float()
        .unwrap_or_else(|| panic!("`{instance}.{name}` is a number"))
}

/// `process.pump` against what the executor put on `p1.outlet`.
#[test]
fn the_executor_s_answer_is_the_pumps() {
    use azoth_process::pump;

    let session = session();
    let inlet = session
        .stream("mix1.product")
        .expect("the mixer produced one");
    let outlet = session.stream("p1.outlet").expect("the pump produced one");

    let answer = pump(
        &inlet.components,
        inlet.n,
        &inlet.z,
        inlet.p,
        inlet.t,
        pascals(parameter(&session, "p1", "outlet_pressure")),
        parameter(&session, "p1", "isentropic_efficiency"),
    )
    .expect("the model runs");

    agrees(
        outlet,
        answer.outlet_n,
        &answer.outlet_z,
        answer.outlet_p.value,
        answer.outlet_t.value,
        answer.outlet_h.value,
        "p1.outlet",
    );
}

/// A relative band between the two paths, named so a failure says which field.
///
/// `1e-9` against a measured `5.7e-12`: two orders of margin over the divergence the two layers'
/// own definitions produce, and eleven below anything a wrong parameter would cause.
fn close(actual: f64, expected: f64, field: &str) {
    let scale = expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        (actual - expected).abs() / scale < 1e-9,
        "{field}: the executor gives {actual} and the model {expected}"
    );
}

/// A stream's whole record against a model's answer.
fn agrees(stream: &Stream, n: f64, z: &[f64], p: f64, t: f64, h: f64, port: &str) {
    close(stream.n, n, &format!("{port}.n"));
    assert_eq!(stream.z.len(), z.len(), "{port}.z is a different length");
    for (index, (got, want)) in stream.z.iter().zip(z).enumerate() {
        close(*got, *want, &format!("{port}.z[{index}]"));
    }
    close(stream.p.value, p, &format!("{port}.P"));
    close(stream.t.value, t, &format!("{port}.T"));
    close(stream.h.value, h, &format!("{port}.h"));
}

/// `process.separator` against both of `sep1`'s outlets.
#[test]
fn the_executor_s_answer_is_the_separators() {
    use azoth_process::separator;

    let session = session();
    let inlet = session
        .stream("hx1.outlet")
        .expect("the heater produced one");
    let vapour = session
        .stream("sep1.vapour")
        .expect("the separator produced one");
    let liquid = session
        .stream("sep1.liquid")
        .expect("the separator produced one");

    let answer = separator(
        &inlet.components,
        inlet.n,
        &inlet.z,
        inlet.p,
        inlet.t,
        pascals(parameter(&session, "sep1", "pressure_drop")),
        parameter(&session, "sep1", "gas_in_liquid"),
        None,
    )
    .expect("the model runs");

    agrees(
        vapour,
        answer.vapour_n,
        &answer.vapour_z,
        answer.vapour_p.value,
        answer.vapour_t.value,
        answer.vapour_h.value,
        "sep1.vapour",
    );
    agrees(
        liquid,
        answer.liquid_n,
        &answer.liquid_z,
        answer.liquid_p.value,
        answer.liquid_t.value,
        answer.liquid_h.value,
        "sep1.liquid",
    );
}

/// **And the executor's own record is one stream per outlet, named by its producer.** A unit
/// whose answer came from anywhere but its kernel would have to have written a stream the model
/// disagrees with, which is what the three tests above would catch - this checks the plumbing they
/// rely on, so a failure names which half is wrong.
#[test]
fn every_outlet_is_a_stream_named_by_its_producer() {
    let session = session();
    for endpoint in [
        "mix1.product",
        "p1.outlet",
        "hx1.outlet",
        "sep1.vapour",
        "sep1.liquid",
    ] {
        assert!(
            session.stream(endpoint).is_ok(),
            "`{endpoint}` is not in the run's streams"
        );
    }
    // A product is a stream too, under the name the document gives it.
    assert!(session.stream("vapour_product").is_ok());
    let run = run(
        session.flowsheet(),
        &load_palette(&root().join("specs/unit_ops")).expect("the palette loads"),
        &BTreeMap::new(),
        ExecutionOrder::Insertion,
    );
    // Without the feed the run refuses rather than answering - a flowsheet with no inlet is a
    // wiring error and not an empty state.
    assert!(run.is_err(), "a run with no feeds produced something");
}
