//! A flowsheet that runs.
//!
//! **The acceptance P12 is held to is that `specs/flowsheets/demo.toml` converges its recycle** -
//! a flowsheet that validates but does not run is a declaration, and every item before this one
//! in the tranche was scaffolding for this loop.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use azoth_core::units::{kelvins, pascals};
use azoth_process::ExecutionOrder;
use azoth_process::executor::{STREAM_FIELDS, Session, run};
use azoth_process::recycle::Acceleration;
use azoth_process::{Flowsheet, Stream, load_palette};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn feed() -> Stream {
    Stream::from_pt(
        ["methane", "n-butane"]
            .iter()
            .map(|name| (*name).to_string())
            .collect(),
        vec![0.9, 0.1],
        1.0,
        pascals(5.0e5),
        kelvins(300.0),
    )
    .expect("the two resolve")
}

/// **The shipped flowsheet, run.** Feed, pump, heat, separate, and the liquid recycled back to
/// the mixer - so the second pass's mixer has two inlets and the tear closes.
#[test]
fn the_shipped_flowsheet_converges_its_recycle() {
    let text = std::fs::read_to_string(root().join("specs/flowsheets/demo.toml"))
        .expect("the shipped flowsheet is there");
    let flowsheet: Flowsheet = toml::from_str(&text).expect("it parses");
    let palette = load_palette(&root().join("specs/unit_ops")).expect("the palette loads");
    let feeds = BTreeMap::from([("feed_1".to_string(), feed())]);

    let report =
        run(&flowsheet, &palette, &feeds, ExecutionOrder::Insertion).expect("the flowsheet runs");

    for name in [
        "feed_1",
        "mix1.product",
        "p1.outlet",
        "sep1.liquid",
        "vapour_product",
    ] {
        assert!(
            report.streams.contains_key(name),
            "`{name}` was never produced"
        );
    }
    assert!(report.converged, "the recycle did not converge");
    assert!(report.iterations > 1, "a recycle takes at least two passes");
    assert!(report.tears[0].solved);
}

/// **A tear that carries something, at two tolerances.**
///
/// The shipped flowsheet's recycle is **empty**: at 20 bar and 320 K that feed is all vapour, so
/// `sep1.liquid` is nil on both passes and `solved()`'s zero-flow floor is what ends the loop.
/// The residuals are exactly zero there because there is nothing to compare, not because a loop
/// closed - so this runs the same graph with the heater cold enough to condense.
///
/// **And the class's default tolerance is loose enough that the loop exits on its second pass
/// anyway.** The tear's flow residual is `4.5e-3` - the difference between a loop with no recycle
/// and one with it - and `flowTolerance` is `1e-2`, so it reads as converged. That is measured
/// rather than argued: the same graph at `1e-8` takes more passes and ends at a genuine fixed
/// point, which is what the second half of this test shows.
#[test]
fn a_tear_that_carries_material_converges() {
    let text = std::fs::read_to_string(root().join("specs/flowsheets/demo.toml"))
        .expect("the shipped flowsheet is there");
    let palette = load_palette(&root().join("specs/unit_ops")).expect("the palette loads");
    let feeds = BTreeMap::from([("feed_1".to_string(), feed())]);

    let condensing = || {
        let mut flowsheet: Flowsheet = toml::from_str(&text).expect("it parses");
        // Below n-butane's dew point at the pump's 20 bar, so the separator has a liquid and the
        // tear carries it.
        flowsheet
            .instances
            .iter_mut()
            .find(|instance| instance.id == "hx1")
            .expect("the heater is declared")
            .parameters
            .insert("outlet_temperature".to_string(), toml::Value::Float(250.0));
        flowsheet
    };

    // The class's own defaults.
    let loose = run(&condensing(), &palette, &feeds, ExecutionOrder::Insertion)
        .expect("the flowsheet runs");
    let loose_residuals = loose.tears[0]
        .residuals
        .expect("a second pass measured them");
    let liquid = loose.streams.get("sep1.liquid").expect("a liquid");
    println!(
        "default: iterations={} converged={} liquid={} residuals={loose_residuals:?}",
        loose.iterations, loose.converged, liquid.n
    );
    assert!(
        liquid.n > 0.0,
        "the tear carries nothing, so this is the shipped test again"
    );
    assert!(loose.converged);
    assert!(
        loose_residuals.flow > 0.0,
        "a first-pass residual of exactly zero means nothing was compared"
    );

    // **And the same graph has no steady state at all.** `demo.toml` declares
    // `vapour_product` as its only product, so the liquid it recycles has **no exit**: the vapour
    // leaves carrying butane and the liquid returns in full, so the loop gains a fixed amount
    // every pass. Measured, `sep1.liquid` runs `0.0862, 0.1725, 0.2587, 0.3450, ...` mol/s - an
    // arithmetic progression whose composition is constant to twelve digits. A tight tolerance
    // therefore never closes, and the residual at pass 100 is the residual at pass 2 to fifteen
    // digits.
    //
    // So the two-pass exit above is the *tolerance*, not the loop: `1e-2` accepts the difference
    // between a flowsheet with no recycle and one with it.
    let mut tight = condensing();
    tight.recycles[0].flow_tolerance = Some(1e-8);
    let closed =
        run(&tight, &palette, &feeds, ExecutionOrder::Insertion).expect("the flowsheet runs");
    let closed_residuals = closed.tears[0]
        .residuals
        .expect("a later pass measured them");
    println!(
        "tight:   iterations={} converged={} residuals={closed_residuals:?}",
        closed.iterations, closed.converged
    );
    assert!(
        !closed.converged && closed.iterations == 100,
        "the accumulating flowsheet converged after {} passes, which the liquid's growth says it \
         cannot",
        closed.iterations
    );
    assert!(
        (closed_residuals.flow - loose_residuals.flow).abs() < 1e-12,
        "the residual moved, so the loop is not the fixed cycle the growth describes"
    );

    // **The purge that would give it one is not expressible today, and that is a gap rather than
    // a decision.** A product connection *reports* a stream and does not consume it - every
    // connection naming a stream reads it - so adding one changes nothing (measured: the residual
    // is unchanged to fifteen digits). What the loop needs is a **split**, and `unit_ops.splitter`
    // declares its outlets as one `many` port named `products`, which a connection cannot address
    // one of. `ExecutionOrder` and this loop are both wired for `many` *inlets* and not for `many`
    // outlets, and closing it means giving a connection a way to name the index.
}

/// The shipped flowsheet as written, with no condensing and therefore no liquid in the loop - the
/// document's own state, and the one the two-pass exit above is measured against.
#[test]
fn the_shipped_flowsheet_is_the_one_with_an_empty_tear() {
    let text = std::fs::read_to_string(root().join("specs/flowsheets/demo.toml"))
        .expect("the shipped flowsheet is there");
    let flowsheet: Flowsheet = toml::from_str(&text).expect("it parses");
    let palette = load_palette(&root().join("specs/unit_ops")).expect("the palette loads");
    let feeds = BTreeMap::from([("feed_1".to_string(), feed())]);
    let report =
        run(&flowsheet, &palette, &feeds, ExecutionOrder::Insertion).expect("the flowsheet runs");
    let liquid = report
        .streams
        .get("sep1.liquid")
        .expect("the outlet exists");
    assert_eq!(liquid.n, 0.0, "the shipped state has no liquid");
    assert_eq!(
        report.tears[0].residuals.expect("measured"),
        azoth_process::recycle::Residuals {
            flow: 0.0,
            absolute_flow_change_kg_per_hr: 0.0,
            composition: 0.0,
            temperature: 0.0,
            pressure: 0.0,
        },
        "an empty tear has nothing to compare, which is the zero-flow floor"
    );
}

/// **The session's named results, and the paths they answer to.**
///
/// The point of the session is that a widget and an agent both point at a *string*, so the paths
/// have to be derivable rather than discovered: `<endpoint>.<field>` where the field is one the
/// port declaration names.
#[test]
fn the_session_addresses_every_value_by_a_stable_path() {
    let text = std::fs::read_to_string(root().join("specs/flowsheets/demo.toml"))
        .expect("the shipped flowsheet is there");
    let flowsheet: Flowsheet = toml::from_str(&text).expect("it parses");
    let palette = load_palette(&root().join("specs/unit_ops")).expect("the palette loads");
    let feeds = BTreeMap::from([("feed_1".to_string(), feed())]);

    let session = Session::run(flowsheet, &palette, &feeds, ExecutionOrder::Insertion)
        .expect("the flowsheet runs");

    // Every endpoint the run produced is addressable, and every field of each.
    for endpoint in [
        "feed_1",
        "mix1.product",
        "p1.outlet",
        "hx1.outlet",
        "sep1.liquid",
    ] {
        assert!(
            session.stream(endpoint).is_ok(),
            "`{endpoint}` is not addressable"
        );
        for field in STREAM_FIELDS {
            let path = format!("{endpoint}.{field}");
            if field == "z" {
                assert!(session.value(&path).is_err(), "`{path}` is a vector");
            } else {
                assert!(session.value(&path).is_ok(), "`{path}` is not readable");
            }
        }
    }

    // The values are the stream's own, and the paths agree with the map.
    let outlet = session.stream("p1.outlet").expect("addressable");
    assert_eq!(session.value("p1.outlet.n").expect("readable"), outlet.n);
    assert_eq!(
        session.value("p1.outlet.T").expect("readable"),
        outlet.t.value
    );
    assert_eq!(
        session.value("p1.outlet.P").expect("readable"),
        outlet.p.value
    );

    // A path nobody wrote is refused with a list rather than a `None`.
    let error = session.value("nosuch.outlet.n").expect_err("it refuses");
    assert!(error.to_string().contains("produced"), "{error}");

    // The tear is addressable by its own name, and it is not a stream.
    let tear = session.tear("recycle_1").expect("named");
    assert_eq!(tear.stream, "recycle_1");
    assert!(session.converged());
    assert_eq!(session.iterations(), 2);
    assert!(session.tear("recycle_9").is_err());

    // And the enumeration a front-end needs.
    let paths = session.paths();
    assert!(paths.contains(&"p1.outlet".to_string()));
    assert!(paths.contains(&"p1.outlet.T".to_string()));
    assert!(paths.contains(&"recycle_1".to_string()));
    assert!(paths.windows(2).all(|pair| pair[0] <= pair[1]), "sorted");
}

/// **The one row NeqSim itself is required for.**
///
/// `validation/neqsim/captures/process_flowsheet.tsv` is the first capture in this repository
/// built by a `ProcessSystem`: every other probe drives a unit standalone, because
/// `SimulationInterface.run()` defaults to `run(UUID.randomUUID())`, but a **recycle** needs the
/// sequential loop around it. So this is the only external measurement of the executor, and it
/// holds the loop's *convergence* - the pass it stopped on and the residual it stopped with - as
/// well as its destination.
///
/// The row is `demo_condensing`: the shipped graph with the heater at 250 K, which is the one
/// whose tear carries material. At 320 K the separator's liquid is nil and NeqSim **deactivates**
/// the recycle outright (see `a_tear_that_carries_material_converges`), so the empty row measures
/// a path this port does not carry.
#[test]
fn the_executor_reproduces_the_captured_convergence() {
    let text = std::fs::read_to_string(root().join("specs/flowsheets/demo.toml"))
        .expect("the shipped flowsheet is there");
    let mut flowsheet: Flowsheet = toml::from_str(&text).expect("it parses");
    let palette = load_palette(&root().join("specs/unit_ops")).expect("the palette loads");
    flowsheet
        .instances
        .iter_mut()
        .find(|instance| instance.id == "hx1")
        .expect("the heater is declared")
        .parameters
        .insert("outlet_temperature".to_string(), toml::Value::Float(250.0));

    let feeds = BTreeMap::from([("feed_1".to_string(), feed())]);
    let session = Session::run(flowsheet, &palette, &feeds, ExecutionOrder::Insertion)
        .expect("the flowsheet runs");

    let capture =
        std::fs::read_to_string(root().join("validation/neqsim/captures/process_flowsheet.tsv"))
            .expect("the capture is there");
    let row = capture_row(&capture, "demo_condensing");

    // **The destination.** The separator's liquid, in mol/s and in composition.
    let liquid = session
        .stream("sep1.liquid")
        .expect("the separator produced one");
    let expected_n = number(&row, "sep1_liquid_n");
    assert!(
        (liquid.n - expected_n).abs() / expected_n < 1e-9,
        "the liquid is {} mol/s and NeqSim's is {expected_n}",
        liquid.n
    );
    let expected_z = composition(row.get("sep1_liquid_z").expect("the composition is there"));
    for (i, want) in expected_z.iter().enumerate() {
        assert!(
            (liquid.z[i] - want).abs() < 1e-12,
            "component {i}: {} against {want}",
            liquid.z[i]
        );
    }

    // **And the convergence, which is what makes this an executor's oracle rather than a
    // kernel's.** A tear that lands in the same place by a different number of passes is a
    // different machine.
    let tear = session.tear("recycle_1").expect("the recycle is named");
    assert_eq!(
        tear.iterations,
        number(&row, "recycle_iterations") as u32,
        "the port stopped on a different pass"
    );
    let residuals = tear.residuals.expect("measured");
    let expected_flow = number(&row, "recycle_error_flow");
    assert!(
        (residuals.flow - expected_flow).abs() / expected_flow < 1e-9,
        "the flow residual is {} and NeqSim's is {expected_flow}",
        residuals.flow
    );
    // The composition residual is a sum of `~1e-14` differences - a rounding-level quantity, so
    // it is bounded rather than compared.
    let expected_composition = number(&row, "recycle_error_composition");
    assert!(
        residuals.composition < expected_composition.max(1e-12),
        "the composition residual is {} against NeqSim's {expected_composition}",
        residuals.composition
    );
    assert!(tear.solved, "the port did not report the tear solved");
    assert!(session.converged());
}

/// One row of the capture's `key=value` lines, selected by its label.
///
/// The rows are separated by a blank line and the first line of each is its label.
fn capture_row(capture: &str, label: &str) -> BTreeMap<String, String> {
    let block = capture
        .split("\n\n")
        .find(|block| block.lines().next() == Some(label))
        .unwrap_or_else(|| panic!("the capture has no `{label}` row"));
    block
        .lines()
        .filter_map(|line| line.split_once('='))
        .map(|(name, value)| (name.to_string(), value.to_string()))
        .collect()
}

/// One field of a row as a number.
fn number(row: &BTreeMap<String, String>, key: &str) -> f64 {
    row.get(key)
        .unwrap_or_else(|| panic!("the row has no `{key}`"))
        .parse()
        .unwrap_or_else(|_| panic!("`{key}` is not a number"))
}

/// A `name:value name:value` composition.
fn composition(field: &str) -> Vec<f64> {
    field
        .split_whitespace()
        .map(|pair| {
            pair.split_once(':')
                .expect("a composition is `name:value`")
                .1
                .parse()
                .expect("a mole fraction is a number")
        })
        .collect()
}

/// **The round trip, over every flowsheet the repository ships.**
///
/// This is `Azoth.Rho`'s `*@P ≅ P` in code: reading a written flowsheet gives the value that was
/// written. Two directions are checked, because one alone passes on a writer that loses
/// something:
///
/// * the value comes back **unchanged** - `from_toml(to_toml(f)) == f`;
/// * writing the re-read value again is **byte-identical** to the first write - `to_toml` is a
///   fixed point of itself. A writer that dropped a field the reader defaulted would satisfy the
///   first and fail this one.
///
/// **The output is not the input's bytes, and the test does not claim it is.** TOML has no
/// comments to round-trip and this schema writes its fields in the struct's order, so the shipped
/// file comes back tidied - `demo.toml`'s blank lines and its `[instances.parameters]` placement
/// both change. What is claimed is the value.
#[test]
fn every_shipped_flowsheet_round_trips() {
    let directory = root().join("specs/flowsheets");
    let mut checked = 0;
    for entry in std::fs::read_dir(&directory).expect("the flowsheets are there") {
        let path = entry.expect("readable").path();
        if path.extension().and_then(|e| e.to_str()) != Some("toml") {
            continue;
        }
        let text = std::fs::read_to_string(&path).expect("readable");
        let original = Flowsheet::from_toml(&text)
            .unwrap_or_else(|error| panic!("{}: {error}", path.display()));

        let written = original
            .to_toml()
            .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
        let reparsed = Flowsheet::from_toml(&written).unwrap_or_else(|error| {
            panic!(
                "{}: the written form does not read back: {error}\n{written}",
                path.display()
            )
        });

        assert_eq!(
            original,
            reparsed,
            "{}: the value changed across the round trip\n{written}",
            path.display()
        );
        assert_eq!(
            written,
            reparsed.to_toml().expect("it writes again"),
            "{}: writing is not a fixed point of itself",
            path.display()
        );
        checked += 1;
    }
    assert!(
        checked > 0,
        "no flowsheet was checked, so this test decided nothing"
    );
}

/// **A declaration that states no convergence parameter writes none of them**, and one that states
/// some writes only those.
///
/// This is the half of the round trip the shipped file cannot exercise: `demo.toml` states no
/// `[[recycles]]` settings at all, so a writer that emitted the defaults would round-trip
/// *semantically* - the reader puts them back - and would still be wrong, because it would have
/// turned a declaration that takes the class's defaults into one that pins them.
#[test]
fn an_unstated_recycle_parameter_is_not_written() {
    let text = std::fs::read_to_string(root().join("specs/flowsheets/demo.toml")).expect("there");
    let flowsheet = Flowsheet::from_toml(&text).expect("it parses");
    let written = flowsheet.to_toml().expect("it writes");
    assert!(
        !written.contains("flow_tolerance"),
        "an unstated tolerance was written out:\n{written}"
    );
    assert!(!written.contains("acceleration_method"), "{written}");

    // And the ones that *are* stated survive - every setting the schema gained, since that is the
    // surface this file exists to hold.
    let mut stated = flowsheet.clone();
    stated.recycles[0].flow_tolerance = Some(1e-6);
    stated.recycles[0].composition_tolerance = Some(2e-6);
    stated.recycles[0].temperature_tolerance = Some(3e-6);
    stated.recycles[0].pressure_tolerance = Some(4e-6);
    stated.recycles[0].max_iterations = Some(42);
    stated.recycles[0].minimum_flow = Some(1e-9);
    stated.recycles[0].acceleration_method = Some("wegstein".to_string());
    let written = stated.to_toml().expect("it writes");
    let round = Flowsheet::from_toml(&written).expect("it reads");
    assert_eq!(
        round, stated,
        "a stated setting did not survive:\n{written}"
    );

    // **And the reader resolves them the way the class does**, so the round trip preserves the
    // machine and not only the document.
    let settings = round.recycles[0].settings().expect("they resolve");
    assert_eq!(settings.flow_tolerance, 1e-6);
    assert_eq!(settings.max_iterations, 42);
    assert_eq!(settings.acceleration, Acceleration::Wegstein);
}
