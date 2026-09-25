//! A flowsheet that runs.
//!
//! **The acceptance P12 is held to is that `specs/flowsheets/demo.toml` converges its recycle** -
//! a flowsheet that validates but does not run is a declaration, and every item before this one
//! in the tranche was scaffolding for this loop.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use azoth_core::units::{kelvins, pascals};
use azoth_process::ExecutionOrder;
use azoth_process::executor::run;
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
