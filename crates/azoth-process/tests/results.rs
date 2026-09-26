//! What a run publishes beside its streams.
//!
//! **A kernel that ran a flash, a column solve or a reactor integration reached numbers that are
//! on no outlet stream**, and until the executor carried them a front end could draw a flowsheet
//! and not read one - which is what the *editor* tranche found: the duty was computed and dropped
//! at the call site. `KernelOutcome` is what carries it, and this file is the gate over it.
//!
//! **Two things are checked, and they are different in kind.** The first is a table against the
//! dispatch table in both directions, so an entry that starts publishing without a row here - or a
//! row whose entry stops - fails rather than passing quietly. The second runs the *shipped*
//! flowsheet and reads the wire, because a table of field names agrees with itself whether or not
//! a number ever crosses.
//!
//! **What is deliberately not here**: one row per unit operation. A unit operation's own numbers
//! are pinned by its case (`specs/models/process/*.toml`) and by the kernel tests beside this
//! file; what this file holds is that the *run* publishes them and names them what the model does.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use azoth_core::CalcResult;
use azoth_process::executor::json::{SessionReport, session_report};
use azoth_process::executor::{DISPATCH, Session};
use azoth_process::models;
use azoth_process::{ExecutionOrder, Flowsheet, load_palette};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// **What each dispatched entry publishes**, as the model's own `FIELDS`.
///
/// `None` is a real answer and not a gap: a mixer's result *is* its mixed stream and a separator's
/// *are* its two, so a run that published them here would be the same numbers under a second key.
/// `python/src/azoth/process/layers.py`'s `NO_INTERIOR` table names five of the same entries that
/// way - that their kernel forms nothing their records do not carry.
fn published() -> Vec<(&'static str, Option<&'static [&'static str]>)> {
    vec![
        (
            "unit_ops.absorption_column",
            Some(models::AbsorptionColumnResult::FIELDS),
        ),
        ("unit_ops.component_splitter", None),
        ("unit_ops.compressor", None),
        ("unit_ops.cooler", Some(models::CoolerResult::FIELDS)),
        (
            "unit_ops.distillation_column",
            Some(models::DistillationColumnResult::FIELDS),
        ),
        ("unit_ops.ejector", None),
        ("unit_ops.expander", None),
        ("unit_ops.filter", Some(models::FilterResult::FIELDS)),
        ("unit_ops.flare", Some(models::FlareResult::FIELDS)),
        ("unit_ops.gas_scrubber", None),
        (
            "unit_ops.gibbs_reactor",
            Some(models::GibbsReactorResult::FIELDS),
        ),
        ("unit_ops.heat_exchanger", None),
        ("unit_ops.heater", Some(models::HeaterResult::FIELDS)),
        ("unit_ops.manifold", None),
        ("unit_ops.mixer", None),
        ("unit_ops.pipe", Some(models::PipeResult::FIELDS)),
        (
            "unit_ops.plug_flow_reactor",
            Some(models::PlugFlowReactorResult::FIELDS),
        ),
        // **The packed column's entry became runnable**, so it appears here for the first time:
        // it dispatches the base column's solve, and its record is not on the wire either - the
        // same statement as its sibling's below.
        ("unit_ops.packed_column", None),
        ("unit_ops.pump", None),
        // **The rate-based packed column reaches an interior and publishes none of it**, which is
        // why its row is `None` rather than a `FIELDS`. What is missing is the step from its
        // kernel's `RateBasedOutcome` to its model's flat record - the model *does* declare the
        // segment profile as `[outputs]`, and builds them from its own flat inputs today. So the
        // row is a statement about the wire rather than about the arithmetic, and the day that
        // constructor lands is the day this row carries a `FIELDS` beside it.
        ("unit_ops.rate_based_packed_column", None),
        ("unit_ops.separator", None),
        (
            "unit_ops.shortcut_distillation_column",
            Some(models::ShortcutDistillationColumnResult::FIELDS),
        ),
        ("unit_ops.splitter", None),
        (
            "unit_ops.stirred_tank_reactor",
            Some(models::StirredTankReactorResult::FIELDS),
        ),
        (
            "unit_ops.stripping_column",
            Some(models::StrippingColumnResult::FIELDS),
        ),
        ("unit_ops.tank", None),
        ("unit_ops.three_phase_separator", None),
        ("unit_ops.throttling_valve", None),
    ]
}

/// **The table above and the dispatch table are the same table.** Both directions, because either
/// one alone leaves a hole: an entry added to `DISPATCH` with no row here, and a row here for an
/// entry that no longer dispatches.
#[test]
fn the_published_table_covers_the_dispatch_table_both_ways() {
    let table = published();
    let named: Vec<&str> = table.iter().map(|(id, _)| *id).collect();
    let dispatched: Vec<&str> = DISPATCH.iter().map(|(id, _)| *id).collect();
    for id in &dispatched {
        assert!(
            named.contains(id),
            "`{id}` dispatches and this file does not say what it publishes"
        );
    }
    for id in &named {
        assert!(
            dispatched.contains(id),
            "this file says what `{id}` publishes and it does not dispatch"
        );
    }
    assert_eq!(table.len(), dispatched.len());
}

/// **The shipped flowsheet, on the wire.**
///
/// The heater's duty is the number no outlet stream carries, and it is the one that made this
/// work necessary: `dispatch.rs` read `Ok(vec![out.outlet])` and the duty went out of scope. The
/// mixer beside it is the counter-case - it publishes nothing - so an implementation that
/// published a restatement of every stream record would fail here.
#[test]
fn the_shipped_demo_publishes_the_duty_no_stream_carries() {
    let text = std::fs::read_to_string(root().join("specs/flowsheets/demo.toml"))
        .expect("the shipped flowsheet is there");
    let flowsheet: Flowsheet = toml::from_str(&text).expect("it parses");
    let palette = load_palette(&root().join("specs/unit_ops")).expect("the palette loads");
    let session = Session::run(
        flowsheet,
        &palette,
        &BTreeMap::new(),
        ExecutionOrder::Insertion,
    )
    .expect("the flowsheet runs");
    let report: SessionReport = session_report(&session).expect("it writes");

    let heater = report
        .results
        .get("hx1")
        .expect("the heater reached a duty and it crosses");
    assert_eq!(heater["outlet_duty"]["unit"], "W");
    let duty = heater["outlet_duty"]["magnitude_si"]
        .as_f64()
        .expect("a magnitude");
    assert!(duty.is_finite() && duty != 0.0, "the duty is {duty}");

    // **Every key is one the model declares**, so the wire cannot grow a spelling of its own.
    let keys: Vec<&str> = heater
        .as_object()
        .expect("an object")
        .keys()
        .map(String::as_str)
        .collect();
    let mut declared: Vec<&str> = models::HeaterResult::FIELDS.to_vec();
    let mut found: Vec<&str> = keys.clone();
    declared.sort_unstable();
    found.sort_unstable();
    assert_eq!(found, declared);

    // A mixer's whole answer is its stream, so it publishes nothing rather than a restatement.
    assert!(
        !report.results.contains_key("mix1"),
        "the mixer published {:?}",
        report.results.get("mix1")
    );

    // **And the record carries what the run knew and the five declared fields do not.** The
    // separator's vapour fraction is the flash's own beta, so a null here would mean the
    // constructor threw the flash's answer away again.
    let vapour = &report.streams["sep1.vapour"];
    assert!(vapour.vapour_fraction.is_some(), "the flash answered");
    assert!(vapour.mass_flow.is_some(), "two resolved components weigh");
    assert_eq!(
        vapour.mass_flow.as_ref().expect("weighed").unit,
        "kg/s",
        "the mass flow crosses in the vocabulary's own unit"
    );
    assert!(vapour.molar_mass.is_some());
}

/// **A result is the last pass's, and it agrees with the stream the same pass published.**
///
/// A result accumulated once and not overwritten would be the *first* pass's, and this flowsheet
/// has a recycle: the mixer is fed the tear, so the first pass runs the heater on a feed the loop
/// has not closed yet. Two numbers from two different passes is the defect, and the only way to
/// see it is to hold the published result against the published stream.
#[test]
fn a_result_is_the_same_passs_as_the_stream_beside_it() {
    let text = std::fs::read_to_string(root().join("specs/flowsheets/demo.toml"))
        .expect("the shipped flowsheet is there");
    let flowsheet: Flowsheet = toml::from_str(&text).expect("it parses");
    let palette = load_palette(&root().join("specs/unit_ops")).expect("the palette loads");
    let session = Session::run(
        flowsheet,
        &palette,
        &BTreeMap::new(),
        ExecutionOrder::Insertion,
    )
    .expect("the flowsheet runs");
    let report: SessionReport = session_report(&session).expect("it writes");
    assert!(report.iterations > 1, "a recycle takes at least two passes");

    let heater = &report.results["hx1"];
    let outlet = &report.streams["hx1.outlet"];
    // **The two shapes differ and both are deliberate.** A *record* writes every scalar as a
    // magnitude and a unit; a model *result* writes the fields its Python dataclass writes, and a
    // flow there is a bare number - which `FIELDS` and the cross-implementation test are what hold
    // the two halves to. So the comparison meets at the magnitude.
    assert_eq!(
        heater["outlet_n"].as_f64(),
        Some(outlet.n.magnitude_si),
        "the heater's flow is not the same pass's stream value"
    );
    for (field, magnitude) in [
        ("outlet_p", outlet.p.magnitude_si),
        ("outlet_t", outlet.temperature.magnitude_si),
    ] {
        assert_eq!(
            heater[field]["magnitude_si"].as_f64(),
            Some(magnitude),
            "the heater's `{field}` is not the same pass's stream value"
        );
    }
}
